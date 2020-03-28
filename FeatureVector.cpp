#include "FeatureVector.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <iostream>

#include <omp.h>
#include <assert.h>
#include <map>
#include <algorithm>
#include "LQSolver.h"

#ifdef DUSE_OPENMP
#define DOMP_END \
}
#else
#define DOMP_END ;
#endif
using namespace std;

double RefMesh::normalScale = 0.3;
const double EPS = 1e-10;

//double RefMesh::normalScale = 0.08;
int RefMesh::root = 0;

RefMesh::~RefMesh() {
    for (int i=0; i<align.size(); i++)
        delete align[i];
}

double cotan(Eigen::Vector3d a, Eigen::Vector3d b) {
    double na = a.norm(), nb = b.norm();
    if (na<Eps || nb<Eps) return 0;
    double cos = a.dot(b) / (na*nb);
    if (cos == 1) return 1;
    return cos / sqrt(1-cos*cos);
}

double tan2(Eigen::Vector3d a, Eigen::Vector3d b) {
    double na = a.norm(), nb = b.norm();
    if (na<Eps || nb<Eps) return 0;
    double cos = a.dot(b) / (na*nb);
    double theta = acos(cos)/2;
    double ans = tan(theta);
    if (ans>=0 && ans<=100) return ans / nb;
    return 1 / nb;
}

double SinValue(const OpenMesh::Vec3d& a, const OpenMesh::Vec3d& b, const OpenMesh::Vec3d& c)
{
	double lab = (b - a).length();
	double lac = (c - a).length();
	return (cross(b-a, c-a)).length()/(lab*lac);
}

double CosValue(const OpenMesh::Vec3d& a, const OpenMesh::Vec3d& b, const OpenMesh::Vec3d& c)
{
	double lab = (b - a).length();
	double lac = (c - a).length();
	double lab2 = (b - a).sqrnorm();
	double lac2 = (c-a).sqrnorm();
	double lbc2 = (b - c).sqrnorm();
	return (lab2+lac2-lbc2)/(2.0*lab*lac);
}

double sexp(double x) {
    if (x <= 0) return exp(x);
    return 1+x;
}

RefMesh::RefMesh()
{
   usec = false, fvNum = 0;
   vvs = 0;
   mesh = NULL;
}

RefMesh::RefMesh(DTriMesh &ms) : mesh(&ms) {
    usec = false, fvNum = 0;
	//// Add vertex normals as default property (ref. previous tutorial)
	ms.request_vertex_normals();
	//// Add face normals as default property
	ms.request_face_normals();

	//ms.update_face_normals();
	//ms.update_vertex_normals();

	ms.update_normals();
	//ms.release_face_normals();
    align.resize(mesh->n_vertices());
    d.resize(mesh->n_vertices());
    rd.resize(mesh->n_vertices());
	this->degree.resize(mesh->n_vertices(),0);
    vvs=0;
    for (int i=0; i<mesh->n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        DTriMesh::Point p = mesh->point(vi);
        std::vector<Eigen::Vector3d> vec;
        d[i] = vvs;
        double lens = 0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            rd[vvi.handle().idx()].push_back(std::make_pair(i, vvs));
            Eigen::Vector3d q = OtoE(mesh->point(vvi) - p);
            vec.push_back(q);
            vvs++;
            lens += q.norm();
			degree[i]++;
        }
        vvs = d[i];
#ifdef chaohui
		OpenMesh::Vec3d currentPosition = mesh->point(vi);
		vector<OpenMesh::Vec3d> neighborPoints; 
		for (DTriMesh::VertexVertexIter vv_it = mesh->vv_iter(vi);  vv_it;  ++vv_it)
		{
			neighborPoints.push_back(mesh->point(vv_it));
		}
		int numNeighbors= (int)neighborPoints.size();
		vector<double> weights(numNeighbors,0);
		if (i==588)
		{
			int a=1;
		}
		//std::cout<<i;
		for (int neighborcounts=0; neighborcounts<numNeighbors; neighborcounts++)
		{
			//precomput pi - pj of reference Mesh (yangjie add)
			OpenMesh::Vec3d edgejk_openmesh = p - neighborPoints[neighborcounts];
			Eigen::Vector3d edgejk(edgejk_openmesh[0], edgejk_openmesh[1], edgejk_openmesh[2]);
			this->edgeLength.push_back(edgejk);
			double w1 = CosValue(neighborPoints[(neighborcounts + numNeighbors -1)%numNeighbors], currentPosition, neighborPoints[neighborcounts])/std::max(SinValue(neighborPoints[(neighborcounts + numNeighbors -1)%numNeighbors], currentPosition, neighborPoints[neighborcounts]), EPS);
			double w2 = CosValue(neighborPoints[(neighborcounts+1)%numNeighbors], currentPosition, neighborPoints[neighborcounts])/std::max(SinValue(neighborPoints[(neighborcounts+1)%numNeighbors], currentPosition, neighborPoints[neighborcounts]), EPS);
			weights[neighborcounts] = 0.5*(w1+w2);
			weights[neighborcounts] = std::max(weights[neighborcounts], EPS);
				//weights[i] = 1.0
			if ((weights[neighborcounts] != weights[neighborcounts])||(weights[neighborcounts]>100))
 				{
 					//cout<<v_it.handle().idx()<<" "<<weights[i]<<endl;
 					//system("pause");
				weights[neighborcounts] = 1;
 				}
#ifdef uniformweight
			weights[neighborcounts]=1;
#endif
			w.push_back(weights[neighborcounts]);
			    vvs++;
		}
		//std::cout<<i;

#else

        double ss = 0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            int j = vvs - d[i];
            int pre = j==0 ? (int)vec.size()-1 : j-1;
            int next = j+1==vec.size() ? 0 : j+1;
/*
            w.push_back(sqrt( std::max(Eps, 0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
*/
#ifndef DUSE_TAN
            w.push_back( sqrt( sexp( 0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
            /*
            if (vec[pre].dot(vec[j]) == vec[pre].norm() * (vec[pre]-vec[j]).norm()) {
                cout << "pl i : " << i << " j : " << j << endl;
            }
            */
#else

            w.push_back(tan2(vec[pre],vec[j]) +
                        tan2(vec[next],vec[j]));
#endif
            //w.push_back(1);

            /*
            w.push_back( sqrt( abs (0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
            */

            ss += w[w.size()-1];
            vvs++;
        }
        for (int j=d[i]; j<w.size(); j++)
            w[j] = w[j] / ss;
#endif
#ifdef NEWWEIGHT
       assert(vec.size()==(w.size()-d[i]));
	   for(int ii = 0; ii < vec.size(); ii++)
	   {
		   vec[ii] = sqrt(w[ii+d[i]])*vec[ii];
			//assert(vec[ii] == vec[ii]);
	   }
#endif
        // add normal
		if (mesh->normal(vi) != mesh->normal(vi)){
			std::cout << mesh->normal(vi) << " " << i << " " << w[i];
		}
		//assert(mesh->normal(vi) == mesh->normal(vi));
        vec.push_back( OtoE(mesh->normal(vi) * (lens / vec.size() * RefMesh::normalScale)) );
        align[i] = new AffineAlign(vec);
    }
    c.resize(w.size(), 0);
    std::vector<bool> visit(mesh->n_vertices(), false);
    int qBeg = 0;
    for (int i=root; i<visit.size(); i++) if (!visit[i]) {
        visit[i] = true;
        bfsq.push_back(make_pair(i,make_pair(-1,-1) ));
        while (qBeg < bfsq.size()) {
            int i = bfsq[qBeg++].first;

            DTriMesh::VertexHandle vi(i);
            int vvs = d[i];

            for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++, vvs++) {
                int j = vvi.handle().idx();
                if (!visit[j]) {
                    visit[j] = true;
                    bfsq.push_back(make_pair( j, make_pair(i,vvs) ));
                }
            }
        }
    }
}

void RefMesh::SetMesh(DTriMesh &ms)
{
    this->mesh = &ms;
    usec = false, fvNum = 0;
	ms.request_vertex_normals();
	//// Add face normals as default property
	ms.request_face_normals();

	//ms.update_face_normals();
	//ms.update_vertex_normals();

	ms.update_normals();
    //ms.update_face_normals();
  //  ms.update_vertex_normals();
    align.resize(mesh->n_vertices());
    d.resize(mesh->n_vertices());
    rd.resize(mesh->n_vertices());
	this->degree.resize(mesh->n_vertices(),0);
    vvs=0;
    for (int i=0; i<mesh->n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        DTriMesh::Point p = mesh->point(vi);
        std::vector<Eigen::Vector3d> vec;
        d[i] = vvs;
        double lens = 0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            rd[vvi.handle().idx()].push_back(std::make_pair(i, vvs));
            Eigen::Vector3d q = OtoE(mesh->point(vvi) - p);
            vec.push_back(q);
            vvs++;
            lens += q.norm();
			degree[i]++;
        }
        vvs = d[i];
#ifdef chaohui
		OpenMesh::Vec3d currentPosition = mesh->point(vi);
		vector<OpenMesh::Vec3d> neighborPoints; 
		for (DTriMesh::VertexVertexIter vv_it = mesh->vv_iter(vi);  vv_it;  ++vv_it)
		{
			neighborPoints.push_back(mesh->point(vv_it));
		}
		int numNeighbors= (int)neighborPoints.size();
		vector<double> weights(numNeighbors,0);
		for (int i = 0; i < numNeighbors; i++)
		{
			//precomput pi - pj of reference Mesh (yangjie add)
			OpenMesh::Vec3d edgejk_openmesh = p - neighborPoints[i];
			Eigen::Vector3d edgejk(edgejk_openmesh[0], edgejk_openmesh[1], edgejk_openmesh[2]);
			this->edgeLength.push_back(edgejk);
			double w1 = CosValue(neighborPoints[(i + numNeighbors - 1) % numNeighbors], currentPosition, neighborPoints[i]) / std::max(SinValue(neighborPoints[(i + numNeighbors - 1) % numNeighbors], currentPosition, neighborPoints[i]), EPS);
			double w2 = CosValue(neighborPoints[(i + 1) % numNeighbors], currentPosition, neighborPoints[i]) / std::max(SinValue(neighborPoints[(i + 1) % numNeighbors], currentPosition, neighborPoints[i]), EPS);
			weights[i] = 0.5*(w1 + w2);
			weights[i] = std::max(weights[i], EPS);
			//weights[i] = 1.0
			if ((weights[i] != weights[i])||(weights[i] > 100))
			{
				//cout<<v_it.handle().idx()<<" "<<weights[i]<<endl;
				//system("pause");
				weights[i] = 1;
			}

#ifdef uniformweight
			weights[i] = 1;
#endif

			w.push_back(weights[i]);
			vvs++;
		}
#else
        double ss = 0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            int j = vvs - d[i];
            int pre = j==0 ? (int)vec.size()-1 : j-1;
            int next = j+1==vec.size() ? 0 : j+1;
/*
            w.push_back(sqrt( std::max(Eps, 0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
*/
#ifndef DUSE_TAN
            w.push_back( sqrt( sexp( 0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
            /*
            if (vec[pre].dot(vec[j]) == vec[pre].norm() * (vec[pre]-vec[j]).norm()) {
                cout << "pl i : " << i << " j : " << j << endl;
            }
            */
#else

            w.push_back(tan2(vec[pre],vec[j]) +
                        tan2(vec[next],vec[j]));
#endif
            //w.push_back(1);

            /*
            w.push_back( sqrt( abs (0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
            */

            ss += w[w.size()-1];
            vvs++;
        }
        for (int j=d[i]; j<w.size(); j++)
            w[j] = w[j] / ss;
#endif

#ifdef NEWWEIGHT
		for (int ii = 0; ii < vec.size(); ii++)
		{
			vec[ii] = vec[ii]*sqrt(w[ii+d[i]]);
		}

#endif
        // add normal
        vec.push_back( OtoE(mesh->normal(vi) * (lens / vec.size() * RefMesh::normalScale)) );
        align[i] = new AffineAlign(vec);
    }
    c.resize(w.size(), 0);
    std::vector<bool> visit(mesh->n_vertices(), false);
    int qBeg = 0;
    for (int i=root; i<visit.size(); i++) if (!visit[i]) {
        visit[i] = true;
        bfsq.push_back(make_pair(i,make_pair(-1,-1) ));
        while (qBeg < bfsq.size()) {
            int i = bfsq[qBeg++].first;

            DTriMesh::VertexHandle vi(i);
            int vvs = d[i];

            for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++, vvs++) {
                int j = vvi.handle().idx();
                if (!visit[j]) {
                    visit[j] = true;
                    bfsq.push_back(make_pair( j, make_pair(i,vvs) ));
                }
            }
        }
    }
}



void RefMesh::getFeature(DTriMesh &ms, FeatureVector &fv) {

    static int callTime = -1;
    callTime ++;
    ms.update_face_normals();
    ms.update_vertex_normals();

    fv.s.resize(align.size());
    fv.r.resize(align.size());
	fv.logr.resize(align.size());
    fv.isConst.resize(align.size(), 0);
    fv.constPoint.resize(align.size());
	fv.isNr.resize(align.size(),0);
	fv.nrPoint.resize(align.size());
	fv.isNrPlane.resize(align.size(),0);
	fv.planePoint.resize(align.size());
/*
    ofstream fout;
    if (callTime == 1) {
        fout.open("color.txt");
        fout << ms.n_vertices() << endl;
    }
    */

    //ofstream cout((string("fvv") + (char)(callTime+'0')).c_str());

	double allres = 0;

    for (int i=0; i<ms.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        DTriMesh::Point p = ms.point(vi);
        std::vector<Eigen::Vector3d> vec;
        double lens = 0;
        for (DTriMesh::VertexVertexIter vvi = ms.vv_iter(vi); vvi; vvi++) {
            Eigen::Vector3d q = OtoE(ms.point(vvi) - p);
            vec.push_back(q);
            lens += q.norm();
        }
#ifdef NEWWEIGHT
		for (int ii = 0; ii < vec.size(); ii++)
		{
			vec[ii] = vec[ii]*sqrt(w[ii+d[i]]);
		}

#endif
        // add normal
        vec.push_back( OtoE(ms.normal(vi) * (lens / vec.size() * RefMesh::normalScale)) );
        if (vec.size() <= 1) {
            fv.s[i] = fv.r[i] = Eigen::Matrix3d::Identity();
            continue;
        }
        Eigen::Matrix3d mat = align[i]->calc(vec);
        polarDec(mat, fv.r[i], fv.s[i]);
		allres+=(align[i]->residualwithoutnormal(mat,vec));
		//yangjie add º∆À„logr
		//for (int i = 0; i < ms.n_vertices(); i++) {
			fv.logr[i] = log(fv.r[i]);
		//}
        /*
        if (i == 12049 && callTime < 2) {
            stringstream ss;
            ss << "tt" << callTime << ".txt";
            ofstream fout(ss.str().c_str());
            fout << "lens : " << lens << endl;
            for (int i=0; i<vec.size(); i++)
                fout << vec[i].transpose() << endl;
        }
        if (fv.s[i].determinant() < 0) {
            fout << 1 << endl;
            cout << "~Warning fv.s " << i << " < 0 id : " << callTime << "\n";
        } else fout << 0 << endl;
        */
    }
	//assert(allres == allres);
	cout<<callTime<<": "<<allres<<std::endl;
//#ifdef _DEBUG
//	fclose(stdout);
//#endif
		
    fv.dr.resize(vvs);
    fv.logdr.resize(vvs);

    if (usec)
        this->fvNum++;

    for (int i=0; i<ms.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int vs = d[i];
        for (DTriMesh::VertexVertexIter vvi = ms.vv_iter(vi); vvi; vvi++) {
            fv.dr[vs] = fv.r[i].transpose() * fv.r[vvi.handle().idx()];
            fv.logdr[vs] = log(fv.dr[vs]);
            if (usec)
                c[vs] = (c[vs] * (this->fvNum-1) + fv.logdr[vs].norm() + (fv.s[i] - Eigen::Matrix3d::Identity()).norm() ) / this->fvNum;
            vs++;
        }
    }
	// get new feature


//	this->checksymm();

}

void RefMesh::getEidmap(DTriMesh &ms, std::vector<std::pair<int,int>>& _eidmap)
{
    for (int i=0; i<ms.n_vertices(); i++)
	{
        DTriMesh::VertexHandle vi(i);
        int vs = d[i];
        for (DTriMesh::VertexVertexIter vvi = ms.vv_iter(vi); vvi; vvi++) 
		{
            std::pair<int,int> tmp = std::pair<int,int>(i,vvi.handle().idx());
            _eidmap.push_back(tmp);            
        }
    }
}

void RefMesh::getEdgeLaplace(DTriMesh &ms, EdgeLaplace& el)
{
	map<int,int> edgemap;
    OpenMesh::EPropHandleT<int> eid;
    ms.add_property(eid);
	int vnum = ms.n_vertices();
	std::vector<std::pair<int,int>> eidmap;
	getEidmap(ms,eidmap);
	el.edgeflipid.resize(eidmap.size(),-1);

	int varnum = 0;
	for(int i = 0; i < eidmap.size(); i++)
	{
		int firid = eidmap[i].first;
		int secid = eidmap[i].second;
		int edgeid = std::min(firid,secid)*vnum+std::max(firid,secid);
        std::map<int,int>::iterator it = edgemap.find(edgeid);
		if(it!=edgemap.end())
		{
			int flipedgeid = it->second;
			el.edgeflipid[i] = flipedgeid;
			el.edgeflipid[flipedgeid] = i;
		}
		else
		{
			edgemap[edgeid] = i;
			varnum++;
		}
	}
    


}




void RefMesh::getFeatureS(DTriMesh &ms, FeatureVector &fv) {

    static int callTime = -1;
    callTime ++;


    ms.update_face_normals();
    ms.update_vertex_normals();

    fv.s.resize(align.size());
    fv.r.resize(align.size());
    fv.isConst.resize(align.size(), 0);
    fv.constPoint.resize(align.size());
/*
    ofstream fout;
    if (callTime == 1) {
        fout.open("color.txt");
        fout << ms.n_vertices() << endl;
    }
    */

    //ofstream cout((string("fvv") + (char)(callTime+'0')).c_str());

    for (int i=0; i<ms.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        DTriMesh::Point p = ms.point(vi);
        std::vector<Eigen::Vector3d> vec;
        double lens = 0;
        for (DTriMesh::VertexVertexIter vvi = ms.vv_iter(vi); vvi; vvi++) {
            Eigen::Vector3d q = OtoE(ms.point(vvi) - p);
            vec.push_back(q);
            lens += q.norm();
        }
        // add normal
        vec.push_back( OtoE(ms.normal(vi) * (lens / vec.size() * RefMesh::normalScale)) );
        if (vec.size() <= 1) {
            fv.s[i] = fv.r[i] = Eigen::Matrix3d::Identity();
            continue;
        }
		// intrinsic (first fundamental form)
		fv.r[i] = Eigen::Matrix3d::Identity();
		// calculate S
		LQSolver ls;
		ls.saveA = false;
		ls.needrs = false;
		int tick = 0;

		vector< pair< pair<int,int>, double > > data;
		vector<double> vb;

		assert( align[i]->p.size()== vec.size());
        for (int _col = 0; _col < vec.size(); _col++)
        {
			// 1 2 3
		    for (int k = 0; k < 3; k++)
		    {
				data.push_back(make_pair(make_pair(tick,k),align[i]->p[_col](k)));
		    }
		    vb.push_back(vec[_col](0));
			tick++;
			// 2 4 5
			data.push_back(make_pair(make_pair(tick,1),align[i]->p[_col](0)));
			data.push_back(make_pair(make_pair(tick,3),align[i]->p[_col](1)));
			data.push_back(make_pair(make_pair(tick,4),align[i]->p[_col](2)));
			vb.push_back(vec[_col](1));
			tick++;
			// 3 5 6
			data.push_back(make_pair(make_pair(tick,2),align[i]->p[_col](0)));
			data.push_back(make_pair(make_pair(tick,4),align[i]->p[_col](1)));
			data.push_back(make_pair(make_pair(tick,5),align[i]->p[_col](2)));
			vb.push_back(vec[_col](2));
			tick++;
        }
		vector<double> res;
		ls.solve(tick, 6, data,vb,res);
		fv.s[i](0,0) = res[0];
		fv.s[i](0,1) = fv.s[i](1,0) = res[1];
		fv.s[i](0,2) = fv.s[i](2,0) = res[2];
		fv.s[i](1,1) = res[3];
		fv.s[i](1,2) = fv.s[i](2,1) = res[4];
		fv.s[i](2,2) = res[5];
		//fv.s[i]
        //Eigen::Matrix3d mat = align[i]->calc(vec);
        //polarDec(mat, fv.r[i], fv.s[i]);


        /*
        if (i == 12049 && callTime < 2) {
            stringstream ss;
            ss << "tt" << callTime << ".txt";
            ofstream fout(ss.str().c_str());
            fout << "lens : " << lens << endl;
            for (int i=0; i<vec.size(); i++)
                fout << vec[i].transpose() << endl;
        }
        if (fv.s[i].determinant() < 0) {
            fout << 1 << endl;
            cout << "~Warning fv.s " << i << " < 0 id : " << callTime << "\n";
        } else fout << 0 << endl;
        */
    }

    fv.dr.resize(vvs);
    fv.logdr.resize(vvs);

    if (usec)
        this->fvNum++;

    for (int i=0; i<ms.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int vs = d[i];
        for (DTriMesh::VertexVertexIter vvi = ms.vv_iter(vi); vvi; vvi++) {
            fv.dr[vs] = fv.r[i].transpose() * fv.r[vvi.handle().idx()];
            fv.logdr[vs] = log(fv.dr[vs]);
            if (usec)
                c[vs] = (c[vs] * (this->fvNum-1) + fv.logdr[vs].norm() + (fv.s[i] - Eigen::Matrix3d::Identity()).norm() ) / this->fvNum;
            vs++;
        }
    }
}


FeatureVector::FeatureVector(std::vector<double> weight, std::vector<FeatureVector> &fvs) {
    s.resize(fvs[0].s.size());
    r.resize(fvs[0].s.size());
    dr.resize(fvs[0].dr.size());
    logdr.resize(fvs[0].dr.size());
    isConst.resize(s.size(), 0);
// 	int snum = s.size();
// 	isConst.resize(snum,false);
// 	isConst.resize(snum);
// 	for (int i = 0; i < snum; i++)
// 	{
// 		isConst[i] = false;
// 	}
    constPoint.resize(s.size());

	isNr.resize(s.size(),false);
    nrPoint.resize(s.size());

	isNrPlane.resize(s.size(),false);
	planePoint.resize(s.size());

    this->blendFrom(weight, fvs);
}


void FeatureVector::IdentityRotation()
{
    for (int i = 0; i < r.size(); i++)
    {
		r[i] = Eigen::Matrix3d::Identity();
    }
}

void FeatureVector::blendFrom(std::vector<double> weight, std::vector<FeatureVector> &fvs) {


#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<s.size(); i++) {
        r[i] = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
        for (int j=0; j<weight.size(); j++) {
            mat += weight[j] * fvs[j].s[i];
            r[i] += weight[j] * log(fvs[j].r[i]);
        }
        r[i] = exp(r[i]);
        s[i] = mat;
    }
DOMP_END;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<dr.size(); i++) {
        Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
        for (int j=0; j<weight.size(); j++)
            mat += weight[j] * fvs[j].logdr[i];
        logdr[i] = mat;
        dr[i] = exp(mat);
    }
DOMP_END;

}


std::ostream& operator<<(std::ostream& cout, FeatureVector &fv) {
    cout << "S : " << std::endl;
    for (int i=0; i<fv.s.size(); i++) {
        cout << i << std::endl << fv.s[i] << std::endl;
        cout << "det : " << fv.s[i].determinant() << endl;
    }
    cout << "dR : " << std::endl;
    for (int i=0; i<fv.dr.size(); i++) {
        cout << i << std::endl << fv.dr[i] << std::endl;
        cout << "det : " << fv.dr[i].determinant() << endl;
    }
    cout << "log(dR) : " << std::endl;
    for (int i=0; i<fv.logdr.size(); i++)
        cout << i << std::endl << fv.logdr[i] << std::endl;
    return cout;
}

void FeatureVector::setConstPoint(int i, Eigen::Vector3d v) {
    isConst[i] = true;
    constPoint[i] = v;
}

void FeatureVector::setNrPoint(int i, Eigen::Vector3d v)
{
	isNr[i] = true;
	nrPoint[i] = v;
}

void FeatureVector::setPlanePoint(int i, Eigen::Vector3d v)
{
	//isNr[i] = true;
	isNrPlane[i] = true;
    planePoint[i] = v;
}


void FeatureVector::loadConstPoint(std::istream& cin) {
    int n;
    cin >> n;
    for (int i=0; i<n; i++) {
        int m;
        cin >> m;
        std::vector<int> ids(m);
        for (int i=0; i<m; i++)
            cin >> ids[i];
        for (int i=0; i<m; i++) {
            double x,y,z;
            cin >> x >> y >> z;
            this->setConstPoint(ids[i], Eigen::Vector3d(x,y,z));
        }
    }
}

void FeatureVector::loadConstPointFixed(std::istream& cin) {
	int n;
	cin >> n;
	for (int i=0; i<n; i++) {
		int m;
		cin >> m;
		std::vector<int> ids(m);
		for (int i=0; i<m; i++)
			cin >> ids[i];
		for (int i=0; i<m; i++) {
			double x,y,z;
			cin >> x >> y >> z;
			this->setConstPoint(ids[i], Eigen::Vector3d(x,y,z));
		}
		if (i==0)
		{
            iddeformed = ids;
		}
		if (i==1)
		{
			idfixed = ids;
		}
	}
}


void FeatureVector::loadHandleRT(std::istream& cin)
{
	double mat[4][4] ={0.0f};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
             cin>>mat[i][j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rhandle(i,j) = mat[i][j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
        thandle(i) = mat[i][3];		
	}

}



void FeatureVector::loadConstPoint(std::istream& cin, DTriMesh& mesh) {
    int n;
    cin >> n;
    for (int i=0; i<n; i++) {
        int m;
        cin >> m;
        std::vector<int> ids(m);
        for (int i=0; i<m; i++)
            cin >> ids[i];
        for (int i=0; i<m; i++) {
            double x,y,z;
            cin >> x >> y >> z;
            this->setConstPoint(ids[i], OtoE(mesh.point(DTriMesh::VertexHandle(ids[i]))));
        }
    }
}

void FeatureVector::loadsoftConstPoint(std::istream& cin) {
	int n;
	cin >> n;
	for (int i = 0; i < n; i++) {
		int m;
		cin >> m;
		std::vector<int> ids(m);
		for (int i = 0; i < m; i++)
			cin >> ids[i];
		for (int i = 0; i < m; i++) {
			double x, y, z;
			cin >> x >> y >> z;
			this->setNrPoint(ids[i], Eigen::Vector3d(x, y, z));
		}
	}
}


void FeatureVector::loadConstPoint(std::vector<int>& ids, std::vector<OpenMesh::Vec3d>& v3ds)
{
	for (int i = 0; i < ids.size(); i++)
	{
        this->setConstPoint(ids[i], OtoE(v3ds[i]));
	}
}

void FeatureVector::loadNrPoint(std::vector<int>& ids, std::vector<OpenMesh::Vec3d>& v3ds)
{
	for (int i = 0; i < ids.size(); i++)
	{
		this->setNrPoint(ids[i], OtoE(v3ds[i]));
	}
}

void FeatureVector::loadNrPlanePoint(std::vector<int>& ids, std::vector<OpenMesh::Vec3d>& v3ds)
{
	for (int i = 0; i < ids.size(); i++)
	{
		this->setPlanePoint(ids[i], OtoE(v3ds[i]));
	}
}

double RefMesh::getWij(int e) {
    //if (fvNum) return exp(-c[e]);// * w[e];
    return w[e];
}

void RefMesh::printCij(std::ofstream &fout) {
    std::vector<double> color;
    for (int i=0; i<mesh->n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int vs = d[i];
        int s=0;
        double w=0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            w += this->getWij(vs);
            s++;
            vs++;
        }
        color.push_back(w/s);
    }
    double mx = color[0], mi = mx;
    for (int i=1; i<color.size(); i++) {
        mx = std::max(mx, color[i]);
        mi = std::min(mi, color[i]);
    }
    fout << color.size() << endl;
    for (int i=0; i<color.size(); i++) {
        fout << (color[i]-mi) / (mx-mi) << '\n';
    }
    std::cout << "Color Max : " << mx << endl;
    std::cout << "Color Min : " << mi << endl;
    std::cout << "Color Max-Min : " << mx - mi << endl;
    fout << endl;
}

void FeatureVector::calcLogRFromR() {
    logr.resize(r.size());
    for (int i=0; i<logr.size(); i++)
        logr[i] = log(r[i]);
}

void FeatureVector::calcRFromLogR() {
    r.resize(logr.size());
    dr.resize(logdr.size());
    for (int i=0; i<logr.size(); i++)
        r[i] = exp(logr[i]);
    for (int i=0; i<logdr.size(); i++)
        dr[i] = exp(logdr[i]);
    isConst.resize(s.size(), false);
	isNr.resize(s.size(),false);
}

void FeatureVector::resize(const FeatureVector &other) {
    s.resize(other.s.size());
    logdr.resize(other.logdr.size());
    logr.resize(other.logr.size());
}

FeatureVector operator +(const FeatureVector &a, const FeatureVector &b) {
    FeatureVector ans;
    ans.resize(a);
    for (int i=0; i<ans.s.size(); i++) {
        ans.s[i] = a.s[i] + b.s[i];
        ans.logr[i] = a.logr[i] + b.logr[i];
    }
    for (int i=0; i<ans.logdr.size(); i++) {
        ans.logdr[i] = a.logdr[i] + b.logdr[i];
    }
    return ans;
}

FeatureVector operator -(const FeatureVector &a, const FeatureVector &b) {
    FeatureVector ans;
    ans.resize(a);
    for (int i=0; i<ans.s.size(); i++) {
        ans.s[i] = a.s[i] - b.s[i];
        ans.logr[i] = a.logr[i] - b.logr[i];
    }
    for (int i=0; i<ans.logdr.size(); i++) {
        ans.logdr[i] = a.logdr[i] - b.logdr[i];
    }
    return ans;
}


double dot(const Eigen::Matrix3d &a, const Eigen::Matrix3d &b) {
    return
            a(0,0) * b(0,0) + a(0,1) * b(0,1) + a(0,2) * b(0,2) +
            a(1,0) * b(1,0) + a(1,1) * b(1,1) + a(1,2) * b(1,2) +
            a(2,0) * b(2,0) + a(2,1) * b(2,1) + a(2,2) * b(2,2);
}

double operator *(const FeatureVector &a, const FeatureVector &b) {
    double ans = 0;
    for (int i=0; i<a.s.size(); i++) {
        ans += dot(a.s[i], b.s[i]);
    }
    for (int i=0; i<a.logdr.size(); i++) {
        ans += dot(a.logdr[i], b.logdr[i]);
    }
    return ans;
}

FeatureVector operator *(const FeatureVector &a, const double &b) {
    FeatureVector ans;
    ans.resize(a);
    for (int i=0; i<ans.s.size(); i++) {
        ans.s[i] = a.s[i] * b;
        ans.logr[i] = a.logr[i] * b;
    }
    for (int i=0; i<ans.logdr.size(); i++) {
        ans.logdr[i] = a.logdr[i] * b;
    }
    return ans;
}

bool RefMesh::checksymm()
{
	map<int,vector<int>> edgesets;
	int nv = mesh->n_vertices();
	for (int i = 0; i < mesh->n_vertices(); i++)
	{
		int ei = d[i];
		DTriMesh::VertexHandle vi(i);
		for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++,ei++) {
			int j = vvi.handle().idx();
			int smalli,smallj;
			if (i<j)
			{
				smalli = i;
				smallj = j;
			}
			else
			{
		        smalli = j;
				smallj = i;
			}
			int id = smalli* nv+smallj;
			std::map<int,vector<int>>::iterator it = edgesets.find(id);

			if (it!= edgesets.end())
			{
				it->second.push_back(ei);
			}
			else
			{
				vector<int> tmp;
				tmp.push_back(ei);
				edgesets.insert(pair<int,vector<int>>(id,tmp));
			}
		}
	}
	vector<double> vd;
	for (map<int,vector<int>>::iterator iter = edgesets.begin(); iter!=edgesets.end(); iter++)
	{
		assert(iter->second.size()==2);
		int edgeid1 = iter->second[0];
		int edgeid2 = iter->second[1];
		double val = abs(getWij(edgeid1)-getWij(edgeid2));
		vd.push_back(val);
		std::cout<<val<<std::endl;
		//assert(abs(getWij(edgeid1)-getWij(edgeid2))<=0.001);
	}
	return true;
}

		//yangjie add

		FeatureVector::FeatureVector(DTriMesh & mesh)
		{
			this->r.resize(mesh.n_vertices(), Eigen::Matrix3d::Identity());
			this->s.resize(mesh.n_vertices(), Eigen::Matrix3d::Identity());
			this->logr.resize(mesh.n_vertices(), Eigen::Matrix3d::Zero());
			this->dr.resize(mesh.n_edges() * 2, Eigen::Matrix3d::Identity());
			this->logdr.resize(mesh.n_edges() * 2, Eigen::Matrix3d::Zero());
			this->isConst.resize(mesh.n_vertices(), false);
			this->constPoint.resize(mesh.n_vertices(), Eigen::Vector3d::Zero());
			this->rots.resize(mesh.n_vertices(),Rot::Rot());
		}

		void FeatureVector::yj_blendFrom(std::vector<double> weight, std::vector<FeatureVector> &fvs, std::string string) {
#ifdef DUSE_OPENMP
#pragma omp parallel
			{
#pragma omp for
#endif
				for (int i = 0; i<s.size(); i++) {
					r[i] = Eigen::Matrix3d::Zero();
					rots[i].r = Eigen::Matrix3d::Zero();
					Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
					for (int j = 0; j<weight.size(); j++) {
						mat += weight[j] * fvs[j].s[i];
						r[i] += weight[j] * fvs[j].logr[i];
						rots[i].r += weight[j] * fvs[j].rots[i].logr;
					}
					//std::cout<<r[i]<<std::endl;
					//std::cout<<rots[i].r;
					r[i] = exp(r[i]);
					rots[i].r = exp(rots[i].r);
					//std::cout<<rots[i].r<<"\n"<<r[i];
					s[i] = mat;
				}
				DOMP_END;

			}

			FeatureVector::FeatureVector(std::vector<double> weight, std::vector<FeatureVector> &fvs, std::string string) {
				s.resize(fvs[0].s.size());
				r.resize(fvs[0].s.size());
				dr.resize(fvs[0].dr.size());
				logdr.resize(fvs[0].dr.size());
				logr.resize(fvs[0].r.size());
				isConst.resize(s.size(), 0);
				rots.resize(fvs[0].rots.size());
				// 	int snum = s.size();
				// 	isConst.resize(snum,false);
				// 	isConst.resize(snum);
				// 	for (int i = 0; i < snum; i++)
				// 	{
				// 		isConst[i] = false;
				// 	}
				constPoint.resize(s.size());

				isNr.resize(s.size(), false);
				nrPoint.resize(s.size());

				isNrPlane.resize(s.size(), false);
				planePoint.resize(s.size());

				this->yj_blendFrom(weight, fvs, string);
			}