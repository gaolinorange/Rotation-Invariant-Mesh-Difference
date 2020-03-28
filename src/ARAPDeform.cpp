#include "ARAPDeform.h"
#include <fstream>
#include <queue>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <sstream>
#include <limits>
#include <algorithm>
#include "coutcolorarap.h"
#include "graph.h"


using namespace std;
using namespace Eigen;

#include <omp.h>

#ifdef DUSE_OPENMP
#define DOMP_END \
}
#else
#define DOMP_END ;
#endif

#define floatTime ((clock()-tt)*1.0 / CLOCKS_PER_SEC)

bool isNan(double fN) {
    return !(fN==fN);
}

bool isNan(Eigen::Matrix3d m) {
    for (int i=0; i<3; i++) for (int j=0; j<3; j++)
        if (isNan(m(i,j))) return true;
    return false;
}


bool isNan(Eigen::Vector3d v) {
    for (int i=0; i<3; i++)
        if (isNan(v(i))) return true;
    return false;
}

ARAPDeform::ARAPDeform(DTriMesh &refMesh, std::vector<DTriMesh*> ms) : ref(refMesh), meshs(ms), lqsolver() {
//ARAPDeform::ARAPDeform(DMEngine &eng,  &refMesh, std::vector<DTriMesh*> ms) : ref(refMesh), meshs(ms), lqsolver(&eng) {
    int cpunum = omp_get_num_procs();
	//cpunum=1;
    omp_set_num_threads(cpunum);
    SRp = 1;

    iterEps = 0.01;
    newtonIterEps = 1e-3;
    needAlign = false;
    maxIterTime = 100;
    newtonIterTime = 5;
    AtAChanged = true;
    needCalcRs = true;
    CalcLocalFirst = false;
    //this->eng = &eng;
    rlimit = make_pair(0,1);
    fvs.resize(ms.size());
	this->getsimfv = false;
	initwithhandle = false;

	//yangjie add 初始化变量
	//this->csrRowIndGTPtr=NULL;
	//this->csrColGTPtr = NULL;
	//this->csrValGTPtr = NULL;
	//this->nnzGT = 0;
	//this->csrRowNumbGT = 0;
	//this->csrColNumbGT = 0;

	//this->csrRowIndGPtr_chol = NULL;
	//this->csrColGPtr_chol = NULL;
	//this->csrValGPtr_chol = NULL;
	//this->nnzG_chol=0;
	//this->csrRowNumbG_chol=0;
	//this->csrColNumbG_chol=0;
	//this->resultXarap = NULL;

#ifndef PARALLEL

//	#ifdef _DEBUG
    //freopen("errorlog.txt","w",stdout);      
//    #endif

    ofstream fout("fv.txt");
    for (int i=0; i<fvs.size(); i++) {
        std::cout << "Calc " << i << " feature vector " << endl;
        ref.getFeature(*meshs[i], fvs[i]);
        fout << "F"<<i<<endl;
        fout << fvs[i] << endl;
		std::cout << getRS(fvs[i], *meshs[i])<<endl;
    }
//#ifdef _DEBUG
//	fclose(stdout);
//#endif
#else
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
	for (int i=0; i<fvs.size(); i++) {
        std::cout << "Calc " << i << " feature vector " << endl;
		if (i==1)
		{
			std::cout<<"Debug"<<endl;
		}
		ref.getFeature(*meshs[i], fvs[i]);
		this->bfscorrot(fvs[i]);
        //fout << "F"<<i<<endl;
        //fout << fvs[i] << endl;
    }
DOMP_END
#endif
	 ref.getEidmap(*meshs[0],_eidmap);
}

void ARAPDeform::GetFeature(DTriMesh &refMesh, std::vector<DTriMesh*> ms)
{
	fvs.resize(ms.size());
	//ofstream fout("fv.txt");
#ifndef PARALLEL
	for (int i=0; i<fvs.size(); i++) {
		std::cout << "Calc " << i << " feature vector " << endl;
		ref.getFeature(*meshs[i], fvs[i]);
// 		fout << "F"<<i<<endl;
// 		fout << fvs[i] << endl;
	}
#else
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
	for (int i=0; i<fvs.size(); i++) {
        std::cout << "Calc " << i << " feature vector " << endl;
        ref.getFeature(*meshs[i], fvs[i]);
        //fout << "F"<<i<<endl;
        //fout << fvs[i] << endl;
    }

DOMP_END
#endif
}

void ARAPDeform::GetFeatureS(DTriMesh &refMesh, std::vector<DTriMesh*> ms)
{
	this->meshs = ms;
	fvs.resize(ms.size());
	//ofstream fout("fv.txt");
	for (int i=0; i<fvs.size(); i++) {
		std::cout << "Calc " << i << " feature vector " << endl;
		ref.getFeatureS(*meshs[i], fvs[i]);
		// 		fout << "F"<<i<<endl;
		// 		fout << fvs[i] << endl;
	}
}

void ARAPDeform::init()
{
	int cpunum = omp_get_num_procs();//cpunum=1;
	omp_set_num_threads(cpunum);
	SRp = 1;
	iterEps = 0.01;
	newtonIterEps = 1e-3;
	needAlign = false;
	maxIterTime = 100;
	newtonIterTime = 5;
	AtAChanged = true;
	needCalcRs = true;
	CalcLocalFirst = false;
	lambda = 1e6;
	lambdaplane = 1;
	//this->eng = &eng;
	rlimit = make_pair(0,1);
	initwithhandle = false;
// 	fvs.resize(ms.size());
// 	ofstream fout("fv.txt");
// 	for (int i=0; i<fvs.size(); i++) {
// 		std::cout << "Calc " << i << " feature vector " << endl;
// 		ref.getFeature(*meshs[i], fvs[i]);
// 		fout << "F"<<i<<endl;
// 		fout << fvs[i] << endl;
// 	}
}

//void ARAPDeform::SetEngine(DMEngine* _eng)
//{
//	this->eng = _eng;
//	lqsolver.SetEngine(_eng);
//}


void ARAPDeform::initMesh(FeatureVector &fv, DTriMesh &mesh) {
    long long tt = clock();
    std::cout << "init feature vector ... ";
    bfsInit2(fv, mesh);
    std::cout << " time : " << floatTime << endl;
    return;
    vector<int> visit(fv.s.size(), 0);
    for (int i=0; i<visit.size(); i++)
        if (!visit[i]) {
            //fv.r[i] = fvs[0].r[i];
            bfsInit(i, visit, fv, mesh);
        }
}

void ARAPDeform::initMeshRotationHandle(FeatureVector &fv, DTriMesh &mesh)
{
    //convert mesh to adjalist
	adjalist adlist;
    ConvertDTriMeshToAdjalist(mesh, adlist);
	assert(mesh.n_vertices()==fv.r.size());
	Eigen::Matrix3d fixed = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d logfixed = log(fixed);
	Eigen::Matrix3d deformed = fv.rhandle;
	Eigen::Matrix3d logdeformed = log(deformed);
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
	for(int i = 0; i < mesh.n_vertices();i++)
	{
		GraPath meshpath;
		Dijkstra(adlist,i, meshpath);
		vector<double> fixeddist;
		vector<double> handledist;
		for (int j = 0; j < fv.idfixed.size(); j++)
		{
			fixeddist.push_back(meshpath.distlist[fv.idfixed[j]]);
		}
		for(int j = 0; j < fv.iddeformed.size(); j++)
		{
			handledist.push_back(meshpath.distlist[fv.iddeformed[j]]);
		}
		double dis2fix = *std::min_element(fixeddist.begin(),fixeddist.end());
		double dis2handle = *std::min_element(handledist.begin(),handledist.end());
		double _dis = dis2fix+dis2handle;
        double wfix,whandle;
		if (abs(_dis)<0.000000001)
		{
			wfix=whandle=0.5;
		}
		wfix = dis2fix/_dis;
		whandle = dis2handle/_dis;
		//Eigen::Matrix3d logr = wfix*logfixed+whandle*logdeformed;
		Eigen::Matrix3d logr = (1-wfix)*logfixed+wfix*logdeformed;
		fv.r[i] = exp(logr);
	}
DOMP_END
}



void ARAPDeform::writeIterMesh(DTriMesh &mesh, string name, int id) {
#ifdef WRITE_ITER_MESH
    DTriMesh ms = mesh;
    RotateAlign::AlignAtoB(ms, *(ref.mesh));
    stringstream ss;
    ss << name << id+1000 << ".obj";
    std::cout << "Write mesh " + ss.str() << endl;
    OpenMesh::IO::write_mesh(ms, ss.str().c_str());
#endif
}

void ARAPDeform::solve(FeatureVector fv, DTriMesh &mesh, bool predec) {
    if (CalcLocalFirst) {
        double presr = this->SRp;
        this->SRp = 0;
        localOptimize(fv, mesh);
        this->SRp = presr;
    } else
	{
		if(!initwithhandle)
		{           		
			initMesh(fv, mesh);
		}
		else
		{
			initMeshRotationHandle(fv,mesh);
		}
	}
        

//	fv.IdentityRotation();


    DTriMesh meshtst;
    //OpenMesh::IO::read_mesh(meshtst,"E:/SIGA2014/recons/barsim2/2.obj");
	//std::cout<<getRS(fv,meshtst)<<endl;

 //   OpenMesh::IO::read_mesh(meshtst,"E:/SIGA2014/recons/barsim2/1.obj");
//	  std::cout<<getRS(fvs[0],meshtst)<<endl;
    if (predec)
        preDecMatrix(fv, mesh);

    double preRes = 1e10;
    double res = globalOptimize(fv, mesh);
	//es.push_back(res);
	std::cout<<"Berfor Iter Res: "<<res<<endl;
    int iter = 0;
    writeIterMesh(mesh, "iter", iter);
//	OpenMesh::IO::write_mesh(mesh,"test1.obj");
    while (abs(preRes - res) > iterEps && iter < maxIterTime) {
        iter ++;
        std::cout << "Iter time : " << iter << endl;
        preRes = res;
        localOptimize(fv, mesh);
		std::cout<<"After Local Optimize: "<<getRS(fv,mesh)<<endl;
//		OpenMesh::IO::write_mesh(mesh,"test2.obj");
        res = globalOptimize(fv, mesh);
		std::cout<<"After Global Optimize: "<<res<<endl;
//		OpenMesh::IO::write_mesh(mesh,"test3.obj");
        writeIterMesh(mesh, "iter", iter);
		es.push_back(res);
        std::cout << "DResidual : " << abs(preRes - res) << endl;
    }

	if (getsimfv)
	{
		this->simfv = fv;
	}

    ckAlign(mesh);
}

void ARAPDeform::gl_solve(FeatureVector fv, DTriMesh &mesh, bool predec) {
	if (CalcLocalFirst) {
		double presr = this->SRp;
		this->SRp = 0;
		localOptimize(fv, mesh);
		this->SRp = presr;
	} else
		initMesh(fv, mesh);

	DTriMesh meshtst;
	//OpenMesh::IO::read_mesh(meshtst,"E:/SIGA2014/recons/barsim2/2.obj");
	//std::cout<<getRS(fv,meshtst)<<endl;

	//   OpenMesh::IO::read_mesh(meshtst,"E:/SIGA2014/recons/barsim2/1.obj");
	//	  std::cout<<getRS(fvs[0],meshtst)<<endl;
	if (predec)
		preDecMatrix(fv, mesh);

	double preRes = 1e10;
	//double res = globalOptimize(fv, mesh);
	double res = gl_globalOptimize_nr(fv,mesh);
	std::cout<<"Berfor Iter Res: "<<res<<endl;
	int iter = 0;
	writeIterMesh(mesh, "iter", iter);
	//	OpenMesh::IO::write_mesh(mesh,"test1.obj");
	while (abs(preRes - res) > iterEps && iter < maxIterTime) {
		iter ++;
		std::cout << "Iter time : " << iter << endl;
		preRes = res;
		localOptimize(fv, mesh);
		std::cout<<coutcmd::yellow<<"After Local Optimize: "<<getRS_nr(fv,mesh)<<coutcmd::white<<endl;
		//		OpenMesh::IO::write_mesh(mesh,"test2.obj");
		res = gl_globalOptimize_nr(fv, mesh);
		std::cout<<"After Global Optimize: "<<coutcmd::yellow<<res<<coutcmd::white<<endl;
		//		OpenMesh::IO::write_mesh(mesh,"test3.obj");
		writeIterMesh(mesh, "iter", iter);
		std::cout << "DResidual : " << abs(preRes - res) << endl;
	}
	//ckAlign(mesh);
}




void ARAPDeform::bfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh) {
    queue<int> q;
    q.push(i);
    visit[i] = true;
    while (!q.empty()) {
        i = q.front(); q.pop();
        Matrix3d &r = fv.r[i];
        DTriMesh::VertexHandle vi(i);
        int vs = ref.d[i];
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++) {
            int j = vvi.handle().idx();
            if (!visit[j]) {
                fv.r[j] = r * fv.dr[vs];
                q.push(j);
                visit[j] = true;
            }
            vs++;
        }
    }
}

void ARAPDeform::bfsInit2(FeatureVector &fv, DTriMesh &mesh) {
    for (int i=0; i<ref.bfsq.size(); i++) {
        int j = ref.bfsq[i].first;
        int fa = ref.bfsq[i].second.first;
        int ei = ref.bfsq[i].second.second;
        if (fa != -1)
            fv.r[j] = fv.r[fa] * fv.dr[ei];
    }
}

void ARAPDeform::bfscorrot(FeatureVector& fv)
{
	fv.rots.resize(fv.r.size());
	vector<bool> tovist(fv.r.size(),false);
    for (int i=0; i<ref.bfsq.size(); i++) {
        int j = ref.bfsq[i].first;
		//if (j==600)
		//{
		//	int a =1;
		//}
        int fa = ref.bfsq[i].second.first;
        int ei = ref.bfsq[i].second.second;
		//std::cout<<j<<" "<<fa<<endl;
		if (j==180)
		{
			//std::cout<<j<<endl;
		}
        if (fa != -1)
		{
			assert(tovist[fa]);
			fv.rots[j] = logrot(fv.r[j],fv.rots[fa]);
			tovist[j] = true;
		}
		else
		{
			fv.rots[j] = logrot(fv.r[j]);
			tovist[j] = true;
		}
//		std::cout<<j<<": "<<fv.rots[j].ToAngle()<<endl;
		Eigen::Matrix3d logr = log(fv.r[j]);
		Eigen::Matrix3d r = exp(fv.rots[j].logr);
		//double _normal = (fv.rots[j].logr-logr).squaredNorm();
		double _normal = (exp(fv.rots[j].logr)-fv.r[j]).squaredNorm();
		if (_normal>0.0001)
		{
			std::cout<<"error"<<endl;
		}
    }
}

void ARAPDeform::dfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh) {
    DTriMesh::VertexHandle vi(i);
    visit[i] = 1;
    Matrix3d &r = fv.r[i];
    int vs = ref.d[i];
    for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++) {
        int j = vvi.handle().idx();
        if (!visit[j]) {
            fv.r[j] = r * fv.dr[vs];
            dfsInit(j, visit, fv, mesh);
        }
        vs++;
    }
}


void ARAPDeform::preDecMatrix(FeatureVector &fv, DTriMesh &mesh) {
    if (!AtAChanged) return;
    AtAChanged = false;
    lqsolver.saveA = true;
    globalOptimize(fv, mesh);
    lqsolver.saveA = false;
}

void ARAPDeform::preDecMatrix_nr(FeatureVector &fv, DTriMesh &mesh) {
    if (!AtAChanged) return;
    AtAChanged = false;
    lqsolver.saveA = true;
    globalOptimize_nr(fv, mesh);
    lqsolver.saveA = false;
}

void ARAPDeform::gl_preDecMatrix_nr(FeatureVector &fv, DTriMesh &mesh) {
    if (!AtAChanged) return;
    AtAChanged = false;
    lqsolver.saveA = true;
    gl_globalOptimize_nr(fv, mesh);
    lqsolver.saveA = false;
}




double ARAPDeform::globalOptimize(FeatureVector &fv, DTriMesh &mesh) {
#ifdef GLNEW
	  return this->gl_globalOptimize(fv,mesh);
#endif
    if (!lqsolver.saveA) return globalOptimizeFast(fv, mesh);

    std::cout << "Global Optimize ...... " << endl;
    long long tt = clock();
    vector< pair< pair<int,int>, double > > data;
    vector<double> vb;
    int n=0;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        if (fv.isConst[i]) {
            for (int dim=0; dim<3; dim++) {
                data.push_back( make_pair( make_pair(n++, i*3+dim), 1 ));
                vb.push_back(0);
            }
        }
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            //double wij = ref.w[ei];
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            DTriMesh::Point qj = ref.mesh->point(vj);
            int ej = ref.d[j];
			double wj = 1/(double)(ref.degree[j]);
            for (DTriMesh::VertexVertexIter vvj = mesh.vv_iter(vj); vvj; vvj++, ej++) {
                int k = vvj.handle().idx();
				double wk = 1/(double)(ref.degree[k]);
                double wjk = ref.w[ej];
				//
                DTriMesh::Point tqjk = ref.mesh->point(DTriMesh::VertexHandle(k)) - qj;
				//p_k-p_j
                Eigen::Vector3d qjk(tqjk[0],tqjk[1],tqjk[2]);
                Eigen::Vector3d c = ri*(rij*(sj*qjk));
                if (fv.isConst[k]) c -= fv.constPoint[k];
                if (fv.isConst[j]) c += fv.constPoint[j];
#ifndef NEWWEIGHT
                 for (int dim=0; dim<3; dim++) {
                    if (!fv.isConst[k])
                        data.push_back(make_pair( make_pair(n, k*3+dim), wij * wjk ));
                    if (!fv.isConst[j])
                        data.push_back(make_pair( make_pair(n, j*3+dim), - wij * wjk ));
                    vb.push_back( wij * wjk * c(dim) );
                    n++;
                }
#else
				 for (int dim=0; dim<3; dim++) {
                    if (!fv.isConst[k])
                        data.push_back(make_pair( make_pair(n, k*3+dim), wj * wjk ));
                    if (!fv.isConst[j])
                        data.push_back(make_pair( make_pair(n, j*3+dim), - wj * wjk ));
                    vb.push_back( wj * wjk * c(dim) );
                    n++;
                }
#endif
            }
        }
    }
    vector<double> result;
    std::cout << "matlab solve" << endl;
    lqsolver.needrs = true;
    double rs = lqsolver.solve(n, fv.s.size()*3, data, vb, result);
    if (lqsolver.saveA) return 0;
    lqsolver.needrs = false;
    std::cout << "copy ans" << endl;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        for (int dim=0; dim<3; dim++)
            mesh.point(vi)[dim] = result[i*3+dim];
    }

    std::cout << "!!!Global Optimize Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << rs*rs << endl;
    return rs*rs;
}

double ARAPDeform::gl_globalOptimize(FeatureVector &fv, DTriMesh &mesh) {
    if (!lqsolver.saveA) return gl_globalOptimizeFast(fv, mesh);	                            
    std::cout << "Global Optimize ...... " << endl;
    long long tt = clock();
    vector< pair< pair<int,int>, double > > data;
	vector<double> vb(ref.mesh->n_vertices()*3,0);
	this->vbconst.clear();
	this->vbconst.resize(ref.mesh->n_vertices()*3,0);
    int n=0;

	//vector<Eigen::Matrix3d> rdrs(ref.mesh->n_vertices(),Eigen::Matrix3d::Zero());
	this->rdrs.clear();
	this->rdrs.resize(ref.mesh->n_vertices(),Eigen::Matrix3d::Zero());

	double* weights = new double[ref.mesh->n_vertices()];
	int * debugtick = new int[ref.mesh->n_vertices()];
    for (int i=0; i<mesh.n_vertices(); i++) {

    	memset(weights,0,ref.mesh->n_vertices()*sizeof(double));
	    memset(debugtick,0,ref.mesh->n_vertices()*sizeof(int));

        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
		double wi = 1/(double)ref.degree[i];
        int ei = ref.d[i];
        if (fv.isConst[i]) 
		{
            for (int dim=0; dim<3; dim++) 
			{
                //data.push_back( make_pair( make_pair(n++, i*3+dim), 1 ));
				data.push_back( make_pair( make_pair(i*3+dim, i*3+dim), 1 ));
                //vb.push_back(0);
				vb[i*3+dim] = fv.constPoint[i][dim];
            }
			continue;
        }

		double sumwij = 0, sumwji = 0;
		for (int k = 0; k < ref.rd[i].size(); k++)
		{
			int j = ref.rd[i][k].first;
			int ej = ref.rd[i][k].second;
			weights[j] = ref.getWij(ej);
			sumwji+=ref.getWij(ej);
			//
			Eigen::Matrix3d &rji = fv.dr[ej];
			Eigen::Matrix3d &rj = fv.r[j];
			this->rdrs[i]+= (rj*rji);
			debugtick[j]++;
		}
		Eigen::Matrix3d &si = fv.s[i];
        rdrs[i] = rdrs[i]*si;
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            //double wij = ref.w[ei];
			double wj = 1/(double)ref.degree[j];
            double wij = ref.getWij(ei);
			sumwij+=wij;
			weights[j] += ref.getWij(ei);
			debugtick[j]++;
			assert(debugtick[j]==2);
			if (!fv.isConst[j])
			{
       			for (int dim = 0; dim < 3; dim++)
     			{
		    		data.push_back(make_pair( make_pair(i*3+dim, j*3+dim), -weights[j]));
	    		}			
			}
			else
			{
				Eigen::Vector3d vec = weights[j]*fv.constPoint[j];
				for (int dim = 0; dim < 3; dim++)
				{
					//vb[j*3+dim]+=vec[dim];
					vbconst[i*3+dim]+=vec[dim];
				}
				//vb[j]+=(weights[j]*fv.constPoint[j]);
			}
            //Eigen::Matrix3d &rij = fv.dr[ei];
            //Eigen::Matrix3d &sj = fv.s[j];
			//rdrs[i] += (ri*rij*sj);
		}
		if(!fv.isConst[i])
		{
		  if(!fv.isNr[i])
		  {
		    for (int dim = 0; dim < 3; dim++)
		    {
			    data.push_back( make_pair( make_pair(i*3+dim, i*3+dim), sumwij+sumwji ));
		    }
		  }
		  else
		  {
			for(int dim = 0; dim < 3; dim++)
			{
				data.push_back(make_pair(make_pair(i*3+dim,i*3+dim),sumwij+sumwji+lambda));
				vbconst[i*3+dim]+=(lambda*fv.nrPoint[i](dim));
			}
		  }
		}
	 }
	 delete[] weights;
	 weights = NULL;
	 delete[] debugtick;
	 debugtick = NULL;

	//assign b
	//parallel

	 this->assignvb(fv,mesh,vb);

    vector<double> result;
    std::cout << "matlab solve" << endl;
    lqsolver.needrs = true;
    double rs = lqsolver.glsolve(fv.s.size()*3, fv.s.size()*3, data, vb, result);
    if (lqsolver.saveA) return 0;
    lqsolver.needrs = false;


    std::cout << "copy ans" << endl;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        for (int dim=0; dim<3; dim++)
            mesh.point(vi)[dim] = result[i*3+dim];
    }
	rs = this->getRS(fv,mesh);
	std::cout << "!!!Global Optimize Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << rs << endl;
	return rs;
    //std::cout << "!!!Global Optimize Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << rs*rs << endl;
    //return rs*rs;
}

void ARAPDeform::assignvb(FeatureVector &fv, DTriMesh &mesh, std::vector<double>& vb)
{
	#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		if(fv.isConst[i])
		{
			for(int dim = 0; dim < 3; dim++)
			{
				vb[i*3+dim] = fv.constPoint[i][dim];
			}
			continue;
		}

		DTriMesh::VertexHandle vi(i);
		int ei = ref.d[i];
		double wi = 1/(double)ref.degree[i];
		DTriMesh::Point qi = ref.mesh->point(DTriMesh::VertexHandle(i)); 
		Eigen::Vector3d vectmp = Eigen::Vector3d::Zero();
		for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++)
		{
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
			double wij = ref.getWij(ei);
			DTriMesh::Point qj = ref.mesh->point(DTriMesh::VertexHandle(j));
            DTriMesh::Point qij = qi-qj;
			Eigen::Vector3d qij1(qij[0],qij[1],qij[2]);
			vectmp += wi*wij*rdrs[i]*qij1;
		}
        
		for (int k = 0; k < ref.rd[i].size(); k++)
		{
			int j = ref.rd[i][k].first;
			int ej = ref.rd[i][k].second;
			double wj = 1/(double)ref.degree[j];
			double wji = ref.getWij(ej);
			DTriMesh::Point qj = ref.mesh->point(DTriMesh::VertexHandle(j));
            DTriMesh::Point qij = qi-qj;
			Eigen::Vector3d qij1(qij[0],qij[1],qij[2]);
			vectmp += wj*wji*rdrs[j]*qij1;
		}

		for(int dim = 0; dim < 3; dim++)
		{
			vb[i*3+dim]+=vectmp[dim];
			vb[i*3+dim]+=vbconst[i*3+dim];

			//vb[i*3+dim]*=0.5;
		}


	}
DOMP_END;
}

double ARAPDeform::getRS(FeatureVector &fv, DTriMesh &mesh)
{
	assert(fv.s.size()==mesh.n_vertices());
	double residual = 0;
    if (needCalcRs) {
		vector<Eigen::Matrix3d> sRs(mesh.n_vertices(),Eigen::Matrix3d::Zero());
		vector<Eigen::Matrix3d> Rs(mesh.n_vertices(),Eigen::Matrix3d::Zero());
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for 
#endif
	for(int j = 0; j < mesh.n_vertices(); j++)
	{
         for(int k = 0; k < ref.rd[j].size(); k++)
		 {
			 int i = ref.rd[j][k].first;
			 int eij = ref.rd[j][k].second;
			 Eigen::Matrix3d mat=fv.r[i]*fv.dr[eij]*fv.s[j];             
			 sRs[j]+=mat.transpose()*mat;
			 Rs[j]+=mat;
		 }
	}
DOMP_END;
  
        vector<double> rss(omp_get_num_procs(), 0);
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for (int i=0; i<mesh.n_vertices(); i++) {
            double &residual = rss[omp_get_thread_num()];
            DTriMesh::VertexHandle vi(i);
            int ei = ref.d[i];
            DTriMesh::Point qi = ref.mesh->point(vi);
            DTriMesh::Point pi = mesh.point(vi);
			double wi = 1/(double)ref.degree[i];
            for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
                int j = vvi.handle().idx();
				double wj = 1/(double)ref.degree[j];
                DTriMesh::VertexHandle vj(j);
                Vector3d qij = OtoE(qi - ref.mesh->point(vj));
				Eigen::Vector3d qij_e (qij[0],qij[1],qij[2]);
                Vector3d pij = OtoE(pi - mesh.point(vj));
                Eigen::Vector3d pij_e(pij[0],pij[1],pij[2]);
                double wij = ref.w[ei];
                //assert(abs(cs[i]) != 0);
                double ds = 0;
                ds += wij*pij.squaredNorm();
                ds -= wi*wij*(Rs[i] * qij_e).dot(pij_e) * 2;
				ds += wi*wij*qij_e.transpose()*sRs[i]*qij_e;
                //ds += (fv.s[i] * qij).squaredNorm() * cs[i];
				//ds += 
//#ifndef NEWWEIGHT
//              residual += ds * wij * wij;
//#else
//				residual += ds * wj * wij;
//#endif
                residual += ds;    
            }
        }
DOMP_END;

        for (int i=0; i<rss.size(); i++) residual += rss[i];

    }
return residual;
}

double ARAPDeform::getRS_nr(FeatureVector &fv, DTriMesh& mesh)
{
	double rs = getRS(fv,mesh);

	vector<double> rss(omp_get_num_procs(), 0);
#ifdef DUSE_OPENMP
#pragma omp parallel
	{
#pragma omp for
#endif
		for (int i = 0; i < mesh.n_vertices(); i++)
		{
			double& residual = rss[omp_get_thread_num()];
			if (fv.isNr[i])
			{
				DTriMesh::VertexHandle vi(i);
				DTriMesh::Point pi = mesh.point(vi);
				residual += lambda*(fv.nrPoint[i]-OtoE(pi)).squaredNorm();
			}
			if (fv.isNrPlane[i])
			{
				DTriMesh::VertexHandle vi(i);
				DTriMesh::Point pi = mesh.point(vi);
				double _val = (fv.nrPoint[i]-OtoE(pi)).dot(fv.planePoint[i]);
				residual += lambdaplane*_val*_val;
			}
		}
	DOMP_END;
	for(int i = 0; i < rss.size(); i++)
		rs += rss[i];
	return rs;
}



double ARAPDeform::globalOptimize_nr(FeatureVector &fv, DTriMesh &mesh) {
    if (!lqsolver.saveA) return globalOptimizeFast_nr(fv, mesh);

    std::cout << "Global Optimize ...... " << endl;
    long long tt = clock();
    vector< pair< pair<int,int>, double > > data;
    vector<double> vb;
    int n=0;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        if (fv.isConst[i]) {
            for (int dim=0; dim<3; dim++) {
                data.push_back( make_pair( make_pair(n++, i*3+dim), 1 ));
                vb.push_back(0);
            }
        }
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            //double wij = ref.w[ei];
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            DTriMesh::Point qj = ref.mesh->point(vj);
            int ej = ref.d[j];
            for (DTriMesh::VertexVertexIter vvj = mesh.vv_iter(vj); vvj; vvj++, ej++) {
                int k = vvj.handle().idx();
                double wjk = ref.w[ej];
                DTriMesh::Point tqjk = ref.mesh->point(DTriMesh::VertexHandle(k)) - qj;
                Eigen::Vector3d qjk(tqjk[0],tqjk[1],tqjk[2]);
                Eigen::Vector3d c = ri*(rij*(sj*qjk));
                if (fv.isConst[k]) c -= fv.constPoint[k];
                if (fv.isConst[j]) c += fv.constPoint[j];
                for (int dim=0; dim<3; dim++) {
                    if (!fv.isConst[k])
                        data.push_back(make_pair( make_pair(n, k*3+dim), wij * wjk ));
                    if (!fv.isConst[j])
                        data.push_back(make_pair( make_pair(n, j*3+dim), - wij * wjk ));
                    vb.push_back( wij * wjk * c(dim) );
                    n++;
                }
            }
        }
    }
	//soft constraints
	for(int i = 0; i < fv.isNr.size(); i++)
	{
		if (fv.isNr[i])
		{
			std::cout << lambda << endl;
	       for(int dim = 0; dim < 3; dim++)
		   {
			   data.push_back(make_pair(make_pair(n,i*3+dim),lambda));
			   vb.push_back(lambda*fv.nrPoint[i](dim));
			   n++;  
		   }          
		}
	}

    vector<double> result;
    std::cout << "matlab solve" << endl;
    lqsolver.needrs = true;
    double rs = lqsolver.solve(n, fv.s.size()*3, data, vb, result);
    if (lqsolver.saveA) return 0;
    lqsolver.needrs = false;
    std::cout << "copy ans" << endl;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        for (int dim=0; dim<3; dim++)
            mesh.point(vi)[dim] = result[i*3+dim];
    }

    std::cout << "!!!Global Optimize Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << rs*rs << endl;
    return rs*rs;
}

double ARAPDeform::gl_globalOptimize_nr(FeatureVector &fv, DTriMesh &mesh) {
    if (!lqsolver.saveA) return gl_globalOptimizeFast_nr(fv, mesh);

	std::cout << "Global Optimize ...... " << endl;
    long long tt = clock();
    vector< pair< pair<int,int>, double > > data;
	vector<double> vb(ref.mesh->n_vertices()*3,0);
	this->vbconst.clear();
	this->vbconst.resize(ref.mesh->n_vertices()*3,0);
    int n=0;

	//vector<Eigen::Matrix3d> rdrs(ref.mesh->n_vertices(),Eigen::Matrix3d::Zero());
	this->rdrs.clear();
	this->rdrs.resize(ref.mesh->n_vertices(),Eigen::Matrix3d::Zero());

	double* weights = new double[ref.mesh->n_vertices()];
	int * debugtick = new int[ref.mesh->n_vertices()];
    for (int i=0; i<mesh.n_vertices(); i++) {

    	memset(weights,0,ref.mesh->n_vertices()*sizeof(double));
	    memset(debugtick,0,ref.mesh->n_vertices()*sizeof(int));

        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
		double wi = 1/(double)ref.degree[i];
        int ei = ref.d[i];
        if (fv.isConst[i]) 
		{
            for (int dim=0; dim<3; dim++)
			{
                //data.push_back( make_pair( make_pair(n++, i*3+dim), 1 ));
				data.push_back( make_pair( make_pair(i*3+dim, i*3+dim), 1 ));
                //vb.push_back(0);
				vb[i*3+dim] = fv.constPoint[i][dim];
            }
        }
		else
		{
		double sumwij = 0, sumwji = 0;
		for (int k = 0; k < ref.rd[i].size(); k++)
		{
			int j = ref.rd[i][k].first;
			int ej = ref.rd[i][k].second;
			weights[j] = ref.getWij(ej);
			sumwji+=ref.getWij(ej);
			//
			Eigen::Matrix3d &rji = fv.dr[ej];
			Eigen::Matrix3d &rj = fv.r[j];
			this->rdrs[i]+= (rj*rji);
			debugtick[j]++;
		}

		Eigen::Matrix3d &si = fv.s[i];
        rdrs[i] = rdrs[i]*si;
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            //double wij = ref.w[ei];
			double wj = 1/(double)ref.degree[j];
            double wij = ref.getWij(ei);
			sumwij+=wij;
			weights[j] += ref.getWij(ei);
			debugtick[j]++;
			assert(debugtick[j]==2);
			if (!fv.isConst[j])
			{
       			for (int dim = 0; dim < 3; dim++)
     			{
		    		data.push_back(make_pair( make_pair(i*3+dim, j*3+dim), -weights[j]));
	    		}			
			}
			else
			{
				Eigen::Vector3d vec = weights[j]*fv.constPoint[j];
				for (int dim = 0; dim < 3; dim++)
				{
					//vb[j*3+dim]+=vec[dim];
					vbconst[i*3+dim]+=vec[dim];
				}
				//vb[j]+=(weights[j]*fv.constPoint[j]);
			}
            //Eigen::Matrix3d &rij = fv.dr[ei];
            //Eigen::Matrix3d &sj = fv.s[j];
			//rdrs[i] += (ri*rij*sj);
		}

		double wii = sumwij+sumwji;
	    if(fv.isNr[i])
		{
			wii+=lambda;
			for(int dim = 0; dim < 3; dim++)
			{
				vbconst[i*3+dim]+=(lambda*fv.nrPoint[i](dim));
			}
		}
		Eigen::Matrix3d nnt;
		if (fv.isNrPlane[i])
		{
			Eigen::Vector3d _normal = fv.planePoint[i];
			_normal.normalize();
			nnt = _normal*_normal.transpose();
 			for (int dim1 = 0; dim1 < 3; dim1++)
 			{
                 for (int dim2 = 0; dim2 < 3; dim2++)
                 {
					  double _val = lambdaplane*nnt(dim1,dim2);
                      if (dim1 == dim2)
                      {
						  _val+=wii;
                      }
					  data.push_back(make_pair(make_pair(i*3+dim1,i*3+dim2),_val));
                 }
 			}     
			Eigen::Vector3d nv = nnt*fv.nrPoint[i];
			for (int dim = 0; dim < 3; dim++)
			{
				vbconst[i*3+dim]+=(lambdaplane*nv(dim));
			}
		}
		else
		{
			for (int dim = 0; dim < 3; dim++)
			{
				//data.push_back( make_pair( make_pair(i*3+dim, i*3+dim), sumwij+sumwji ));
				data.push_back( make_pair( make_pair(i*3+dim, i*3+dim), wii ));
			}
		}
		}
	 }
	 delete[] weights;
	 weights = NULL;
	 delete[] debugtick;
	 debugtick = NULL;
	 
	//assign b
	//parallel

    this->assignvb(fv,mesh,vb);

	//soft constraints
	/*
	for(int i = 0; i < fv.isNr.size(); i++)
	{
		if (fv.isNr[i])
		{
	       for(int dim = 0; dim < 3; dim++)
		   {
			   vb[3*i+dim]+=(lambda*fv.nrPoint[i](dim));
		   }          
		}
	}
	*/
    vector<double> result;
    std::cout << "matlab solve" << endl;
    lqsolver.needrs = true;
    double rs = lqsolver.glsolve(fv.s.size()*3, fv.s.size()*3, data, vb, result);
    if (lqsolver.saveA) return 0;
    lqsolver.needrs = false;


    std::cout << "copy ans" << endl;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        for (int dim=0; dim<3; dim++)
            mesh.point(vi)[dim] = result[i*3+dim];
    }
	rs = this->getRS_nr(fv,mesh);
	std::cout << "!!!Global Optimize Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << rs << endl;
	return rs;
}



double ARAPDeform::gl_globalOptimizeFast(FeatureVector &fv, DTriMesh &mesh)
{
	std::cout << "Global Optimize Fast ...... " << endl;
    long long tt = clock();
	rdrs.clear();
    this->rdrs.resize(ref.mesh->n_vertices(),Eigen::Matrix3d::Zero());
	vector<double> vb(ref.mesh->n_vertices()*3,0);
	tic
#ifdef DUSE_OPENMP
#pragma omp parallel
	{
#pragma omp for
#endif
    for (int i = 0; i < mesh.n_vertices(); i++)
    {
		 DTriMesh::VertexHandle vi(i);
		 int ei = ref.d[i];
		 Eigen::Matrix3d &ri = fv.r[i];
//		 for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) 
//		 {
//			 int j = vvi.handle().idx();
//			 Eigen::Matrix3d &rij = fv.dr[ei];
//			 Eigen::Matrix3d &sj = fv.s[j];
//			 rdrs[i] += (ri*rij*sj);
//		 }
		 for (int k = 0; k < ref.rd[i].size(); k++)
		 {
			 int j = ref.rd[i][k].first;
			 int eji = ref.rd[i][k].second;
			 Eigen::Matrix3d &rji = fv.dr[eji];			 
			 Eigen::Matrix3d &rj = fv.r[j];
			 rdrs[i]+=(rj*rji);
		 }
		 Eigen::Matrix3d &si = fv.s[i];
		 rdrs[i] = rdrs[i]*si;
    }
 	DOMP_END;
	std::cout<<"global time difference: "<<endl;
	tocp
	//assign rdrs
	this->assignvb(fv,mesh,vb);
	vector<double> dataB(mesh.n_vertices()*3), result;
/*
#ifdef DUSE_OPENMP
#pragma omp parallel
	{
#pragma omp for
#endif
    for (int i = 0; i < mesh.n_vertices(); i++)
    {

		 for(int dim = 0; dim < 3; dim++)
		 {
			 vb[3*i+dim]+=vbconst[3*i+dim];
             //dataB[3*i+dim] = vb[i][dim];
		 }
    }
DOMP_END;
*/
	dataB = vb;
#ifdef _DEBUG
//	DTriMesh testmesh;
//	OpenMesh::IO::read_mesh(testmesh,"E:/SIGA2014/recons/barsim2/2.obj");
//	vector<double> vbmesh(mesh.n_vertices()*3);
//	for(int i = 0; i < mesh.n_vertices(); i++)
//	{
//		for(int dim = 0; dim < 3; dim++)
//		{
//			vbmesh[3*i+dim] = ref.mesh->point(DTriMesh::VertexHandle(i))[dim];
//		}
//	}
//	DMMatrix ref(*eng, "refb", vbmesh.size(), 1, &vbmesh[0]);
#endif

	lqsolver.glsolve(dataB.size(),dataB,result);
	for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        mesh.point(vi)[0] = result[i*3+0];
        mesh.point(vi)[1] = result[i*3+1];
        mesh.point(vi)[2] = result[i*3+2];
    }

	double residual = 0;
    if (needCalcRs) 
	{
		residual = getRS(fv,mesh);
	}
	std::cout << "Global Optimize Fast Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << residual << endl;
    return residual;
}


double ARAPDeform::gl_globalOptimizeFast_nr(FeatureVector &fv, DTriMesh &mesh)
{
	std::cout << "Global Optimize Fast ...... " << endl;
    long long tt = clock();
	rdrs.clear();
    this->rdrs.resize(ref.mesh->n_vertices(),Eigen::Matrix3d::Zero());
	vector<double> vb(ref.mesh->n_vertices()*3,0);
#ifdef DUSE_OPENMP
#pragma omp parallel
	{
#pragma omp for
#endif
    for (int i = 0; i < mesh.n_vertices(); i++)
    {
		 DTriMesh::VertexHandle vi(i);
		 int ei = ref.d[i];
		 Eigen::Matrix3d &ri = fv.r[i];
		 for (int k = 0; k < ref.rd[i].size(); k++)
		 {
			 int j = ref.rd[i][k].first;
			 int eji = ref.rd[i][k].second;
			 Eigen::Matrix3d &rji = fv.dr[eji];			 
			 Eigen::Matrix3d &rj = fv.r[j];
			 rdrs[i]+=(rj*rji);
		 }
		 Eigen::Matrix3d &si = fv.s[i];
		 rdrs[i] = rdrs[i]*si;
    }
 	DOMP_END;
	//assign rdrs
	this->assignvb(fv,mesh,vb);
	vector<double> dataB(mesh.n_vertices()*3), result;
	/*
#ifdef DUSE_OPENMP
#pragma omp parallel
	{
#pragma omp for
#endif
    for (int i = 0; i < mesh.n_vertices(); i++)
    {

		 for(int dim = 0; dim < 3; dim++)
		 {
			 vb[3*i+dim]+=vbconst[3*i+dim];
		 }
		 //if(fv.isNr[i])
		 //{
		 //	 for(int dim = 0; dim < 3; dim++)
		 //	 vb[3*i+dim]+=lambda*(fv.nrPoint[i](dim));
		 //}
    }
DOMP_END;
*/
	dataB = vb;
	lqsolver.glsolve(dataB.size(),dataB,result);
	for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        mesh.point(vi)[0] = result[i*3+0];
        mesh.point(vi)[1] = result[i*3+1];
        mesh.point(vi)[2] = result[i*3+2];
    }

	double residual = 0;
    if (needCalcRs) 
	{
		residual = getRS_nr(fv,mesh);
		std::cout<<coutcmd::green<<residual<<coutcmd::white<<endl;
	}
	std::cout << "Global Optimize Fast Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << residual << endl;
    return residual;

}





double ARAPDeform::globalOptimizeFast(FeatureVector &fv, DTriMesh &mesh) {
#ifdef GLNEW
	return gl_globalOptimizeFast(fv,mesh);
#endif


    std::cout << "Global Optimize Fast ...... " << endl;
    long long tt = clock();
    vector<Vector3d> cv(fv.s.size()), rcv(fv.logdr.size(), Vector3d::Zero());
    vector<double> cs(fv.s.size());
    vector<Matrix3d> cm(fv.s.size());
    auto &rd = ref.rd;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int j=0; j<mesh.n_vertices(); j++) {
        Eigen::Matrix3d &cmj = cm[j];
        double &csj = cs[j];

        cv[j] = Vector3d::Zero();
        cm[j] = Matrix3d::Zero();
        cs[j] = 0;
		double wj = 1/(double)(ref.degree[j]);
        for (int k=0; k<rd[j].size(); k++) {
            int i = rd[j][k].first;
            int ei = rd[j][k].second;

            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
#ifdef NEWWEIGHT
			cmj += (wj*wij) * fv.r[i] * rij * sj;
			csj += (wj*wij);
#else
            cmj += (wij*wij) * fv.r[i] * rij * sj;
            csj += (wij*wij);
#endif

        }
    }
/*
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            cm[j] += (wij*wij) * ri * rij * sj;
            cs[j] += (wij*wij);
        }
    }
*/
DOMP_END;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int ei = ref.d[i];
		double wi = 1/(double)ref.degree[i];
        DTriMesh::Point qi = ref.mesh->point(vi);
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            Vector3d qij = OtoE(ref.mesh->point(vj) - qi);
            double wij = ref.w[ei];
			double wj = 1/(double)ref.degree[j];
#ifndef NEWWEIGHT
			Eigen::Vector3d c = cm[i] * (qij * (wij*wij));
            if (fv.isConst[i]) c += fv.constPoint[i] * (cs[i] * (wij*wij));
            if (fv.isConst[j]) c -= fv.constPoint[j] * (cs[i] * (wij*wij));
#else
			Eigen::Vector3d c = cm[i] * (qij*(wi * wij));
			if (fv.isConst[i]) c += fv.constPoint[i] * (cs[i] * (wi * wij));
            if (fv.isConst[j]) c -= fv.constPoint[j] * (cs[i] * (wi * wij));
#endif
            
            if (!fv.isConst[i]) cv[i] -= c;
            //if (!fv.isConst[j]) cv[j] += c;
            if (!fv.isConst[j]) rcv[ei] += c;
        }
    }

DOMP_END;

    vector<double> dataB(cv.size()*3), result;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int j=0; j<mesh.n_vertices(); j++) {
        for (int k=0; k<rd[j].size(); k++) {
            int ei = rd[j][k].second;
            cv[j] += rcv[ei];
        }
        dataB[j*3] = cv[j](0);
        dataB[j*3+1] = cv[j](1);
        dataB[j*3+2] = cv[j](2);
    }

DOMP_END;


    lqsolver.solve(dataB.size(), dataB, result);

    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        mesh.point(vi)[0] = result[i*3+0];
        mesh.point(vi)[1] = result[i*3+1];
        mesh.point(vi)[2] = result[i*3+2];
    }

    double residual = 0;
    if (needCalcRs) {
        vector<double> rss(omp_get_num_procs(), 0);

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

        for (int i=0; i<mesh.n_vertices(); i++) {
            double &residual = rss[omp_get_thread_num()];
            DTriMesh::VertexHandle vi(i);
            int ei = ref.d[i];
            DTriMesh::Point qi = ref.mesh->point(vi);
            DTriMesh::Point pi = mesh.point(vi);
            for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
                int j = vvi.handle().idx();
				double wj = 1/(double)ref.degree[j];
                DTriMesh::VertexHandle vj(j);
                Vector3d qij = OtoE(ref.mesh->point(vj) - qi);
                Vector3d pij = OtoE(mesh.point(vj) - pi);
                double wij = ref.w[ei];
                assert(abs(cs[i]) != 0);
                double ds = 0;
                ds += pij.squaredNorm() * cs[i];
                ds -= (cm[i] * qij).dot(pij) * 2;
                ds += (fv.s[i] * qij).squaredNorm() * cs[i];
#ifndef NEWWEIGHT
                residual += ds * wij * wij;
#else
				residual += ds * wj * wij;
#endif

            }
        }
DOMP_END;

        for (int i=0; i<rss.size(); i++) residual += rss[i];

    }

    std::cout << "Global Optimize Fast Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << residual << endl;
    return residual;
}

double ARAPDeform::globalOptimizeFast_nr(FeatureVector &fv, DTriMesh &mesh) {
    std::cout << "Global Optimize Fast ...... " << endl;
    long long tt = clock();
    vector<Vector3d> cv(fv.s.size()), rcv(fv.logdr.size(), Vector3d::Zero());
    vector<double> cs(fv.s.size());
    vector<Matrix3d> cm(fv.s.size());
    auto &rd = ref.rd;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int j=0; j<mesh.n_vertices(); j++) {
        Eigen::Matrix3d &cmj = cm[j];
        double &csj = cs[j];

        cv[j] = Vector3d::Zero();
        cm[j] = Matrix3d::Zero();
        cs[j] = 0;

        for (int k=0; k<rd[j].size(); k++) {
            int i = rd[j][k].first;
            int ei = rd[j][k].second;

            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];

            cmj += (wij*wij) * fv.r[i] * rij * sj;
            csj += (wij*wij);
        }
    }
/*
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            cm[j] += (wij*wij) * ri * rij * sj;
            cs[j] += (wij*wij);
        }
    }
*/
DOMP_END;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int ei = ref.d[i];
        DTriMesh::Point qi = ref.mesh->point(vi);
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            Vector3d qij = OtoE(ref.mesh->point(vj) - qi);
            double wij = ref.w[ei];
            Eigen::Vector3d c = cm[i] * (qij * (wij*wij));
            if (fv.isConst[i]) c += fv.constPoint[i] * (cs[i] * (wij*wij));
            if (fv.isConst[j]) c -= fv.constPoint[j] * (cs[i] * (wij*wij));
            if (!fv.isConst[i]) cv[i] -= c;
            //if (!fv.isConst[j]) cv[j] += c;
            if (!fv.isConst[j]) rcv[ei] += c;
        }
    }

DOMP_END;

    vector<double> dataB(cv.size()*3), result;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int j=0; j<mesh.n_vertices(); j++) {
        for (int k=0; k<rd[j].size(); k++) {
            int ei = rd[j][k].second;
            cv[j] += rcv[ei];
        }
        dataB[j*3] = cv[j](0);
        dataB[j*3+1] = cv[j](1);
        dataB[j*3+2] = cv[j](2);
    }

DOMP_END;

    for(int i = 0; i < mesh.n_vertices(); i++)
	{
		if(fv.isNr[i])
		{
			for(int dim = 0; dim < 3; dim++)
			{
				//dataB.push_back(lambda*fv.nrPoint[i](dim));
                dataB[i*3+dim] += lambda*lambda*fv.nrPoint[i](dim);
			}
		}
	}

    lqsolver.solve(dataB.size(), dataB, result);

    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        mesh.point(vi)[0] = result[i*3+0];
        mesh.point(vi)[1] = result[i*3+1];
        mesh.point(vi)[2] = result[i*3+2];
    }

    double residual = 0;
    if (needCalcRs) {
        vector<double> rss(omp_get_num_procs(), 0);

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

        for (int i=0; i<mesh.n_vertices(); i++) {
            double &residual = rss[omp_get_thread_num()];
            DTriMesh::VertexHandle vi(i);
            int ei = ref.d[i];
            DTriMesh::Point qi = ref.mesh->point(vi);
            DTriMesh::Point pi = mesh.point(vi);
            for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
                int j = vvi.handle().idx();
                DTriMesh::VertexHandle vj(j);
                Vector3d qij = OtoE(ref.mesh->point(vj) - qi);
                Vector3d pij = OtoE(mesh.point(vj) - pi);
                double wij = ref.w[ei];
                assert(abs(cs[i]) != 0);
                double ds = 0;
                ds += pij.squaredNorm() * cs[i];
                ds -= (cm[i] * qij).dot(pij) * 2;
                ds += (fv.s[i] * qij).squaredNorm() * cs[i];
                residual += ds * wij * wij;
				if (fv.isNr[i])
				{
				    residual += lambda*(OtoE(pi)-fv.nrPoint[i]).squaredNorm();
				}
            }
        }
DOMP_END;

        for (int i=0; i<rss.size(); i++) residual += rss[i];

    }

    std::cout << "Global Optimize Fast Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << residual << endl;
    return residual;
}



double ARAPDeform::localOptimize(FeatureVector &fv, DTriMesh &mesh) {
    localOptimizeFast(fv, mesh);
    return 0;

    std::cout << "Local Optimize ...... " << endl;
    long long tt = clock();
    double residual = 0;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        vector<Eigen::Vector3d> vq, vp;
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
			double wj = 1/(double)ref.degree[j];
            DTriMesh::VertexHandle vj(j);
            //double wij = ref.w[ei];
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            DTriMesh::Point qj = ref.mesh->point(vj);
            DTriMesh::Point pj = mesh.point(vj);
            int ej = ref.d[j];
            for (DTriMesh::VertexVertexIter vvj = mesh.vv_iter(vj); vvj; vvj++, ej++) {
                int k = vvj.handle().idx();
                double wjk = ref.w[ej];
                DTriMesh::Point tqjk = ref.mesh->point(DTriMesh::VertexHandle(k)) - qj;
                DTriMesh::Point tpjk = mesh.point(DTriMesh::VertexHandle(k)) - pj;
                Eigen::Vector3d qjk(tqjk[0],tqjk[1],tqjk[2]);
                Eigen::Vector3d pjk(tpjk[0],tpjk[1],tpjk[2]);
                qjk = (rij*(sj*qjk));
#ifndef NEWWEIGHT
                vq.push_back(qjk * (wij*wjk));
                vp.push_back(pjk * (wij*wjk));
#else
			    vq.push_back(qjk * (wj*wjk));
                vp.push_back(pjk * (wj*wjk));
#endif
            }
        }
        RotateAlign align(vp, vq);
        ri = align.calc();
        residual += align.res;
    }
    std::cout << "!!!Local Optimize Done , residual : " << residual << " Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << endl;
    return residual;
}

void ARAPDeform::localOptimizeFast(FeatureVector &fv, DTriMesh &mesh) {
    std::cout << "Local Optimize Fast...... " << endl;
    long long tt = clock();
    double residual = 0;
    vector<Eigen::Matrix3d> mats(fv.s.size());
tic
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int ei = ref.d[i];
        Eigen::Matrix3d &mat = mats[i];
        mat = Matrix3d::Zero();
        DTriMesh::Point qi = ref.mesh->point(vi);
        DTriMesh::Point pi = mesh.point(vi);
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            double wij = ref.w[ei];
            Vector3d tqij = OtoE(ref.mesh->point(vj) - qi);
            Vector3d tpij = OtoE(mesh.point(vj) - pi);
#ifdef NEWWEIGHT
			mat += (wij * tqij) * tpij.transpose();
#else
			mat += ((wij * wij) * tqij) * tpij.transpose();
#endif            
        }
    }
DOMP_END;
std::cout<<"Time local difference: "<<endl;
tocp
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        Matrix3d mat(Matrix3d::Zero());
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
			double wj = 1/(double)ref.degree[j];
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
#ifdef NEWWEIGHT
			mat += wj* rij*sj*mats[j];
#else
			mat += rij * sj * mats[j] * (wij * wij);
#endif
            
        }
        polarDec(mat, ri);
        ri.transposeInPlace();
    }
DOMP_END;

    std::cout << "Local Optimize Fast Done , " << " Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << endl;
}


double ARAPDeform::getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh) {
    std::cout << "GetWR : ";
    for (int i=0; i<weight.size(); i++)
        std::cout << weight[i] << " ";
    std::cout << endl;
    fv.blendFrom(weight, this->fvs);
    initMesh(fv, mesh);
    return globalOptimize(fv, mesh);
}


double ARAPDeform::getWeightResidual_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh) {
    std::cout << "GetWR : ";
    for (int i=0; i<weight.size(); i++)
        std::cout << weight[i] << " ";
    std::cout << endl;
    fv.blendFrom(weight, this->fvs);
    initMesh(fv, mesh);
    return globalOptimize_nr(fv, mesh);
}

double ARAPDeform::gl_getWeightResidual_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh) {
    std::cout << "GetWR : ";
    for (int i=0; i<weight.size(); i++)
        std::cout << weight[i] << " ";
    std::cout << endl;
    fv.blendFrom(weight, this->fvs);
    initMesh(fv, mesh);
    return gl_globalOptimize_nr(fv, mesh);
}

std::vector<double> ARAPDeform::getGradient(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult) {
    vector<double> g, w=weight;
    double rs = getWeightResidual(fv, weight, mesh);
    zeroResult = rs;
    double eps = 1e-6;
    for (int i=1; i<weight.size(); i++) {
        w[0] -= eps;
        w[i] += eps;

        double rs2 = getWeightResidual(fv, w, mesh);
        g.push_back(- (rs2 - rs) / eps);

        w[0] += eps;
        w[i] -= eps;
    }
    return g;
}



std::vector<double> ARAPDeform::getGradient_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult) {
    vector<double> g, w=weight;
    double rs = getWeightResidual_nr(fv, weight, mesh);
    zeroResult = rs;
    double eps = 1e-6;
    for (int i=1; i<weight.size(); i++) {
        w[0] -= eps;
        w[i] += eps;

        double rs2 = getWeightResidual_nr(fv, w, mesh);
        g.push_back(- (rs2 - rs) / eps);

        w[0] += eps;
        w[i] -= eps;
    }
    return g;
}

std::vector<double> ARAPDeform::gl_getGradient_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult) {
    vector<double> g, w=weight;
    double rs = gl_getWeightResidual_nr(fv, weight, mesh);
    zeroResult = rs;
    double eps = 1e-6;
    for (int i=1; i<weight.size(); i++) {
        w[0] -= eps;
        w[i] += eps;

        double rs2 = gl_getWeightResidual_nr(fv, w, mesh);
        g.push_back(- (rs2 - rs) / eps);

        w[0] += eps;
        w[i] -= eps;
    }
    return g;
}




double unifromGradient(std::vector<double> &weight) {
    double ans=0,s=0;
    for (int i=0; i<weight.size(); i++) {
        ans += weight[i] * weight[i];
        s += weight[i];
    }
    weight.insert(weight.begin(), -s);
    s = 1/sqrt(s*s + ans);
    for (int i=0; i<weight.size(); i++) weight[i] *= s;
    return sqrt(ans);
}

void ARAPDeform::solve(std::vector<double> weight, DTriMesh &mesh) {
    std::cout << "blend weight" <<endl;
    FeatureVector fv(weight, fvs);
    std::cout << "solve fv" << endl;
    solve(fv, mesh);
}

void ARAPDeform::ckAlign(DTriMesh &mesh) {
    if (needAlign)
        RotateAlign::AlignAtoB(mesh, *(ref.mesh));
}

void ARAPDeform::solve(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    initWeight(weight);
    fv.blendFrom(weight, fvs);
    initMesh(fv, mesh);

    globalOptimize(fv, mesh);

    double preRes = 1e10;
    double res = 1e9;
    int iter = 0;
    writeIterMesh(mesh, "iter", iter);

    while (abs(preRes - res) > iterEps && iter < maxIterTime) {
        iter ++;
        std::cout << "Iter time : " << iter << endl;
        preRes = res;

        solve(mesh, weight);
        fv.blendFrom(weight, fvs);

        res = localOptimize(fv, mesh);

        globalOptimize(fv, mesh);

        writeIterMesh(mesh, "iter", iter);
        std::cout << "DResidual : " << abs(preRes - res) << endl;
    }

    ckAlign(mesh);
}

std::vector<double> muladd(double x, std::vector<double> a, std::vector<double> b) {
    vector<double> ans(a.size());
    for (int i=0; i<a.size(); i++)
        ans[i] = x*a[i] + b[i];
    return ans;
}

void ARAPDeform::solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    preDecMatrix(fv, mesh);
    int t=0;
    double teps = this->newtonIterEps;
    double best = sqrt(weight.size()*1.), num;
    TriNewtonIterSolver cvs;
    std::vector<double> g;
    cvs.eps = teps;
    cvs.convexFunction = [&](double x) -> double {
        return getWeightResidual(fv, muladd(x,g,weight), mesh);
    };
    if (weight.size() > 1)
    while (t++ < newtonIterTime) {
        double l=-abs(best)/sqrt(weight.size()*1.)/2, r=abs(best);
        g = getGradient(fv, weight, mesh, cvs.zeroResult);
        double wr = unifromGradient(g);
        std::cout << "Gradient Norm : " << wr << endl;
        std::cout << "Uni-Gradient Vector : ";
        for (int i=0; i<g.size(); i++)
            std::cout << g[i] << " ";
        std::cout << endl;
        best = 0;
        cvs.solve(l, r, best, num);
        std::cout << "Best : " << best << endl;
        weight = muladd(best,g,weight);
        std::cout << "Weight : ";
        for (int i=0; i<weight.size(); i++)
            std::cout << weight[i] << " ";
        std::cout << endl;
        writeIterMesh(mesh, "sf", t);
        if (abs(best) < teps) break;
    }
    fv.blendFrom(weight, fvs);
    solve(fv, mesh, false);
}

void ARAPDeform::solve3(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
	preDecMatrix(fv, mesh);
	int t=0;
	double teps = this->newtonIterEps;
	double best = sqrt(weight.size()*1.), num;
	TriNewtonIterSolver cvs;
	std::vector<double> g;
	cvs.eps = teps;
	cvs.convexFunction = [&](double x) -> double {
		return getWeightResidual(fv, muladd(x,g,weight), mesh);
	};
	if (weight.size() > 1)
		while (t++ < newtonIterTime) {
			double l=-abs(best)/sqrt(weight.size()*1.)/2, r=abs(best);
			g = getGradient(fv, weight, mesh, cvs.zeroResult);
			double wr = unifromGradient(g);
			std::cout << "Gradient Norm : " << wr << endl;
			std::cout << "Uni-Gradient Vector : ";
			for (int i=0; i<g.size(); i++)
				std::cout << g[i] << " ";
			std::cout << endl;
			best = 0;
			cvs.solve(l, r, best, num);
			std::cout << "Best : " << best << endl;
			weight = muladd(best,g,weight);
			std::cout << "Weight : ";
			for (int i=0; i<weight.size(); i++)
				std::cout << weight[i] << " ";
			std::cout << endl;
			writeIterMesh(mesh, "sf", t);
			if (abs(best) < teps) break;
		}
		fv.blendFrom(weight, fvs);
		int cnum = fv.isConst.size();
		fv.isConst.resize(cnum,0);
		fv.constPoint.resize(cnum);
		solve(fv, mesh, false);
}

void ARAPDeform::solve_nr(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
	DTriMesh meshtmp = mesh;
 	preDecMatrix_nr(fv, mesh);
 	int t=0;
 	double teps = this->newtonIterEps;
 	double best = sqrt(weight.size()*1.), num;
 	TriNewtonIterSolver cvs;
 	std::vector<double> g;
 	cvs.eps = teps;
 	cvs.convexFunction = [&](double x) -> double {
 		return getWeightResidual_nr(fv, muladd(x,g,weight), mesh);
 	};
	if (weight.size() > 1)
		while (t++ < newtonIterTime) {
			double l = -abs(best) / sqrt(weight.size()*1.) / 2, r = abs(best);
			g = getGradient_nr(fv, weight, mesh, cvs.zeroResult);
			double wr = unifromGradient(g);
			std::cout << "Gradient Norm : " << wr << endl;
			std::cout << "Uni-Gradient Vector : ";
			for (int i = 0; i < g.size(); i++)
				std::cout << g[i] << " ";
			std::cout << endl;
			best = 0;
			cvs.solve(l, r, best, num);
			std::cout << "Best : " << best << endl;
			weight = muladd(best, g, weight);
			std::cout << "Weight : ";
			for (int i = 0; i < weight.size(); i++)
				std::cout << weight[i] << " ";
			std::cout << endl;
			writeIterMesh(mesh, "sf", t);
			if (abs(best) < teps) break;
		}
	fv.blendFrom(weight, fvs);
	int cnum = fv.isConst.size();
	//fv.isConst.resize(cnum,0);
	//fv.constPoint.resize(cnum);
//		solve(fv, meshtmp, false);
	this->AtAChanged = true;
	solve(fv, meshtmp, true);
	mesh = meshtmp;
}

void ARAPDeform::gl_solve_nr(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
	DTriMesh meshtmp = mesh;
	gl_preDecMatrix_nr(fv, mesh);
	int t=0;
	double teps = this->newtonIterEps;
	double best = sqrt(weight.size()*1.), num;
	TriNewtonIterSolver cvs;
	std::vector<double> g;
	cvs.eps = teps;
	cvs.convexFunction = [&](double x) -> double {
		return gl_getWeightResidual_nr(fv, muladd(x,g,weight), mesh);
	};

	//newtonIterTime = 0;

	if (weight.size() > 1)
		while (t++ < newtonIterTime) {
			double l=-abs(best)/sqrt(weight.size()*1.)/2, r=abs(best);
			g = gl_getGradient_nr(fv, weight, mesh, cvs.zeroResult);
			double wr = unifromGradient(g);
			std::cout << "Gradient Norm : " << wr << endl;
			std::cout << "Uni-Gradient Vector : ";
			for (int i=0; i<g.size(); i++)
				std::cout << g[i] << " ";
			std::cout << endl;
			best = 0;
			cvs.solve(l, r, best, num);
			std::cout << "Best : " << best << endl;
			weight = muladd(best,g,weight);
			std::cout << "Weight : ";
			for (int i=0; i<weight.size(); i++)
				std::cout << weight[i] << " ";
			std::cout << endl;
			writeIterMesh(mesh, "sf", t);

			FeatureVector fvtst = fv;
			fvtst.blendFrom(weight, fvs);
		    //gl_solve(fvtst, meshtmp, true);
			gl_globalOptimize_nr(fvtst,meshtmp);
			double rstst = this->getRS_nr(fvtst,meshtmp);
			std::cout<<coutcmd::yellow<<"Energy:"<<rstst<<coutcmd::white<<endl;

			if (abs(best) < teps) break;
		}
		fv.blendFrom(weight, fvs);
		int cnum = fv.isConst.size();
		//fv.isConst.resize(cnum,0);
		//fv.constPoint.resize(cnum);
		//fv.isNr.clear();
		//fv.isNr.resize(cnum,0);
		//		solve(fv, meshtmp, false);
		//this->AtAChanged = true;
		//solve(fv, meshtmp, true);
         
		//this->AtAChanged = true;
		//solve(fv,meshtmp,true);
		//this->initMesh(fv,meshtmp);
		//gl_globalOptimize_nr(fv,meshtmp);
		gl_solve(fv,meshtmp,false);
		mesh = meshtmp;
}




void ARAPDeform::solve(FeatureVector &fv, std::vector<double> &weight) {
    initWeight(weight);

    vector< pair< pair<int,int>, double > > data;
    vector<double> vb;
    int n=0;

    for (int i=0; i<fv.s.size(); i++) {
        for (int x=0; x<3; x++) for (int y=0; y<3; y++) {
            vb.push_back(fv.s[i](x,y));
            for (int j=0; j<fvs.size(); j++)
                data.push_back( make_pair( make_pair(n, j), fvs[j].s[i](x,y) ) );
            n++;
        }
    }

    for (int i=0; i<fv.logdr.size(); i++) {
        for (int x=0; x<3; x++) for (int y=0; y<3; y++) {
            vb.push_back(fv.logdr[i](x,y));
            for (int j=0; j<fvs.size(); j++)
                data.push_back( make_pair( make_pair(n, j), fvs[j].logdr[i](x,y) ) );
            n++;
        }
    }

    lqsolver.solve(n, weight.size(), data, vb, weight, true);

    std::cout << "Weights : ";
    double s = 0;
    for (int i=0; i<weight.size(); i++) {
        s += weight[i];
        std::cout << weight[i] << " ";
    }
    std::cout << endl << "S : " << s << endl;
}

void ARAPDeform::leastsquare(FeatureVector &fv, std::vector<double> & weight, DTriMesh &mesh, bool sumone)
{
	initWeight(weight);
    vector< pair< pair<int,int>, double > > data;
    vector<double> vb;
    int n=0;

    for (int i=0; i<fv.s.size(); i++) {
        for (int x=0; x<3; x++) for (int y=0; y<3; y++) {
            vb.push_back(fv.s[i](x,y));
            for (int j=0; j<fvs.size(); j++)
                data.push_back( make_pair( make_pair(n, j), fvs[j].s[i](x,y) ) );
            n++;
        }
    }

    for (int i=0; i<fv.logdr.size(); i++) {
        for (int x=0; x<3; x++) for (int y=0; y<3; y++) {
            vb.push_back(fv.logdr[i](x,y));
            for (int j=0; j<fvs.size(); j++)
                data.push_back( make_pair( make_pair(n, j), fvs[j].logdr[i](x,y) ) );
            n++;
        }
    }

    lqsolver.solve(n, weight.size(), data, vb, weight, sumone);

    std::cout << "Weights : ";
    double s = 0;
    for (int i=0; i<weight.size(); i++) {
        s += weight[i];
        std::cout << weight[i] << " ";
    }
    std::cout << endl << "S : " << s << endl;
	solve(fv, mesh, true);
}


void ARAPDeform::solve(DTriMesh &mesh, std::vector<double> &weight) {
    initWeight(weight);
    FeatureVector fv;
    ref.getFeature(mesh, fv);
    solve(fv, weight);
}

void ARAPDeform::initWeight(std::vector<double> &weight) {
    if (weight.size() != fvs.size()) {
        weight.resize(fvs.size(), 0);
        weight[0] = 1;
    }
}

double T_ARAPDeform::getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh) {
    double prers = ARAPDeform::getWeightResidual(fv, weight, mesh);
    Vector3d s1(Vector3d::Zero()), s2(s1);
    int s=0;
    for (int i=0; i<ev.size(); i++)
        if (ev[i]) {
            s++;
            s1 += fv.constPoint[i];
            s2 += OtoE(mesh.point(DTriMesh::VertexHandle(i)));
        }
    s1 = (s1-s2) / s;

    double rs = 0;
    for (int i=0; i<ev.size(); i++)
        if (ev[i]) {
            std::cout << "fvCP : " << fv.constPoint[i].transpose();
            std::cout << " meshP : " << OtoE(mesh.point(DTriMesh::VertexHandle(i))).transpose();
            std::cout << " ds : " <<(fv.constPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1).transpose();
            std::cout << endl;
            rs += (fv.constPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1).squaredNorm();
        }
    std::cout << "VpRS : " << rs << endl;
    double ws = 0;
    for (int i=0; i<weight.size(); i++) {
        ws += abs(weight[i]);
    }
    std::cout << "ws : " << ws << endl;
    rs = rs * ws + std::max(0.0, rs - 1.2);
    std::cout << "WSRS : " << rs << endl;
    return rs;
}

double T_ARAPDeform::getWeightResidual_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh) {
	double prers = ARAPDeform::getWeightResidual_nr(fv, weight, mesh);
	Vector3d s1(Vector3d::Zero()), s2(s1);
	int s=0;
	for (int i=0; i<ev.size(); i++)
		if (ev[i]) {
			s++;
			s1 += fv.constPoint[i];
			s2 += OtoE(mesh.point(DTriMesh::VertexHandle(i)));
		}

		if(s!=0)
		{
			s1 = (s1-s2) / s;
		}
		
		double rs = 0;
		for (int i=0; i<ev.size(); i++)
			if (ev[i]) {
				std::cout << "fvCP : " << fv.constPoint[i].transpose();
				std::cout << " meshP : " << OtoE(mesh.point(DTriMesh::VertexHandle(i))).transpose();
				std::cout << " ds : " <<(fv.constPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1).transpose();
				std::cout << endl;
				rs += (fv.constPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1).squaredNorm();
			}
			std::cout << "VpRS : " << rs << endl;
			double ws = 0;
			for (int i=0; i<weight.size(); i++) {
				ws += abs(weight[i]);
			}
			std::cout << "ws : " << ws << endl;
			rs = rs * ws + std::max(0.0, rs - 1.2);
			std::cout << "WSRS : " << rs << endl;

			Vector3d s1nr(Vector3d::Zero()), s2nr(s1nr);
			int snr=0;
			for (int i=0; i<fv.isNr.size(); i++)
			if(fv.isNr[i]){
					snr++;
					s1nr += fv.nrPoint[i];
					s2nr += OtoE(mesh.point(DTriMesh::VertexHandle(i)));
			}
			if (snr!=0)
			{
				s1nr = (s1nr-s2nr) / snr;
			}		
			double nrrs = 0;
			for (int i=0; i<fv.isNr.size(); i++)
				if (fv.isNr[i]) 
				{
					//std::cout << "fvCP : " << fv.nrPoint[i].transpose();
					//std::cout << " meshP : " << OtoE(mesh.point(DTriMesh::VertexHandle(fv.isNr[i]))).transpose();
					//std::cout << " ds : " <<(fv.constPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1).transpose();
					//std::cout << endl;
					nrrs += (fv.nrPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1nr).squaredNorm();
				}
			rs+=(lambda*nrrs);
			return rs;
}

void T_ARAPDeform::solve_nr(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh)
{
	ev.resize(fv.isConst.size(), false);
    fv.isConst.swap(ev);
	//ev = fv.isConst;
	DTriMesh meshtmp = mesh;
	preDecMatrix_nr(fv, mesh);
	int t=0;
	double teps = this->newtonIterEps;
	double best = sqrt(weight.size()*1.), num;
	TriNewtonIterSolver cvs;
	std::vector<double> g;
	cvs.eps = teps;
	cvs.convexFunction = [&](double x) -> double {
		return getWeightResidual_nr(fv, muladd(x,g,weight), mesh);
	};
	if (weight.size() > 1)
		while (t++ < newtonIterTime) {
			double l=-abs(best)/sqrt(weight.size()*1.)/2, r=abs(best);
			g = getGradient_nr(fv, weight, mesh, cvs.zeroResult);
			double wr = unifromGradient(g);
			std::cout << "Gradient Norm : " << wr << endl;
			std::cout << "Uni-Gradient Vector : ";
			for (int i=0; i<g.size(); i++)
				std::cout << g[i] << " ";
			std::cout << endl;
			best = 0;
			cvs.solve(l, r, best, num);
			std::cout << "Best : " << best << endl;
			weight = muladd(best,g,weight);
			std::cout << "Weight : ";
			for (int i=0; i<weight.size(); i++)
				std::cout << weight[i] << " ";
			std::cout << endl;
			writeIterMesh(mesh, "sf", t);
			if (abs(best) < teps) break;
		}
		fv.blendFrom(weight, fvs);
		int cnum = fv.isConst.size();
		//fv.isConst.resize(cnum,0);
		//fv.constPoint.resize(cnum);
		//		solve(fv, meshtmp, false);
		fv.isConst.swap(ev);
		this->AtAChanged = true;
		solve(fv, meshtmp, true);
		mesh = meshtmp;
}

void T_ARAPDeform::solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    ev.resize(fv.isConst.size(), false);
    fv.isConst.swap(ev);

    org->AtAChanged = this->AtAChanged;
    org->maxIterTime = this->maxIterTime;
    org->newtonIterEps = this->newtonIterEps;
    org->newtonIterTime = this->newtonIterTime;
    org->iterEps = this->iterEps;

    std::cout << "TSolve2" << endl;

    preDecMatrix(fv, *(this->ref.mesh));
    int t=0;
    double teps = this->newtonIterEps;
    double best = sqrt(sqrt(weight.size()*1.)), num;
    TriNewtonIterSolver cvs;
    std::vector<double> g;
    cvs.eps = teps;
    cvs.convexFunction = [&](double x) -> double {
        return getWeightResidual(fv, muladd(x,g,weight), mesh);
    };
    while (t++ < newtonIterTime) {
        double l=-abs(best)/sqrt(weight.size()*1.)/2, r=abs(best);
        g = getGradient(fv, weight, mesh, cvs.zeroResult);
        double wr = unifromGradient(g);
        std::cout << "Gradient Norm : " << wr << endl;
        std::cout << "Uni-Gradient Vector : ";
        for (int i=0; i<g.size(); i++)
            std::cout << g[i] << " ";
        std::cout << endl;
        best = 0;
        cvs.solve(l, r, best, num);
        std::cout << "Best : " << best << endl;
        weight = muladd(best,g,weight);
        std::cout << "Weight : ";
        for (int i=0; i<weight.size(); i++)
            std::cout << weight[i] << " ";
        std::cout << endl;
        writeIterMesh(mesh, "sf", t);
        if (abs(best) < teps) break;
    }
    fv.isConst.swap(ev);
    org->solve(fv, mesh);
}

//T_ARAPDeform::
//T_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms)
//    : ARAPDeform(eng, refMesh, ms) {
//    org = new ARAPDeform(eng2, refMesh, ms);
//}

T_ARAPDeform::
T_ARAPDeform(DTriMesh &refMesh, std::vector<DTriMesh*> ms)
	: ARAPDeform(refMesh, ms) {
	org = new ARAPDeform(refMesh, ms);
}

// T2 code

void T2_ARAPDeform::solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    ev.resize(fv.isConst.size(), false);
    fv.isConst.swap(ev);

    org->AtAChanged = this->AtAChanged;
    org->maxIterTime = this->maxIterTime;
    org->newtonIterEps = this->newtonIterEps;
    org->newtonIterTime = this->newtonIterTime;
    org->iterEps = this->iterEps;
    this->needCalcRs = false;

    std::cout << "T2Solve2" << endl;

    preDecMatrix(fv, *(this->ref.mesh));
    double teps = this->newtonIterEps;

    LQSolverWithLimit::newtonIterEps = this->newtonIterEps;
    LQSolverWithLimit::newtonMaxIterTime = this->newtonIterTime;
    //LQSolverWithLimit::iterEps = this->iterEps;

    vector<int> constList;
    vector<Vector3d> constPList;
    for (int i=0; i<ev.size(); i++)
        if (ev[i]) {
            constList.push_back(i);
            constPList.push_back(fv.constPoint[i]);
        }

    auto convert = [&](vector<Vector3d> &v) -> VectorXd {
        VectorXd res(v.size() * 3);
        Vector3d s(Vector3d::Zero());
        for (int i=0; i<v.size(); i++) s+=v[i];
        s /= v.size();
        for (int i=0; i<v.size(); i++) {
            auto x = v[i] - s;
            res(i*3) = x(0);
            res(i*3+1) = x(1);
            res(i*3+2) = x(2);
        };
        return res;
    };

    VectorXd targetB = convert(constPList);

    auto getX = [&](void) -> Eigen::VectorXd {
        ARAPDeform::getWeightResidual(fv, weight, mesh);
        Eigen::VectorXd res(constList.size() * 3);
        vector<Vector3d> v(constPList.size());
        for (int i=0; i<constList.size(); i++)
            v[i] = OtoE(mesh.point(DTriMesh::VertexHandle(constList[i])));
        return convert(v);
    };

    LQSolverWithLimit::GetJacobiFunction setJacobi =
            [&](
            std::vector<Eigen::VectorXd> &A,
            Eigen::VectorXd &b,
            std::vector< std::pair<double, double> > &limit
            ) {
        limit.resize(fvs.size()-1);
        for (int i=1; i<fvs.size(); i++)
            limit[i-1] = make_pair(rlimit.first - weight[i], rlimit.second - weight[i]);
        b = getX();
        A.resize(fvs.size()-1);
        double eps2 = 1e-6;
        for (int i=1; i<fvs.size(); i++) {
            weight[0] -= eps2;
            weight[i] += eps2;

            A[i-1] = (getX() - b) / eps2;

            weight[0] += eps2;
            weight[i] -= eps2;
        }
        b = targetB - b;
    };

    LQSolverWithLimit::ReturnRsFunction getRs =
            [&](std::vector<double> &x) {
        std::cout << "Get Delta : ";
        for (int i=0; i<x.size(); i++) {
            weight[i+1] += x[i];
            weight[0] -=  x[i];
            std::cout << weight[i+1] << "(" << x[i] << ") ";
        }
        std::cout << endl;
    };

    if (weight.size() > 1)
    LQSolverWithLimit::NewtonIter(setJacobi, getRs);

    fv.isConst.swap(ev);
// 
//     for (int i = 0; i < fv.s.size(); i++)
//     {
// 		fv.s[i] = fv.s[i]*0.8;
//     }

    org->solve(fv, mesh);
}

//T2_ARAPDeform::
//T2_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms)
//    : ARAPDeform(eng, refMesh, ms) {
//    org = new ARAPDeform(eng2, refMesh, ms);
//}
T2_ARAPDeform::
T2_ARAPDeform(DTriMesh &refMesh, std::vector<DTriMesh*> ms)
	: ARAPDeform(refMesh, ms) {
	org = new ARAPDeform(refMesh, ms);
}
//multi scale operation

bool ARAPDeform::loaddensemesh(const char* filename)
{
	//DTriMesh temp;
	DTriMesh* tmp = new DTriMesh;
	bool _status = OpenMesh::IO::read_mesh(*tmp,filename);
	densemesh.SetMesh(*tmp);
	return _status;
}


bool ARAPDeform::loadd2smap(const char* filename)
{
	d2smap.clear();	
	assert(densemesh.mesh->n_vertices()>0);
	ifstream inputfile(filename);
	if (!inputfile)
	{
		inputfile.close();
		return false;
	}
	d2smap.resize(densemesh.mesh->n_vertices());	
	string strtmp;
	int vtick = 0;
	while(!inputfile.eof())
	{
		getline(inputfile,strtmp);
		vector<string>  splitstr;
		coutcmd::multistringsplit(strtmp,splitstr);					
		if (splitstr.size() == 2)
		{
			int orinid = atoi(splitstr[0].c_str());
			int destid = atoi(splitstr[1].c_str());
			assert(orinid>=0 && orinid < densemesh.mesh->n_vertices());
			assert(destid>=0 && destid < meshs[0]->n_vertices());
			d2smap[orinid] = destid;
			vtick++;
		}
	}
	assert(vtick == densemesh.mesh->n_vertices());
	inputfile.close();
	return true;
}


bool ARAPDeform::setFeatureVector(FeatureVector& simfv, FeatureVector& densefv, bool findconst)
{
	simfv.t.resize(simfv.s.size());
#ifdef DUSE_OPENMP
#pragma omp parallel
	{
#pragma omp for
#endif
		for (int i = 0; i < simfv.t.size(); i++)
		{
			simfv.t[i] = simfv.r[i]*simfv.s[i];
		}
		DOMP_END;
		assert(d2smap.size()==densemesh.mesh->n_vertices());
		densefv.t.resize(densemesh.mesh->n_vertices());
#ifdef DUSE_OPENMP
#pragma omp parallel
		{
#pragma omp for
#endif
			for (int i = 0; i < d2smap.size(); i++)
			{
				densefv.t[i] = simfv.t[d2smap[i]];
			}
			DOMP_END;
			if (findconst)
			{
				vector<double> simfvdis(simfv.t.size(),std::numeric_limits<double>::max());
				vector<double> densedis(densefv.t.size(),std::numeric_limits<double>::max());
				int densen = densefv.t.size();
				densefv.isConst.resize(densen,false);
				densefv.constPoint.resize(densen,Eigen::Vector3d::Zero());
				densefv.isNr.resize(densen,false);
				densefv.nrPoint.resize(densen,Eigen::Vector3d::Zero());
				for (int i = 0; i < densen; i++)
				{
					int simid = d2smap[i];
					if(simfv.isConst[simid])
					{
						OpenMesh::Vec3d simvec = ref.mesh->point(OpenMesh::VertexHandle(simid));
						//OpenMesh::Vec3d densevec = ref.mesh->point(OpenMesh::VertexHandle(i));
						OpenMesh::Vec3d densevec = densemesh.mesh->point(OpenMesh::VertexHandle(i));
						double dis = (simvec-densevec).length();
						densedis[i] = dis;
						if (dis<simfvdis[simid])
						{
							simfvdis[simid] = dis; 
						}
					}
				}
				for (int i = 0; i < densen; i++)
				{
					int simid = d2smap[i];
					if (simfv.isConst[simid] && (densedis[i] <= simfvdis[simid]))
					{
						densefv.isConst[i] = true;
						densefv.constPoint[i] = simfv.constPoint[simid];
					}
				}
			}
			return true;
		}


		void ARAPDeform::preSetMatrix(FeatureVector& fv)
		{
			std::cout << "Global Optimize ...... " << endl;
			long long tt = clock();
			vector< pair< pair<int,int>, double > > data;
			this->densevbconst.clear();
			this->densevbconst.resize(densemesh.mesh->n_vertices()*3,0);
			int n=0;
			double* weights = new double[densemesh.mesh->n_vertices()];
			int * debugtick = new int[densemesh.mesh->n_vertices()];
			for (int i=0; i<densemesh.mesh->n_vertices(); i++) {
				memset(weights,0,densemesh.mesh->n_vertices()*sizeof(double));
				memset(debugtick,0,densemesh.mesh->n_vertices()*sizeof(int));
				DTriMesh::VertexHandle vi(i);
				int ei = densemesh.d[i];
				if (fv.isConst[i]) 
				{
					for (int dim=0; dim<3; dim++)
					{
						data.push_back( make_pair( make_pair(i*3+dim, i*3+dim), 1 ));
					}
					continue;
				}

				double sumwij = 0, sumwji = 0;
				for (int k = 0; k < densemesh.rd[i].size(); k++)
				{
					int j = densemesh.rd[i][k].first;
					int ej = densemesh.rd[i][k].second;
					weights[j] = densemesh.getWij(ej);
					sumwji+=densemesh.getWij(ej);
					debugtick[j]++;
				}
				for (DTriMesh::VertexVertexIter vvi = densemesh.mesh->vv_iter(vi); vvi; vvi++, ei++) {
					int j = vvi.handle().idx();
					DTriMesh::VertexHandle vj(j);
					double wij = densemesh.getWij(ei);
					sumwij+=wij;
					weights[j] += densemesh.getWij(ei);
					debugtick[j]++;
					assert(debugtick[j]==2);
					if (!fv.isConst[j])
					{
						for (int dim = 0; dim < 3; dim++)
						{
							data.push_back(make_pair( make_pair(i*3+dim, j*3+dim), -weights[j]));
						}			
					}
					else
					{
						Eigen::Vector3d vec = weights[j]*fv.constPoint[j];
						for (int dim = 0; dim < 3; dim++)
						{
							densevbconst[i*3+dim]+=vec[dim];
						}
					}
				}
				if(!fv.isConst[i])
				{
					if(!fv.isNr[i])
					{
						for (int dim = 0; dim < 3; dim++)
						{
							data.push_back( make_pair( make_pair(i*3+dim, i*3+dim), sumwij+sumwji ));
						}
					}
					else
					{
						for(int dim = 0; dim < 3; dim++)
						{
							data.push_back(make_pair(make_pair(i*3+dim,i*3+dim),sumwij+sumwji+lambda));
							densevbconst[i*3+dim]+=(lambda*fv.nrPoint[i](dim));
						}
					}
				}
			}
			delete[] weights;
			weights = NULL;
			delete[] debugtick;
			debugtick = NULL;

			this->A.resize(fv.t.size() * 3, fv.t.size() * 3);
			std::vector<T> coefficients;
			coefficients.resize(data.size(), T(0, 0, 0));
			for (int i = 0; i < data.size(); i++) { //遍历行
				coefficients[i] = T(double(data[i].first.second), double(data[i].first.first), data[i].second);
			}
			this->A.setFromTriplets(coefficients.begin(), coefficients.end());
			this->Alu.compute(this->A);
			//this->A.compute

			//DMSpMatrix_gl A(*eng, "Adense", fv.t.size()*3, fv.t.size()*3, data, true);

			std::cout << "Save Matrix for pre-dec" << endl;
			std::cout << "Calc LU decompose ... ";
			//eng->Eval("[ldense,udense,pdense,qdense]=lu(Adense);");
			std::cout << "done" << endl;
			return ;
		}

		void ARAPDeform::solvefast(FeatureVector& fv, DTriMesh& mesh)
		{
			//std::cout << "Solving the dense mesh ......" <<endl;
			//long long tt = clock();
			vector<double> vb(densemesh.mesh->n_vertices()*3,0);
#ifdef DUSE_OPENMP
#pragma omp parallel
			{
#pragma omp for
#endif
				for (int i = 0; i < densemesh.mesh->n_vertices(); i++)
				{
					if (fv.isConst[i])
					{
						for (int dim = 0; dim < 3; dim++)
						{
							vb[i*3+dim] = fv.constPoint[i](dim);
						}
						continue;
					}
					DTriMesh::VertexHandle vi(i);
					int ei = densemesh.d[i];
					OpenMesh::Vec3d pi = densemesh.mesh->point(vi);
					for (DTriMesh::VertexVertexIter vvi = densemesh.mesh->vv_iter(vi);vvi;vvi++,ei++)
					{
						OpenMesh::Vec3d pj = densemesh.mesh->point(vvi);
						Eigen::Vector3d pij = OtoE(pi-pj);
						int j = vvi.handle().idx();
						double wij = densemesh.getWij(ei);
						assert(wij == wij);
						Eigen::Vector3d tmp = wij*fv.t[i]*pij;
						for (int dim = 0; dim < 3; dim++)
						{
							vb[i*3+dim]+=tmp(dim);
							vb[j*3+dim]-=tmp(dim);
						}
					}
					for (int dim = 0; dim < 3; dim++)
					{
						//vb[i*3+dim] += vbconst[i*3+dim];
						vb[i*3+dim] += densevbconst[i*3+dim];
					}
				}
				DOMP_END;
				int n = densemesh.mesh->n_vertices()*3;
				std::vector<double> result(densemesh.mesh->n_vertices()*3);
				//long long tt = clock();
				//tt = clock();
				//std::cout << "solve2 copy data to matlab" << endl;
				//DMMatrix b(*eng, "bbdense", n, 1, &vb[0]);
				Eigen::VectorXd x(n);
				double * ptr = &vb[0];
				Eigen::Map<Eigen::VectorXd> b(ptr, vb.size());
				//std::cout << "matlab compute" << endl;
				//eng->Eval("xdense=qdense*(udense\\(ldense\\(pdense*bbdense)));");
				this->Alu.solve(b);
				//DMMatrix x(*eng, "xdense", n, 1);
				result.resize(n);
				for (int i=0; i<n; i++)
					result[i] = x(i);

				for (int i=0; i<mesh.n_vertices(); i++) {
					DTriMesh::VertexHandle vi(i);

					if (fv.isConst[i]) {
						mesh.point(vi) = EtoO(fv.constPoint[i]);
						continue;
					}

					mesh.point(vi)[0] = result[i*3+0];
					mesh.point(vi)[1] = result[i*3+1];
					mesh.point(vi)[2] = result[i*3+2];
				}
				return;
			}

			void ARAPDeform::preprocess(const char* modelname)
			{
				string meshname(modelname);
				int npos = meshname.rfind('.');
				string prefix = meshname.substr(0,npos+1);
				string corrname = prefix+string("txt");
				//load obj
				//load cor
				this->loaddensemesh(modelname);
				this->loadd2smap(corrname.c_str());            
			}


void ARAPTest() {
    //DMEngine eng;
    DTriMesh mesh1, mesh2, mesh3, result;

    OpenMesh::IO::read_mesh(mesh1, "E:/project/SIGA2014/dataset/scape_smooth/49.obj");
    OpenMesh::IO::read_mesh(mesh2, "E:/project/SIGA2014/dataset/scape_smooth/50.obj");
    OpenMesh::IO::read_mesh(mesh3, "tmp.obj");

    vector<DTriMesh*> meshs;
    meshs.push_back(&mesh1);
    meshs.push_back(&mesh2);

    T_ARAPDeform deform(mesh1, meshs);

    FeatureVector fv = deform.fvs[0];
    fv.loadConstPoint(ifstream("E:/project/SIGA2014/dataset/scape_smooth/handles/c1.txt"), mesh3);

    vector<double> weight;
    weight.resize(meshs.size(),0);
    weight[0] = 1;

    result = mesh1;

    deform.solve2(fv, weight, result);

    OpenMesh::IO::write_mesh(result, "result.obj");
    system("pause");
}


void ARAPDeformFun(const char* foldername, int filenum, const char* constname, const char* outname)
{
	//DMEngine eng;
	vector<DTriMesh> meshes;	
	for (int i = 0; i < filenum; i++)
	{
		DTriMesh tmp;
        stringstream ss;
		ss<<(i+1);
		string filename = string(foldername)+"\\"+ss.str()+".obj";
		OpenMesh::IO::read_mesh(tmp,filename.c_str());
		meshes.push_back(tmp);
	}  
	vector<DTriMesh*> meshs;
	for (int i = 0; i < meshes.size(); i++)
	{
        meshs.push_back(&meshes[i]);        
	}
	T2_ARAPDeform deform(meshes[0], meshs);
	FeatureVector fv = deform.fvs[0];
	fv.loadConstPoint(ifstream(constname));
	vector<double> weight;
	weight.resize(meshs.size(),0);

	weight[0] = 0.5;
	weight[1] = 1-weight[0];
	//weight[1] = 0;

    DTriMesh result;
	result = meshes[0];
	deform.iterEps = 0.01;
	deform.maxIterTime = 5;
	deform.solve2(fv, weight, result);
	
	std::cout<<coutcmd::green<<"Weights: ";
	for (int i = 0; i < weight.size();i++)
	{
         std::cout<<weight[i]<<" ";
	}
	std::cout<<endl;
	OpenMesh::IO::write_mesh(result, outname);
	system("pause");
}


void BarDeform(const double w)
{
	const char* foldername = "E:\\SIGGRAPH\\SIGA2015\\bar\\shape";
	const int filenum = 2;
	const char* constname = "E:\\SIGGRAPH\\SIGA2015\\bar\\90.txt";
	//const char* outfolder = "E:\\SIGGRAPH\\SIGA2015\\bar\\result\\deform180";
	const char* outfolder = "E:\\SIGGRAPH\\SIGA2015\\bar\\result\\deform90";
	//DMEngine eng;
	vector<DTriMesh> meshes;	
	for (int i = 0; i < filenum; i++)
	{
		DTriMesh tmp;
		stringstream ss;
		ss<<(i+1);
		string filename = string(foldername)+"\\"+ss.str()+".obj";
		OpenMesh::IO::read_mesh(tmp,filename.c_str());
		meshes.push_back(tmp);
	}  
	vector<DTriMesh*> meshs;
	for (int i = 0; i < meshes.size(); i++)
	{
		meshs.push_back(&meshes[i]);        
	}
	T2_ARAPDeform deform(meshes[0], meshs);
	FeatureVector fv = deform.fvs[0];
	fv.loadConstPoint(ifstream(constname));
	vector<double> weight;
	weight.resize(meshs.size(),0);

	weight[0] = w;
	weight[1] = 1-w;

	DTriMesh result;
	result = meshes[0];
	deform.iterEps = 0.00000;
	deform.maxIterTime = 10000;
	deform.solve2(fv, weight, result);
	stringstream ss;
	ss<<w;
	string outname = string(outfolder)+ss.str()+string(".obj");
	std::cout<<coutcmd::green<<outname.c_str()<<coutcmd::white<<std::endl;
	OpenMesh::IO::write_mesh(result, outname);
}



// void ARAPDeformFun(const char* foldername, int filenum, const char* constname, const char* outname)
// {
// 	DMEngine eng;
// 	vector<DTriMesh> meshes;	
// 	for (int i = 0; i < filenum; i++)
// 	{
// 		DTriMesh tmp;
// 		stringstream ss;
// 		ss<<(i+1);
// 		string filename = string(foldername)+"\\"+ss.str()+".obj";
// 		OpenMesh::IO::read_mesh(tmp,filename.c_str());
// 		meshes.push_back(tmp);
// 	}  
// 	vector<DTriMesh*> meshs;
// 	for (int i = 0; i < meshes.size(); i++)
// 	{
// 		meshs.push_back(&meshes[i]);        
// 	}
// 	T2_ARAPDeform deform(eng, meshes[0], meshs);
// 	FeatureVector fv = deform.fvs[0];
// 	fv.loadConstPoint(ifstream(constname));
// 	vector<double> weight;
// 	weight.resize(meshs.size(),0);
// 	weight[0] = 1;
// 	DTriMesh result;
// 	result = meshes[0];
// 	deform.org->iterEps = 0.0001;
// 	deform.org->maxIterTime = 600;
// 	deform.solve2(fv, weight, result);
// 	OpenMesh::IO::write_mesh(result, outname);
// 	system("pause");
// }

void ARAPDeformHandle(const char* filename, const char* constname, const char* outname,const int iternum)
{
	//DMEngine eng;
	vector<DTriMesh> meshes;	
	DTriMesh tmp;
	OpenMesh::IO::read_mesh(tmp,filename);
	meshes.push_back(tmp);
//	for (int i = 0; i < filenum; i++)
//	{
//		DTriMesh tmp;
//		stringstream ss;
//		ss<<(i+1);
//		string filename = string(foldername)+"\\"+ss.str()+".obj";
//		OpenMesh::IO::read_mesh(tmp,filename.c_str());
//		meshes.push_back(tmp);
//	}  
	vector<DTriMesh*> meshs;
	for (int i = 0; i < meshes.size(); i++)
	{
		meshs.push_back(&meshes[i]);        
	}
	T2_ARAPDeform deform(meshes[0], meshs);
	deform.init();
	FeatureVector fv = deform.fvs[0];
	fv.loadConstPoint(ifstream(constname));
	vector<double> weight;
	weight.resize(meshs.size(),0);
	weight[0] = 1;
	DTriMesh result;
	result = meshes[0];
	//deform.org->iterEps = 0.0000;
	//deform.org->maxIterTime = iternum;
	deform.iterEps = 0; 
	deform.maxIterTime = iternum;
	deform.solve2(fv, weight, result);
	OpenMesh::IO::write_mesh(result, outname);
}

void ARAPDeformHandleRot(const char* filename, const char* constname, const char* rotname, const char* outname,const int iternum)
{
	//DMEngine eng;
	vector<DTriMesh> meshes;	
	DTriMesh tmp;
	OpenMesh::IO::read_mesh(tmp,filename);
	meshes.push_back(tmp);
	vector<DTriMesh*> meshs;
	for (int i = 0; i < meshes.size(); i++)
	{
		meshs.push_back(&meshes[i]);        
	}
	T2_ARAPDeform deform(meshes[0], meshs);
	deform.initwithhandle = true;
	FeatureVector fv = deform.fvs[0];
	//fv.loadConstPoint(ifstream(constname));
	fv.loadConstPointFixed(ifstream(constname));
	ifstream inputfile(rotname);
	if(!inputfile)
	{
		return ;
	}
	fv.loadHandleRT(ifstream(rotname));
	vector<double> weight;
	weight.resize(meshs.size(),0);
	weight[0] = 1;
	DTriMesh result;
	result = meshes[0];
	//deform.org->iterEps = 0.0000;
	//deform.org->maxIterTime = iternum;
	deform.iterEps = 0; 
	deform.maxIterTime = iternum;
	deform.org->initwithhandle = true;
	deform.solve2(fv, weight, result);
	OpenMesh::IO::write_mesh(result, outname);
}



void ARAPDeformInter(const char* srcname, const char* tarname, const int filenum, const char* outfolder)
{
	//DMEngine eng;
	vector<DTriMesh> meshes;	

	DTriMesh srctmp;
	DTriMesh tartmp;
//	stringstream ss;
//	ss<<(i+1);
//	string filename = string(foldername)+"\\"+ss.str()+".obj";
	OpenMesh::IO::read_mesh(srctmp,srcname);
	OpenMesh::IO::read_mesh(tartmp,tarname);
	meshes.push_back(srctmp);
	meshes.push_back(tartmp);

//	for (int i = 0; i < filenum; i++)
//	{
//		DTriMesh tmp;
//		stringstream ss;
//		ss<<(i+1);
//		string filename = string(foldername)+"\\"+ss.str()+".obj";
//		OpenMesh::IO::read_mesh(tmp,filename.c_str());
//		meshes.push_back(tmp);
//	}  
	vector<DTriMesh*> meshs;
	for (int i = 0; i < meshes.size(); i++)
	{
		meshs.push_back(&meshes[i]);        
	}

	//T2_ARAPDeform deform(eng, meshes[0], meshs);
	ARAPDeform deform(meshes[0],meshs);
	for (int i=0; i<filenum; i++) 
	{
	    vector<double> weight;
		double time = (double)i/(double)(filenum-1);
		weight.resize(meshs.size(),0);
		//time = 0.5;
		weight[0] = 1-time;
		weight[1] = time;
    	DTriMesh result;
		FeatureVector fv = deform.fvs[0];
		result = meshes[0];
		//deform.org->iterEps = 0.0001;
		//deform.org->maxIterTime = 600;
		//deform.solve2(fv, weight, result);
		deform.iterEps = 0.0001;
		deform.maxIterTime = 100;
		deform.solve(weight,result);
		stringstream ss;
		//ss<<outputname;
		ss << outfolder << "/" << i+1 << ".obj";
		std::cout << "solve " << ss.str() << endl;
		//lle.blendWeight(ids, idw, ss.str());
	    OpenMesh::IO::write_mesh(result, ss.str());
		//
//		deform.es.Save((string(outfolder)+string("/es.txt")).c_str());
//		deform.es.Save((string(outfolder)+string("/esbfs.txt")).c_str());
	}

	//FeatureVector fv = deform.fvs[0];
	//fv.loadConstPoint(ifstream(constname));
	//vector<double> weight;
	//weight.resize(meshs.size(),0);
	//weight[0] = 1;
	//DTriMesh result;
	//result = meshes[0];
	//deform.org->iterEps = 0.0001;
	//deform.org->maxIterTime = 600;
	//deform.solve2(fv, weight, result);
}



void ARAPDeformMultiScale(const char* foldername, int filenum, const char* constname, const char* densename, const char* outfolder)
{
	//DMEngine eng;
	vector<DTriMesh> meshes;	
	for (int i = 0; i < filenum; i++)
	{
		DTriMesh tmp;
		stringstream ss;
		ss<<(i+1);
		string filename = string(foldername)+"\\"+ss.str()+".obj";
		OpenMesh::IO::read_mesh(tmp,filename.c_str());
		meshes.push_back(tmp);
	}
	vector<DTriMesh*> meshs;
	for (int i = 0; i < meshes.size(); i++)
	{
		meshs.push_back(&meshes[i]);        
	}	
	T2_ARAPDeform deform(meshes[0], meshs);
	FeatureVector fv = deform.fvs[0];
	fv.loadConstPoint(ifstream(constname));
	vector<double> weight;
	weight.resize(meshs.size(),0);
	weight[0] = 1;
	//weight[1] = 1;
	DTriMesh result;
	DTriMesh denseresult;
	result = meshes[0];
	//load dense mesh and map
	deform.preprocess(densename);  
	deform.org->getsimfv = true;
	deform.solve2(fv, weight, result);
	deform.org->getsimfv = false;
	FeatureVector densefv;
	deform.setFeatureVector(deform.org->simfv,densefv,true);
	deform.preSetMatrix(densefv);
	denseresult = *deform.densemesh.mesh;
	long long t2 = clock();
	deform.solvefast(densefv,denseresult);
	std::cout<<"MultiScaleTime:"<<clock()-t2<<endl;
	system("pause");
	string simdeformname = string(outfolder)+"\\simdeform.obj";
	string densedeformname = string(outfolder)+"\\densedeform.obj";
	OpenMesh::IO::write_mesh(result, simdeformname.c_str());
	OpenMesh::IO::write_mesh(denseresult,densedeformname.c_str());
}
