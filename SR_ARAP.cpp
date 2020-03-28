#include "ARAPDeform.h"

using namespace std;
using namespace Eigen;


double SR_ARAPDeform::globalOptimize(FeatureVector &fv, DTriMesh &mesh) {
    if (lqsolver.saveA) return 0;
    cout << "SR_ARAP Global Optimize ...... " << endl;
    long long tt = clock();
    vector< pair< pair<int,int>, double > > data;
    vector<double> vb;
    int n=0;
    for (int j=0; j<mesh.n_vertices(); j++) {
            DTriMesh::VertexHandle vj(j);
            Eigen::Matrix3d &rj = fv.r[j];
            Eigen::Matrix3d &sj = fv.s[j];
            DTriMesh::Point qj = ref.mesh->point(vj);
            int ej = ref.d[j];
            for (DTriMesh::VertexVertexIter vvj = mesh.vv_iter(vj); vvj; vvj++, ej++) {
                int k = vvj.handle().idx();
                double wjk = ref.w[ej];
                DTriMesh::Point tqjk = ref.mesh->point(DTriMesh::VertexHandle(k)) - qj;
                Eigen::Vector3d qjk(tqjk[0],tqjk[1],tqjk[2]);
                Eigen::Vector3d c = rj*(sj*qjk);
                if (fv.isConst[k]) c -= fv.constPoint[k];
                if (fv.isConst[j]) c += fv.constPoint[j];
                for (int dim=0; dim<3; dim++) {
                    if (!fv.isConst[k])
                        data.push_back(make_pair( make_pair(n, k*3+dim), wjk ));
                    if (!fv.isConst[j])
                        data.push_back(make_pair( make_pair(n, j*3+dim), - wjk ));
                    vb.push_back( wjk * c(dim) );
                    n++;
                }
            }
    }
    vector<double> result;
    cout << "matlab solve" << endl;
    lqsolver.needrs = true;
    double rs = lqsolver.solve(n, fv.s.size()*3, data, vb, result);
    lqsolver.needrs = false;
    cout << "copy ans" << endl;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        for (int dim=0; dim<3; dim++)
            mesh.point(vi)[dim] = result[i*3+dim];
    }

    cout << "SR_ARAP Global Optimize Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << rs*rs << endl;
    return rs*rs;
}

double SR_ARAPDeform::localOptimize(FeatureVector &fv, DTriMesh &mesh) {
    cout << "SR_ARAP Local Optimize ...... " << endl;
    long long tt = clock();
    for (int j=0; j<mesh.n_vertices(); j++) {
        vector<Eigen::Vector3d> vq2;
        Matrix3d mat1(Matrix3d::Zero()), mat2(Matrix3d::Zero());

            DTriMesh::VertexHandle vj(j);
            Eigen::Matrix3d &rj = fv.r[j];
            Eigen::Matrix3d &sj = fv.s[j];
            DTriMesh::Point qj = ref.mesh->point(vj);
            DTriMesh::Point pj = mesh.point(vj);
            int ej = ref.d[j];
            for (DTriMesh::VertexVertexIter vvj = mesh.vv_iter(vj); vvj; vvj++, ej++) {
                Eigen::Matrix3d &rjk = fv.dr[ej];
                int k = vvj.handle().idx();
                double wjk = ref.w[ej];
                DTriMesh::Point tqjk = ref.mesh->point(DTriMesh::VertexHandle(k)) - qj;
                DTriMesh::Point tpjk = mesh.point(DTriMesh::VertexHandle(k)) - pj;
                Eigen::Vector3d qjk(tqjk[0],tqjk[1],tqjk[2]);
                Eigen::Vector3d pjk(tpjk[0],tpjk[1],tpjk[2]);
                vq2.push_back(qjk);
                qjk = (rjk*(sj*qjk));
                qjk *= wjk;
                pjk *= wjk;
                mat1 += qjk * pjk.transpose();
                mat2 += (fv.r[k] * rjk.transpose()).transpose() * wjk;
            }
            vq2.push_back(vq2[0]);
            double area = 0;
            for (int i=0; i+1<vq2.size(); i++)
                area += vq2[i].cross(vq2[i+1]).norm();

        mat1 += SRp * area * mat2;
        polarDec(mat1, rj);
        rj.transposeInPlace();
    }
    cout << "SR_ARAP Local Optimize Done Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << endl;
    return 0;
}
