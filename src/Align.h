#ifndef DALIGN_H
#define DALIGN_H

#include <Eigen/Eigen>
#include <vector>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#define PI 3.141592653589793
#define Eps 1e-10

extern int flag;

struct DTraits : public OpenMesh::DefaultTraits {
    typedef OpenMesh::Vec3d Point; // use double-values points
    typedef OpenMesh::Vec3d Normal; // use double-values points
    VertexAttributes( OpenMesh::Attributes::Normal );
    FaceAttributes( OpenMesh::Attributes::Normal );
};

typedef OpenMesh::TriMesh_ArrayKernelT<DTraits> DTriMesh;

class AffineAlign {
public :
    std::vector<Eigen::Vector3d> p;
    Eigen::Matrix3d AtA;
    AffineAlign(std::vector<Eigen::Vector3d> &v);
    Eigen::Matrix3d calc(const std::vector<Eigen::Vector3d> &v);
    double residual(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
	double residualwithoutnormal(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
    bool ckt(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
};

class RotateAlign {
public :
    std::vector<Eigen::Vector3d> *p, *q;
    double res;

    RotateAlign(std::vector<Eigen::Vector3d>& v1, std::vector<Eigen::Vector3d>& v2);
    Eigen::Matrix3d calc();
	Eigen::Matrix3d svd1();
    double residual(const Eigen::Matrix3d &m);
    bool ckt(Eigen::Matrix3d m);

    static void AlignAtoB(DTriMesh &A, DTriMesh &B);
    static void AlignAtoBCenter(DTriMesh &A, DTriMesh &B);
};

double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r, Eigen::Matrix3d &s);
double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r);

OpenMesh::Vec3d EtoO(const Eigen::Vector3d &v);
Eigen::Vector3d OtoE(const OpenMesh::Vec3d &v);

class Rot
{
public:
	Eigen::Vector3d axis;
	double theta;
	double circlek;
    Eigen::Matrix3d logr;
	Eigen::Matrix3d r;
	Rot() {circlek = theta = 0; axis=Eigen::Vector3d::Zero(); axis(0)=1;}
	void ToLogR();
	double ToAngle();
};

Eigen::Matrix3d exp(Eigen::Matrix3d);
Eigen::Matrix3d log(Eigen::Matrix3d);

Rot logrot(Eigen::Matrix3d, Rot);
Rot logrot(Eigen::Matrix3d);
#endif
