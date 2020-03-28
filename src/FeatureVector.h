#ifndef DFEATURE_VECTOR_H
#define DFEATURE_VECTOR_H

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <vector>
#include <Eigen/Eigen>

#include <fstream>
#include "Align.h"
//#include "DMEngine.h"

class FeatureVector;

class LaplaceMatrix
{
public:
	
};

class EdgeLaplace
{
public:
	int varnum;
	std::vector<int> edgeflipid;
	std::vector<int> edgevarid;
	std::vector<int> faceid;
	EdgeLaplace() {varnum = 0;}
};



class RefMesh {
public :
	int vvs;
	DTriMesh *mesh;
	std::vector<AffineAlign*> align;
	std::vector<int> d;
	std::vector< std::vector< std::pair<int,int> > > rd;
	std::vector<double> w,c;
	std::vector< std::pair<int, std::pair<int,int> > > bfsq;
	std::vector<Eigen::Vector3d> edgeLength;
	std::vector<int> degree;
	bool usec;
	int fvNum;
	~RefMesh();

	RefMesh();
	RefMesh(DTriMesh &ms);
	void SetMesh(DTriMesh &ms);

	void getFeature(DTriMesh &ms, FeatureVector &fv);
	void getFeatureS(DTriMesh &ms, FeatureVector &fv);
	void getEidmap(DTriMesh &ms, std::vector<std::pair<int,int>>& _eidmap);
	void getEdgeLaplace(DTriMesh &ms, EdgeLaplace& el);

	void printCij(std::ofstream &fout);
	double getWij(int e);
	bool checksymm();

	static double normalScale;
	static int root;
};

class FeatureVector {
public :
	std::vector<Eigen::Matrix3d> s,r,dr,logdr,logr,t;
	std::vector<Rot> rots;
	Eigen::Matrix3d rhandle;
	//Eigen::Matrix3d thandle;
	Eigen::Vector3d thandle;
	std::vector<bool> isConst;
	std::vector<int> idfixed;
	std::vector<int> iddeformed;
	std::vector<Eigen::Vector3d> constPoint;
	std::vector<bool> isNr;
	std::vector<bool> isNrPlane;
	std::vector<Eigen::Vector3d> nrPoint;
	std::vector<Eigen::Vector3d> planePoint; 
	FeatureVector() {}
	FeatureVector(std::vector<double> weight, std::vector<FeatureVector> &fvs);
	void setConstPoint(int i, Eigen::Vector3d v);
	void setNrPoint(int i, Eigen::Vector3d v);
	void setPlanePoint(int i, Eigen::Vector3d v);
	void loadConstPoint(std::istream& cin);
	void loadConstPointFixed(std::istream& cin);
	void loadsoftConstPoint(std::istream& cin);
	void loadConstPoint(std::istream& cin, DTriMesh& mesh);
	void loadHandleRT(std::istream& cin);
	void loadConstPoint(std::vector<int>& ids, std::vector<OpenMesh::Vec3d>& v3ds);
	void loadNrPoint(std::vector<int>& ids, std::vector<OpenMesh::Vec3d>& v3ds);
	void loadNrPlanePoint(std::vector<int>& ids, std::vector<OpenMesh::Vec3d>& v3ds);
	void blendFrom(std::vector<double> weight, std::vector<FeatureVector> &fvs);
	void calcLogRFromR();
	void calcRFromLogR();
	void resize(const FeatureVector &other);
	void IdentityRotation();

	//yangjie add
	FeatureVector(DTriMesh & mesh);
	FeatureVector(std::vector<double> weight, std::vector<FeatureVector> &fvs, std::string string);
	void yj_blendFrom(std::vector<double> weight, std::vector<FeatureVector> &fvs, std::string string) ;

};

std::ostream& operator<<(std::ostream& cout, FeatureVector &fv);

FeatureVector operator +(const FeatureVector &a, const FeatureVector &b);
FeatureVector operator -(const FeatureVector &a, const FeatureVector &b);
double operator *(const FeatureVector &a, const FeatureVector &b);
FeatureVector operator *(const FeatureVector &a, const double &b);
#endif
