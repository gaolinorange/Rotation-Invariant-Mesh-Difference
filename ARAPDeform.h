#ifndef ARAP_DEFORM_H
#define ARAP_DEFORM_H


//#include "MatEngine.h"
//#include "DMEngine.h"
#include "FeatureVector.h"
#include "LQSolver.h"
#include "tools.h"
#include "coutcolorarap.h"
//#include "yj_deformpara.h"
//#include "mymatrix.h"

#define tic long long timetick = clock();
#define tocv clock()-timetick
#define tocp std::cout<<coutcmd::blue<<tocv<<coutcmd::white<<std::endl;

class ARAPDeform {
public :
    //DMEngine *eng;
    RefMesh ref;
    std::vector<DTriMesh*> meshs;
    std::vector<FeatureVector> fvs;
	std::vector<double> vbconst;
	std::vector<Eigen::Matrix3d> rdrs;
    LQSolver lqsolver;
    bool needAlign;
    bool needCalcRs;
    int maxIterTime;
    int newtonIterTime;
    double iterEps;
    double newtonIterEps;
    std::pair<double,double> rlimit;
    double SRp;
    bool CalcLocalFirst;

    bool AtAChanged;
	//new parameter
	double lambda;
	double lambdaplane;

	EnergyStatistic es;

	bool initwithhandle;
	std::vector<std::pair<int,int>> _eidmap;

	ARAPDeform() {} 

    //ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms);
	ARAPDeform(DTriMesh &refMesh, std::vector<DTriMesh*> ms);

	void init();

	//void SetEngine(DMEngine* _eng);

	void GetFeature(DTriMesh &refMesh, std::vector<DTriMesh*> ms);

	void GetFeatureS(DTriMesh &refMesh, std::vector<DTriMesh*> ms);

    // given feature vector, solve mesh
    void solve(FeatureVector fv, DTriMesh &mesh, bool predec = true);
	//
	void gl_solve(FeatureVector fv, DTriMesh &mesh, bool predec = true);

    // given feature vector weights, solve mesh
    void solve(std::vector<double> weight, DTriMesh &mesh);

    // given some const point, try to blend fvs to  it
    void solve(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
    virtual void solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
	void solve3(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);

    // given feature vector,  try to blend fvs to  it
    void solve(FeatureVector &fv, std::vector<double> &weight);

    // given mesh, try to blend fvs to it
    void solve(DTriMesh &mesh, std::vector<double> &weight);
	void leastsquare(FeatureVector &fv, std::vector<double> & weight, DTriMesh &mesh, bool sumone);		
    void preDecMatrix(FeatureVector &fv, DTriMesh &mesh);
    virtual double globalOptimize(FeatureVector &fv, DTriMesh &mesh);
    double globalOptimizeFast(FeatureVector &fv, DTriMesh &mesh);
	virtual double gl_globalOptimize(FeatureVector &fv, DTriMesh &mesh);
	void assignvb(FeatureVector &fv, DTriMesh &mesh, std::vector<double>& vb);
	double getRS(FeatureVector &fv, DTriMesh &mesh);
	double gl_globalOptimizeFast(FeatureVector &fv, DTriMesh &mesh);
    virtual double localOptimize(FeatureVector &fv, DTriMesh &mesh);
    void localOptimizeFast(FeatureVector &fv, DTriMesh &mesh);
    virtual double getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh);
    std::vector<double> getGradient(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult);
	//nr function
	virtual void solve_nr(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
	void preDecMatrix_nr(FeatureVector &fv, DTriMesh &mesh);
	virtual double globalOptimize_nr(FeatureVector &fv, DTriMesh &mesh);
	double globalOptimizeFast_nr(FeatureVector &fv, DTriMesh &mesh);
	virtual double getWeightResidual_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh);
    std::vector<double> getGradient_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult);

	virtual void gl_solve_nr(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
	void gl_preDecMatrix_nr(FeatureVector &fv, DTriMesh &mesh);
	virtual double gl_globalOptimize_nr(FeatureVector &fv, DTriMesh &mesh);
	double gl_globalOptimizeFast_nr(FeatureVector &fv, DTriMesh &mesh);
	virtual double gl_getWeightResidual_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh);
	std::vector<double> gl_getGradient_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult);
    double getRS_nr(FeatureVector &fv, DTriMesh& mesh);


    void dfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh);
    void bfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh);
    void bfsInit2(FeatureVector &fv, DTriMesh &mesh);
	void bfscorrot(FeatureVector &fv);
    void writeIterMesh(DTriMesh &mesh, std::string name, int id);
    void ckAlign(DTriMesh &mesh);
    void initWeight(std::vector<double> &weight);
    void initMesh(FeatureVector &fv, DTriMesh &mesh);
	void initMeshRotationHandle(FeatureVector &fv, DTriMesh &mesh);
	//multi scale operation
	bool getsimfv;
	FeatureVector simfv;

	std::vector<int> d2smap;
	bool loadd2smap(const char* filename);
	//DTriMesh densemesh;
	RefMesh densemesh;
	std::vector<double> densevbconst;
	bool setFeatureVector(FeatureVector& simfv, FeatureVector& densefv, bool findconst = false);
	bool loaddensemesh(const char* filename);
    void preSetMatrix(FeatureVector& fv);
	void solvefast(FeatureVector&fv, DTriMesh& mesh);
	//
	void preprocess(const char* modelname);
	void reconstruct(FeatureVector& simfv, DTriMesh& densemesh);
	//void solvemesh();
	SpMat A;
	Eigen::SparseLU<SpMat> Alu;

	//yangjie add

	//bool deformationEnable;	//energy item enable used by gaussion newton method
	//bool dataDrivenEnable;
	//bool itershape; //whether output intermediate result
	//bool weightNormalizationEnable;//feature vector weight normization
	//bool globalRotationEnable;
	//bool deformationdeltar;
	//bool cvx;
	//bool deltaw;
	//bool CUDAEnable;//GPU speed up
	//bool controalPointEnable;//whether use control point
	//int CUDAMaxIterNumb;//CUDA solve Ax=B's max iteration time
	//double lambdaDeformation;//deformation constrain item coefficient
	//double lambdaGlobalRotation;//ICP global rotation constrain item coefficient
	//double scalelambdaDeformation;
	//double scalelambdaGlobalRotation;
	//std::string interfolder;
	//Utility::MatEngine matEngine;
	//Utility::MySparseMatrix * matrixAMatlab;

	////chol factor matrix
	//int * csrRowIndGTPtr;
	//int * csrColGTPtr;
	//float * csrValGTPtr;
	//int nnzGT;
	//int csrRowNumbGT;
	//int csrColNumbGT;

	//int * csrRowIndGPtr_chol;
	//int * csrColGPtr_chol;
	//float * csrValGPtr_chol;
	//int nnzG_chol;
	//int csrRowNumbG_chol;
	//int csrColNumbG_chol;
	//float *resultXarap;

	//void yj_meshik_pre(DTriMesh &trackingMesh, FeatureVector & trackingFeatureVector);
	//void yj_leastsquarefv(FeatureVector &fv, std::vector<double> & weight,bool sumone = true);
	//void yj_set_meshik_deformpara(const DeformOption dopt);
	//void yj_meshik_deform(const std::vector<FeatureVector> &orthogonalBasis, std::vector<double> &featureVectorWeight, DTriMesh &trackingMesh, FeatureVector & trackingFeatureVector);


};

class SR_ARAPDeform : public ARAPDeform {
public:
    //SR_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms)
	SR_ARAPDeform(DTriMesh &refMesh, std::vector<DTriMesh*> ms)
        //: ARAPDeform(eng, refMesh, ms) {}
		: ARAPDeform(refMesh, ms) {}
    double globalOptimize(FeatureVector &fv, DTriMesh &mesh);
    double localOptimize(FeatureVector &fv, DTriMesh &mesh);
};

class T_ARAPDeform : public ARAPDeform {
public:
    ARAPDeform *org;
    //DMEngine eng2;
    std::vector<bool> ev;

	T_ARAPDeform() {}
    //T_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms);
	T_ARAPDeform(DTriMesh &refMesh, std::vector<DTriMesh*> ms);
    ~T_ARAPDeform() {delete org;}
    double getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh);
	double getWeightResidual_nr(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh);
    void solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
	void solve_nr(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);


};

class T2_ARAPDeform : public ARAPDeform {
public:
    ARAPDeform *org;
    //DMEngine eng2;
    std::vector<bool> ev;


    //T2_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms);
	T2_ARAPDeform(DTriMesh &refMesh, std::vector<DTriMesh*> ms);
    ~T2_ARAPDeform() {delete org;}
    void solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
};


void ARAPTest();
void ARAPDeformFun(const char* foldername, int filenum, const char* constname, const char* outname);
void ARAPDeformMultiScale(const char* foldername, int filenum, const char* constname, const char* densename, const char* outfolder);
void ARAPDeformHandle(const char* filename, const char* constname, const char* outname, const int iternum = 600);
void ARAPDeformHandleRot(const char* filename, const char* constname, const char* rotname, const char* outname,const int iternum = 600);
void ARAPDeformInter(const char* srcname, const char* tarname, const int num, const char* outfolder);
void BarDeform(const double w);


#endif
