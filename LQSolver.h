#ifndef DLQ_SOLVER_H
#define DLQ_SOLVER_H

//#include "DMEngine.h"
#include <vector>
#include <functional>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/SparseCholesky>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

typedef Eigen::SparseMatrix<double> SpMat; // 声明一个列优先的双精度稀疏矩阵类型
typedef Eigen::Triplet<double> T; //三元组（行，列，值）

class LQSolver {
public :
    //DMEngine *eng;
    bool needrs;
    bool saveA;
	SpMat A;
	Eigen::SimplicialCholesky<SpMat> Achol;
	Eigen::SparseLU<SpMat> Alu;
	Eigen::FullPivLU<Eigen::MatrixXd> Alu_nons;

	//LQSolver();

    //LQSolver(DMEngine *eng, bool needrs = false);
	LQSolver(bool needrs = false);

	//void SetEngine(DMEngine *eng, bool needrs = false);


    double solve(int n, int m,
               std::vector< std::pair< std::pair<int,int>, double> > &matData,
               std::vector<double> &bData,
               std::vector<double> &result,
               bool fixedOne = false);

    void solve(int n,
               std::vector<double> &bData,
               std::vector<double> &result
               );

	double glsolve(int n, int m,
               std::vector< std::pair< std::pair<int,int>, double> > &matData,
               std::vector<double> &bData,
               std::vector<double> &result,
               bool fixedOne = false);

    void glsolve(int n,
               std::vector<double> &bData,
               std::vector<double> &result
               );

};



class ConvexSolver {
public:
    double zeroResult;
    double eps;
    int iterTime;
    std::function<double(double)> convexFunction;
    virtual void solve(double l, double r, double &best, double &number) {};
    ConvexSolver() : iterTime(0) {}
    double ConvexFunction(double x) {iterTime++; return convexFunction(x);}
};

class TriSecSolver : public ConvexSolver {
public:
    void solve(double l, double r, double &best, double &number);
};

class NewtonIterSolver : public ConvexSolver {
public:
    void solve(double l, double r, double &best, double &number);
};

class TriNewtonIterSolver : public ConvexSolver {
public:
    void solve(double l, double r, double &best, double &number);
};

class LQSolverWithLimit {
public:
    static double iterEps;
    static double newtonIterEps;
    static double newtonMaxIterTime;
    static void solve(std::vector<Eigen::VectorXd> A,
            std::vector<double> &x,
            Eigen::VectorXd b,
            std::vector< std::pair<double, double> > limit);

    static void ckLQWL();
    typedef std::function<void(
            std::vector<Eigen::VectorXd> &,
            Eigen::VectorXd &,
            std::vector< std::pair<double, double> > &
    )> GetJacobiFunction;

    typedef std::function<void(
            std::vector<double> &
    )> ReturnRsFunction;
    static void NewtonIter(
            GetJacobiFunction getJacobi,
            ReturnRsFunction returnRs);
};

#endif
