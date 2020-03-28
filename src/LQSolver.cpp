#include "LQSolver.h"
#include <Eigen/Eigen>
#include <fstream>
#include <iomanip>
#include <iostream>


using namespace std;

//LQSolver::LQSolver() {
//	//this->eng = NULL;
//	this->needrs = false;
//	saveA = false;
//}

//LQSolver::LQSolver(DMEngine *eng, bool needrs) {
//    this->eng = eng;
//    this->needrs = needrs;
//    saveA = false;
//}
LQSolver::LQSolver(bool needrs) {
	//this->eng = eng;
	this->A.resize(1, 1);
	this->needrs = needrs;
	saveA = false;
}


//void LQSolver::SetEngine(DMEngine *eng, bool needrs /* = false */)
//{
//	this->eng = eng;
//	this->needrs = needrs;
//	saveA = false;
//}


//double LQSolver::
//solve(int n, int m,
//        std::vector< std::pair< std::pair<int,int>, double> > &matData,
//        std::vector<double> &bData,
//        std::vector<double> &result, bool fixedOne) {
//
//    //std::cout << "copy data to matlab" << endl;
//    DMSpMatrix A(*eng, "A", n, m, matData, true);
//    DMMatrix b(*eng, "b", n, 1, &bData[0]);
//    //std::cout << "matlab compute" << endl;
//    eng->Eval("ab=A'*b;");
//    eng->Eval("AtA=A'*A;");
//    if (saveA) {
//        std::cout << "Save Matrix for pre-dec" << endl;
//        std::cout << "Calc LU decompose ... ";
//        eng->Eval("[l,u,p,q]=lu(AtA);");
//        std::cout << "done" << endl;
//        return 0;
//    }
//    if (fixedOne)
//        eng->Eval("x=LinearSolveSumFixed(AtA,ab,1)");
//    else
//        eng->Eval("x=AtA\\ab;");
//    //std::cout << "get x " << endl;
//    DMMatrix x(*eng, "x", m, 1);
//    result.resize(m);
//    for (int i=0; i<m; i++)
//        result[i] = x.GetNum(i,0);
//    if (needrs) {
//        eng->Eval("rs=norm(A*x-b);");
//        DMMatrix rs(*eng, "rs", 1, 1);
//        return rs.GetNum(0,0);
//    } else return 0;
//}

double LQSolver::
solve(int n, int m,
	std::vector< std::pair< std::pair<int, int>, double> > &matData,
	std::vector<double> &bData,
	std::vector<double> &result, bool fixedOne) {

	//std::cout << "copy data to matlab" << endl;
	//DMSpMatrix A(*eng, "A", n, m, matData, true);
	//DMMatrix b(*eng, "b", n, 1, &bData[0]);

	//SpMat A(n,m);
	this->A.resize(n, m);
	std::vector<T> coefficients;
	coefficients.resize(matData.size(), T(0, 0, 0));
	for (int i = 0; i < matData.size(); i++) { //遍历行
		coefficients[i] = T(double(matData[i].first.second), double(matData[i].first.first), matData[i].second);
	}
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	// fill A
	Eigen::VectorXd x(n);
	double * ptr = &bData[0];
	Eigen::Map<Eigen::VectorXd> b(ptr, bData.size());
	Eigen::VectorXd Ab = A.transpose()*b;
	Eigen::SparseMatrix<double> AtA = A.transpose()*A;

	//Eigen::SimplicialCholesky<SpMat> chol(A);
	// fill b
	//std::cout << "matlab compute" << endl;
	if (saveA) {
		std::cout << "Save Matrix for pre-dec" << endl;
		std::cout << "Calc LU decompose ... ";
		this->Achol.compute(AtA);
		//eng->Eval("[l,u,p,q]=lu(AtA);");
		std::cout << "done" << endl;
		return 0;
	}

	//eng->Eval("x=AtA\\ab;");
	x = this->Achol.solve(Ab);
	//std::cout << "get x " << endl;
	//DMMatrix x(*eng, "x", m, 1);
	result.resize(m);
	for (int i = 0; i < m; i++)
		result[i] = x(i);
	if (needrs) {
		//eng->Eval("rs=norm(A*x-b);");
		Eigen::VectorXd rs = A * x - b;
		//DMMatrix rs(*eng, "rs", 1, 1);
		//return rs.GetNum(0, 0);
		return rs.norm();
	}
	else return 0;
}



//double LQSolver::
//glsolve(int n, int m,
//        std::vector< std::pair< std::pair<int,int>, double> > &matData,
//        std::vector<double> &bData,
//        std::vector<double> &result, bool fixedOne) 
//{
//    //std::cout << "copy data to matlab" << endl;
//    //DMSpMatrix A(*eng, "A", n, m, matData, true);
//	DMSpMatrix_gl A(*eng, "A", n, m, matData, true);
//    DMMatrix b(*eng, "b", n, 1, &bData[0]);
//    //std::cout << "matlab compute" << endl;
//    //eng->Eval("ab=A'*b;");
//    //eng->Eval("AtA=A'*A;");
//    if (saveA) {
//        std::cout << "Save Matrix for pre-dec" << endl;
//        std::cout << "Calc LU decompose ... ";
//        eng->Eval("[l,u,p,q]=lu(A);");
//        std::cout << "done" << endl;
//        return 0;
//    }
//    if (fixedOne)
//        eng->Eval("x=LinearSolveSumFixed(A,b,1)");
//    else
//        eng->Eval("x=A\\b;");
//    //std::cout << "get x " << endl;
//    DMMatrix x(*eng, "x", m, 1);
//    result.resize(m);
//    for (int i=0; i<m; i++)
//        result[i] = x.GetNum(i,0);
//	return 0;
//    //if (needrs) {
//    //   eng->Eval("rs=norm(A*x-b);");
//    //   DMMatrix rs(*eng, "rs", 1, 1);
//    //   return rs.GetNum(0,0);
//    //} 
//    //else return 0;
//}


double LQSolver::
glsolve(int n, int m,
	std::vector< std::pair< std::pair<int, int>, double> > &matData,
	std::vector<double> &bData,
	std::vector<double> &result, bool fixedOne) {

	//std::cout << "copy data to matlab" << endl;
	//DMSpMatrix A(*eng, "A", n, m, matData, true);
	//DMMatrix b(*eng, "b", n, 1, &bData[0]);

	//SpMat A(n, m);
	this->A.resize(n, m);
	std::vector<T> coefficients;
	coefficients.resize(matData.size(), T(0, 0, 0));
	//int a1 = 0;
	//int a2 = 0;

	for (int i = 0; i < matData.size(); i++) { //遍历行
		coefficients[i] = T(double(matData[i].first.second), double(matData[i].first.first), matData[i].second);
		//if (matData[i].first.first > a1)
		//	a1 = matData[i].first.first;
		//if (matData[i].first.second > a2)
		//	a2 = matData[i].first.second;
	}
	this->A.setFromTriplets(coefficients.begin(), coefficients.end());
	// fill A
	Eigen::VectorXd x(n);
	double * ptr = &bData[0];
	Eigen::Map<Eigen::VectorXd> b(ptr, bData.size());
	//Eigen::VectorXd Ab = A*b;
	//Eigen::SparseMatrix<double> AtA = A.transpose()*A;
	std::cout << this->A.rows() << " " << this->A.cols() << std::endl;
	assert(this->A.rows() == this->A.cols());

	//Eigen::SimplicialCholesky<SpMat> chol(A);
	// fill b
	//std::cout << "matlab compute" << endl;
	if (saveA) {
		std::cout << "Save Matrix for pre-dec" << endl;
		std::cout << "Calc LU decompose ... ";
		//Eigen::MatrixXd denseA = Eigen::MatrixXd(this->A);
		this->Alu.compute(this->A);

		//this->Alu_nons.compute(denseA);
		//Eigen::MatrixXd B;
		//B.full
		//eng->Eval("[l,u,p,q]=lu(AtA);");
		std::cout << "done" << endl;
		return 0;
	}

	//eng->Eval("x=AtA\\ab;");
	x = this->Alu.solve(b);
	//x = this->Alu_nons.solve(b);
	//std::cout << "get x " << endl;
	//DMMatrix x(*eng, "x", m, 1);
	result.resize(m);
	for (int i = 0; i < m; i++)
		result[i] = x(i);
	return 0;
	//if (needrs) {
	//	//eng->Eval("rs=norm(A*x-b);");
	//	Eigen::VectorXd rs = A * x - b;
	//	//DMMatrix rs(*eng, "rs", 1, 1);
	//	//return rs.GetNum(0, 0);
	//	return rs.norm();
	//}
	//else return 0;
}



//void LQSolver::solve(int n,
//           std::vector<double> &bData,
//           std::vector<double> &result
//           ) {
//    long long tt = clock();
//    std::cout << "solve2 copy data to matlab" << endl;
//    DMMatrix b(*eng, "bb", n, 1, &bData[0]);
//    std::cout << "matlab compute" << endl;
//    eng->Eval("x=q*(u\\(l\\(p*bb)));");
//    DMMatrix x(*eng, "x", n, 1);
//    result.resize(n);
//    for (int i=0; i<n; i++)
//        result[i] = x.GetNum(i,0);
//
//    std::cout << "Linear solve done, time : " << (clock()-tt)*1./CLOCKS_PER_SEC << endl;
//}

void LQSolver::solve(int n,
	std::vector<double> &bData,
	std::vector<double> &result
) {
	long long tt = clock();
	std::cout << "solve2 copy data to matlab" << endl;
	//DMMatrix b(*eng, "bb", n, 1, &bData[0]);
	Eigen::VectorXd x(n);
	double * ptr = &bData[0];
	Eigen::Map<Eigen::VectorXd> b(ptr, bData.size());
	Eigen::VectorXd Ab = this->A.transpose()* b;
	//std::cout << "matlab compute" << endl;
	x = this->Achol.solve(Ab);
	//eng->Eval("x=q*(u\\(l\\(p*bb)));");
	//DMMatrix x(*eng, "x", n, 1);
	result.resize(n);
	for (int i = 0; i < n; i++)
		result[i] = x(i);

	std::cout << "Linear solve done, time : " << (clock() - tt)*1. / CLOCKS_PER_SEC << endl;
}

//
//void LQSolver::glsolve(int n,
//           std::vector<double> &bData,
//           std::vector<double> &result
//           ) {
//    long long tt = clock();
//    std::cout << "solve2 copy data to matlab" << endl;
//    DMMatrix b(*eng, "bb", n, 1, &bData[0]);
//    std::cout << "matlab compute" << endl;
//    eng->Eval("x=q*(u\\(l\\(p*bb)));");
//    DMMatrix x(*eng, "x", n, 1);
//    result.resize(n);
//    for (int i=0; i<n; i++)
//        result[i] = x.GetNum(i,0);
//
//    std::cout << "Linear solve done, time : " << (clock()-tt)*1./CLOCKS_PER_SEC << endl;
//}


void LQSolver::glsolve(int n,
	std::vector<double> &bData,
	std::vector<double> &result
) {
	long long tt = clock();
	std::cout << "solve2 copy data to matlab" << endl;
	//DMMatrix b(*eng, "bb", n, 1, &bData[0]);
	Eigen::VectorXd x(n);
	double * ptr = &bData[0];
	Eigen::Map<Eigen::VectorXd> b(ptr, bData.size());
	//Eigen::VectorXd Ab = this->A* b;
	std::cout << "matlab compute" << endl;
	//eng->Eval("x=q*(u\\(l\\(p*bb)));");
	x = this->Alu.solve(b);
	//DMMatrix x(*eng, "x", n, 1);
	result.resize(n);
	for (int i = 0; i < n; i++)
		result[i] = x(i);

	std::cout << "Linear solve done, time : " << (clock() - tt)*1. / CLOCKS_PER_SEC << endl;
}


void TriSecSolver::solve(double l, double r, double &best, double &number) {
     while (r-l>eps) {
         double ml = l + (r-l) * 0.4;
         double mr = l + (r-l) * 0.6;
         double sl = ConvexFunction(ml);
         double sr = ConvexFunction(mr);
         std::cout << "Q : " << l << " " << r << endl;
         std::cout << "M : " << ml << " " << mr << endl;
         std::cout << "S : " << sl << " " << sr << endl;
         if (sl < sr)
             r = mr, best = ml, number = sl;
         else
             l = ml, best = mr, number = sr;
     }
}

void addPoint(double x, double y, double k, Eigen::Matrix3d &AtA, Eigen::Vector3d &Atb) {
    double x2 = x*x, x3 = x2*x, x4 = x3*x;
    AtA *= k;
    Atb *= k;

    AtA(2,2) += 1;
    AtA(2,1) += x, AtA(1,2) += x;
    AtA(2,0) += x2, AtA(1,1) += x2, AtA(0,2) += x2;
    AtA(1,0) += x3, AtA(0,1) += x3;
    AtA(0,0) += x4;

    Atb += Eigen::Vector3d(y*x2, y*x, y);
}

double getMin(Eigen::Matrix3d &AtA, Eigen::Vector3d &Atb) {
    Eigen::Vector3d v = AtA.colPivHouseholderQr().solve(Atb);
    if (abs(v(0)) == 0)
        assert(1);
    return v(1) / (-2 * v(0));
}

void NewtonIterSolver::solve(double l, double r, double &best, double &number) {
    Eigen::Matrix3d AtA(Eigen::Matrix3d::Zero());
    Eigen::Vector3d Atb(Eigen::Vector3d::Zero());
    addPoint(l, ConvexFunction(l), 1, AtA, Atb);
    addPoint(r, ConvexFunction(r), 1, AtA, Atb);
    addPoint((l+r)/2, ConvexFunction((l+r)/2), 1, AtA, Atb);
    double pre = (l+r)/2, now = getMin(AtA, Atb);
    while (abs(pre-now) > eps) {
        std::cout << "Pre : " << pre << " Now : " << now << endl;
        addPoint(now, number = ConvexFunction(now), 1./10, AtA, Atb);
        pre = now, best = now = getMin(AtA, Atb);
    }
}

vector< pair<double,double> > findConvex(vector<pair<double,double>> v, pair<double,double> &best, pair<double,double> &best2, bool &isConvex) {
    sort(v.begin(), v.end(), [](pair<double,double> a, pair<double,double> b) -> bool {return a.second < b.second;});
    int fs=3;
    isConvex = true;
    if (v[0].first < v[1].first && v[1].first < v[2].first) {
        while (fs < v.size() && !(v[fs].first < v[0].first)) fs++;
        if (fs == v.size()) {
            std::cout << "Warning : function not convex" << endl;
            fs = 3;
            isConvex = false;
        }
    } else
    if (v[0].first > v[1].first && v[1].first > v[2].first) {
        while (fs < v.size() && !(v[fs].first > v[0].first)) fs++;
        if (fs == v.size()) {
            std::cout << "Warning : function not convex" << endl;
            fs = 3;
            isConvex = false;
        }
    }
    best = v[0], best2 = v[1];
    v[3] = v[fs];
    sort(v.begin(), v.begin() + 4);
    if (!(v[0].second > v[1].second && v[2].second < v[3].second)) {
        std::cout << "Warning : function not convex" << endl;
        isConvex = false;
    }
    v.erase(v.begin()+4, v.end());
    return v;
}

double getMin(vector< pair<double,double> > &cv) {
    double k0 = (cv[1].second - cv[0].second) / (cv[1].first - cv[0].first);
    if (cv[1].first == cv[0].first) k0 = -1;
    double k1 = (cv[3].second - cv[2].second) / (cv[3].first - cv[2].first);
    if (cv[3].first == cv[2].first) k0 = 1;
    double m0 = (cv[1].first + cv[0].first) * 0.5;
    double m1 = (cv[3].first + cv[2].first) * 0.5;
    if (k0 == k1) return (cv[1].first + cv[2].first) * 0.5;
    double kk = (m1 - m0) / (k1 - k0);
    return m0 - kk * k0;
}

void TriNewtonIterSolver::solve(double l, double r, double &best, double &number) {
    l = l, r = r;
    int calltime=0;
    Eigen::Matrix3d AtA(Eigen::Matrix3d::Zero());
    Eigen::Vector3d Atb(Eigen::Vector3d::Zero());
    vector< pair<double,double> > res;
    pair<double,double> tmp,bs;
    bool isConvex;

    tmp = make_pair(l, ConvexFunction(l));
    addPoint(tmp.first, tmp.second, 1, AtA, Atb);
    res.push_back(tmp);

    tmp = make_pair(r, ConvexFunction(r));
    addPoint(tmp.first, tmp.second, 1, AtA, Atb);
    res.push_back(tmp);

    tmp = make_pair((l+r)/2, ConvexFunction((l+r)/2));
    addPoint(tmp.first, tmp.second, 1, AtA, Atb);
    res.push_back(tmp);

    tmp.first = 0, tmp.second = this->zeroResult;
    res.push_back(tmp);

    double prego = 1e30, ts = 1;

    while (1) {
        vector< pair<double,double> > cv = findConvex(res, bs, tmp, isConvex);
        if (isConvex) {
            best = getMin(cv);
            if (prego == best) {
                std::cout << "# Error same value 1 inc ts :" << ts << endl;
                isConvex = false;
                ts *= 2;
            }
        }
        trick :
        if (!isConvex) {
            best = (bs.first * (ts+1) - tmp.first) / ts;
            std::cout << "# Fixed Convex function , best : ("<< bs.first << ", " << bs.second << ") ";
            std::cout << " second : (" << tmp.first << ", " << tmp.second << ")" << endl;
        }
        if (prego == best) {
            std::cout << "# Error same value 2 inc ts :" << ts << endl;
            ts *= 2;
            goto trick;
        }
        prego = best;
        number = bs.second;
        std::cout << "Pre : " << bs.first << " Now : " << best << " ds : " << abs(bs.first - best) << endl;
        if (abs(bs.first - best) < eps) return;
        number = ConvexFunction(best);
        std::cout << "Calculate " << best << " get rs :" << number << endl;
        res.push_back(make_pair(best, number));
    }
}

#define Eps 1e-10

double LQSolverWithLimit::iterEps = 1e-9;

double LQSolverWithLimit::newtonIterEps = 1e-3;
double LQSolverWithLimit::newtonMaxIterTime = 6;

void LQSolverWithLimit::solve(
        std::vector< Eigen::VectorXd > A,
        std::vector<double> &x,
        Eigen::VectorXd b,
        std::vector< std::pair<double, double> > limit
) {
    auto tt = clock();
    std::cout << "LQSolverWithLimit ..." << endl;
    vector< Eigen::VectorXd > blend;
    blend.resize(x.size(), Eigen::VectorXd::Zero(x.size()));
    for (int i=0; i<x.size(); i++)
        blend[i](i) = 1;
    vector<int> o(x.size(), 0);

    auto getLimit = [&](int i, double &l, double &r) {
        l = limit[i].first;
        r = limit[i].second;
        return;
        l = -1e300, r = -l;
        Eigen::VectorXd &bl = blend[i];
        for (int j=0; j<x.size(); j++) {
            if (abs(bl(j)) > Eps) {
                if (bl(j) < 0) {
                    l = max(l,limit[j].second / bl(j));
                    r = min(r,limit[j].first / bl(j));
                } else {
                    l = max(l,limit[j].first / bl(j));
                    r = min(r,limit[j].second / bl(j));
                }
            }
        }
        if (l > r) {
            std::cout << "!!Error l > r i : " << i <<
                    " l : " << l <<
                    " r : " << r << endl;
            for (int j=0; j<x.size(); j++) {
                if (abs(bl(j)) > Eps) {
                    std::cout << "j : " << j <<
                            " limit : " << limit[j].first << ", " << limit[j].second <<
                            " bj : " << bl(j) << endl;
                    l = max(l,limit[j].first / bl(j));
                    r = min(r,limit[j].second / bl(j));
                }
            }
        }
    };

    auto calcRs = [&](int i, double x) -> double {
        return (b - A[i] * x).squaredNorm();
    };

    auto getBest = [&](int i, double &rx, double &rs) {
        Eigen::VectorXd &a = A[i];
        rx = a.dot(b) / a.squaredNorm();
        double l,r;
        getLimit(i,l,r);
        // std::cout << "rg : " << l << ", " << r << " rx : " << rx << endl;
        if (rx < l) rx = l; else
            if (rx > r) rx = r;
        rs = calcRs(i,rx);
    };

    while (1) {
        double prers = calcRs(0,0);
        double bestrs = 1e300;
        double bestx;
        int besti = -1;
        for (int i=0; i<x.size(); i++) if (o[i] != 1) {
            if (A[i].norm() < Eps) {
                o[i] = 1;
                continue;
            }
            double brs, bx;
            getBest(i, bx, brs);

            //std::cout << "i : " << i <<
            //        " bx : " << bx <<
            //        " brs : " << brs << endl;

            if (bx < 0) brs = prers - (prers - brs) * 64;

            if (brs < bestrs) {
                bestrs = brs;
                bestx = bx;
                besti = i;
            }
        }
        //std::cout << "besti : " << besti <<
        //        " bestRs : " << bestrs <<
        //        " bestX : " << bestx << endl;
        if (besti == -1) break;
        if (abs(bestx) < iterEps) break;
        b -= bestx * A[besti];
        //std::cout << "realRs : " << b.squaredNorm() << endl;
        o[besti] = 2;
        for (int i=0; i<x.size(); i++) {
            double dx = bestx * blend[besti](i);
            x[i] += dx;
            limit[i].first -= dx;
            limit[i].second -= dx;
            if (o[besti] == 0)
            if (o[i] == 0 && i != besti) {
                double rx = A[i].dot(A[besti]) / A[besti].squaredNorm();
                A[i] -= A[besti] * rx;
                blend[i] -= blend[besti] * rx;
            }
        }
    }
    std::cout << "LQSolverWithLimit Done time : " << (clock()-tt)*1. / CLOCKS_PER_SEC << endl;
}

void LQSolverWithLimit::ckLQWL() {
    auto sand = time(0);
    //sand = 1400086731ll;
    srand(sand);
    int n = 200;
    int m = 160;
    vector< Eigen::VectorXd > A;
    vector< double > x;
    vector< pair<double,double> > limit;

    for (int i=0; i<m; i++) {
        A.push_back(Eigen::VectorXd::Random(n));
        x.push_back(0);
        double rd = rand() * 1. / RAND_MAX;
        limit.push_back( make_pair(0-rd,1-rd) );
    }

    Eigen::VectorXd b = Eigen::VectorXd::Random(n) * m;

    solve(A, x, b, limit);

    for (int i=0; i<x.size(); i++)
        std::cout << x[i] << endl;

    auto calcRs = [&](void) -> double {
        Eigen::VectorXd bb = b;
        for (int i=0; i<x.size(); i++) {
            if (x[i] < limit[i].first || x[i] > limit[i].second)
                return 1e30;
            bb -= x[i] * A[i];
        }
        return bb.squaredNorm();
    };

    double delta = 1e-3;
    //x[0] = 0.5;
    double prers = calcRs();

    std::cout << "X : " << endl;
    for (int i=0; i<x.size(); i++) std::cout << x[i] << endl;

    for (int i=0; i<x.size(); i++) for (int f=-1; f<=1; f+=2) {
        x[i] += delta * f;
        double nowrs = calcRs();
        if (nowrs < prers)
            std::cout << "!Error i : " << i <<
                    " f : " << f <<
                    " rs : " << nowrs << " , " << prers << endl;
        x[i] -= delta * f;
    }

    std::cout << "Rs : " << calcRs() << endl;

    std::cout << b.transpose() << endl;
    std::cout << "sand : " << sand << endl;
}


void LQSolverWithLimit::NewtonIter(
        GetJacobiFunction getJacobi,
        ReturnRsFunction returnRs) {
    int iterTime = 0;
    while (1) {
        std::cout << "NewtonIter " << iterTime+1 << endl;
        vector< Eigen::VectorXd > A;
        vector< double > x;
        vector< pair<double,double> > limit;
        Eigen::VectorXd b;

        getJacobi(A, b, limit);
        x.resize(A.size(), 0);
        solve(A, x, b, limit);
        returnRs(x);
        double rsx = 0;
        for (int i=0; i<x.size(); i++) rsx += x[i]*x[i];
        rsx = sqrt(rsx);
        if (rsx < newtonIterEps) return;
        if (iterTime++ >= newtonMaxIterTime) return;
    }
}
