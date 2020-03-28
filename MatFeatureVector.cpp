#include "MatFeatureVector.h"
#include "ARAPDeform.h"
#include <fstream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

#ifdef DUSE_OPENMP
#define DOMP_END \
}
#else
#define DOMP_END ;
#endif

template<typename T>
T load_txt1(std::string path)
{
	std::ifstream in;
	in.open(path);
	std::string line;
	std::vector<double> values;
	UINT rows = 0;
	while (std::getline(in, line))
	{
		std::stringstream ss(line);
		std::string cell;
		while (std::getline(ss, cell, ','))
		{
			double val = std::stod(cell);
			values.push_back(val);
		}
		++rows;
	}
	return Eigen::Map<const Eigen::Matrix<
		typename T::Scalar,
		T::RowsAtCompileTime,
		T::ColsAtCompileTime,
		Eigen::RowMajor>>(values.data(), rows, values.size() / rows);
}


namespace MatFV
{
	  bool MatFeatureVector::LoadMesh(const char* filelist)
	  {
		  ifstream inputfile(filelist);
		  if (!filelist)
		  {
			cout << "Load File List Error" << endl;
			  inputfile.close();
			  return false;
		  }
		  vector<string> vs;
		int num = 0;
		  //inputfile>>num;
		  string tmp;
		getline(inputfile, tmp);
		  num = atoi(tmp.c_str());
		  for (int i = 0; i < num; i++)
		  {
			getline(inputfile, tmp);
				vs.push_back(tmp);
		  }
		  inputfile.close();
		  this->LoadMesh(vs);
		  return true;
	  }

	  bool MatFeatureVector::LoadMesh(const char* foldername, int num)
	  {
		  string fn(foldername);
		  vector<string> vs;
		  for (int i = 1; i <= num; i++)
		  {
			    char buffer[512];
			memset(buffer, 0, 512 * sizeof(char));
			sprintf_s(buffer, "%d.obj", i);
				string iname(buffer);
			vs.push_back(fn + "\\" + iname);
		  }
		  this->LoadMesh(vs);
          return true;
	  }


	  bool MatFeatureVector::LoadMeshWithR(const char* foldername, int num)
	  {
		  string fn(foldername);
		  vector<string> vs;
		  for (int i = 1; i <= num; i++)
		  {
			    char buffer[512];
			memset(buffer, 0, 512 * sizeof(char));
			sprintf_s(buffer, "%d.obj", i);
				string iname(buffer);
			vs.push_back(fn + "\\" + iname);
		  }
		  this->LoadMesh(vs);
          return true;
	  }




	  bool MatFeatureVector::LoadMeshS(const char* foldername, int num)
	  {
		  string fn(foldername);
		  vector<string> vs;
		  for (int i = 1; i <= num; i++)
		  {
			  char buffer[512];
			memset(buffer, 0, 512 * sizeof(char));
			sprintf_s(buffer, "%d.obj", i);
			  string iname(buffer);
			vs.push_back(fn + "\\" + iname);
		  }
		  this->LoadMeshS(vs);
		  return true;
	  }

//	  bool MatFeatureVector::PutFVS(DMEngine &eng, std::vector<FeatureVector>& fvs/*,  const char* matname*/)
//	  {
//		int m = fvs[0].s.size() * 9;
//		  int n = fvs.size();
//		  //Utility::MyMatrix smat(fvs.size(), m);
//		mxArray* p = mxCreateDoubleMatrix(n, m, mxREAL);
//		  for (int k = 0; k < fvs.size(); k++)
//		  {
//			  int tick = 0;
//			  for (int s = 0; s < fvs[k].s.size(); s++)
//			  {
//				  for (int i = 0; i < 3; i++)
//				  {
//					  for (int j = 0; j < 3; j++)
//					  {
//						  //smat.SetElement(k,tick,fvs[k].s[s].coeff(i,j));
//						mxGetPr(p)[k + tick*n] = fvs[k].s[s].coeff(i, j);
//						  tick++;
//					  } 
//				  }
//			  }
//		  }
//		  engPutVariable(eng.GetEngine(), "S", p);
//		  //DMMatrix(eng, "S",n,m,smat.data);
//		  // mat.PutVariable()
//                   
//// 		  m = fvs[0].r.size()*9;
//// 		  n = fvs.size();
//// 		  Utility::MyMatrix rmat(fvs.size(),m);
//// 
//// 		  for (int k = 0; k < fvs.size(); k++)
//// 		  {
//// 			  int tick = 0;
//// 			  for (int s = 0; s < fvs[k].r.size(); s++)
//// 			  {
//// 				  for (int i = 0; i < 3; i++)
//// 				  {
//// 					  for (int j = 0; j < 3; j++)
//// 					  {
//// 						  rmat.SetElement(k,tick,fvs[k].r[s].coeff(i,j));
//// 						  tick++;
//// 					  } 
//// 				  }
//// 			  }
//// 		  }
//// 		  DMMatrix(eng, "R",n, m, rmat.data);
//// 
//// 		  //dr
//// 		  m = fvs[0].dr.size()*9;
//// 		  n = fvs.size();
//// 		  Utility::MyMatrix drmat(fvs.size(),m);
//// 
//// 		  for (int k = 0; k < fvs.size(); k++)
//// 		  {
//// 			  int tick = 0;
//// 			  for (int s = 0; s < fvs[k].dr.size(); s++)
//// 			  {
//// 				  for (int i = 0; i < 3; i++)
//// 				  {
//// 					  for (int j = 0; j < 3; j++)
//// 					  {
//// 						  drmat.SetElement(k,tick,fvs[k].dr[s].coeff(i,j));
//// 						  tick++;
//// 					  } 
//// 				  }
//// 			  }
//// 		  }
//// 		  DMMatrix(eng, "DR",n,m,drmat.data);
//
//		  //logdr
//		m = fvs[0].logdr.size() * 9;
//		  n = fvs.size();
//		  //Utility::MyMatrix logdrmat(fvs.size(),m);
//		mxArray* plogdr = mxCreateDoubleMatrix(n, m, mxREAL);
//		  for (int k = 0; k < fvs.size(); k++)
//		  {
//			  int tick = 0;
//			  for (int s = 0; s < fvs[k].logdr.size(); s++)
//			  {
//				  for (int i = 0; i < 3; i++)
//				  {
//					  for (int j = 0; j < 3; j++)
//					  {
//						  //logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
//						mxGetPr(plogdr)[k + tick*n] = fvs[k].logdr[s].coeff(i, j);
//						  tick++;
//					  } 
//				  }
//			  }
//		  }
//		  //DMMatrix(eng, "LOGDR",n ,m, logdrmat.data);
//		  engPutVariable(eng.GetEngine(), "LOGDR", plogdr);
//		  eng.Eval("save('E:/SIGA2014/workspace/fv.mat');");
////		  eng.Eval("save('F:/SIGA2014/workspace/fv.mat');");
////		  system("pause");
//		  mxDestroyArray(p);
//		  mxDestroyArray(plogdr);
//		  // save mat
//		  //eng.Eval("save('E:/SIGA2014/workspace/fv.mat');");
////		  system("pause");
//		  return true;
//	  }


	  bool MatFeatureVector::PutFVS(std::vector<FeatureVector>& fvs/*,  const char* matname*/)
	  {
		  int m = fvs[0].s.size() * 9;
		  int n = fvs.size();
		  //Utility::MyMatrix smat(fvs.size(), m);
		  //mxArray* p = mxCreateDoubleMatrix(n, m, mxREAL);
		  std::vector<Eigen::MatrixXd> R_Ptr(fvs.size());
		  std::vector<Eigen::MatrixXd> S_Ptr(fvs.size());
		  for (int k = 0; k < fvs.size(); k++)
		  {
			  int tick = 0;
			  std::vector<double> S_Ptr_(fvs[k].s.size() * 9, 0);
			  for (int s = 0; s < fvs[k].s.size(); s++)
			  {				  
				  for (int i = 0; i < 3; i++)
				  {
					  for (int j = 0; j < 3; j++)
					  {
						  //smat.SetElement(k,tick,fvs[k].s[s].coeff(i,j));
						  S_Ptr_[tick++] = fvs[k].s[s].coeff(i, j);
						  //tick++;
					  }
				  }
			  }
			  Eigen::MatrixXd SPtr = Eigen::VectorXd::Map(&S_Ptr_[0], S_Ptr_.size());
			  S_Ptr[k] = SPtr.transpose();
		  }
		  //engPutVariable(eng.GetEngine(), "S", p);
		  //DMMatrix(eng, "S",n,m,smat.data);
		  // mat.PutVariable()

// 		  m = fvs[0].r.size()*9;
// 		  n = fvs.size();
// 		  Utility::MyMatrix rmat(fvs.size(),m);
// 
// 		  for (int k = 0; k < fvs.size(); k++)
// 		  {
// 			  int tick = 0;
// 			  for (int s = 0; s < fvs[k].r.size(); s++)
// 			  {
// 				  for (int i = 0; i < 3; i++)
// 				  {
// 					  for (int j = 0; j < 3; j++)
// 					  {
// 						  rmat.SetElement(k,tick,fvs[k].r[s].coeff(i,j));
// 						  tick++;
// 					  } 
// 				  }
// 			  }
// 		  }
// 		  DMMatrix(eng, "R",n, m, rmat.data);
// 
// 		  //dr
// 		  m = fvs[0].dr.size()*9;
// 		  n = fvs.size();
// 		  Utility::MyMatrix drmat(fvs.size(),m);
// 
// 		  for (int k = 0; k < fvs.size(); k++)
// 		  {
// 			  int tick = 0;
// 			  for (int s = 0; s < fvs[k].dr.size(); s++)
// 			  {
// 				  for (int i = 0; i < 3; i++)
// 				  {
// 					  for (int j = 0; j < 3; j++)
// 					  {
// 						  drmat.SetElement(k,tick,fvs[k].dr[s].coeff(i,j));
// 						  tick++;
// 					  } 
// 				  }
// 			  }
// 		  }
// 		  DMMatrix(eng, "DR",n,m,drmat.data);

		  //logdr
		  m = fvs[0].logdr.size() * 9;
		  n = fvs.size();
		  //Utility::MyMatrix logdrmat(fvs.size(),m);
		  //mxArray* plogdr = mxCreateDoubleMatrix(n, m, mxREAL);
		  for (int k = 0; k < fvs.size(); k++)
		  {
			  int tick = 0;
			  std::vector<double> R_Ptr_(fvs[k].logdr.size() * 9, 0);
			  for (int s = 0; s < fvs[k].logdr.size(); s++)
			  {
				  for (int i = 0; i < 3; i++)
				  {
					  for (int j = 0; j < 3; j++)
					  {
						  //logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
						  R_Ptr_[tick++] = fvs[k].logdr[s].coeff(i, j);
						  //tick++;
					  }
				  }
			  }
			  Eigen::MatrixXd RPtr = Eigen::VectorXd::Map(&R_Ptr_[0], R_Ptr_.size());
			  R_Ptr[k] = RPtr.transpose();
		  }
		  //DMMatrix(eng, "LOGDR",n ,m, logdrmat.data);
		  //engPutVariable(eng.GetEngine(), "LOGDR", plogdr);
		  //eng.Eval("save('E:/SIGA2014/workspace/fv.mat');");
		  std::ofstream LOGRNEW;
		  std::ofstream S;
		  LOGRNEW.open("D:/LOGRNEW.txt", std::ofstream::out);
		  S.open("D:/S.txt", std::ofstream::out);
		  std::copy(R_Ptr.begin(), R_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(LOGRNEW, "\n"));
		  std::copy(S_Ptr.begin(), S_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(S, "\n"));
		  LOGRNEW.close();
		  S.close();
		  //		  eng.Eval("save('F:/SIGA2014/workspace/fv.mat');");
		  //		  system("pause");
		  //mxDestroyArray(p);
		  //mxDestroyArray(plogdr);
		  // save mat
		  //eng.Eval("save('E:/SIGA2014/workspace/fv.mat');");
//		  system("pause");
		  return true;
	  }

//	  bool MatFeatureVector::OutputR = false;

	 // bool MatFeatureVector::PutFVSWithR(DMEngine &eng, std::vector<FeatureVector>& fvs)
	 // {
		//int m = fvs[0].s.size() * 9;
		//  int n = fvs.size();
		//  //Utility::MyMatrix smat(fvs.size(), m);
		//mxArray* p = mxCreateDoubleMatrix(n, m, mxREAL);
		//  for (int k = 0; k < fvs.size(); k++)
		//  {
		//	  int tick = 0;
		//	  for (int s = 0; s < fvs[k].s.size(); s++)
		//	  {
		//		  for (int i = 0; i < 3; i++)
		//		  {
		//			  for (int j = 0; j < 3; j++)
		//			  {
		//				  //smat.SetElement(k,tick,fvs[k].s[s].coeff(i,j));
		//				mxGetPr(p)[k + tick*n] = fvs[k].s[s].coeff(i, j);
		//				  tick++;
		//			  } 
		//		  }
		//	  }
		//  }
		//  engPutVariable(eng.GetEngine(), "S", p);
		//  //R
		//m = fvs[0].r.size() * 9;
		//  n = fvs.size();
		//  //Utility::MyMatrix smat(fvs.size(), m);
		//mxArray* pr = mxCreateDoubleMatrix(n, m, mxREAL);
		//  for (int k = 0; k < fvs.size(); k++)
		//  {
		//	  int tick = 0;
		//	  for (int s = 0; s < fvs[k].r.size(); s++)
		//	  {
		//		  Eigen::Matrix3d logr = log(fvs[k].r[s]);
  //                
		//		  for (int i = 0; i < 3; i++)
		//		  {
		//			  for (int j = 0; j < 3; j++)
		//			  {
		//				  //smat.SetElement(k,tick,fvs[k].s[s].coeff(i,j));
		//				mxGetPr(pr)[k + tick*n] = logr.coeff(i, j);
		//				  tick++;
		//			  } 
		//		  }
		//	  }
		//  }
		//  engPutVariable(eng.GetEngine(), "LOGR", pr);
		// 
		//  //LOGR NEW
		//m = fvs[0].rots.size() * 9;
		//  n = fvs.size();
		//  //Utility::MyMatrix smat(fvs.size(), m);
		//mxArray* prnew = mxCreateDoubleMatrix(n, m, mxREAL);
		//  for (int k = 0; k < fvs.size(); k++)
		//  {
		//	  int tick = 0;
		//	  for (int s = 0; s < fvs[k].rots.size(); s++)
		//	  {
		//		  Eigen::Matrix3d logr = fvs[k].rots[s].logr;  
		//		  Eigen::Matrix3d logr1 = log(fvs[k].r[s]);
		//		double _norm = (logr - logr1).squaredNorm();
		//		if (_norm > 0.0001)
		//		  {
		//			  //cout<<"error out"<<endl;
		//		  }
		//		  for (int i = 0; i < 3; i++)
		//		  {
		//			  for (int j = 0; j < 3; j++)
		//			  {
		//				  //smat.SetElement(k,tick,fvs[k].s[s].coeff(i,j));
		//				mxGetPr(prnew)[k + tick*n] = logr.coeff(i, j);
		//				  tick++;
		//			  } 
		//		  }
		//	  }
		//  }
		//  engPutVariable(eng.GetEngine(), "LOGRNEW", prnew);
		//  //theta
		//  m = fvs[0].rots.size();
		//  n = fvs.size();
		//  //Utility::MyMatrix smat(fvs.size(), m);
		//mxArray* theta = mxCreateDoubleMatrix(n, m, mxREAL);
		//  for (int k = 0; k < fvs.size(); k++)
		//  {
		//	  int tick = 0;
		//	  for (int s = 0; s < fvs[k].rots.size(); s++)
		//	  {
		//		  double _theta = fvs[k].rots[s].ToAngle();
		//		mxGetPr(theta)[k + tick*n] = _theta;
		//		  tick++;
		//	  }
		//  }
		//  engPutVariable(eng.GetEngine(), "THETA", theta);
		////axis
		//mxArray* axis = mxCreateDoubleMatrix(n * 3, m, mxREAL);
		//for (int k = 0; k < fvs.size(); k++)
		//{
		//	int tick = 0;
		//	for (int s = 0; s < fvs[k].rots.size(); s++)
		//	{
		//		// double _theta = fvs[k].rots[s].ToAngle();
		//		mxGetPr(axis)[k * 3 + tick*n * 3 + 0] = fvs[k].rots[s].axis(0);
		//		mxGetPr(axis)[k * 3 + tick*n * 3 + 1] = fvs[k].rots[s].axis(1);
		//		mxGetPr(axis)[k * 3 + tick*n * 3 + 2] = fvs[k].rots[s].axis(2);
		//		tick++;
		//	}
		//}
		//engPutVariable(eng.GetEngine(), "axis", axis);
		//  //logdr
		//m = fvs[0].logdr.size() * 9;
		//  n = fvs.size();
		//  //Utility::MyMatrix logdrmat(fvs.size(),m);
		//mxArray* plogdr = mxCreateDoubleMatrix(n, m, mxREAL);
		//  for (int k = 0; k < fvs.size(); k++)
		//  {
		//	  int tick = 0;
		//	  for (int s = 0; s < fvs[k].logdr.size(); s++)
		//	  {
		//		  for (int i = 0; i < 3; i++)
		//		  {
		//			  for (int j = 0; j < 3; j++)
		//			  {
		//				  //logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
		//				mxGetPr(plogdr)[k + tick*n] = fvs[k].logdr[s].coeff(i, j);
		//				  tick++;
		//			  } 
		//		  }
		//	  }
		//  }
		//  //DMMatrix(eng, "LOGDR",n ,m, logdrmat.data);
		//  engPutVariable(eng.GetEngine(), "LOGDR", plogdr);

		//  eng.Eval("save('E:/SIGA2014/workspace/fv.mat');");
		//  //		  eng.Eval("save('F:/SIGA2014/workspace/fv.mat');");
		//  //		  system("pause");
		//  mxDestroyArray(p);
		//  mxDestroyArray(plogdr);
		//  mxDestroyArray(pr);
		//  mxDestroyArray(prnew);
		//  mxDestroyArray(theta);
		//  return true;
	 // }


bool MatFeatureVector::PutFVSWithR(std::vector<FeatureVector>& fvs)
{
	int m = fvs[0].s.size() * 9;
	int n = fvs.size();
	//Utility::MyMatrix smat(fvs.size(), m);
	std::vector<Eigen::MatrixXd> R_Ptr(fvs.size());
	std::vector<Eigen::MatrixXd> S_Ptr(fvs.size());
	std::vector<Eigen::MatrixXd> RN_Ptr(fvs.size());
	std::vector<Eigen::MatrixXd> Theta_Ptr(fvs.size());
	std::vector<Eigen::MatrixXd> Axis_Ptr(fvs.size());
	std::vector<Eigen::MatrixXd> DR_Ptr(fvs.size());
	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> S_Ptr_(fvs[k].s.size() * 9, 0);
		for (int s = 0; s < fvs[k].s.size(); s++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					//smat.SetElement(k,tick,fvs[k].s[s].coeff(i,j));
					S_Ptr_[tick++] = fvs[k].s[s].coeff(i, j);
					//tick++;
				}
			}
		}
		Eigen::MatrixXd SPtr = Eigen::VectorXd::Map(&S_Ptr_[0], S_Ptr_.size());
		S_Ptr[k] = SPtr.transpose();
	}
	
	//R
	m = fvs[0].r.size() * 9;
	n = fvs.size();
	//Utility::MyMatrix smat(fvs.size(), m);
	//mxArray* plogdr = mxCreateDoubleMatrix(n, m, mxREAL);
	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> R_Ptr_(fvs[k].r.size() * 9, 0);
		
		for (int s = 0; s < fvs[k].r.size(); s++)
		{
			Eigen::Matrix3d logr = log(fvs[k].r[s]);
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					//logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
					R_Ptr_[tick++] = logr.coeff(i, j);
					//tick++;
				}
			}
		}
		Eigen::MatrixXd RPtr = Eigen::VectorXd::Map(&R_Ptr_[0], R_Ptr_.size());
		R_Ptr[k] = RPtr.transpose();
	}

	//LOGR NEW
	m = fvs[0].rots.size() * 9;
	n = fvs.size();
	//Utility::MyMatrix smat(fvs.size(), m);
	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> RN_Ptr_(fvs[k].rots.size() * 9, 0);
		for (int s = 0; s < fvs[k].rots.size(); s++)
		{
			Eigen::Matrix3d logr = fvs[k].rots[s].logr;
			Eigen::Matrix3d logr1 = log(fvs[k].r[s]);
			double _norm = (logr - logr1).squaredNorm();
			if (_norm > 0.0001)
			{
				//cout<<"error out"<<endl;
			}
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					//logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
					RN_Ptr_[tick++] = logr.coeff(i, j);
					//tick++;
				}
			}
		}
		Eigen::MatrixXd RNPtr = Eigen::VectorXd::Map(&RN_Ptr_[0], RN_Ptr_.size());
		RN_Ptr[k] = RNPtr.transpose();
	}

	//theta
	m = fvs[0].rots.size();
	n = fvs.size();
	//Utility::MyMatrix smat(fvs.size(), m);
	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> Theta_Ptr_(fvs[k].rots.size(), 0);
		for (int s = 0; s < fvs[k].rots.size(); s++)
		{
			double _theta = fvs[k].rots[s].ToAngle();
			Theta_Ptr_[tick++] = _theta;
		}
		Eigen::MatrixXd ThetaPtr = Eigen::VectorXd::Map(&Theta_Ptr_[0], Theta_Ptr_.size());
		Theta_Ptr[k] = ThetaPtr.transpose();
	}

	//axis
	//mxArray* axis = mxCreateDoubleMatrix(n * 3, m, mxREAL);
	//for (int k = 0; k < fvs.size(); k++)
	//{
	//	int tick = 0;
	//	for (int s = 0; s < fvs[k].rots.size(); s++)
	//	{
	//		// double _theta = fvs[k].rots[s].ToAngle();
	//		mxGetPr(axis)[k * 3 + tick * n * 3 + 0] = fvs[k].rots[s].axis(0);
	//		mxGetPr(axis)[k * 3 + tick * n * 3 + 1] = fvs[k].rots[s].axis(1);
	//		mxGetPr(axis)[k * 3 + tick * n * 3 + 2] = fvs[k].rots[s].axis(2);
	//		tick++;
	//	}
	//}


	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> Axis_Ptr_(fvs[k].rots.size()*3, 0);
		for (int s = 0; s < fvs[k].rots.size(); s++)
		{
			//double _theta = fvs[k].rots[s].ToAngle();
			Axis_Ptr_[tick++] = fvs[k].rots[s].axis(0);
			Axis_Ptr_[tick++] = fvs[k].rots[s].axis(1);
			Axis_Ptr_[tick++] = fvs[k].rots[s].axis(2);
			//tick++;
			//Axis_Ptr_[tick++] = _theta;
		}
		Eigen::MatrixXd AxisPtr = Eigen::VectorXd::Map(&Axis_Ptr_[0], Axis_Ptr_.size());
		Axis_Ptr[k] = AxisPtr.transpose();
	}
	//logdr
	m = fvs[0].logdr.size() * 9;
	n = fvs.size();
	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> DR_Ptr_(fvs[k].logdr.size() * 9, 0);
		for (int s = 0; s < fvs[k].logdr.size(); s++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					//logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
					DR_Ptr_[tick++] = fvs[k].logdr[s].coeff(i, j);
					//tick++;
				}
			}
		}
		Eigen::MatrixXd DRPtr = Eigen::VectorXd::Map(&DR_Ptr_[0], DR_Ptr_.size());
		DR_Ptr[k] = DRPtr.transpose();
	}

	//eng.Eval("save('E:/SIGA2014/workspace/fv.mat');");
	std::ofstream LOGR;
	std::ofstream S;
	std::ofstream LOGRNEW;
	std::ofstream Theta;
	std::ofstream Axis;
	std::ofstream LOGDR;

	LOGR.open("D:/LOGR.txt", std::ofstream::out);
	S.open("D:/S.txt", std::ofstream::out);
	LOGRNEW.open("D:/LOGRNEW.txt", std::ofstream::out);
	Theta.open("D:/Theta.txt", std::ofstream::out);
	Axis.open("D:/Axis.txt", std::ofstream::out);
	LOGDR.open("D:/LOGDR.txt", std::ofstream::out);
	std::cout << "D:/LOGR.txt" << std::endl;
	std::copy(R_Ptr.begin(), R_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(LOGR, "\n"));
	std::cout << "D:/S.txt" << std::endl;
	std::copy(S_Ptr.begin(), S_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(S, "\n"));
	std::cout << "D:/LOGRNEW.txt" << std::endl;
	std::copy(RN_Ptr.begin(), RN_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(LOGRNEW, "\n"));
	std::cout << "D:/Theta.txt" << std::endl;
	std::copy(Theta_Ptr.begin(), Theta_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(Theta, "\n"));
	std::cout << "D:/Axis.txt" << std::endl;
	std::copy(Axis_Ptr.begin(), Axis_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(Axis, "\n"));
	std::cout << "D:/LOGDR.txt" << std::endl;
	std::copy(DR_Ptr.begin(), DR_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(LOGDR, "\n"));
	LOGR.close();
	S.close();
	LOGRNEW.close();
	Theta.close();
	Axis.close();
	LOGDR.close();
	//		  eng.Eval("save('F:/SIGA2014/workspace/fv.mat');");
	//		  system("pause");


	return true;
}



	//bool MatFeatureVector::PutFVS(DMEngine &eng, std::vector<FeatureVector>& fvs, std::vector<std::pair<int, int>>& eidmap/*,const char* matname = "Tmp"*/)
	//  {
	//	int m = fvs[0].s.size() * 9;
	//	  int n = fvs.size();
	//	mxArray* p = mxCreateDoubleMatrix(n, m, mxREAL);
	//	  for (int k = 0; k < fvs.size(); k++)
	//	  {
	//		  int tick = 0;
	//		  for (int s = 0; s < fvs[k].s.size(); s++)
	//		  {
	//			  for (int i = 0; i < 3; i++)
	//			  {
	//				  for (int j = 0; j < 3; j++)
	//				  {
	//					mxGetPr(p)[k + tick*n] = fvs[k].s[s].coeff(i, j);
	//					  tick++;
	//				  } 
	//			  }
	//		  }
	//	  }
	//	  engPutVariable(eng.GetEngine(), "S", p);
	//	  //logdr
	//	m = fvs[0].logdr.size() * 9;
	//	  n = fvs.size();
	//	  //Utility::MyMatrix logdrmat(fvs.size(),m);
	//	mxArray* plogdr = mxCreateDoubleMatrix(n, m, mxREAL);
	//	  for (int k = 0; k < fvs.size(); k++)
	//	  {
	//		  int tick = 0;
	//		  for (int s = 0; s < fvs[k].logdr.size(); s++)
	//		  {
	//			  for (int i = 0; i < 3; i++)
	//			  {
	//				  for (int j = 0; j < 3; j++)
	//				  {
	//					  //logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
	//					mxGetPr(plogdr)[k + tick*n] = fvs[k].logdr[s].coeff(i, j);
	//					  tick++;
	//				  } 
	//			  }
	//		  }
	//	  }
	//	  engPutVariable(eng.GetEngine(), "LOGDR", plogdr);
	//	  //eidmap
	//	  m = 2;
	//	  n = eidmap.size();
	//	  //Utility::MyMatrix logdrmat(fvs.size(),m);
	//	mxArray* peidmap = mxCreateDoubleMatrix(n, m, mxREAL);
	//	  for (int k = 0; k < eidmap.size(); k++)
	//	  {
	//		  int tick = 0;
	//		mxGetPr(peidmap)[k + tick*n] = eidmap[k].first;
	//		  tick++;
	//		mxGetPr(peidmap)[k + tick*n] = eidmap[k].second;
	//		  tick++;
	//	  }
	//	  engPutVariable(eng.GetEngine(), "eidmap", peidmap);

	//	  eng.Eval("save('E:/SIGA2014/workspace/fv.mat');");
	//	  mxDestroyArray(p);
	//	  mxDestroyArray(plogdr);
	//	  mxDestroyArray(peidmap);
	//	  return true;
	//  }

bool MatFeatureVector::PutFVS(std::vector<FeatureVector>& fvs, std::vector<std::pair<int, int>>& eidmap/*,const char* matname = "Tmp"*/)
{
	int m = fvs[0].s.size() * 9;
	int n = fvs.size();

	std::vector<Eigen::MatrixXd> S_Ptr(fvs.size());
	std::vector<Eigen::MatrixXd> DR_Ptr(fvs.size());
	std::vector<Eigen::MatrixXd> eid_Ptr(fvs.size());
	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> S_Ptr_(fvs[k].s.size() * 9, 0);
		for (int s = 0; s < fvs[k].s.size(); s++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					//smat.SetElement(k,tick,fvs[k].s[s].coeff(i,j));
					S_Ptr_[tick++] = fvs[k].s[s].coeff(i, j);
					//tick++;
				}
			}
		}
		Eigen::MatrixXd SPtr = Eigen::VectorXd::Map(&S_Ptr_[0], S_Ptr_.size());
		S_Ptr[k] = SPtr.transpose();
	}
	//logdr
	m = fvs[0].logdr.size() * 9;
	n = fvs.size();
	//Utility::MyMatrix logdrmat(fvs.size(),m);
	for (int k = 0; k < fvs.size(); k++)
	{
		int tick = 0;
		std::vector<double> DR_Ptr_(fvs[k].logdr.size() * 9, 0);
		for (int s = 0; s < fvs[k].logdr.size(); s++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					//logdrmat.SetElement(k,tick, fvs[k].logdr[s].coeff(i,j));
					DR_Ptr_[tick++] = fvs[k].logdr[s].coeff(i, j);
					//tick++;
				}
			}
		}
		Eigen::MatrixXd DRPtr = Eigen::VectorXd::Map(&DR_Ptr_[0], DR_Ptr_.size());
		DR_Ptr[k] = DRPtr.transpose();
	}
	//eidmap
	m = 2;
	n = eidmap.size();

	for (int k = 0; k < eidmap.size(); k++)
	{

		std::vector<double> eid_Ptr_(2, 0);

		int tick = 0;
		eid_Ptr_[tick++] = eidmap[k].first;

		eid_Ptr_[tick++] = eidmap[k].second;
		
		Eigen::MatrixXd eidPtr = Eigen::VectorXd::Map(&eid_Ptr_[0], eid_Ptr_.size());
		eid_Ptr[k] = eidPtr.transpose();
	}

	std::ofstream LOGDR;
	std::ofstream S;
	std::ofstream Eid;


	LOGDR.open("D:/LOGDR.txt", std::ofstream::out);
	S.open("D:/S.txt", std::ofstream::out);
	Eid.open("D:/Eidmap.txt", std::ofstream::out);

	std::cout << "D:/LOGR.txt" << std::endl;
	std::copy(DR_Ptr.begin(), DR_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(LOGDR, "\n"));
	std::cout << "D:/S.txt" << std::endl;
	std::copy(S_Ptr.begin(), S_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(S, "\n"));
	std::cout << "D:/LOGRNEW.txt" << std::endl;
	std::copy(eid_Ptr.begin(), eid_Ptr.end(), std::ostream_iterator<Eigen::MatrixXd>(Eid, "\n"));

	LOGDR.close();
	S.close();
	Eid.close();


	return true;
}

	  bool MatFeatureVector::LoadMesh(const vector<string>& meshname)
	  {

 		  //DMEngine eng;
 		  //DTriMesh mesh1, mesh2;
 		  vector<DTriMesh> vdms(meshname.size());
		  vector<DTriMesh*> msp(meshname.size());
 		  cerr << "load mesh" << endl;
		for (int i = 0; i < meshname.size(); i++)
		  {
			cout << i << endl;
			OpenMesh::IO::read_mesh(vdms[i], meshname[i]);
			  msp[i] = &vdms[i];
		  }
// 		  //OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bars/2.obj");
// 		  //OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bars/4.obj");
// 		  OpenMesh::IO::read_mesh(mesh1, "E:/SIGA2014/dataset/bar/bars/2.obj");
// 		  OpenMesh::IO::read_mesh(mesh2, "E:/SIGA2014/dataset/bar/bars/3.obj");
// 		  vector<DTriMesh*> ms;
// 		  ms.push_back(&mesh1);
// 		  ms.push_back(&mesh2);
		  RefMesh::root = _root;
		assert(RefMesh::root < vdms[0].n_vertices());
 		  cerr << "calc feature vector" << endl;
 		  ARAPDeform deform(vdms[0], msp);
//		  deform.ref.root = _root;

		if (this->OutputR)
		  {
			this->PutFVSWithR(deform.fvs);
			  return true;
		  }

		  if (deform._eidmap.empty())
		  {
			this->PutFVS(deform.fvs);
		  }
		  else
		  {
			this->PutFVS(deform.fvs, deform._eidmap);
		  }
		//eng.Eval("exit");
		  return true;
	  }

      bool MatFeatureVector::LoadMeshWithR(const vector<string>& meshname)
	  {

 		  //DMEngine eng;
 		  //DTriMesh mesh1, mesh2;
 		  vector<DTriMesh> vdms(meshname.size());
		  vector<DTriMesh*> msp(meshname.size());
 		  cerr << "load mesh" << endl;
		for (int i = 0; i < meshname.size(); i++)
		  {
			cout << i << endl;
			OpenMesh::IO::read_mesh(vdms[i], meshname[i]);
			  msp[i] = &vdms[i];
		  }
 		  cerr << "calc feature vector" << endl;
 		  ARAPDeform deform(vdms[0], msp);
		this->PutFVSWithR(deform.fvs);
/*
		  if (deform._eidmap.empty())
		  {
			 this->PutFVS(eng,deform.fvs);
		  }
		  else
		  {
			this->PutFVS(eng,deform.fvs,deform._eidmap);
		  }
*/		  
		  return true;
	  }



	  bool MatFeatureVector::LoadMeshS(const vector<string>& meshname)
	  {
		  //DMEngine eng;
		  //DTriMesh mesh1, mesh2;
		  vector<DTriMesh> vdms(meshname.size());
		  vector<DTriMesh*> msp(meshname.size());
		  cerr << "load mesh" << endl;
		for (int i = 0; i < meshname.size(); i++)
		  {
			cout << i << endl;
			OpenMesh::IO::read_mesh(vdms[i], meshname[i]);
			  msp[i] = &vdms[i];
		  }
		  // 		  //OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bars/2.obj");
		  // 		  //OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bars/4.obj");
		  // 		  OpenMesh::IO::read_mesh(mesh1, "E:/SIGA2014/dataset/bar/bars/2.obj");
		  // 		  OpenMesh::IO::read_mesh(mesh2, "E:/SIGA2014/dataset/bar/bars/3.obj");
		  // 		  vector<DTriMesh*> ms;
		  // 		  ms.push_back(&mesh1);
		  // 		  ms.push_back(&mesh2);
		  cerr << "calc feature vector" << endl;
		  //ARAPDeform deform(eng, vdms[0], msp);
		  ARAPDeform deform;
		  //DMEngine::firstopen = true;
		  //DMEngine eng(false);
		  //deform.SetEngine(&eng);
      	  deform.ref.SetMesh(vdms[0]);
	  	  deform.meshs = msp;
		  deform.init();
		  deform.needAlign = true;
		  deform.maxIterTime = 100;
		  //deform.init();
		deform.GetFeatureS(vdms[0], msp);
		  //arap.SetEngine(dmeng);
		  //arap.GetFeature(*((DTriMesh*)&shapes[0]),ms);
		this->PutFVS(deform.fvs);
		  return true;
	  }



	  



	  bool MatFeatureVector::LoadFV(const char* logrmat, const char* smat, const char* meshname, const char* savemesh, int iternum)
	  {

#ifdef _DEBUG
         // freopen("errorlog.txt","w",stdout);      
#endif

		  //DMEngine eng(true);

		  DTriMesh mesh1;
		OpenMesh::IO::read_mesh(mesh1, meshname);
#ifdef  ALIGN
		  DTriMesh basemesh(mesh1);
#endif
		  vector<DTriMesh*> msp;
		  msp.push_back(&mesh1);
		  DTriMesh basemesh2(mesh1);
		ARAPDeform ad(basemesh2, msp);

		  //string sfvmat(fvmat);


		  Eigen::MatrixXd logr = load_txt1<Eigen::MatrixXd>(std::string(logrmat));
		  Eigen::MatrixXd s = load_txt1<Eigen::MatrixXd>(std::string(smat));


		  

		//string loadcmd = string("load('") + sfvmat + string("');");
		  //eng.Eval(loadcmd);
		  //Utility::MatEngine me;
		  //me.SetEngine(eng.GetEngine());
		  //Utility::MyMatrix NS;
		  //Utility::MyMatrix NLOGDR;
		//NS.GetVariable(me, "NS");
		//NLOGDR.GetVariable(me, "NLOGDR");
		assert(logr.rows() == s.rows());
		  if (s.rows()==1)
		  {
			  FeatureVector fv;
			  fv = ad.fvs[0];
			  int sm = 1;
			  int sn = s.cols() / 9;
			  fv.s.resize(sn);
			  int tick = 0;
			  for (int k = 0; k < sn; k++)
			  {
				  for (int i = 0; i < 3; i++)
				  {
					  for (int j = 0; j < 3; j++)
					  {
						  double tmp = 0;
						//NS.GetElement(0, tick, tmp);
						  
						fv.s[k](i, j) = s(0,tick);
						  tick++;
					  }
				  }
			  }
			  int logdrm = 1;
			  int logdrn = logr.cols() / 9;
			  fv.logdr.resize(logdrn);
			  tick = 0;
			  for (int k = 0; k < logdrn; k++)
			  {
				  for (int i = 0; i < 3; i++)
				  {
					  for (int j = 0; j < 3; j++)
					  {
						  double tmp = 0;
						//NLOGDR.GetElement(0, tick, tmp);
						  
						fv.logdr[k](i, j) = logr(0,tick);
						  tick++;
					  }
				  }
			  }

			  fv.dr.resize(logdrn);
			  for (int i = 0; i < logdrn; i++) {
				  Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
				fv.dr[i] = fv.logdr[i].exp();
			  }
//			  ad.maxIterTime =1;
			  ad.maxIterTime = iternum;
			  ad.iterEps = 0;
			ad.solve(fv, mesh1);
#ifdef       ALIGN
			RotateAlign::AlignAtoBCenter(mesh1, basemesh);
#endif
			OpenMesh::IO::write_mesh(mesh1, savemesh);
		  }
		  else
		  {
			  for (int iter = 0; iter < s.rows(); iter++)
			  {
				  FeatureVector fv;
				  fv = ad.fvs[0];
				  int sm = 1;
				  int sn = s.rows() / 9;
				  fv.s.resize(sn);
				  int tick = 0;
				  for (int k = 0; k < sn; k++)
				  {
					  for (int i = 0; i < 3; i++)
					  {
						  for (int j = 0; j < 3; j++)
						  {
							  double tmp = 0;
							//NS.GetElement(iter, tick, tmp);
							fv.s[k](i, j) = s(iter,tick);
							  tick++;
						  }
					  }
				  }
				  int logdrm = 1;
				int logdrn = logr.cols() / 9;
				  fv.logdr.resize(logdrn);
				  tick = 0;
				  for (int k = 0; k < logdrn; k++)
				  {
					  for (int i = 0; i < 3; i++)
					  {
						  for (int j = 0; j < 3; j++)
						  {
							  double tmp = 0;
							//NLOGDR.GetElement(iter, tick, tmp);
							fv.logdr[k](i, j) = logr(iter,tick);
							  tick++;
						  }
					  }
				  }

				  fv.dr.resize(logdrn);
				for (int i = 0; i < logdrn; i++) {
					  Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
					fv.dr[i] = fv.logdr[i].exp();
				  }
				ad.maxIterTime = 1;
				ad.solve(fv, mesh1);
				  string savefolder(savemesh);
				  char buffer[256];
				memset(buffer, 0, 256 * sizeof(int));
				sprintf_s(buffer, "\\%d.obj", iter + 1);
				  string iname(buffer);
				string itermesh = savefolder + iname;
#ifdef     ALIGN
				RotateAlign::AlignAtoBCenter(mesh1, basemesh);
#endif
				OpenMesh::IO::write_mesh(mesh1, itermesh.c_str());
			  }
		  }
		  return true;
	  }

	  bool MatFeatureVector::LoadFVWithR(const char* logrmat, const char* smat, const char* meshname, const char* savemesh, int iternum/* = 1*/)
	  {

#ifdef _DEBUG
		  // freopen("errorlog.txt","w",stdout);      
#endif
		 // DMEngine eng(true);

		  DTriMesh mesh1;
		OpenMesh::IO::read_mesh(mesh1, meshname);
#ifdef  ALIGN
		  DTriMesh basemesh(mesh1);
#endif
		  vector<DTriMesh*> msp;
		  msp.push_back(&mesh1);
		  DTriMesh basemesh2(mesh1);
		ARAPDeform ad(basemesh2, msp);
		  ad.densemesh.SetMesh(basemesh2);

		  Eigen::MatrixXd NR = load_txt1<Eigen::MatrixXd>(std::string(logrmat).c_str());
		  Eigen::MatrixXd NS = load_txt1<Eigen::MatrixXd>(std::string(smat));
		 // string sfvmat(fvmat);
		//string loadcmd = string("load('") + sfvmat + string("');");
		//  eng.Eval(loadcmd);
		  //Utility::MatEngine me;
		 // me.SetEngine(eng.GetEngine());
		 // Utility::MyMatrix NS;
		 // Utility::MyMatrix NR;
		//NS.GetVariable(me, "NS");
		//NR.GetVariable(me, "NLOGR");
		assert(NS.rows() == NR.rows());
		  if (NS.rows() == 1)
		  {
			  FeatureVector fv;
			  fv = ad.fvs[0];
			  int sm = 1;
			int sn = NS.cols() / 9;
			  fv.s.resize(sn);
			  int tick = 0;
			  for (int k = 0; k < sn; k++)
			  {
				  for (int i = 0; i < 3; i++)
				  {
					  for (int j = 0; j < 3; j++)
					  {
						  double tmp = 0;
						//NS.GetElement(0, tick, tmp);
						fv.s[k](i, j) = NS(0,tick);
						  tick++;
					  }
				  }
			  }
			  int rm = 1;
			int rn = NR.cols() / 9;
			  fv.r.resize(rn);
			  tick = 0;
			  for (int k = 0; k < rn; k++)
			  {
				  for (int i = 0; i < 3; i++)
				  {
					  for (int j = 0; j < 3; j++)
					  {
						  double tmp = 0;
						//NR.GetElement(0, tick, tmp);
						fv.r[k](i, j) = NR(0,tick);
						  tick++;
					  }
				  }
				  fv.r[k] = exp(fv.r[k]);
			  }
			  fv.t.resize(fv.s.size());
#ifdef DUSE_OPENMP
#pragma omp parallel
			  {
#pragma omp for
#endif
				  for (int i = 0; i < fv.t.size(); i++)
				  {
					fv.t[i] = fv.r[i] * fv.s[i];
				  }
DOMP_END;
			  ad.maxIterTime = iternum;
			  ad.iterEps = 0;
	          ad.preSetMatrix(fv);
				ad.solvefast(fv, mesh1);
#ifdef     ALIGN
				RotateAlign::AlignAtoBCenter(mesh1, basemesh);
#endif
				OpenMesh::IO::write_mesh(mesh1, savemesh);
		  }
		  else
		  {
			  for (int iter = 0; iter < NS.rows(); iter++)
			  {
				  FeatureVector fv;
				  fv = ad.fvs[0];
				  int sm = 1;
				int sn = NS.cols() / 9;
				  fv.s.resize(sn);
				  int tick = 0;
				  for (int k = 0; k < sn; k++)
				  {
					  for (int i = 0; i < 3; i++)
					  {
						  for (int j = 0; j < 3; j++)
						  {
							  double tmp = 0;
							//NS.GetElement(iter, tick, tmp);
							fv.s[k](i, j) = NS(iter,tick);
							  tick++;
						  }
					  }
				  }
				  int rm = 1;
				int rn = NR.cols() / 9;
				  fv.r.resize(rn);
				  tick = 0;
				  for (int k = 0; k < rn; k++)
				  {
					  for (int i = 0; i < 3; i++)
					  {
						  for (int j = 0; j < 3; j++)
						  {
							  double tmp = 0;
							//NR.GetElement(iter, tick, tmp);
							fv.r[k](i, j) = NR(iter,tick);
							  tick++;
						  }
					  }
				  }

				  fv.r.resize(rn);

//				  for (int i=0; i<rn; i++) {
//					  Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
//					  fv.dr[i]=fv.logdr[i].exp();
//				  }
				  #ifdef DUSE_OPENMP
#pragma omp parallel
			  {
#pragma omp for
#endif
				  for (int i = 0; i < fv.t.size(); i++)
				  {
						fv.t[i] = fv.r[i] * fv.s[i];
				  }
DOMP_END;
			      ad.maxIterTime = iternum;
			      ad.iterEps = 0;
	              ad.preSetMatrix(fv);
					ad.solvefast(fv, mesh1);
//				  ad.maxIterTime =1;
//				  ad.solve(fv,mesh1);
				  string savefolder(savemesh);
				  char buffer[256];
					memset(buffer, 0, 256 * sizeof(int));
					sprintf_s(buffer, "\\%d.obj", iter + 1);
				  string iname(buffer);
					string itermesh = savefolder + iname;
#ifdef     ALIGN
					RotateAlign::AlignAtoBCenter(mesh1, basemesh);
#endif
					OpenMesh::IO::write_mesh(mesh1, itermesh.c_str());
			  }
		  }
			//eng.Eval("exit");
		  return true;

	  }


	  bool MatFeatureVector::LoadFVS(const char* logrmat, const char* smat, const char* meshname, const char* savemesh)
	  {
		 // DMEngine eng;
		  DTriMesh mesh1;
			OpenMesh::IO::read_mesh(mesh1, meshname);
#ifdef  ALIGN
		  DTriMesh basemesh(mesh1);
#endif
		  vector<DTriMesh*> msp;
		  msp.push_back(&mesh1);
			ARAPDeform ad(mesh1, msp);

			Eigen::MatrixXd NLOGDR = load_txt1<Eigen::MatrixXd>(std::string(logrmat).c_str());
			Eigen::MatrixXd NS = load_txt1<Eigen::MatrixXd>(std::string(smat));
		  
		  //string sfvmat(fvmat);
		//	string loadcmd = string("load('") + sfvmat + string("');");
		  //eng.Eval(loadcmd);
		  //Utility::MatEngine me;
		  //me.SetEngine(eng.GetEngine());
		  //Utility::MyMatrix NS;
		  //Utility::MyMatrix NLOGDR;
		//	NS.GetVariable(me, "NS");
			//NLOGDR.GetVariable(me, "NLOGDR");
			assert(NS.rows() == NLOGDR.rows());
		  for (int iter = 0; iter < NS.rows(); iter++)
		  {
			    // eng.Eval("clear");
				  FeatureVector fv;
				  fv = ad.fvs[0];
				  int sm = 1;
				int sn = NS.cols() / 9;
				  fv.s.resize(sn);
				  int tick = 0;
				  for (int k = 0; k < sn; k++)
				  {
					  for (int i = 0; i < 3; i++)
					  {
						  for (int j = 0; j < 3; j++)
						  {
							  double tmp = 0;
							//NS.GetElement(iter, tick, tmp);
							fv.s[k](i, j) = NS(iter,tick);
							  tick++;
						  }
					  }
				  }
				  int logdrm = 1;
				int logdrn = NLOGDR.cols() / 9;
				  fv.logdr.resize(logdrn);
				  tick = 0;
				  for (int k = 0; k < logdrn; k++)
				  {
					  for (int i = 0; i < 3; i++)
					  {
						  for (int j = 0; j < 3; j++)
						  {
							  double tmp = 0;
							//NLOGDR.GetElement(iter, tick, tmp);
							fv.logdr[k](i, j) = NLOGDR(iter, tick);
							  tick++;
						  }
					  }
				  }

				  fv.dr.resize(logdrn);
				for (int i = 0; i < logdrn; i++) {
					  Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
					fv.dr[i] = fv.logdr[i].exp();
				  }
				ad.maxIterTime = 1;
				ad.solve(fv, mesh1);
				  string savefolder(savemesh);
				  char buffer[256];
				memset(buffer, 0, 256 * sizeof(char));
				sprintf_s(buffer, "\\%d.obj", iter + 1);
				  string iname(buffer);
				string itermesh = savefolder + iname;
#ifdef     ALIGN
				RotateAlign::AlignAtoB(mesh1, basemesh);
// 				  ALIGNMENT::ModelAlignment ma(eng.GetEngine());
// 				  ma.AlignAtoB(mesh1,)
#endif
				OpenMesh::IO::write_mesh(mesh1, itermesh.c_str());
			  }		  
		  return true;
	  }





}