#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <ctime>
#include <iostream>

#include "ARAPDeform.h"
#include "MatFeatureVector.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <direct.h>
//#include "iniparser\INIHelper.h"
#include "INIHelper.h"
#include "FVAnalysis.h"
//#include "deformsf.h"
#include "FeatureTool.h"


using namespace std;
using namespace Eigen;

int flag = 0;
void arapdeformm_(char * doptFile);
void ckAlign() {
	srand(time(0));
	Matrix3d t = Matrix3d::Random();
	Matrix3d r,s;
	polarDec(t,r,s);
	vector<Vector3d> v1, v2;
	for (int i=0; i<2; i++) {
		Vector3d v(Vector3d::Random());
		v1.push_back(v);
		//v2.push_back(r*v);
		v2.push_back(Vector3d::Random());
	}
	RotateAlign align(v1, v2);
	Matrix3d rr = align.calc();

	cout << rr << endl;
	cout << r << endl;
	cout << (r-rr.transpose()).norm() << endl;

	align.ckt(rr);
}




int test0()
{
	//DMEngine eng(true);
	DTriMesh mesh1, mesh2;
	cerr << "load mesh" << endl;
	//OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bars/2.obj");
	//OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bars/4.obj");
	OpenMesh::IO::read_mesh(mesh1, "F:/workspace/card/0.obj");
	OpenMesh::IO::read_mesh(mesh2, "F:/workspace/card/30.obj");
	vector<DTriMesh*> ms;
	ms.push_back(&mesh1);
	ms.push_back(&mesh2);
	cerr << "calc feature vector" << endl;
	ARAPDeform deform(mesh1, ms);
	deform.needAlign = true;
	deform.maxIterTime = 100;
	vector<double> weight;
	double t = 0;
	weight.push_back(1-t);
	weight.push_back(t);
	DTriMesh result(mesh1);
	cerr << "solve" << endl;
	FeatureVector fv(weight, deform.fvs);
	fv.loadConstPoint(ifstream("F:/workspace/card/cardhandle2.txt"));

	//deform.solve(fv, result);
	deform.solve2(fv, weight, result);
	cout<<"weight"<<endl;
	for (int i = 0; i < weight.size(); i++) cout<<"weight["<<i<<"]="<<weight[i];
	cerr << "write mesh" << endl;
	OpenMesh::IO::write_mesh(result, "F:/workspace/card/result_solve2.obj");

	ARAPDeform deform1(mesh1, ms);
	deform1.needAlign = true;
	deform1.maxIterTime = 100;
	FeatureVector fvweight(weight, deform1.fvs);
	deform1.solve(fvweight,result);
	OpenMesh::IO::write_mesh(result, "F:/workspace/card/result_weight1.obj");
	system("pause");
	return 1;
}

int main1() {

	//DMEngine eng;
	DTriMesh mesh1, mesh2;
	cerr << "load mesh" << endl;
	//OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bars/2.obj");
	//OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bars/4.obj");
	OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/box3D/box3D.obj");
	OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/box3D/box3D_90.obj");
	vector<DTriMesh*> ms;
	ms.push_back(&mesh1);
	ms.push_back(&mesh2);
	cerr << "calc feature vector" << endl;
	ARAPDeform deform(mesh1, ms);
	deform.needAlign = true;
	deform.maxIterTime = 100;
	vector<double> weight;
	double t = 0;
	weight.push_back(1-t);
	weight.push_back(t);
	DTriMesh result(mesh1);
	cerr << "solve" << endl;
	FeatureVector fv(weight, deform.fvs);
	fv.loadConstPoint(ifstream("C:/Users/CJLD/Desktop/Lab/box3D/constPoint/ts2.txt"));

	deform.solve(fv, result);
	//deform.solve(fv, weight, result);

	cerr << "write mesh" << endl;
	OpenMesh::IO::write_mesh(result, "result.obj");

	return 0;
}

int test5(const char* foldername, int modelnum, const char* csname, const char* saveobj)
{
	//DMEngine eng;
	//    DTriMesh mesh1, mesh2;
	vector<DTriMesh> meshs(modelnum);
	vector<DTriMesh*> ms;
	cerr << "load mesh" << endl;
	//OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bars/2.obj");
	//OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bars/4.obj");
	string _folder(foldername);

	for(int i = 0; i < modelnum; i++)
	{
		char buffer[256];
		memset((char*)buffer, 0, 256*sizeof(char));
		sprintf((char*)buffer, "\\%d.obj", i+1);
		std::string iname = _folder+std::string(buffer);
		OpenMesh::IO::read_mesh(meshs[i],iname.c_str());
		ms.push_back(&meshs[i]);
		//OpenMesh::IO::read_mesh(meshs[i],"");
	}
	/*
	OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/box3D/box3D.obj");
	OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/box3D/box3D_90.obj");

	ms.push_back(&mesh1);
	ms.push_back(&mesh2);
	*/
	cerr << "calc feature vector" << endl;
	//    ARAPDeform deform(eng, mesh1, ms);
	ARAPDeform deform( meshs[0], ms);
	deform.needAlign = true;
	deform.maxIterTime = 100;
	vector<double> weight(modelnum,0);
	double t = 0;
	weight[0] = 1;
	/*
	weight.push_back(1-t);
	weight.push_back(t);
	*/
	DTriMesh result(meshs[0]);
	cerr << "solve" << endl;
	FeatureVector fv(weight, deform.fvs);
	//  fv.loadConstPoint(ifstream("C:/Users/CJLD/Desktop/Lab/box3D/constPoint/ts2.txt"));
	fv.loadConstPoint(ifstream(csname));
	//deform.solve(fv, result);
	//deform.solve(fv, weight, result);
	deform.solve2(fv,weight,result);
	cerr << "write mesh" << endl;
	//OpenMesh::IO::write_mesh(result, "result.obj");
	OpenMesh::IO::write_mesh(result,saveobj);
	return 0;
}

int test1(const char* listname)
{
	//DMEngine eng;
	// 	DTriMesh mesh1, mesh2;
	// 	cerr << "load mesh" << endl;
	// 	//OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bars/2.obj");
	// 	//OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bars/4.obj");
	// 	OpenMesh::IO::read_mesh(mesh1, "E:/SIGA2014/dataset/bar/bars/2.obj");
	// 	OpenMesh::IO::read_mesh(mesh2, "E:/SIGA2014/dataset/bar/bars/3.obj");
	// 
	// 	vector<DTriMesh*> ms;
	// 	ms.push_back(&mesh1);
	// 	ms.push_back(&mesh2);
	//const char* listname = "E:/SIGA2014/workspace/file.txt";
	// 	vector<DTriMesh> meshs;
	// 	vector<DTriMesh*> pms;
	MatFV::MatFeatureVector mfv;
	mfv.LoadMesh(listname);
	// 	cerr << "calc feature vector" << endl;
	// 	ARAPDeform deform(eng, mesh1, ms);

	// 	deform.needAlign = true;
	// 	deform.maxIterTime = 100;
	// 	vector<double> weight;
	// 	double t = 1;
	// 	weight.push_back(1-t);
	// 	weight.push_back(t);
	// 	DTriMesh result(mesh1);
	// 	cerr << "solve" << endl;
	// 	FeatureVector fv(weight, deform.fvs);
	// 	//fv.loadConstPoint(ifstream("C:/Users/CJLD/Desktop/Lab/box3D/constPoint/ts2.txt"));
	// 	deform.solve(fv, result);
	// 	//deform.solve(fv, weight, result);
	// 	cerr << "write mesh" << endl;
	// 	OpenMesh::IO::write_mesh(result, "E:/SIGA2014/dataset/bar/bars/result/result.obj");

	return 1;
}

int test2(const char* ininame)
{

	// 	const char* fvmat = "E:/SIGA2014/workspace/nfv.mat";
	// 	const char* meshname = "E:/SIGA2014/dataset/scape/1.obj";
	// 	const char* savename =  "E:/SIGA2014/workspace/1_49_5.obj";

	Base::Utility::MyINIParser ini;
	if (!ini.LoadINIFile(ininame))
	{
		cout<<"Load File Error"<<endl;
		return 0;
	}

	cout<<ini.GetNumSections()<<endl;
	const char* logrmat = ini.GetString("default:logrmat",NULL);
	const char* smat = ini.GetString("default:smat", NULL);
	const char* meshname = ini.GetString("default:meshname",NULL);
	const char* savename = ini.GetString("default:savename",NULL);
	int iternum = ini.GetInt("default:iternum",-1);
	MatFV::MatFeatureVector mfv;
	if (iternum==-1)
	{
		mfv.LoadFV(logrmat, smat, meshname, savename);
	}
	else
	{
		assert(iternum>=0);
		mfv.LoadFV(logrmat, smat, meshname,savename,iternum);
	}
	return 1;
}

int GenShapeWithR(const char* ininame)
{
	Base::Utility::MyINIParser ini;
	if (!ini.LoadINIFile(ininame))
	{
		cout<<"Load File Error"<<endl;
		return 0;
	}

	cout<<ini.GetNumSections()<<endl;
	const char* logrmat = ini.GetString("default:logrmat", NULL);
	const char* smat = ini.GetString("default:smat", NULL);
	const char* meshname = ini.GetString("default:meshname",NULL);
	const char* savename = ini.GetString("default:savename",NULL);
	int iternum = ini.GetInt("default:iternum",-1);
	MatFV::MatFeatureVector mfv;
	if (iternum==-1)
	{
		//mfv.LoadFV(fvmat,meshname, savename);
		mfv.LoadFVWithR(logrmat, smat,meshname, savename);
	}
	else
	{
		assert(iternum>=0);
		//mfv.LoadFV(fvmat, meshname,savename,iternum);
		mfv.LoadFVWithR(logrmat,smat,meshname,savename,iternum);
	}
	return 1;
}

int test3(const char* foldername, int objnum)
{
	// 	const char* foldername = "E:/SIGA2014/dataset/scape";
	// 	int objnum = 71;
	MatFV::MatFeatureVector mfv;
	mfv.LoadMesh(foldername, objnum);
	return 1;
}

int GetFetureWithR(const char* foldername, int objnum, int _root = 0)
{
	MatFV::MatFeatureVector mfv;
	mfv._root = _root;
	mfv.LoadMeshWithR(foldername, objnum);
	return 1;
}

int test8(const char* foldername, int objnum)
{
	return 1;
}

int test9(const char* foldername, int filenum, const char* constname, const char* outname)
{
	ARAPDeformFun(foldername, filenum, constname, outname);
	return 1;
}

int test6(const char* foldername, int objnum)
{
	MatFV::MatFeatureVector mfv;
	mfv.LoadMeshS(foldername, objnum);
	return 1;
}

int test4(const char* ininame)
{
	Base::Utility::MyINIParser ini;
	if (!ini.LoadINIFile(ininame))
	{
		cout<<"Load File Error"<<endl;
		return 0;
	}
	cout<<ini.GetNumSections()<<endl;
	const char* logrmat = ini.GetString("default:logrmat",NULL);
	const char* smat = ini.GetString("default:smat", NULL);
	const char* meshname = ini.GetString("default:meshname",NULL);
	const char* savename = ini.GetString("default:savename",NULL);
	MatFV::MatFeatureVector mfv;
	mfv.LoadFVS(logrmat, smat,meshname, savename);
	return 1;
}

void testFVA() {
	//DMEngine eng;
	vector<DTriMesh*> ms;
	string path = "E:/project/SIGA2014/dataset/horse/";
	int num = 49;
	cout << "Input Path : " << endl;
	cin >> path;
	cout << "Input Number : " << endl;
	cin >> num;
	for (int i=1; i <= num; i++) {
		DTriMesh *mesh = new DTriMesh();
		stringstream ss;
		ss << path << "/" << i << ".obj";
		cout << "Read Mesh : " << ss.str() << endl;
		OpenMesh::IO::read_mesh(*mesh, ss.str().c_str());
		ms.push_back(mesh);
	}
	cout << "Ref Id : " << endl;
	int x;
	cin >> x;
	swap(ms[0], ms[x]);
	ARAPDeform deform(*ms[0], ms);
	cout << "ARAP construct done" << endl;
	FVAnalysis fva(&deform);
	cout << "FVA construct done" << endl;
	fva.work();
}

void FVA(const char* pathname, int num, int refid = 0)
{
	//DMEngine eng;
	vector<DTriMesh*> ms;
	string path = string(pathname); //"E:/project/SIGA2014/dataset/horse/";
	string pathpca = path+string("\\pca\\");
	mkdir(pathpca.c_str());
	//int num = 49;
	cout << "Input Path : " << endl;
	//cin >> path;
	cout << "Input Number : " << endl;
	//cin >> num;
	for (int i=1; i <= num; i++) {
		DTriMesh *mesh = new DTriMesh();
		stringstream ss;
		ss << path << "/" << i << ".obj";
		cout << "Read Mesh : " << ss.str() << endl;
		OpenMesh::IO::read_mesh(*mesh, ss.str().c_str());
		ms.push_back(mesh);
	}
	cout << "Ref Id : " << endl;
	int x = refid;
	//cin >> x;
	swap(ms[0], ms[x]);
	ARAPDeform deform(*ms[0], ms);
	cout << "ARAP construct done" << endl;
	FVAnalysis fva(&deform);
	fva.pathpca = pathpca;
	cout << "FVA construct done" << endl;
	fva.work();
}


void FVALOG(const char* pathname, int num, int refid = 0)
{
	//DMEngine eng;
	vector<DTriMesh*> ms;
	string path = string(pathname); //"E:/project/SIGA2014/dataset/horse/";
	string pathpca = path+string("\\pca\\");
	mkdir(pathpca.c_str());
	//int num = 49;
	cout << "Input Path : " << endl;
	//cin >> path;
	cout << "Input Number : " << endl;
	//cin >> num;
	for (int i=1; i <= num; i++) {
		DTriMesh *mesh = new DTriMesh();
		stringstream ss;
		ss << path << "/" << i << ".obj";
		cout << "Read Mesh : " << ss.str() << endl;
		OpenMesh::IO::read_mesh(*mesh, ss.str().c_str());
		ms.push_back(mesh);
	}
	cout << "Ref Id : " << endl;
	int x = refid;
	//cin >> x;
	swap(ms[0], ms[x]);
	ARAPDeform deform(*ms[0], ms);
	cout << "ARAP construct done" << endl;
	FVAnalysis fva(&deform);
	fva.pathpca = pathpca;
	cout << "FVA construct done" << endl;
	fva.worklog();
}


// int test2(const char* ininame)
// {
// 
// 	// 	const char* fvmat = "E:/SIGA2014/workspace/nfv.mat";
// 	// 	const char* meshname = "E:/SIGA2014/dataset/scape/1.obj";
// 	// 	const char* savename =  "E:/SIGA2014/workspace/1_49_5.obj";
// 
// 	Base::Utility::MyINIParser ini;
// 	if (!ini.LoadINIFile(ininame))
// 	{
// 		cout<<"Load File Error"<<endl;
// 		return 0;
// 	}
// 
// 	cout<<ini.GetNumSections()<<endl;
// 	const char* fvmat = ini.GetString("default:fvmat",NULL);
// 	const char* meshname = ini.GetString("default:meshname",NULL);
// 	const char* savename = ini.GetString("default:savename",NULL);
// 
// 	MatFV::MatFeatureVector mfv;
// 	mfv.LoadFV(fvmat,meshname, savename);
// 	return 1;
// }
#ifdef EXE

int main(int argc, char* argv[])
{
	// 	test0();
	//testsf();
	//return 1;
	//BarDeform(atof(argv[1]));
	//return 1;
	int cmd = atoi(argv[1]);
	switch (cmd)
	{
	case 1:
		test1(argv[2]);
		break;
	case 2:
		test2(argv[2]);
		break;
	case 3:
		test3(argv[2], atoi(argv[3]));
		break;
	case 4:
		test4(argv[2]);
		break;
	case 5:
		test5(argv[2],atoi(argv[3]), argv[4], argv[5]);
		break;
	case 6:
		test6(argv[2], atoi(argv[3]));
		break;
	case 7:
		if (argc==4)
		{
			FVA(argv[2], atoi(argv[3]));
		}
		else if (argc == 5)
		{
			FVA(argv[2], atoi(argv[3]), atoi(argv[4]));
		}
		break;
	case 8:
		//test8();
		break;
	case 9:
		ARAPDeformFun(argv[2],atoi(argv[3]),argv[4],argv[5]);
		break;
	case 10:
		ARAPDeformMultiScale(argv[2],atoi(argv[3]),argv[4], argv[5], argv[6]);
		break;
	case 11:
		ARAPDeformHandle(argv[2],argv[3], argv[4],atoi(argv[5]));
		break;
	case 12:
		ARAPDeformInter(argv[2],argv[3],atoi(argv[4]),argv[5]);
		break;
	case 14:
		ARAPDeformHandleRot(argv[2], argv[3], argv[4], argv[5],atoi(argv[6]));
		break;
	case 15: 
		if (argc==4)
		{
			FVALOG(argv[2], atoi(argv[3]));
		}
		else if (argc == 5)
		{
			FVALOG(argv[2], atoi(argv[3]), atoi(argv[4]));
		}
		break;
	case 16:
		//get feature with R
		if(argc==4)
		{
			GetFeatureWithR(argv[2], atoi(argv[3]));
		}
		if(argc==5)
		{
			GetFeatureWithR(argv[2], atoi(argv[3]),atoi(argv[4]));
		}
		break;
	case 17:
		//reconstruct with R
		ReconstructWithR(argv[2]);
		break;
	case 18:
		//char * inifile = argv[2];
		arapdeformm_(argv[2]);
		//system("pause");
		break;
	case 19:

		break;
	}
	//	system("pause");
	return 1;
}

#endif
// 	if ( cmd ==0 )
// 	{
// 		test1();
// 	}
// 	else
// 	{
// 		test2();
// 	}
//	return 1;
//}

int main(int argc, char* argv[])
{
	// 	test0();
	//testsf();
	//return 1;
	//BarDeform(atof(argv[1]));
	//return 1;
	int cmd = atoi(argv[1]);
	ARAPDeformHandle(argv[2], argv[3], argv[4], atoi(argv[5]));
}