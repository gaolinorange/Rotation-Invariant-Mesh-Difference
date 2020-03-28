#pragma once

#include "FeatureVector.h"
//#include "mymatrix.h"
#include <string>
#include <vector>
//#include "DMEngine.h"

#define ALIGN

namespace MatFV
{
	using namespace std;
	class MatFeatureVector
	{
	public:
		//bool PutFVS(DMEngine &eng, std::vector<FeatureVector>& fvs /*,const char* matname = "Tmp"*/);
		//bool PutFVS(DMEngine &eng, std::vector<FeatureVector>& fvs,  std::vector<std::pair<int,int>>& eidmap/*,const char* matname = "Tmp"*/);
		//bool PutFVSWithR(DMEngine &eng, std::vector<FeatureVector>& fvs);
		bool PutFVS(std::vector<FeatureVector>& fvs /*,const char* matname = "Tmp"*/);
		bool PutFVS(std::vector<FeatureVector>& fvs, std::vector<std::pair<int, int>>& eidmap/*,const char* matname = "Tmp"*/);
		bool PutFVSWithR(std::vector<FeatureVector>& fvs);
 		bool LoadMesh(const char* filelist);
 		bool LoadMesh(const char* foldername, int num);
		bool LoadMeshWithR(const char* foldername, int num);
		bool LoadMeshS(const char* foldername, int num);
 		bool LoadMesh(const vector<string>& meshname);
		bool LoadMeshS(const vector<string>& meshname);
		bool LoadMeshWithR(const vector<string>& meshname);
		//bool LoadFV(const char* fvmat, const char* meshname, const char* savemesh, int iternum = 1);
		//bool LoadFVS(const char* fvmat, const char* meshname, const char* savemesh);
		//bool LoadFVWithR(const char* fvmat, const char* meshname, const char* savemesh, int iternum = 1);
		bool LoadFV(const char* logrmat, const char* smat, const char* meshname, const char* savemesh, int iternum = 1);
		bool LoadFVS(const char* logrmat, const char* smat, const char* meshname, const char* savemesh);
		bool LoadFVWithR(const char* logrmat, const char* smat, const char* meshname, const char* savemesh, int iternum = 1);
	    bool OutputR;
		int _root;
		MatFeatureVector() {OutputR = false;_root = 0;}
	};
}