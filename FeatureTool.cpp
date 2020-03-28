#include "FeatureTool.h"

int GetFeatureWithR(const char* foldername, int objnum, int _root)
{
	MatFV::MatFeatureVector mfv;
	mfv.OutputR = true;
	mfv._root = _root;
	mfv.LoadMesh(foldername, objnum);
	return 1;
}

int ReconstructWithR(const char* ininame)
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
	int iternum = ini.GetInt("default:iternum",-1);
	MatFV::MatFeatureVector mfv;
	if (iternum==-1)
	{
		mfv.LoadFVWithR(logrmat,smat,meshname, savename);
	}
	else
	{
		assert(iternum>=0);
		mfv.LoadFVWithR(logrmat, smat, meshname,savename,iternum);
	}
	return 1;
}