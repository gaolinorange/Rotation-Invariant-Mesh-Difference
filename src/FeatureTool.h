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

using namespace std;
using namespace Eigen;


int GetFeatureWithR(const char* foldername, int objnum, int _root = 0);
int ReconstructWithR(const char* ininame);
