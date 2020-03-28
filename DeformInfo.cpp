#include "DeformInfo.h"

using namespace std;
using namespace Eigen;

const char DeformInfo::name1[] = "MInfo.txt";
const char DeformInfo::name2[] = "HInfo.txt";

std::string ItoS(int i) {
    std::stringstream ss;
    ss << i;
    return ss.str();
}

std::string operator +(std::string s, int i) {return s + ItoS(i);}
std::string operator +(int i, std::string s) {return ItoS(i) + s;}
