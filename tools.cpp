#include "tools.h"
#include <fstream>
using namespace std;
bool EnergyStatistic::Save(const char* filename)
{
	ofstream outfile(filename);
	for (int i = 0;i < envals.size(); i++)
	{
		outfile<<envals[i]<<endl;
	}
	outfile.close();
	return true;
}