#pragma once
#include <iostream>
#include <vector>

class EnergyStatistic
{
public:
	std::vector<double> envals;
	inline void push_back(double& val) {envals.push_back(val);}
	bool Save(const char* filename);
};