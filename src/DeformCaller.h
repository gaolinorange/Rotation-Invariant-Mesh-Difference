#ifndef DEFORMCALLER_H
#define DEFORMCALLER_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include "DeformInfo.h"

class DeformCaller {
public:
    DeformInfo info;

    virtual void calc() {};
};

class DVCaller : public DeformCaller {
public:
    void calc();
};


void CallMain();
#endif // DEFORMCALLER_H
