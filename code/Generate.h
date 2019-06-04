//
// Created by WANGJade on 10/12/2018.
//

#ifndef SSP1_GENERATE_H
#define SSP1_GENERATE_H

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>
#include <vector>
#include <list>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include "Constants.h"
#include "rnorm.h"
#include "runif.h"
#include "runifVec.h"
#include "strToVec.h"

using namespace std;
using namespace Eigen;

class Generate
{
private:
    // covariates cn * cp
    MatrixXd x;
    // response col vec len=cn
    VectorXd y;
    // manifest variale in CFA
    MatrixXd v;
    // missing indicator of v, =1 miss
    MatrixXd r;
    // cov in missing data model
    MatrixXd xm1; // n*r, no "ones" for intercept
    MatrixXd xm; // n*(r+1)

public:
    // gene data
    int generateX(Constants constants); // generate x
    int generateY(Constants constants); // generate y
    int generateV(Constants constants); // generate v, manifest variable in CFA
    int generateR(Constants constants); // missing indicator of v, =1 miss
    int generateXm(Constants constants); // cov for missing n*(r+1)
    // get data
    MatrixXd getX(void);
    VectorXd getY(void);
    MatrixXd getV(void);
    MatrixXd getR(void);
    MatrixXd getXm(void);
    // write data
    int writeX(void);
    int writeY(void);
    int writeV(void);
    int writeR(void);
    int writeXm(void);

    Generate(Constants constants);
};


#endif //SSP1_GENERATE_H
