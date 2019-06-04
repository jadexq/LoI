//
// Created by WANGJade on 10/12/2018.
//

#ifndef SSP1_INITIAL_H
#define SSP1_INITIAL_H

#include <iostream>
#include <Eigen/Dense>
#include "rnorm.h"
#include "runif.h"
#include "runifVec.h"
#include "rbeta.h"
#include "Constants.h"
#include "Generate.h"
using namespace std;
using namespace Eigen;

class Initial {
public:
    // initial values of parameters
    VectorXd beta; // regression coefficient
    VectorXd tau; // ajust lam0 for beta when gam = 0
    VectorXd sdx; // sd of x
    double psd; // var (psd = sigma^2)
    // initial value of theta
    double theta;
    // initial values of indicator gamma_j
    VectorXd gamma;
    // clambda, factor loading in CFA
    MatrixXd clambda;
    // psx, residual var in CFA
    VectorXd psx;
    // initial value of latent variable, y
    VectorXd y;
    // if there is missing, initial value of manifest variable
    MatrixXd v;
    // MAR or MNAR, missing data model coefficient
    MatrixXd betam;

    // set initial value of y according to y = x*beta + e
    int setY(Constants constants, Generate generate);
    // write initial value of y (latent variable) into file
    int writeY();

    // set initial value of v miss, actually set all values in v
    int setV(Constants constants, Generate generate);
    // write v initial into file
    int writeV();

    // set initial value of tau, sd(x)
    int setTauSdx(Constants constants, Generate generate);

    Initial(Constants constants);
};


#endif //SSP1_INITIAL_H
