//
// Created by WANGJade on 10/12/2018.
//

#ifndef SSP1_MCMC_H
#define SSP1_MCMC_H

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include "Constants.h"
#include "Generate.h"
#include "Initial.h"
#include "runif.h"
#include "runifVec.h"
#include "rnorm.h"
#include "mvnorm.h"
#include "rgamma.h"
#include "normalPdf.h"
#include "normalPdfVec.h"
#include "rbeta.h"

using namespace std;
using namespace Eigen;


class Mcmc
{
private:
    // para to be updated //
    VectorXd beta;
    double psd = 0;
    // indicator, gamma
    VectorXd gamma;
    // tau
    VectorXd tau;
    VectorXd a4; // IG(a4,b4) prior for tau, b4 was specified in Constants
    // mixture prob, theta
    double theta;
    // clambda, factor loading in CFA
    MatrixXd clambda;
    // psx, residual error in CFA
    VectorXd psx;
    // data //
    MatrixXd x;
    MatrixXd v; // manifest variable in CFA
    // latent variable //
    VectorXd y;
    // missing //
    MatrixXd r; // missing indicator of v
    MatrixXd xm; // cov in missing data model
    MatrixXd betam; // para in missing data model (r+1)*q

    // calculate the log full conditional dist of Y, log(Y|.)
    double logPostBeta(VectorXd beta, Constants constants);


public:
    // count #ac of each element of beta
    VectorXd ac_beta;
    // count #ac of each element of beta
    MatrixXd ac_betam;
    // count #ac of v_miss
    MatrixXd ac_v;
    // count #ac of tau
    VectorXd ac_tau;

    // read data, x, v, r
    int readData(Constants constants, Generate generate);
    // set initial value
    int setInitial(Constants constants, Initial initial);

    // update para
    int updateBeta(Constants constants);
    // update psd
    int updatePsd(Constants constants);
    // update gamma
    int updateGamma(Constants constants);
    // update tau
    int updateTau(Constants constants);
    // update theta, the prob of mixture prior
    int updateTheta(Constants constants);
    // update clambda, the factor loading in CFA, and update psx, residual var in CFA
    int updateClambdaPsx(Constants constants);

    // update y, latent variable
    int updateY(Constants constants);

    // if there is missing in v, update v
    int updateV(Constants constants);

    // if MAR or MNAR update betam
    int updateBetam(Constants constants);


    // get current (updated) value of para
    VectorXd getBeta(void);
    double getPsd(void);
    VectorXd getGamma(void);
    VectorXd getTau(void);
    double getTheta(void);
    VectorXd getClambda(void); // clambda should be transposed to col vec first
    VectorXd getPsx(void);
    VectorXd getBetam(void); // betam should be transposed to col vec first
    // get latent variable y
    VectorXd getY(void);
    // get updated v, manifest variable
    MatrixXd getV(void);

    // constructor
    Mcmc(Constants constants);
};

#endif //SSP1_MCMC_H
