//
// Created by WANGJade on 10/12/2018.
//

#ifndef SSP1_CONSTANTS_H
#define SSP1_CONSTANTS_H

#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Constants
{
public:
    // model type //
    // model setup: 0 ~ simulation2, 1 ~ simulation1 (also can be used for real data analysis)
    int setting = 1;
    // missing mechanism: 0 ~ fully observed, 1 ~ MCAR, 2 ~ MAR, 3 ~ MNAR
    int model_type = 3;
    // SSL prior choice:  0 ~ SSL, 1 ~ adjusted SSL
    int ssl = 1;

    // mcmc choice //
    int n_mcmc = 15000; // #mcmc iterations
    int n_burn = 5000; // #mcmc burn-in
    int n_rep = 10; // #replication

    // data //
    int cn = 800; // N sample size
    int cp = 300; // P dim of beta or X
    int cq = 4; // Q dim of manifest variables in CFA

    // parameter true values declaration, refer to .cpp file to change the value //
    VectorXd beta; // regression coefficient
    VectorXd ind_beta; // sparsity of beta; non-zero element, beta_ind = 1
    double psd; // var (psd = sigma^2)
    MatrixXd clambda; // factor loading
    MatrixXd ind_clambda; // factor loading indicator, 1 = free
    VectorXd psx; // var of CFA model

    // hyper parameters //
    // beta_k in LoI, SSL
    double lambda1 = 1; // nu_1 parameter in \varphi_1
    double lambda0 = 30; // nu_0 parameter in \varphi_0(.)
    // psi_delta in LoI, prior invert gamma (a1,b1)
    double a1 = 6;
    double b1 = 10;
    // pi, prior beta(a2,b2)
    double a2 = 1;
    double b2 = 1 * cp;
    // clambda0 normal prior for factor loading
    MatrixXd clambda0; // prior mean, refer to cpp to change the value
    double sig0_clambda = 1000; // prior var
    // psi_epsilon Ivert-Gamma(a3,b3) prior
    double a3 = 6;
    double b3 = 10;
    // tau IG(a4,b4) prior, for adjusted SSL, a4 will be specified in Mcmc.cpp
    double b4 = 1;

    // tunings //
    // MH tuning for beta_k in LoI
    VectorXd tun_beta;
    // MH tuning for tau_k in adjusted SSL
    VectorXd tun_tau;

    // missing  //
    double p_miss = 0.4; // missing completely at random prob
    int cr = 2; // number of covariates in missing data model, NOT include intercept (thus total dim = cr+cq+1)
    MatrixXd betam; // para in missing data model, NOTE dim = cr+cq+1, under MAR set ind for v = 0
    MatrixXd ind_betam; // 1~free, 0~fixed
    // MH tuning for betam
    MatrixXd tun_betam;
    // normal prior, refer to .cpp to change its value
    MatrixXd betam0;
    double sig0_betam;
    // MH tuning for v_miss, length=q vector
    VectorXd tun_v;

    // set indicator (one element is fixed or free) of betam according to model_type, and ajust betam accordingly
    int setIndBetam(void);

    Constants();
};


#endif //SSP1_CONSTANTS_H
