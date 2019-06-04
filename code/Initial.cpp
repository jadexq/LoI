//
// Created by WANGJade on 10/12/2018.
//

// steps:
// 0. Initial initial(constants)
// 1. initial.setY(constants, generate)
//    initial.writeY()
// 2. initial.setY(constants) if there is missing
//    initial.writeV()

#include "Initial.h"

Initial::Initial(Constants constants)
        : beta(constants.cp),
          gamma(constants.cp),
          tau(constants.cp),
          sdx(constants.cp),
          clambda(2, constants.cq),
          psx(constants.cq),
          y(constants.cn),
          betam(constants.cr+constants.cq+1, constants.cq)
{
    // parameter initial values
    // beta, col vec
    beta = VectorXd::Constant(constants.cp, 0); // ini value
    // beta = constants.beta; // true value
    // psd, scalar
    psd = 0.5; // ini value
    // psd = constants.psd; // true value
    // initial value of theta
    theta = rbeta(1, constants.a2, constants.b2)(0); // rbeta returns a VectorXd
    // theta = double(constants.ind_beta.sum())/constants.cp; // "true" value
    // initial values of indicator gamma_j, gam_j=1 with prob theta
    gamma = (runifVec(constants.cp).array()<theta).cast<double>(); // ini value
    // gamma = constants.ind_beta; // true value
    // clambda, factor loading in CFA
    clambda = MatrixXd::Constant(2, constants.cq, 1); // ini value
    clambda(1,0) = 1; // fixed at 1 for indentification
    //clambda = constants.clambda; // true value
    // psx, residual var in CFA
    psx = VectorXd::Constant(constants.cq, 0.5); // ini value
    // psx = constants.psx; // true value
    // create y, latent variable
    y = VectorXd::Constant(constants.cn, 0);
    // betam
    betam = MatrixXd::Constant(constants.cr+constants.cq+1, constants.cq, 0); // ini value
    //betam = constants.betam; // true value
    // set tau=0 currently, later use setTau function
    tau = VectorXd::Constant(constants.cp, 0);
    sdx = VectorXd::Constant(constants.cp, 0);
}

// def initializeY
int Initial::setY(Constants constants, Generate generate)
{
    // ini value
    MatrixXd x = generate.getX();
    y = x*beta + rnorm(constants.cn)*sqrt(psd);
    // true value
    //y = generate.getY();

    return 0;
}

// def initialize v
int Initial::setV(Constants constants, Generate generate)
{
    // ini
    MatrixXd v_true = generate.getV(); // true v
    MatrixXd r = generate.getR(); // miss indicator
    // generate initial
    VectorXd rnorm_vec = rnorm(constants.cn * constants.cq);
    Map<MatrixXd> rnorm_mat(rnorm_vec.data(), constants.cn, constants.cq);
    MatrixXd v_ini = (rnorm_mat.array().rowwise()*sqrt(psx.transpose().array()) + (y*clambda.row(1)).array()).matrix();
    VectorXd ones = VectorXd::Constant(constants.cn, 1);
    v_ini = (v_ini.array() + (ones * clambda.row(0)).array()).matrix(); // intercept
    // assign
    v = (r.array() * v_ini.array() + (1-r.array()) * v_true.array()).matrix();

    // true
    // v = v_true;

    return 0;
}


// def initialize tau as sd(x)
int Initial::setTauSdx(Constants constants, Generate generate)
{
    // get x
    MatrixXd x = generate.getX();
    sdx = x.array().pow(2).colwise().sum().transpose();
    sdx = (sdx.array()/constants.cn).sqrt().matrix();
    tau = sdx;
    // debug

    if (constants.ssl == 0)
    {
        tau = ArrayXd::Constant(constants.cp, 1);
    }

    return 0;
}

int Initial::writeY(void)
{
    ofstream file;
    file.open("./out/y_ini.txt");
    file << y << endl;
    file.close();

    return 0;
}

int Initial::writeV(void)
{
    ofstream file;
    file.open("./out/v_ini.txt");
    file << v << endl;
    file.close();

    return 0;
}
