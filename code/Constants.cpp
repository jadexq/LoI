//
// Created by WANGJade on 10/12/2018.
//

#include "Constants.h"

Constants::Constants() // use initializer list to construct the matrices
        : beta(cp),
          tun_beta(cp),
          clambda(2, cq),
          ind_clambda(2, cq),
          clambda0(2, cq),
          betam(cr+cq+1, cq),
          ind_betam(cr+cq+1,cq),
          betam0(cr+cq+1,cq),
          tun_betam(cr+cq+1,cq),
          tun_v(cq),
          tun_tau(cp)
{
    // parameter true values //
    // beta, col vec, non-zero value
    beta = VectorXd::Constant(cp, 1);
    // the sparsity of beta
    ind_beta = VectorXd::Constant(cp, 0);
    ind_beta(0) = 0.8;
    ind_beta(1) = 0.8;
    ind_beta(2) = 0.8;
    ind_beta(49) = 0.8;
    ind_beta(99) = 0.8;
    ind_beta(149) = 0.8;
    ind_beta(199) = 0.8;
    // ind_beta(299) = 0.8;
    // sparse beta
    beta = (beta.array() * ind_beta.array()).matrix();
    // psd
    psd = 0.3;
    // tunings //
    // tuning for beta
    tun_beta = VectorXd::Constant(cp, 5); // proposal 2
    // tuning for tau
    tun_tau = VectorXd::LinSpaced(cp, 1.2, 0.6);

    // CFA //
    // clambda, factor loading
    clambda = MatrixXd::Constant(2, cq, 0); // value of intercept = 0
    clambda.row(1) = RowVectorXd::Constant(cq, 0.8); // true value of slop = 0.8
    clambda(1,0) = 1; // fixed at 1
    ind_clambda = MatrixXd::Constant(2, cq, 1);
    ind_clambda(1,0) = 0; // fixed element
    clambda0 = MatrixXd::Constant(2, cq, 0);
    // psx
    psx = VectorXd::Constant(cq, 0.3);

    // missing data model //
    betam = MatrixXd::Constant(cr+cq+1, cq, 0.5); // para in missing data model, NOTE dim = cr+cq+1
    betam.row(0) = RowVectorXd::Constant(cq, -0.5); // intercept
    ind_betam = MatrixXd::Constant(cr+cq+1, cq, 1); // 1~free, 0~fixed at zero
    // tuning for betam
    tun_betam = MatrixXd::Constant(cr+cq+1, cq, 5); // 5 ajust pvar to aqrt(pvar) elementwise; 45 elementwise; 1.3 rowwise
    // prior
    betam0 = MatrixXd::Constant(cr+cq+1, cq, 0);
    sig0_betam = 1000;
    // tuning v
    tun_v = VectorXd::Constant(cq, 1); // 1
}

// set ind of betam and ajust betam accordingly
int Constants::setIndBetam()
{
    // nig = 1 only when model_type=3, MNAR
    int nig = model_type/3;
    // ind for v, NOTE assume it is diagonal!
    VectorXd diag_vec = VectorXd::Constant(cq, nig);
    MatrixXd diag = MatrixXd::Constant(cq, cq, 0);
    diag.diagonal() = diag_vec;
    // set ind_betam
    ind_betam.bottomRows(cq) = diag;
    // ajust betam according to ind_betam
    betam = (betam.array() * ind_betam.array()).matrix();

//    cout << "ind_betam" << ind_betam << endl;
//    cout << "betam" << betam << endl;

    return 0;
}
