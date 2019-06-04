//
// Created by WANGJade on 28/10/2018.
//

#ifndef CLASS_MCMC_NROMALPDFVEC_H
#define CLASS_MCMC_NROMALPDFVEC_H

#include <iostream>
#include <Eigen/Dense>
#include "normalPdf.h"
using namespace std;
using namespace Eigen;

VectorXd normalPdfVec(VectorXd X, VectorXd M, VectorXd S);

#endif //CLASS_MCMC_NROMALPDFVEC_H
