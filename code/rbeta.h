//
// Created by WANGJade on 11/12/2018.
//

#ifndef SSP2_RBETA_H
#define SSP2_RBETA_H

// generate beta random number from gamma random number
/*
 * x ~ gam(a, theta)
 * y ~ gam(b, theta), let theta = 1
 * x/(x+y) beta(a,b)
*/

#include <iostream>
#include <Eigen/Dense> // Eigen
#include <random> // random number generator
#include "rgamma.h"
using namespace Eigen;

VectorXd rbeta(int n, double a, double b);

#endif //SSP2_RBETA_H
