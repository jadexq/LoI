//
// Created by WANGJade on 11/12/2018.
//

// generate beta random number from gamma random number
/*
 * x ~ gam(a, theta)
 * y ~ gam(b, theta), let theta = 1
 * x/(x+y) beta(a,b)
*/

#include "rbeta.h"
using namespace Eigen;

VectorXd rbeta(int n, double a, double b)
{
    VectorXd x = VectorXd::Zero(n);
    VectorXd y = VectorXd::Zero(n);
    VectorXd z = VectorXd::Zero(n); // for store random beta(a,b)
    x = rgamma(n, a, 1); // gam(a,1)
    y = rgamma(n, b, 1); // gam(b,1)
    // z = x / (x + y)
    z = x.array()/(x.array() + y.array());
    return z;
}

