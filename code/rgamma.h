//
// Created by WANGJade on 17/10/2018.
//

// here gamma(k, theta)
// E = k * theta
// Var = k * theta^2

// however in R: R_gamma (a, b), E = a/b, Var = a/b^2

// R_gamma(a, b) is same as gamma(a, 1/b)

#ifndef TEST_RGAMMA_H
#define TEST_RGAMMA_H

#include <iostream>
#include <string>
#include <random> // random number generator
#include <Eigen/Dense> // Eigen
using namespace Eigen;

VectorXd rgamma(int n, double a, double b);

#endif //TEST_RGAMMA_H
