//
// Created by WANGJade on 17/10/2018.
//

#ifndef TEST_MVRNORM_H
#define TEST_MVRNORM_H

#include "rnorm.h"

#include <iostream>
#include <string>
#include <random> // random number generator ??? do I need it here ???
#include <Eigen/Dense> // Eigen
using namespace Eigen;

MatrixXd mvnorm(int n, VectorXd mean, MatrixXd cov);

#endif //TEST_MVRNORM_H
