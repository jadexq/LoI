//
// Created by WANGJade on 16/10/2018.
//

// rnorm(int n)
// return: Eigen::VectorXd (column vec)

// description:
// for generating standard normal random numbers

// arguments:
// n - int, number of observations, > 1

// begin:
#include "rnorm.h"
using namespace Eigen;

VectorXd rnorm(int n)
{
    VectorXd result = VectorXd::Zero(n); // store the results
    std::random_device rd;
    static std::default_random_engine generator;
    generator.seed(rd());
    static std::normal_distribution<double> distribution(0,1);
    //std::default_random_engine generator;
    //std::normal_distribution<double> distribution(0,1);
    for (int i = 0; i < n; i++)
    {
        double number = distribution(generator);
        result(i) = number;
    }

    return result;
}


