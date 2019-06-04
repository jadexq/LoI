//
// Created by WANGJade on 17/10/2018.
//

// here gamma(k, theta)
// E = k * theta
// Var = k * theta^2

// however in R: R_gamma (a, b), E = a/b, Var = a/b^2

// R_gamma(a, b) is same as gamma(a, 1/b)

#include "rgamma.h"
using namespace Eigen;

VectorXd rgamma(int n, double a, double b)
{
    VectorXd result = VectorXd::Zero(n); // store the results
    std::random_device rd;
    //static std::default_random_engine generator;
    std::default_random_engine generator;
    generator.seed(rd());
    //static std::gamma_distribution<double> distribution(a, b);
    std::gamma_distribution<double> distribution(a, b);
    for (int i = 0; i < n; i++)
    {
        double number = distribution(generator);
        result(i) = number;
    }

    return result;
}
