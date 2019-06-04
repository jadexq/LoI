//
// Created by WANGJade on 16/10/2018.
//

// mvnorm(n, mean, cov)
// returen: MatrixXd (dim, n)
// use rnorm function

// description:
// for generating multivariate normal random numbers
// use Cholesky decomposition

// arguments:
// n - int, number of observations, >= 1
// mean - double, vector of means
// cov - double, matrix of covariance

// begin:

#include "mvnorm.h"
using namespace Eigen;

MatrixXd mvnorm(int n, VectorXd mean, MatrixXd cov)
{
    // get the dim of cov
    int dim = 0;
    dim = cov.cols();
    // compute Cholesky decomposition (matrix L) for positive definite matrix
    LLT<MatrixXd> lltOfcov(cov);
    MatrixXd L = lltOfcov.matrixL();
    // random mvnorm:
    VectorXd rnorm_vec = rnorm(dim*n);
    Map<MatrixXd> result(rnorm_vec.data(),dim,n);

    result = L * result ; // adjust the cov
    result.colwise() += mean; // add the mean vec

    return result;
}

