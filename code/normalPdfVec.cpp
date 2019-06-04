//
// Created by WANGJade on 28/10/2018.
//

// for compute the normal pdf of a vector
// arguments: X, M(mean), S(std); vector of the same length
// output vector of normal prob

#include "normalPdfVec.h"

VectorXd normalPdfVec(VectorXd X, VectorXd M, VectorXd S)
{
    int length = X.rows();
    VectorXd result(length);

    for (int i=0; i<length; i++)
    {
        result(i) = normalPdf(X(i), M(i), S(i));
    }

    return result;
}