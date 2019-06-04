//
// Created by WANGJade on 28/10/2018.
//

// generate uniform random numbers
// arguments n vector length
// output VectorXd of random numbers of length n

#include "runifVec.h"


using namespace std;
using namespace Eigen;

VectorXd runifVec(int n)
{
    VectorXd result(n);
    for (int i=0; i<n; i++)
    {
        result(i) = runif();
    }

    return result;
}

