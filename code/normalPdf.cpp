//
// Created by WANGJade on 28/10/2018.
//

// for compute normal density
// arguments: x, m, s; all float
// output normal prob; float

#include "normalPdf.h"

double normalPdf(double x, double m, double s)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;
    // return the value of pdf
    return inv_sqrt_2pi / s * exp(-0.5f * a * a);
    // return the value of log pdf
//    return log(inv_sqrt_2pi) - log(s) - 0.5f * a * a;
}