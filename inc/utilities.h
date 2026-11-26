#ifndef UTILITIES_H
#define UTILITIES_H

namespace utilities
{
    // to find the range [x1, x2] in which func(x) = 0
    void brackit(double (*func)(double), double& x1, double& x2);
    // to find a root of func(x) within range [x1, x2]
    double bisect(double (*func)(double), double x1, double x2, double epsi);
    // to determine if a is approximately equal to b, within tolerance
    bool lmb_eq(double a, double b, double tol);
    // to determine if a is approximately equal to b, within tolerance, and with a < b
    bool lmb_lt(double a, double b, double tol);
    // to perform an adaptative Gaussian integration 
    // author T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
    // reference THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
    double dgadap(double a0, double b0, double (*f)(double), double eps);
}

#endif 