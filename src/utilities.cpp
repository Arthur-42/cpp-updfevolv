#include "utilities.h"
#include <cmath>
#include <algorithm>

using std::abs, std::max, std::min, std::log;

// to find the range [x1, x2] in which func(x) = 0
void utilities::brackit(double (*func)(double), double& x1, double& x2)
/*
    *func : pointer to the function of interest
    x1    : (passed by reference) lower bound of the interval
    x2    : (passed by reference) upper bound of the interval
*/
{
    const double factor = 1.6;
    const int ntry = 50;
    double f1, f2;

    f1 = func(x1);
    f2 = func(x2);
    for(int j=1; j <= ntry; j++){
        if(f1*f2 < 0.) return; // if func flips sign, it has a root
        // expand the bounds to continue search
        if(abs(f1) < abs(f2)){
            x1 = x1 + factor*(x1 - x2);
            f1 = func(x1);
        }
        else{
            x2 = x2 + factor*(x2 - x1);
            f2 = func(x2);
        }
    }

    return;
}

// to find a root of func(x) within range [x1, x2]
double utilities::bisect(double (*func)(double), double x1, double x2, double epsi)
/*
    *func : pointer to the function of interest
    x1    : lower bound of the interval
    x2    : upper bound of the interval
    epsi  : tolerance
*/
{
    const int jmax = 40;
    double fmid, f, output, dx, xmid;

    // initialisation
    fmid = func(x2);
    f = func(x1);
    // we want to approach the root from the negative part of the function
    if(f < 0.){
        output = x1;
        dx = x2 - x1; // the xmid will be increased
    }
    else{
        output = x2;
        dx = x1 - x2; // the xmid will be decreased
    }

    // find the root
    for(int j = 1, j <= jmax; j++){
        dx = 0.5 * dx;      // use a finer increment
        xmid = output + dx; // increment xmid
        fmid = func(xmid);
        if(fmid <= 0.) output = xmid; // we don't update if f(xmid) becomes positive
        if((abs(dx) < epsi) || (fmid == 0)) return output; // we got sufficiently close to zero
    }

    return 1.E11; // if none found (must be amended)
}

// to determine if a is approximately equal to b, within tolerance
bool utilities::lmb_eq(double a, double b, double tol)
/*
    a,b : values to compare
    tol : tolerance (if positive, absolute; if negative, relative)
*/
{
    double eps;
    // absolute or relative, that is the question
    if(tol>0) eps = tol;
    else eps = abs(tol)*max(abs(a),abs(b));

    return abs(a-b) < eps;
}

// to determine if a is approximately equal to b, within tolerance, and with a < b
bool utilities::lmb_lt(double a, double b, double tol)
/*
    a,b : values to compare
    tol : tolerance (if positive, absolute; if negative, relative)
*/
{
    double eps;
    // absolute or relative, that is the question
    if(tol>0) eps = tol;
    else eps = abs(tol)*max(abs(a),abs(b));

    return (abs(a-b) > eps) && ((a-b)<0);
}