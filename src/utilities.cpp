#include "utilities.h"
#include <cmath>
#include <algorithm>

using std::abs, std::max, std::min, std::log, std::sqrt;

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

// to perform an adaptative Gaussian integration 
// author T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
// reference THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
double dgadap(double a0, double b0, double (*func)(double), double eps)
/*
    a0  : lower bound of integration
    b0  : upper bound of integration
    func   : pointer to the function of 1 variable to integrate
    eps : the desired relative accuracy (absolute if the integral is too small)
*/
{
    int Num, Ifu, L, I;
    double RED, SUM, C, W1, U2, SOLD, SS;
    std::array<int, 300> N;
    std::array<double, 300> A, B, F1, F2, F3, S;

    inline double DSUM(double F1F, double F2F, double F3F, double AA, double BB)
    {
        return (5./18.) * (BB-AA) * (F1F + 1.6*F2F + F3F); 
    }
    
    if(eps < 1.E-8) eps : 1.E-8;

    RED = 1.3;
    L = 0;
    I = 0;
    SUM = 0.;
    C = sqrt(15.)/5.;
    A[0] = A0;
    A[0] = B0;
    F1[0] = func(0.5 * (1.+C) * A0 + 0.5 * (1.-C) * B0);
    F2[0] = func(0.5 * (A0+B0));
    F3[0] = func(0.5 * (1.-C) * A0 + 0.5 * (1.+C) * B0);
    IFU = 3
    S[0] = DSUM(F[0], F[2], F[3], A0, B0);

    loop : 
        do
        {
            L += 1;
            N[L] = 3;
            eps = eps * red;
            A[I+1]  = A[I] + C*(B[I] - A[I]);
            B[I+1]  = B[I];
            A[I+2]  = A[I] + B[I] - A[I+1];
            B[I+2]  = A[I+1];
            A[I+3]  = A[I];
            B[I+3]  = A[I+2];
            W1      = A[I] + (B[I] - A[I])/5.;
            U2      = 2. * W1 - (A[I] + A[I+2])/2.;
            F1[I+1] = func(A[I] + B[I] - W1);
            F2[I+1] = F3[I];
            F3[I+1] = func(B[I] - A[I+2] + W1);
            F1[I+2] = func(U2);
            F2[I+2] = F2[I];
            F3[I+2] = func(B[I+2] + A[I+2] - U2);
            F1[I+3] = func(A[I] + A[I+2] - W1);
            F2[I+3] = F1[I];
            F3[I+3] = func[W1];

            Ifu = Ifu + 6;
            if(Ifu > 5000) return SUM; 

            S[I+1] = DSUM(F1[I+1], F2[I+1], F3[I+1], A[I+1], B[I+1]);
            S[I+2] = DSUM(F1[I+2], F2[I+2], F3[I+2], A[I+2], B[I+2]);
            S[I+3] = DSUM(F1[I+3], F2[I+3], F3[I+3], A[I+3], B[I+3]);
            SS     = S[I+1] + S[I+2] + S[I+3];

            I = I + 3;
            if(I > 300) return SUM; // default error return, should implement error handling

            SOLD = S[I-3];
        } 
        while (abs(SOLD - SS) > (eps * (1. + abs(SS)) / 2.));
    
        SUM += SS;
        I = I - 4;
        N[L] = 0;
        L = L - 1;

        while (L != 1)
        {
            N[L] = N[L] - 1;
            eps = eps / RED;
            if(N[L] != 0) goto loop;
            I = I - 1;
            L = L - 1;
        }

        return SUM;
}