#include "nlo_splittings.h"
#include <cmath>
#include <TMath.h>

using std::log;

/*
g->gg      : GG1SFUNC
q->gq      : FG1SFUNC
q_i->q_ig  : PP1SFUNC, PM1SFUNC, FF1SFUNC
g->qq      : GF1SFUNC
q_i->aq_ig : done
q_i->q_jg  : done
q_i->aq_jg : done
*/

static const double PI = 3.141592653589793;
static const double CPI2S3 = PI*PI/3.;
static const double C11S6 = 11./6.;
static const double C20S3 = 20./3.;
static const double C2S3 = 2./3.;
static const double C5S3 = 5./3.;
static const double C4S3 = 4./3.;
static const double C16S9 = 16./9.;
static const double C4S9 = 4./9.;
static const double CPIA = 67./18. - CPI2S3/2.;
static const double CPIE = 5. - CPI2S3;
static const double C182S9 = 182./9.;
static const double C14S9 = 14./9.;
static const double C40S9 = 40./9.;
static const double C136S3 = 136./3.;
static const double C38S3 = 38./3.;
static const double C44S3 = 44./3.;
static const double CPIF = CPI2S3 - 218./9.;
static const double C14S3 = 14./3.;
static const double C16S3 = 16./3.;
static const double C10S9 = 10./9.;
static const double C40S3 = 40./3.;
static const double C112S9 = 112./9.;

// WIP !!!

// 
double GG1SFUNC(double z, int nflav)
{
    double AAA, BBB, CCC;
    double CX2    = z*z;
    double C1PX   = 1. + z;
    double C1MX   = 1. - z;
    double CLX    = log(z);
    double CLX2   = CLX*CLX;
    double CL1MX  = log(C1MX);
    double CL1PX  = log(C1PX);
    double CL1MX2 = CL1MX*CL1MX;
    double CS3X   = -TMath::DiLog(-z);
    double CS2X   = 0.5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX);
    double CPGG   = 1./C1MX + 1./z -2. + z - CX2;
    double CMPGG  = 1./C1PX - 1./z -2. - z - CX2;

    AAA   = -16.+ 8.*z+ C20S3*CX2 + C4S3/z + (-6.-10.*z)*CLX + (-2.)*C1PX*CLX2;

    BBB   = 2.* C1MX +  26./9.*(CX2-1./z) - C4S3*C1PX*CLX - 20./9.*CPGG;

    CCC   = 27./2.*C1MX + 67./9.*(CX2-1./z)+(-25./3.+11./3.*z-44./3.*CX2)*CLX+4.*C1PX*CLX2+
            (67./9.-4.*CLX*CL1MX + CLX2-CPI2S3)*CPGG + 2.*CMPGG*CS2X;

    return C2S3*nflav*AAA + 1.5*nflav*BBB + 9.* CCC;
}

//
double FG1SFUNC(double z, int nflav)
{
    double AAA, BBB, CCC;
    double CX2    = z*z;
    double C1PX   = 1. + z;
    double C1MX   = 1. - z;
    double CLX    = log(z);
    double CLX2   = CLX*CLX;
    double CL1MX  = log(C1MX);
    double CL1PX  = log(C1PX);
    double CL1MX2 = CL1MX*CL1MX;
    double CPFGX  = (1.+C1MX*C1MX) / z;
    double CPFGMX = - (1.+C1PX*C1PX) / z;
    double CS3X   = -TMath::DiLog(-z);
    double CS2X   = 0.5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX);

    AAA   = -5./2.- 7./2.*z+(2.+7./2.*z)*CLX+(-1.+0.5*z)*CLX2 - 2.*z*CL1MX+ (-3.*CL1MX-CL1MX2)*CPFGX;

    BBB   = 28./9.+65./18.*z+44./9.*CX2+(-12.-5.*z-8./3.*CX2)*CLX+(4.+z)*CLX2+2.*z*CL1MX
            + (-2.*CLX*CL1MX+0.5*CLX2 + 11./3.*CL1MX+CL1MX2-0.5*CPI2S3+0.5)*CPFGX+CPFGMX*CS2X;

    CCC   = -C4S3*z- (20./9.+C4S3*CL1MX)*CPFGX;

    return C16S9*AAA + 4.*BBB + 2./3.*nflav*CCC;
}

//
double PP1SFUNC(double z, int nflav)
{
    double AAA, BBB, CCC, PQQ, PQQB;
    double CX2    = z*z;
    double C1PX   = 1. + z;
    double C1MX   = 1. - z;
    double CLX    = log(z);
    double CLX2   = CLX*CLX;
    double CL1MX  = log(C1MX);
    double CL1PX  = log(C1PX);
    double CPFFX  = (1. + CX2) / C1MX;
    double CPFFMX = (1. + CX2) / C1PX;
    double CS3X   = -TMath::DiLog(-z);
    double CS2X   = 0.5*(CLX2 - CPI2S3) + 2.*(CS3X - CLX * CL1PX);
 
    // P_F(x) in FP eq. (4.52)
    AAA = - CPFFX*2.*CLX*CL1MX - (2.*z+3./C1MX)*CLX - .5*C1PX*CLX2 - 5.*C1MX;
    // 0.5*P_G(x) in FP eq. (4.53)     
    BBB = CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX;
    // P_N(x) in FP eq. (4.54)   
    CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX;
    // P_QQ(x) in FP eq. (4.50) 
    PQQ  = C16S9*AAA + 4.*BBB + C2S3*nflav*CCC;
    // 0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
    PQQB = - C4S9 * (CPFFMX*CS2X + C1PX*CLX + 2.*C1MX);
 
    return PQQ + PQQB;
}

//
double PM1SFUNC(double z, int nflav)
{
    double AAA, BBB, CCC, PQQ, PQQB;
    double CX2    = z*z;
    double C1PX   = 1. + z;
    double C1MX   = 1. - z;
    double CLX    = log(z);
    double CLX2   = CLX*CLX;
    double CL1MX  = log(C1MX);
    double CL1PX  = log(C1PX);
    double CPFFX  = (1. + CX2) / C1MX;
    double CPFFMX = (1. + CX2) / C1PX;
    double CS3X   = -TMath::DiLog(-z);
    double CS2X   = 0.5*(CLX2 - CPI2S3) + 2.*(CS3X -CLX * CL1PX);

    // P_F(x) in FP eq. (4.52) 
    AAA = - CPFFX*2.*CLX*CL1MX - (2.*z+3./C1MX)*CLX - .5*C1PX*CLX2 - 5.*C1MX;
    // 0.5*P_G(x) in FP eq. (4.53)     
    BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX;
    // P_N(x) in FP eq. (4.54)      
    CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX;
    // P_QQ(x) in FP eq. (4.50) 
    PQQ  =   C16S9*AAA + 4.*BBB + C2S3*NF*CCC;
    // 0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
    PQQB = - C4S9 * (CPFFMX*CS2X + C1PX*CLX + 2.*C1MX);
 
    return PQQ - PQQB;
}

//
double FF1SFUNC(double z, int nflav)
{
    double AAA, BBB, CCC;
    double CX2    = z*z;
    double C1PX   = 1. + z;
    double C1MX   = 1. - z;
    double CLX    = log(z);
    double CLX2   = CLX*CLX;
    double CL1MX  = log(C1MX);
    double CL1PX  = log(C1PX);
    double CPFFX  = (1.+CX2) / C1MX;
    double CPFFMX = (1.+CX2) / C1PX;
    double CS3X   = -TMath::DiLog(-z);
    double CS2X   = 0.5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX);

    AAA = - CPFFX*CLX*(1.5+2.*CL1MX) + 2.*CPFFMX*CS2X - 1. + z + (.5-1.5*z)*CLX - .5*C1PX*CLX2;

    BBB =   CPFFX*(C11S6*CLX+.5*CLX2+CPIA) - CPFFMX*CS2X + C14S3*C1MX;

    CCC = - CPFFX*(C10S9+C2S3*CLX) + C40S9/z - 2.*C1PX*CLX2 - C16S3 + C40S3*z 
            + (10.*z+C16S3*CX2+2.)*CLX - C112S9*CX2;

    return C16S9*AAA + 4.*BBB + C2S3*nflav*CCC;
}

//
double GF1SFUNC(double z, int nflav)
{
    double AAA, DDD;
    double CX2    = z*z;
    double C1PX   = 1. + X;
    double C1MX   = 1. - z;
    double CLX    = log(z);
    double CLX2   = CLX*CLX;
    double CL1MX  = log(C1MX);
    double CL1PX  = log(C1PX);
    double CL1MX2 = CL1MX*CL1MX;
    double CPGFX  = CX2 + C1MX*C1MX;
    double CPGFMX = CX2 + C1PX*C1PX;
    double CS3X   = -TMath::DiLog(-z);
    double CS2X   = 0.5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX);

    AAA =   4. - 9.*z + (4.*z-1.)*CLX + (2.*z-1.)*CLX2 + 4.*CL1MX
            + (2.*CLX-2.*CLX*CL1MX+CLX2-2.*CL1MX+CL1MX2+CPIE) * 2. * CPGFX;

    DDD =   C182S9 + C14S9*X + C40S9/z + (C136S3*z-C38S3)*CLX - 4.*CL1MX - (2.+8.*z)*CLX2
            + 2.*CS2X*CPGFMX + (C44S3*CLX-CLX2-2.*CL1MX2+4.*CL1MX+CPIF) * CPGFX;

    return C2S3*nflav*AAA + 1.5*nflav*DDD;
}
