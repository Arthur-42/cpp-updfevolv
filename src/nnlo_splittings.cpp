#include "nnlo_splittings.h"
#include <cmath>
#include "specfunc.h"
#include <TMath.h>

using std::log;

/*
g->gg      : P2GGA, P2GGB, A2GG, SOFTG
q->gq      : P2GQA, A2GQ 
q_i->q_ig  : P2NSB, P2NSPA, P2NSMA, P2PSA, P2NSSA, A2QQNS, SOFTQ2
g->qq      : P2QGA, A2QG
q_i->aq_ig : done
q_i->q_jg  : A2QQ
q_i->aq_jg : done
*/

// WIP !!! Needs proofreading !!!

// regular piece of the gluon-gluon splitting function
double nnlo_splittings::P2GGA(double z, int nflav)
{
    double P2GGA0, P2GGA1, P2GGA2;
    double DL  = log(z);
    double DL1 = log(1. - z);

    P2GGA0 = 2675.8 * DL/z + 14214./ z - 144. * (DL*DL*DL*DL) + 72. * (DL*DL*DL)
             - 7471. * (DL*DL) + 274.4 * DL + 3589. * DL1 - 20852. 
             + 3968.* z - 3363. * (z*z) + 4848. * (z*z*z) 
             + DL*DL1 * ( 7305. + 8757. * DL );

    P2GGA1 = 157.27 * DL/z + 182.96 / z + 512./27. * (DL*DL*DL*DL)
             + 832./9. * (DL*DL*DL) + 491.3 * (DL*DL) + 1541. * DL
             - 320.0 * DL1 - 350.2 + 755.7 * z - 713.8 * (z*z) 
             + 559.3 * (z*z*z) + DL*DL1 * ( 26.15 - 808.7 * DL );

    P2GGA2 = - 680./(243. * z) - 32./27. * (DL*DL*DL) + 9.680 * (DL*DL)
             - 3.422 * DL - 13.878 + 153.4 * z - 187.7 * (z*z) 
             + 52.75 * (z*z*z) - DL*DL1 * (115.6 - 85.25* z + 63.23* DL);

    return P2GGA0 + nflav * ( P2GGA1 + nflav * P2GGA2 );
}

// singular piece of the gluon-gluon splitting function
double nnlo_splittings::P2GGB(double z, int nflav)
{
    return (2643.521 - nflav * 412.172 - nflav*nflav * 16./9.) / (1. - z);
}

// from ome.f assumed FS2=HM2=OR=10
// this turns off the heavy-mass corrections
double nnlo_splittings::A2GG(double z)
{
    double A0, B0; 
    double DLZ = log(z);
    double DLZ2 = DLZ*DLZ;
    double DLZ3 = DLZ2*DLZ;
    double DLM = log(1. - z);
    // DL = 0, DL2 = 0

    A0 = 4.*(1. + z)*DLZ3/3. + (6.+10.*z)*DLZ2
         + (32. + 48.*z)*DLZ
         - 8./z + 80. - 24.*z*z-48.*z;

    B0 = 4.*(1. + z)*DLZ2/3. + (52. + 88.*z)*DLZ/9.
         - 4.*z*DLM/3.
         + (556./z - 628. + 548.*z - 700.*z*z)/27.;

    return (2.*A0/3. + 1.5*B0);
}

//
double nnlo_splittings::SOFTG(double z)
{
    double B0;
    double DM = 1./(1. - z);
    
    B0 = 224./27.

    return 1.5*B0*DM;
}

// quark->gluon splitting function
double nnlo_splittings::P2GQA(double z, int nflav)
{
    double P2GQ0, P2GQ1, P2GQ2;
    double DL  = log(z);
    double DL1 = log(1. - z);

    P2GQ0 =   1189.3 * DL/z + 6163.1 / z - 4288./81. * (DL*DL*DL*DL)
            + 1568./9. * (DL*DL*DL) - 1794. * (DL*DL) + 4033. * DL
            + 400./81. * (DL1*DL1*DL1*DL1) + 2200./27. * (DL1*DL1*DL1)
            + 606.3 * (DL1*DL1) + 2193.* DL1 
            - 4307. + 489.3 * z + 1452.* (z*z) + 146.0 * (z*z*z)
            - 447.3 * (DL*DL)*DL1 - 972.9 * Y*(DL*DL);

    P2GQ1 =   71.082 * DL/z  - 46.41 / z + 128./27. * (DL*DL*DL*DL)
            + 704/81. * (DL*DL*DL) + 20.39  * (DL*DL) + 174.8 * DL
            - 400./81. * (DL1*DL1*DL1) - 68.069 * (DL1*DL1) - 296.7 * DL1
            - 183.8 + 33.35 * z - 277.9 * (z*z) + 108.6 * Y*(DL*DL)
            - 49.68 * DL*DL1;

    P2GQ2 = ( 64. * ( - 1./z + 1. + 2.* z)
            + 320.* DL1 * ( 1./z - 1. + 0.8 * z)
            + 96.* (DL1*DL1) * ( 1./z - 1. + 0.5 * z) ) / 27.;

    return (P2GQ0 + nflav * (P2GQ1 + nflav * P2GQ2));
}

//
double nnlo_splittings::A2GQ(double z)
{
    double A0;
    double DLM = log(1. - z);
    double DLM2 = DLM*DLM;

    A0 = 4.*(2./z - 2. + z)*DLM2/3. + 8.*(10./z - 10. + 8.*z)*DLM/9.
         + (448./z - 448. + 344.*z)/27.;

    return 2.*A0/3.;
}

//
double nnlo_splittings::P2NSB(double z, int nflav)
{
    return (1174.898 - nflav*183.187 - nflav*nflav*64./81.) / (1. - z);
}

//
double nnlo_splittings::P2NSPA(double z, int nflav)
{
    double P2NSPA;
    double DL  = log(z);
    double DL1 = log(1. - z);
    const double D81 = 1./81.;

    P2NSPA = 1641.1 - 3135.* z + 243.6 * (z*z) - 522.1 * (z*z*z)
             + 128.*D81 * (DL*DL*DL*DL) + 2400.*D81 * (DL*DL*DL)
             + 294.9 * (DL*DL) + 1258.* DL
             + 714.1 * DL1 + DL*DL1 * (563.9 + 256.8 * DL)
             + nflav * (-197.0 + 381.1 * z + 72.94 * (z*z) + 44.79 * (z*z*z)
                     - 192.*D81 * (DL*DL*DL)  - 2608.*D81 * (DL*DL) - 152.6 * DL
                     - 5120.*D81 * DL1 - 56.66 * DL*DL1 - 1.497 * z*(DL*DL*DL) )
             + (nflav*nflav) * (32.* z*DL/(1.-z) * (3.* DL + 10.) + 64.
                        + (48.* (DL*DL) + 352.* DL + 384.) * (1.-z)) * D81;

    return P2NSPA; 
}

//
double nnlo_splittings::P2NSMA(double z, int nflav)
{
    double P2NSMA;
    double DL  = log(z);
    double DL1 = log(1. - z);
    const double D81 = 1./81.;

    P2NSMA = 1860.2 - 3505.* z + 297.0 * (z*z) - 433.2 * (z*z*z)
             + 116.*D81 * (DL*DL*DL*DL) + 2880.*D81 * (DL*DL*DL) 
             + 399.2 * (DL*DL) + 1465.2 * DL
             + 714.1 * DL1 + DL*DL1 * (684.0 + 251.2 * DL)
             + nflav * (-216.62 + 406.5 * z + 77.89 * (z*z) + 34.76 * (z*z*z)
                     - 256.*D81 * (DL*DL*DL)  - 3216.*D81 * (DL*DL) - 172.69 * DL 
                     - 5120.*D81 * DL1 - 65.43 * DL*DL1 - 1.136 * Y*(DL*DL*DL))
             + NF**2 * (32.* z*DL/(1.-z) * (3.* DL + 10.) + 64.
                        + (48.* (DL*DL) + 352.* DL + 384.) * (1.-z)) * D81;

    return P2NSMA;
}

// regular pure-singlet splitting function
double nnlo_splittings::P2PSA(double z, int nflav)
{
    double P2PS1, P2PS2;
    double DL  = log(z);
    double DL1 = log(1. - z);

    P2PS1 = - 3584./(27.*z) * DL - 506./ z + 160./27. * (DL*DL*DL*DL)
            - 400./9. * (DL*DL*DL) + 131.4 * (DL*DL) - 661.6 * DL
            - 5.926  * (DL1*DL1*DL1) - 9.751 * (DL1*DL1) - 72.11 * DL1
            + 177.4 + 392.9 * z - 101.4 * (z*z) - 57.04 * DL*DL1;

    P2PS2 =   256./(81.*z) + 32./27. * (DL*DL*DL) + 17.89 * (DL*DL)
            + 61.75 * DL + 1.778 * (DL1*DL1) + 5.944 * DL1 + 100.1
            - 125.2 * z + 49.26 * (z*z) - 12.59 * (z*z*z) 
            - 1.889 * DL*DL1; 

    return (1.-z) * nflav * (P2PS1 + nflav * P2PS2);
}

//
double nnlo_splittings::P2NSSA(double z, int nflav)
{
    double P2NSSA;
    const double D27 = 1./27.;
    double DL  = log(z);
    double Y1  = 1. - z;
    double DL1 = log(Y1);

    P2NSSA = Y1* ( 151.49 + 44.51 * z - 43.12 * (z*z) + 4.820 * (z*z*z) )
             + 40.*D27 * (DL*DL*DL*DL) - 80.*D27 * (DL*DL*DL) + 6.892 * (DL*DL)
             + 178.04 * DL + DL*DL1 * ( - 173.1 + 46.18 * DL )
             + Y1*DL1 * ( - 163.9 / z - 7.208 * z ); 

    return nflav * P2NSSA; 
}

//
double nnlo_splittings::A2QQNS(double z)
{
    double A0;
    double DLZ = log(z);
    double DLZ2 = DLZ*DLZ;

    A0 = (1. + z*z)*(2.*DLZ2/3. + 20.*DLZ/9.)/(1.-Z) 
         + 8.*(1. - Z)*DLZ/3. + 44./27. - 268.*z/27.;

    return 2.*A0/3.;
}

//
double nnlo_splittings::SOFTQ2(double z)
{
    double A0;
    double DM = 1./(1. - z);

    A0 = 224./27.;

    return 2.*A0/3.;
}

// gluon->quark splitting function
double nnlo_splittings::P2QGA(double z, int nflav)
{
    double P2QG1, P2QG2;
    double DL  = log(z);
    double DL1 = log(1. - z);

    P2QG1 = - 896./(3.*z) * DL - 1268.3 / z + 536./27. * (DL*DL*DL*DL) 
            - 44./3. * (DL*DL*DL) + 881.5 * (DL*DL) + 424.9 * DL 
            + 100./27. * (DL1*DL1*DL1*DL1) - 70./9. * (DL1*DL1*DL1) 
            - 120.5 * (DL1*DL1) + 104.42 * DL1
            + 2522. - 3316.* z + 2126.* (z*z)
            + DL*DL1 * (1823. - 25.22 * DL) - 252.5 * z*(DL*DL*DL);  

    P2QG2 =   1112./(243.*z) - 16./9. * (DL*DL*DL*DL) 
            - 376./27. * (DL*DL*DL) - 90.8 * (DL*DL) - 254.0 * DL   
            + 20./27. * (DL1*DL1*DL1) + 200./27. * (DL1*DL1) - 5.496 * DL1
            - 252.0  + 158.0 * z + 145.4 * (z*z) - 139.28 * (z*z*z)
            - DL*DL1 * (53.09  + 80.616 * DL) - 98.07 * z*(DL*DL)
            + 11.70 * z*(DL*DL*DL);

    return nflav * (P2QG1 + nflav * P2QG2); 
}

// from xa2hgp.f instead of the one with polylogs IS IT CORRECT ?
double nnlo_splittings::A2QG(double z)
{
    double A2HGA;
    double DL = log(z);
    double DL1 = log(1. - z);

    A2HGA = - 24.89 / z - 187.8 + 249.6 * z - 146.8 * (DL*DL) * DL1  
            - 1.556 * (DL*DL*DL)  - 3.292 * (DL*DL)  - 93.68 * DL
            - 1.111 * (DL1*DL1*DL1) - 0.400 * (DL1*DL1) - 2.770 * DL1;

    return A2HGA;
}

//
double nnlo_splittings::A2QQ(double z, int nflav)
{
    double A0;
    // S_1,2(x) = deb_2(-log(1-x))/2 
    double S121MZ = specfunc::debye_2(-log(z))/2.; // here x = 1-z
    // S_1,1(x) = dilog(x)
    double S111MZ = TMath::DiLog(1.-z); // here x = 1-z
    double ZETA2 = 1.6449340668;
    double DLZ = log(z);
    double DLZ2 = DLZ*DLZ;
    double DLZ3 = DLZ2*DLZ;

    A0 =   (1. + z)*(32.*S121MZ + 16.*DLZ*S111MZ - 16.*ZETA2*DLZ - 4.*DLZ3/3.) 
         + (32./z/3. + 8. - 8.*z - 32.*z*z/3.)*S111MZ 
         + (-32./z/3. - 8. + 8.*z + 32.*z*z/3.)*ZETA2
         + (2. + 10.*z + 16.*z*z/3.)*DLZ2
         - (56./3. + 88.*z/3. + 448.*z*z/9.)*DLZ
         - 448./z/27. - 4./3. - 124.*z/3. + 1600.*z*z/27.;

    return 2.*A0/3.;
}
