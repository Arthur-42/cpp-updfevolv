#include "SudakovTables.h"
#include "lo_splittings.h"
#include "nlo_splittings.h"
#include "nnlo_splittings.h"
#include "AlphaStrong.h"
#include "utilities.h"
#include <cmath>
#include <algorithm>
#include <array>

using std::sqrt, std::log, std::abs, std::max, std::pow, std::exp;
namespace lo = lo_splittings;
namespace nlo = nlo_splittings;
namespace nnlo = nnlo_splittings;
namespace util = utilities;

static const double pi = 3.14159265358979;
static const double aqed = 0.0072973525693; // QED coupling (at mZ)
static const double sin2w = 0.23;  // sin^2(theta_W)
static const double cos2w = 0.77;  // cos^2(theta_W)
static const double e2u = 4./9.;   // up quark QED charge squared
static const double e2d = 1./9.;   // down quark QED charge squared
static const double t3u = 0.5;     // up quark weak isospin
static const double t3d = -0.5;   // down quark weak isospin

/* branchtypes and itype reference
0  : g      -> g g
1  : q      -> g q
2  : q_i    -> q_i g
3  : g      -> q q
4  : q_i    -> qbar_i g
5  : q_i    -> q_j g  (i neq j)
6  : q_i    -> qbar_j g  (i neq j)
7  : uq     -> photon uq
8  : uq     -> uq photon 
9  : photon -> uq uq
10 : dq     -> photon dq
11 : dq     -> dq photon
12 : photon -> dq dq
13 : uq     -> Z uq
14 : uq     -> uq Z
15 : Z      -> uq uq
16 : dq     -> Z dq
17 : dq     -> dq Z
18 : Z      -> dq dq
19 : q      -> W+ q
20 : q      -> q W+
21 : W+     -> q q
22 : q      -> W- q
23 : q      -> q W-
24 : W-     -> q q
25 : q      -> H q
26 : q      -> q H
27 : H      -> q q

itype only : 
28 : g      -> (= 0+3)
29 : uq     -> (= 1+2+4+5+6+7+8+13+14+19+20+22+23+25+26)
30 : dq     -> (= 1+2+4+5+6+10+11+16+17+19+20+22+23+25+26)
31 : photon -> (= 9+12)
32 : Z      -> (= 15+18)

Note that the remaining initial states are given in 21, 24, 27
*/

// constructor
SudakovTables::SudakovTables(double qmin, double qmax, double eps, AlphaStrong& alphas, 
                             int order, int quordering, bool smbosons, double charmmass =1.27, 
                             double bottommass =4.19, double topmass =172.6, double Wmass =80.36, double Zmass =91.188, 
                             double Hmass =125., double q0ord =1., double zmaxfixed =0.99)
/*
    qmin       : lower scale of the tables
    qmax       : upper scale of the tables
    eps        : precision of the integrations
    alphas     : instance of AlphaStrong (the running coupling)
    order      : the order of the computations 1: LO; 2: NLO; 3: NNLO
    quordering : the ordering scheme 0: angular ordering; 1: virtual ordering; 2: no z dependance in a_s(q) 
    smbosons   : do we add the non-QCD bosons contributions to the branchings
    [...]mass  : mass of the [...] quark/boson
    q0ord      : reference scale used for angular or virtual ordering
    zmaxfixed  : fixed zmax used if quordering = 2
*/
{
    m_order = order;
    m_qordering = quordering;
    m_smbosons = smbosons;
    m_charmmass = charmmass;
    m_bottommass = bottommass;
    m_topmass = topmass;
    m_Wmass = Wmass;
    m_Zmass = Zmass;
    m_Hmass = Hmass;
    m_q0ord = q0ord;
    m_zmaxfixed = zmaxfixed;

    filltables(qmin, qmax, eps, alphas);
}

// to fill the sudakov tables
void SudakovTables::filltables(double qmin, double qmax, double eps, AlphaStrong& alphas)
/*
    qmin   : lower scale of the table
    qmax   : upper scale of the table
    eps    : precision of the integrations
    alphas : instance of AlphaStrong (the running coupling) to be passed to the splitting function
*/
{
    double q2min, q2max;
    double x3lmin, x3lmax, x3ldif;
    double qnow, temp_sudz, temp_sudint;
    int nf_light, n_boson;
    //double power = 1. / (static_cast<double>(tablelength - 1.));

    x3lmin = log(qmin);
    x3lmax = log(qmax);
    x3ldif = (x3lmax - x3lmin) / (static_cast<double>(tablelength));

    for(int itype = 0; itype < 33; itype++)
    {
        q2min = qmin * qmin;

        m_sudint[itype][0] = 0.;
        if(itype == 0) m_scales[0] = q2min; // we only note the scale once
        m_sudmin[itype] = 0.;
        m_sudmax[itype] = -9999999.; // default fake value

        if(itype <= 27){
            q2max = q2min;
            m_sudz[itype][0] = integrate_splitting(q2max, q2min, itype, eps, alphas);
        }

        for(int inode = 1; inode < tablelength; inode++)
        {
            qnow = exp(x3lmin + x3ldif * inode);
            // force first eligible node to heavy quark or boson mass
            nf_light = 3;
            if((qnow >= m_charmmass) && (nf_light != 4)){
                nf_light = 4;
                qnow = m_charmmass;
            }
            if((qnow >= m_bottommass) && (nf_light != 5)){
                nf_light = 5;
                qnow = m_bottommass;
            }
            if((qnow >= m_topmass) && (nf_light != 6)){
                nf_light = 6;
                qnow = m_topmass;
            }
            n_boson = 1;
            if((qnow >= m_Wmass) && (n_boson != 3)){
                n_boson = 3;
                qnow = m_Wmass;
            }
            if((qnow >= m_Zmass) && (n_boson != 4)){
                n_boson = 4;
                qnow = m_Zmass;
            }
            if((qnow >= m_Hmass) && (n_boson != 5)){
                n_boson = 5;
                qnow = m_Hmass;
            }

            q2max = qnow * qnow;

            if(itype <= 27){
                temp_sudz = integrate_splitting(q2max, q2min, itype, eps, alphas);
                m_sudz[itype][inode] = temp_sudz;
                // approximate q integral (a_s fixed after z integration)
                temp_sudint = m_sudint[itype][inode-1] + temp_sudz * log(q2max/q2min);
                m_sudint[itype][inode] = temp_sudint;

                q2min = q2max;

                if(temp_sudint > m_sudmax[itype]) m_sudmax[itype] = temp_sudint;
            }
            // gluon initial states
            else if(itype == 28){
                m_sudz[itype][inode] = m_sudz[0][inode] + m_sudz[3][inode]; // is it useful ?
                temp_sudint = m_sudint[0][inode] + m_sudint[3][inode];
                m_sudint[itype][inode] = temp_sudint;
                if(temp_sudint < m_sudmin[itype]) m_sudmin[itype] = temp_sudint;
                if(temp_sudint > m_sudmax[itype]) m_sudmax[itype] = temp_sudint;
            }
            // up type quark initial states
            else if(itype == 29){
                m_sudz[itype][inode] = m_sudz[1][inode] + m_sudz[2][inode] + m_sudz[4][inode] + m_sudz[5][inode]
                                     + m_sudz[6][inode] + m_sudz[7][inode] + m_sudz[8][inode] + m_sudz[13][inode]
                                     + m_sudz[14][inode] + m_sudz[19][inode] + m_sudz[20][inode] + m_sudz[22][inode]
                                     + m_sudz[23][inode] + m_sudz[25][inode] + m_sudz[26][inode]; // is it useful ?
                temp_sudint = m_sudint[1][inode] + m_sudint[2][inode] + m_sudint[4][inode] + m_sudint[5][inode]
                            + m_sudint[6][inode] + m_sudint[7][inode] + m_sudint[8][inode] + m_sudint[13][inode]
                            + m_sudint[14][inode] + m_sudint[19][inode] + m_sudint[20][inode] + m_sudint[22][inode]
                            + m_sudint[23][inode] + m_sudint[25][inode] + m_sudint[26][inode];
                m_sudint[itype][inode] = temp_sudint;
                if(temp_sudint < m_sudmin[itype]) m_sudmin[itype] = temp_sudint;
                if(temp_sudint > m_sudmax[itype]) m_sudmax[itype] = temp_sudint;
            }
            // down type quark initial states
            else if(itype == 30){
                m_sudz[itype][inode] = m_sudz[1][inode] + m_sudz[2][inode] + m_sudz[4][inode] + m_sudz[5][inode]
                                     + m_sudz[6][inode] + m_sudz[10][inode] + m_sudz[11][inode] + m_sudz[16][inode]
                                     + m_sudz[17][inode] + m_sudz[19][inode] + m_sudz[20][inode] + m_sudz[22][inode]
                                     + m_sudz[23][inode] + m_sudz[25][inode] + m_sudz[26][inode]; // is it useful ?
                temp_sudint = m_sudint[1][inode] + m_sudint[2][inode] + m_sudint[4][inode] + m_sudint[5][inode]
                            + m_sudint[6][inode] + m_sudint[10][inode] + m_sudint[11][inode] + m_sudint[16][inode]
                            + m_sudint[17][inode] + m_sudint[19][inode] + m_sudint[20][inode] + m_sudint[22][inode]
                            + m_sudint[23][inode] + m_sudint[25][inode] + m_sudint[26][inode];
                m_sudint[itype][inode] = temp_sudint;
                if(temp_sudint < m_sudmin[itype]) m_sudmin[itype] = temp_sudint;
                if(temp_sudint > m_sudmax[itype]) m_sudmax[itype] = temp_sudint;
            }
            // photon initial states
            else if(itype == 31){
                m_sudz[itype][inode] = m_sudz[9][inode] + m_sudz[12][inode]; // is it useful ?
                temp_sudint = m_sudint[9][inode] + m_sudint[12][inode];
                m_sudint[itype][inode] = temp_sudint;
                if(temp_sudint < m_sudmin[itype]) m_sudmin[itype] = temp_sudint;
                if(temp_sudint > m_sudmax[itype]) m_sudmax[itype] = temp_sudint;
            }
            // Z boson initial states
            else if(itype == 32){
                m_sudz[itype][inode] = m_sudz[15][inode] + m_sudz[18][inode]; // is it useful ?
                temp_sudint = m_sudint[15][inode] + m_sudint[18][inode];
                m_sudint[itype][inode] = temp_sudint;
                if(temp_sudint < m_sudmin[itype]) m_sudmin[itype] = temp_sudint;
                if(temp_sudint > m_sudmax[itype]) m_sudmax[itype] = temp_sudint;
            }

            if(itype == 0) m_scales[inode] = q2max; // we only note the scale once
        }   
    }
}

// z times the coupling times the splitting function
double SudakovTables::splitting_integrand(double z, double q, double aijnorm, int branchtype, AlphaStrong& alphas)
/*
    z          : the splitting variable
    q          : the scale
    aijnorm    : the normalisation for the near quark mass threshold term (expected 10/log(q2max/q2min))
    branchtype : the type of splitting, see below
    alphas     : instance of AlphaStrong passed by reference
*/
{
    const double tolerance = 0.005; // tolerance for neartresh 
    double integrand = 0.;  
    double psplit = 0.;
    bool neartresh = false; // are we close to a mass threshold ?
    bool iabs = true;
    int nflav = 3;
    double aijnorm = 1.;
    double as1;

    if(q >= m_charmmass){
        nflav = 4;
        if(util::lmb_eq(q, m_charmmass, -tolerance)) neartresh = true; // relative tolerance, slightly different from the original
    }
    if(q >= m_bottommass){
        nflav = 5;
        if(util::lmb_eq(q, m_bottommass, -tolerance)) neartresh = true;
    }
    if(q >= m_topmass){
        nflav = 6;
        if(util::lmb_eq(q, m_topmass, -tolerance)) neartresh = true;
    }

    // computation of a_s
    switch (m_qordering)
    {
        // angular ordering, dynamical 
        case 0:   
            as1 = alphas.running(q*q*(1.-z)*(1.-z));
            break;
        
        // virtuality ordering, dynamical !!! error in the original
        case 1:   
            as1 = alphas.running(q*q*(1.-z));
            break;

        // angular ordering, fixed in alpha_s
        case 0:   
            as1 = alphas.running(q*q);
            break;
    }

    // computation of the splitting function
    switch (branchtype)
    {
        // g -> gg
        case 0:
            // LO
            psplit = lo::dqcP0GGS(z) * lo::dqcP0GGR(z) + lo::dqcP0GGA(z);
            // NLO
            if(m_order > 1){
                psplit += as1 * nlo::GG1SFUNC(z, nflav);
            }
            // NNLO
            if(m_order > 2){
                psplit += as1 * as1 * (nnlo::P2GGA(z, nflav) + nnlo::P2GGB(z, nflav)) / 8.;
                // near mass threshold
                if(neartresh){
                    psplit += aijnorm * as1 * as1 * (nnlo::A2GG(z) + nnlo::SOFTG(z)) / 4.;
                    if(iabs) psplit = abs(psplit);
                }
            }
            break;

        // q -> gq
        case 1:
            // LO
            psplit = lo::dqcP0GFA(z);
            // NLO
            if(m_order > 1){
                psplit += as1 * nlo::FG1SFUNC(z, nflav);
            }
            // NNLO
            if(m_order > 2){
                psplit += as1 * as1 * nnlo::P2GQA(z, nflav) / 8.;
                // near mass threshold
                if(neartresh){
                    psplit += aijnorm * as1 * as1 * nnlo::A2GQ(z) / 4.;
                    if(iabs) psplit = abs(psplit);
                }
            }
            break;
    
        // q_i -> q_i g
        case 2:
            // LO
            psplit = lo::dqcP0FFR(z) * lo::dqcP0FFS(z);
            // NLO
            if(m_order > 1){
                psplit += as1 * ((nlo::PP1SFUNC(z, nflav) + nlo::PM1SFUNC(z, nflav))/2. + 
                                 (nlo::FF1SFUNC(z, nflav) - nlo::PP1SFUNC(z, nflav))/(2.*nflav));
            }
            // NNLO
            if(m_order > 2){
                psplit += as1 * as1 * (nnlo::P2NSB(z, nflav) + (nnlo::P2NSPA(z, nflav) + nnlo::P2NSMA(z, nflav))/2. 
                                       + (nnlo::P2PSA(z, nflav) + nnlo::P2NSSA(z,nflav))/(2.*nflav)) / 8.;
                // near mass threshold
                if(neartresh){
                    psplit += aijnorm * as1 * as1 * (nnlo::A2QQNS(z) + nnlo::SOFTQ2(z)) / 4.;
                    if(iabs) psplit = abs(psplit);
                }
            }
            break;

        // g -> qq
        case 3:
            // LO
            psplit = 2 * nflav * lo::dqcP0FGA(z);
            // NLO
            if(m_order > 1){
                psplit += as1 * nlo::GF1SFUNC(z, nflav);
            }
            // NNLO
            if(m_order > 2){
                psplit += as1 * as1 * nnlo::P2GQA(z, nflav) / 8.;
                // near mass threshold
                if(neartresh){
                    psplit += aijnorm * as1 * as1 * nnlo::A2QG(z) / 4.;
                    if(iabs) psplit = abs(psplit);
                }
            }
            break;

        // q_i -> qbar_i g
        case 4:
            // LO
            psplit = 0.;
            // NLO
            if(m_order > 1){
                psplit += as1 * ((nlo::PP1SFUNC(z, 1) - nlo::PM1SFUNC(z, 1))/2. // nflav = 1 for those ?
                                 + (nlo::FF1SFUNC(z, nflav) - nlo::PP1SFUNC(z, nflav))/(2.*nflav));  
            }
            // NNLO
            if(m_order > 2){
                psplit += as1 * as1 * ((nnlo::P2NSPA(z, nflav) - nnlo::P2NSMA(z, nflav))/2. 
                                       + (nnlo::P2PSA(z, nflav) - nnlo::P2NSSA(z, nflav))/(2.*nflav)) / 8.;
                // no near mass treshold ?   
            }
            break;

        // q_i -> q_j g
        case 5:
            // LO
            psplit = 0.;
            // NLO
            if(m_order > 1){
                psplit += as1 * (nflav-1) * (nlo::FF1SFUNC(z, nflav) - nlo::PP1SFUNC(z, nflav)) / (2.*nflav);
            }
            // NNLO
            if(m_order > 2){
                psplit += as1 * as1 * (nflav-1) * (nnlo::P2PSA(z, nflav) + nnlo::P2NSSA(z, nflav)) / (2.*8.*nflav);
                // near mass threshold
                if(neartresh){
                    psplit += aijnorm * as1 * as1 * nnlo::A2QQ(z) / (4.*2.); // divide by 2 to separate q and qbar
                    if(iabs) psplit = abs(psplit);
                }
            }
            break;

        // q_i -> qbar_j g
        case 6:
            // LO
            psplit = 0.;
            // NLO
            if(m_order > 1){
                psplit += as1 * (nflav-1) * (nlo::FF1SFUNC(z, nflav) - nlo::PP1SFUNC(z, nflav)) / (2.*nflav);
            }
            // NNLO
            if(m_order > 2){
                psplit += as1 * as1 * (nflav-1) * (nnlo::P2PSA(z, nflav) - nnlo::P2NSSA(z, nflav)) / (2.*8.*nflav);
                // near mass threshold
                if(neartresh){
                    psplit += aijnorm * as1 * as1 * nnlo::A2QQ(z) / (4.*2.);
                    if(iabs) psplit = abs(psplit);
                }
            }
            break;

        // q_i -> boson + q
        // 7/10: photon & u/d ; 13/16: Z & u/d ; 19: W+ ; 22: W- ; 25: H
        case 7: case 10: case 13: case 16: case 19: case 22: case 25:
            if(m_smbosons) psplit = lo::dqedP0GFA(z);
            break;

        // q_i -> q + boson
        // 8/11: photon & u/d ; 14/17: Z & u/d ; 20: W+ ; 23: W- ; 26: H
        case 8: case 11: case 14: case 17: case 20: case 23: case 26:
            if(m_smbosons) psplit = lo::dqedP0FFR(z) * lo::dqcP0FFS(z);
            break;

        // boson -> q q
        // 9/12: photon & u/d ; 15/18: Z & u/d ; 21: W+ ; 24: W- ; 27: H
        case 9: case 12: case 15: case 18: case 21: case 24: case 27:
            if(m_smbosons) psplit = lo::dqedP0FGA(z);
            break;
    }

    // Ievol = 1 in the original
    psplit = psplit * z;

    // putting the correct prefactors 
    if(branchtype <= 6){
        integrand = as1 * psplit;
    }
    // photon & up type quark
    else if((branchtype >= 7) && (branchtype <= 9)){
        integrand = aqed * e2u * psplit / (2*pi); // is it correct ?
    }
    // photon & down type quark
    else if((branchtype >= 10) && (branchtype <= 12)){
        integrand = aqed * e2d * psplit / (2*pi); // is it correct ?
    }
    // Z & up type quark
    else if((branchtype >= 12) && (branchtype <= 15)){
        integrand = (aqed / (4.*sin2w*cos2w)) * (psplit / (2*pi)) * ((t3u - 2.*sin2w)*(t3u - 2.*sin2w) + t3u*t3u);
    }
    // Z & down type quark
    else if((branchtype >= 16) && (branchtype <= 18)){
        integrand = (aqed / (4.*sin2w*cos2w)) * (psplit / (2*pi)) * ((t3d - 2.*sin2w)*(t3d - 2.*sin2w) + t3d*t3d);
    }
    // W+/-
    else if((branchtype >= 19) && (branchtype <= 24)){
        integrand = (aqed / (4.*sin2w)) * (psplit / (2*pi));
    }
    // H
    else if((branchtype >= 25) && (branchtype <= 27)){
        integrand = aqed * psplit / (2*pi);
    }

    return integrand;
}

// wrapper to sample splitting_indegrand following a 1/(1-z) distribution, for the numerical integration
double SudakovTables::sample_integrand_z(double y,double q, double aijnorm, int branchtype, AlphaStrong& alphas)
/*
    y          : the sampling variable to generate z
### parameters for splitting integrand ###
    q          : the scale
    aijnorm    : the normalisation for the near quark mass threshold term (expected 10/log(q2max/q2min))
    branchtype : the type of splitting, see below
    alphas     : instance of AlphaStrong passed by reference
*/
{
    const double epsilon = 1.E-15;
    double z1min = 0.1;
    double z1max, z1, z, z1weight;
    // TBA : the possibility to specify a zmin

    z1max = 1. - epsilon;
    if(m_qordering == 0) z1min = max(epsilon, (m_q0ord/sqrt(q2max)));
    else if(m_qordering == 1) z1min = max(epsilon, (m_q0ord*m_q0ord / q2max));
    else if(m_qordering == 2) z1min = 1. - m_zmaxfixed;

    z1 = z1min * pow((z1max / z1min), y);
    z1weight = z1 * log(z1max / z1min);

    z = 1. - z1;

    return z1weight * splitting_integrand(z, q, aijnorm, branchtype, alphas);
}

// to integrate a splitting function along z, at set q, using sample_integrand_z
double integrate_splitting(double q2max, double q2min, int branchtype, double eps, AlphaStrong& alphas)
/*
    q2max      : upper bound of scale squared (in the grid), also the scale at which to evaluate
    q2min      : lower bound of scale squared (in the grid), needed to calculate aijnorm
    branchtype : the type of branching
    eps        : the desired precision for the integration
    alphas     : reference to an instance of AlphaStrong (the running coupling)
*/
{
    double q;

    q = sqrt(q2max);
    aijnorm = 0;
    if(q2max != q2min) aijnorm = 10./log(q2max/q2min); 

    // lambda expression to keep only the y dependance, for the integration
    double splitting_lmbd = [q, aijnorm, branchtype, alphas](double y){
        return sample_integrand_z(y, q, aijnorm, branchtype, alphas);
    };

    return util::dgadap(0., 1., &splitting_lmbd, eps);
}