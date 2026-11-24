#include "SudakovTables.h"
#include "lo_splittings.h"
#include "nlo_splittings.h"
#include "nnlo_splittings.h"
#include "AlphaStrong.h"
#include "utilities.h"
#include <cmath>

using std::sqrt, std::log, std::abs;
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


// 
double SudakovTables::splitting_integrand(double z, double q2min, double q2max, int branchtype, bool iabs, AlphaStrong& alphas)
/*

*/
{
    const double tolerance = 0.005; // tolerance for neartresh 
    double integrand = 0.;  
    double psplit = 0.;
    bool neartresh = false; // are we close to a mass threshold ?
    int nflav = 3;
    double eq2sum = 5./18.; // default average of the sum of the quark masses squared ?
    double aijnorm = 1.;
    double as1;

    eq2sum = 2*e2d+e2u;
    if(q2max >= m_charmmass2){
        nflav = 4;
        eq2sum = 2*e2d+2*e2u;
        if(util::lmb_eq(sqrt(q2max), sqrt(m_charmmass2), -tolerance)) neartresh = true; // relative tolerance, slightly different from the original
    }
    if(q2max >= m_bottommass2){
        nflav = 5;
        eq2sum = 3*e2d+2*e2u;
        if(util::lmb_eq(sqrt(q2max), sqrt(m_bottommass2), -tolerance)) neartresh = true;
    }
    if(q2max >= m_topmass2){
        nflav = 6;
        eq2sum = 3*e2d+3*e2u;
        if(util::lmb_eq(sqrt(q2max), sqrt(m_topmass2), -tolerance)) neartresh = true;
    }
    aijnorm = 10./log(q2max/q2min);

    // computation of a_s
    switch (m_qordering)
    {
        // angular ordering, dynamical 
        case 0:   
            as1 = alphas.running(q2max*(1-z)*(1-z));
            break;
        
        // virtuality ordering, dynamical !!! error in the original
        case 1:   
            as1 = alphas.running(q2max*(1-z));
            break;

        // angular ordering, fixed in alpha_s
        case 0:   
            as1 = alphas.running(q2max);
            break;
    }

    // computation of the splitting function
    switch (branchtype)
    {
        // g -> gg
        case 1:
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
        case 2:
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
        case 3:
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
    }
    

}