#include "SudakovTables.h"
#include "lo_splittings.h"
#include "nlo_splittings.h"
#include "nnlo_splittings.h"
#include "AlphaStrong.h"

namespace lo = lo_splittings;
namespace nlo = nlo_splittings;
namespace nnlo = nnlo_splittings;

static const double pi = 3.14159265358979;
static const double aqed = 0.0072973525693; // QED coupling (at mZ)
static const double sin2w = 0.23;  // sin^2(theta_W)
static const double cos2w = 0.77;  // cos^2(theta_W)
static const double e2u = 4./9.;   // up quark QED charge squared
static const double e2d = 1./9.;   // down quark QED charge squared


// 
double SudakovTables::splitting_integrand(double z, double q2min, double q2max, int branchtype, AlphaStrong& as)
/*

*/
{
    double aijnorm;
    double integrand = 0.;
    double psplit = 0.;
    bool neartresh = false;
    int nflav = 3;
    double eq2sum = 5./18.; // default average of the sum of the quark masses squared ?

    eq2sum = 2*e2d+e2u;
    if(q2max >= m_charmmass2){
        nflav = 4;
        eq2sum = 2*e2d+2*e2u;
    }
    if(q2max >= m_bottommass2){
        nflav = 5;
        eq2sum = 3*e2d+2*e2u;
    }
    if(q2max >= m_topmass2){
        nflav = 6;
        eq2sum = 3*e2d+3*e2u;
    }

}