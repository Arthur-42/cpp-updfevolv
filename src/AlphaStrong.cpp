#include "AlphaStrong.h"
#include "utilities.h"
#include <cmath>
#include <algorithm>


using std::abs, std::max, std::min, std::log;
namespace util = utilities;

static const double pi = 3.14159265358979;

// TBA : error handling if out of perturbative regime

/* Adapted from QCDNUM 18-00-00 srcAlfa.f */

// constructor
AlphaStrong::AlphaStrong(int order, double mu20, double as0, double cmass, double bmass, double tmass)
/*
    order : order of all computations 1 = L0, 2 = NLO, 3 = NNLO
    mu20  : reference scale (squared)
    as20  : referance value of as at mu20
    cmass : charm quark mass
    bmass : bottom quark mass
    tmass : top quark mass
*/
{
    // initialise member variables
    m_order = order;
    m_mu20 = mu20;
    m_as0 = as0;
    m_quarkmasses2 = {cmass*cmass, bmass*bmass, tmass*tmass};
    m_tol = 1.E-9 // default QCDNUM value
    m_regions = 0;
    m_starts = {0., 0., 0., 0.};
    m_stops = {0., 0., 0., 0.};
    m_flavs = {0, 0, 0, 0};
    m_direct = 0;
    m_temp_asplus = 0.1;
    m_temp_mu2 = 10.;
    m_coeffs = {
        {0., 2./3., 0.},
        {14./3., 38./3., 4./9.}
    };
    m_errorflag = false;
}

// to get a_s / 2pi at the desired scale
double AlphaStrong::running(double mu2)
/*
    mu2 : the desired scale (squared)
*/
{
    double asend, rend, rstart;
    bool ihit = false;
    int nfdum, nfjump;

    // get evolution limits
    get_limits(m_mu20, mu2);

    // nothing to evolve 
    if(m_direct == 0){
        return m_as0;
    }

    // upward evolution
    if(m_direct == 1){
        asend = m_as0;
        for(int i = 0; i < m_regions; i++){
            rend = m_stops[i];
            asend = integrate_beta(m_starts[i], rend, asend, m_flavs[i]);
            // is the endpoint at a flavour threshold ?
            nfdum = get_nflav(rend, ihit);
            // add discontinuity if it's the case
            if(ihit){
                nfjump = 1;
                asend = do_jump(asend, rend, nfjump);
            }
        }
    }
    // downward evolution
    else{
        asend = m_as0;
        for(int j = 0; j < m_regions; j++){
            rstart = m_starts[i];
            nfdum = get_nflav(rstart, ihit);
            // substract discontinuity if start point at flavour threshold
            if(ihit){
                nfjump = -1;
                asend = do_jump(asend, rstart, nfjump);
            }
            asend = integrate_beta(rstart, m_stops[i], asend, m_flavs[i]);
        }
    }

    // unlike the original, we return a_s / 2pi
    return asend/(2.*pi);
}

// to indicate which flavour thresholds are in the way
void AlphaStrong::get_limits(double mu2ini, double mu2fin)
/*
    mu2ini : intial scale (squared)
    mu2fin : final scale (squared)
    -> fills m_starts, m_stops, m_flavs and m_regions
*/
{
    bool ihit1 = false;
    bool ihit2 = false;
    int nf1, nf2, n, m;

    // no evolution required
    if(util::lmb_eq(mu2ini, mu2fin, -m_tol)){
        m_direct = 0; 
        return;
    }

    // specify the evolution direction
    if(util::lmb_lt(mu2ini, mu2fin, -m_tol)) m_direct = 1;
    else m_direct = -1;

    // number of flavours at the start and at the end
    nf1 = get_nflav(mu2ini, ihit1);
    nf2 = get_nflav(mu2fin, ihit2);

    if(m_direct == 1){
    // upward evolution
    // if the end point is at a flavour threshold, nf2 -> nf2-1
        if(ihit2) nf2 = nf2-1;
        n = 0;
        for(int i = nf1; i <= nf2; i++){
            if(i == 3){
                m_starts[n] = mu2ini;
                m_stops[n] = min(mu2fin, m_quarkmasses2[0]);
                m_flavs[n] = 3;
            }
            else if(i == 4 || i == 5){
                m_starts[n] = min(max(mu2ini, m_quarkmasses2[i-4]), m_quarkmasses2[i-3]); 
                m_stops[n] = min(max(mu2fin, m_quarkmasses2[i-4]), m_quarkmasses2[i-3]);
                m_flavs[n] = i;
            }
            else if(i == 6){
                m_starts[n] = max(mu2ini, m_quarkmasses2[2]);
                m_stops[n] = mu2fin;
                m_flavs[n] = 6;
            }
            n += 1;
        }
        m_regions = n;
    }
    else{
    // downward evolution
    // is the start poin is at a flavour threshold, nf1 -> nf1-1
        if(ihit1) nf1 = nf1-1;
        m = 0;
        for(int i = nf1; i >= nf2; i--){ 
            if(i == 6){
                m_starts[m] = mu2ini;
                m_stops[m] = max(mu2fin, m_quarkmasses2[6]);
                m_flavs[m] = 3;
            }
            else if(i == 4 || i == 5){
                m_starts[m] = min(max(mu2ini, m_quarkmasses2[i-4]), m_quarkmasses2[i-3]); 
                m_stops[m] = min(max(mu2fin, m_quarkmasses2[i-4]), m_quarkmasses2[i-3]);
                m_flavs[m] = i;
            }
            else if(i == 3){
                m_starts[m] = min(mu2ini, m_quarkmasses2[0]);
                m_stops[m] = mu2fin;
                m_flavs[m] = 3;
            }
            m += 1;
        }
        m_regions = m;
    }
    return;
}

// to integrate the beta function via Runge Kutta 4 method
double AlphaStrong::integrate_beta(double mu2ini, double mu2fin, double asini, int nflav)
/*
    mu2ini : initial scale (squared)
    mu2fin : final scale (squared)
    asini  : value of a_s at the initial scale
    nflav  : number of active flavours
*/
{
    const int nsteps = 10;
    const double dadrlim = 1.E4;
    const double C2S3 = 2./3.;
    const double C38S3 = 38./3.;
    const double C5033S18 = 5033./18.;
    const double C325S54 = 325./54.;
    double dlq, as, xk0, xk1, xk2, xk3, bet;

    // beta function in powers of a = as/(4*pi)
    inline double fbeta(a){
        double value = 0;
        if(m_order > 0) value += - (11. - C2S3*nflav);
        if(m_order > 1) value += - a * (102. - C38S3.*nflav);
        if(m_order > 2) value += - (a*a) * (1428.5 - C5033S18*nflav + C325S54*nflav*nflav);
        return (a*a) * value;
    }

    // RK4 method
    dlq = log(mu2fin/mu2ini)/nsteps;
    as = asini/(4*pi);
    for(int k1=1; k1 <= nsteps; k1++){
        bet = fbeta(as);
        if(bet < -dardlim) return 0.; // error
        xk0 = dlq * bet;
        xk1 = dlq * fbeta(as+0.5*xk0);
        xk2 = dlq * fbeta(as+0.5*xk1);
        xk3 = dlq * fbeta(as+xk2);
        as = as + 1./6. * (xk0 + 2.*xk1 + 2.xk2 + xk3);
    }

    return as*4.*pi;
}

// to get the number of flavours active at the current scale
int AlphaStrong::get_nflav(double mu2, bool& ihit)
/*
    mu2 : current scale (squared)
    ihit : (passed by reference) are we close to a treshold ?
*/
{
    int nf = 3;

    ihit = false;
    for(int i=0; <3, i++){
        if(mu2 > m_quarkmasses2[i]) nf = 3+i;
        if(util::lmb_eq(mu2, m_quarkmasses2[i], -m_tol)){
            nf = 3+i;
            ihit = true;
        }
    }

    return nf;
}

// to recalculate a_s at a flavour treshold
double AlphaStrong::do_jump(double asold, double mu2, int nfjump)
/*
    asold  : initial value of a_s
    mu2    : current scale (squared)
    nfjump : magnitude of the flavour jump (+1 or -1)
*/
{
    double asnew, minint, maxint;

    // no jump at LO
    if(m_order == 1) return asold;

    // jump for upward evolution from nf to nf+1
    if(nfjump == 1) asnew = asold + get_jump(asold, mu2);
    // jump for downward evolution from nf to nf-1
    else if(nfjump == -1){
        // we have to solve as(nf) = as(nf-1) + asjump[as(nf-1)]
        // thus solve_jump(asnew) = asold - asnew - get_jump(asnew) = 0
        m_temp_asplus = asold;
        m_temp_mu2 = mu2;
        // initial guess of the interval of the root
        minint = 0.95 * asold;
        maxint = 1.05 * asold;
        // find the correct bracket if needed
        util::brackit(&solve_jump, minint, maxint);
        // find root within tolerance 10^-6
        asnew = util::bisect(&solve_jump, minint, maxint, 1.E-6);
    }

    return asnew;
}

// to compute the jump in a_s (at a flavour threshold)
double AlphaStrong::get_jump(double asold, double mu2)
/*
    asold : initial value of a_s
    mu2   : current scale (squared)
*/
{
    double asini, ajump, asfin, cterm, rterm;

    // no jump at LO
    if(m_order <= 1) return 0.;

    // formulas for as/4pi
    asini = asold/(4*pi);
    ajump = 0.;
    asfin = asini;
    for(int n=0; n < (m_order-1); n++){
        cterm = 0.
        for(int j=0; j <= n; j++){
            cterm = cterm + m_coeffs[n][j]*rterm;
            rterm = 0.;
        }
        asfin = asfin*asini;
        ajump = ajump + asfin*cterm;
    }

    return ajump*4.*pi;
}

// function to be solved to get downward jump in flavour
double AlphaStrong::solve_jump(double asmin)
/*
    asmin : a_s before the flavour threshold
*/
{
    return m_temp_asplus - asmin - get_jump(asmin, m_temp_mu2);
}