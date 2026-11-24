#include "lo_splittings.h"

/*
from dqcP0GGS to dqcP0FGA : taken from QCDNUM 18-00-00, file pij_lo.f, 
"written by Ouraou and Virchaud, 1988"
from dqedP0GFA to dqed P0FGA : taken from uPDFevolv-2.5.03, file Photon_Splitting.f
*/

/*
g->gg      : dqcP0GGS, dqcP0GGR, dqcP0GGA
q->gq      : dqcP0GFA
q_i->q_ig  : dqcP0FFR, dqcP0FFS
g->qq      : dqcP0FGA
q_i->aq_ig : /
q_i->q_jg  : /
q_i->aq_jg : /
q->bq      : dqedPOGFA
q->qb      : dqedP0FFR, dqcP0FFS
b->qq      : dqedP0FGA
*/

// utility : to avoid compiler warning for unused int variable
//void lo_splittings::lo_ignore(const int&){} // does nothing, shouldn't even get compiled

// 
double lo_splittings::dqcP0GGS(double z)
{
    //lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 1. / (1. - z);
}

// 
double lo_splittings::dqcP0GGR(double z)
{
    return 6. * z;
}

// 
double lo_splittings::dqcP0GGA(double z)
{
    return 6. * ((1. - z) / z + z * (1. - z));
}

// 
double lo_splittings::dqcP0GFA(double z)
{
    return 4. * (1. + (1. - z)*(1. - z)) / (3. * z);
}

// 
double lo_splittings::dqcP0FFR(double z)
{
    return 4. * (1. + z*z) / 3.;
}

// 
double lo_splittings::dqcP0FFS(double z)
{
    return 1. / (1. - z);
}

// 
double lo_splittings::dqcP0FGA(double z)
{
    return 0.5 * (z*z + (1. - z)*(1. - z));
}

// 
double lo_splittings::dqedP0GFA(double z)
{
    return (1. + (1. - z)*(1. - z)) / z;
}

// 
double lo_splittings::dqedP0FFR(double z)
{
    return (1. + z*z);
}

// 
double lo_splittings::dqedP0FGA(double z)
{
    return 3. * (z*z + (1. - z)*(1. - z));
}
