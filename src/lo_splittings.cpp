#include "lo_splittings.h"

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
void lo_splittings::lo_ignore(const int&){} // does nothing, shouldn't even get compiled

// 
double lo_splittings::dqcP0GGS(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 1. / (1. - z);
}

// 
double lo_splittings::dqcP0GGR(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 6. * z;
}

// 
double lo_splittings::dqcP0GGA(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 6. * ((1. - z) / z + z * (1. - z));
}

// 
double lo_splittings::dqcP0GFA(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 4. * (1. + (1. - z)*(1. - z)) / (3. * z);
}

// 
double lo_splittings::dqcP0FFR(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 4. * (1. + z*z) / 3.;
}

// 
double lo_splittings::dqcP0FFS(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 1. / (1. - z);
}

// 
double lo_splittings::dqcP0FGA(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 0.5 * (z*z + (1. - z)*(1. - z));
}

// 
double lo_splittings::dqedP0GFA(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return (1. + (1. - z)*(1. - z)) / z;
}

// 
double lo_splittings::dqedP0FFR(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return (1. + z*z);
}

// 
double lo_splittings::dqedP0FGA(double z, int nflav)
{
    lo_ignore(nflav); // avoid compiler warning for unused nflav
    return 3. * (z*z + (1. - z)*(1. - z));
}
