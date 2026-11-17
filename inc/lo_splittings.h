#ifndef LO_SPLITTINGS_H
#define LO_SPLITTINGS_H

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

namespace lo_splittings
{
    // 
    double dqcP0GGS(double z, int nflav);
    // 
    double dqcP0GGR(double z, int nflav);
    // 
    double dqcP0GGA(double z, int nflav);
    // 
    double dqcP0GFA(double z, int nflav);
    // 
    double dqcP0FFR(double z, int nflav);
    // 
    double dqcP0FFS(double z, int nflav);
    // 
    double dqcP0FGA(double z, int nflav);
    // 
    double dqedP0GFA(double z, int nflav);
    // 
    double dqedP0FFR(double z, int nflav);
    // 
    double dqedP0FGA(double z, int nflav);

    // utility : to avoid compiler warning
    void lo_ignore(const int&);
}

#endif