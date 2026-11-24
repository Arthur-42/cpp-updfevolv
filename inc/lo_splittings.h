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
    double dqcP0GGS(double z);
    // 
    double dqcP0GGR(double z);
    // 
    double dqcP0GGA(double z);
    // 
    double dqcP0GFA(double z);
    // 
    double dqcP0FFR(double z);
    // 
    double dqcP0FFS(double z);
    // 
    double dqcP0FGA(double z);
    // 
    double dqedP0GFA(double z);
    // 
    double dqedP0FFR(double z);
    // 
    double dqedP0FGA(double z);

    // utility : to avoid compiler warning
    //void lo_ignore(const int&);
}

#endif