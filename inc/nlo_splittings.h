#ifndef NLO_SPLITTINGS_H
#define NLO_SPLITTINGS_H

/*
g->gg      : GG1SFUNC
q->gq      : FG1SFUNC
q_i->q_ig  : PP1SFUNC, PM1SFUNC, FF1SFUNC
g->qq      : GF1SFUNC
q_i->aq_ig : done
q_i->q_jg  : done
q_i->aq_jg : done
*/

namespace nlo_splittings
{
    //
    double GG1SFUNC(double z, int nflav);
    //
    double FG1SFUNC(double z, int nflav);
    //
    double PP1SFUNC(double z, int nflav);
    //
    double PM1SFUNC(double z, int nflav);
    //
    double FF1SFUNC(double z, int nflav);
    //
    double GF1SFUNC(double z, int nflav); 
}

#endif