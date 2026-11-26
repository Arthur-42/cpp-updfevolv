#ifndef NNLO_SPLITTINGS_H
#define NNLO_SPLITTINGS_H

/*
g->gg      : P2GGA, P2GGB, A2GG, SOFTG
q->gq      : P2GQA, A2GQ 
q_i->q_ig  : P2NSB, P2NSPA, P2NSMA, P2PSA, P2NSSA, A2QQNS, SOFTQ2
g->qq      : P2QGA, A2QG
q_i->aq_ig : done
q_i->q_jg  : A2QQ
q_i->aq_jg : done
*/

namespace nnlo_splittings
{
    //
    double P2GGA(double z, int nflav);
    //
    double P2GGB(double z, int nflav);
    //
    double A2GG(double z);
    //
    double SOFTG(double z);
    //
    double P2GQA(double z, int nflav);
    //
    double A2GQ(double z);
    //
    double P2NSB(double z, int nflav);
    //
    double P2NSPA(double z, int nflav);
    //
    double P2NSMA(double z, int nflav);
    //
    double P2PSA(double z, int nflav);
    //
    double P2NSSA(double z, int nflav);
    //
    double A2QQNS(double z);
    //
    double SOFTQ2(double z);
    //
    double P2QGA(double z, int nflav);
    //
    double A2QG(double z);
    //
    double A2QQ(double z);
}

#endif