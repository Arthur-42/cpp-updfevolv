#ifndef SUDAKOVTABLES_H
#define SUDAKOVTABLES_H

class SudakovTables
{
    private:
    /* data */
    // scale ordering 
    int m_qordering;
    // heavy quarks masses squared
    double m_charmmass2;
    double m_bottommass2;
    double m_topmass2;

    /* methods */
    // typically a_s P(z)
    double splitting_integrand(double z,double q2min, double q2max, int branchtype, AlphaStrong& as);

    public:
};

#endif