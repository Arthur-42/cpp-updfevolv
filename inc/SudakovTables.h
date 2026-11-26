#ifndef SUDAKOVTABLES_H
#define SUDAKOVTABLES_H

// number of grid points
#define tablelength 64

class SudakovTables
{
    private:
    /* data */
    // order of computations
    int m_order;
    // scale ordering 
    int m_qordering;
    // do we include photon, Z, W+/-, H branchings ?
    bool m_smbosons;
    // heavy quarks and bosons masses
    double m_charmmass;
    double m_bottommass;
    double m_topmass;
    double m_Wmass;
    double m_Zmass;
    double m_Hmass;
    // reference scale for dynamical z_max
    double m_q0ord;
    // z_max if fixed
    double m_zmaxfixed;

    // list of the scales at which the grid is evaluated
    double m_scales[tablelength];
    // splitting functions integrated along z
    double m_sudz[33][tablelength];
    // full integral in the sudakov form factors
    double m_sudint[33][tablelength];
    // lists of min and max values per type
    double m_sudmin[33];
    double m_sudmax[33];

    /* methods */
    // to fill the sudakov tables
    void filltables(double qmin, double qmax, double eps, AlphaStrong& alphas);
    //
    double splitting_integrand(double z, double q, double aijnorm, int branchtype, AlphaStrong& alphas);
    // wrapper to sample splitting_integrand along z following a 1/(1-x) distribution, for the numerical integration
    double sample_integrand_z(double y,double q, double aijnorm, int branchtype, AlphaStrong& alphas);
    // to integrate a splitting function along z, at set q, using sample_integrand_z
    double integrate_splitting(double q2max, double q2min, int branchtype, double eps, AlphaStrong& alphas);

    public:
};

#endif