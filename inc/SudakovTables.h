#ifndef SUDAKOVTABLES_H
#define SUDAKOVTABLES_H
#include <array>

// number of grid points
#define tablelength 64

// alias for the 2d arrays
template <typename T, std::size_t Row, std::size_t Col>
using array2d = std::array<std::array<T, Col>, Row>;

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
    std::array<double, tablelength> m_scales;
    // splitting functions integrated along z
    array2d<double, 33, tablelength> m_sudz;
    // full integral in the sudakov form factors
    array2d<double, 33, tablelength> m_sudint;
    // lists of min and max values per type
    std::array<double, 33> m_sudmin;
    std::array<double, 33> m_sudmax;

    /* methods */
    // to fill the sudakov tables
    void filltables(double qmin, double qmax, double eps, AlphaStrong& alphas);
    // z times the coupling times the splitting function
    double splitting_integrand(double z, double q, double aijnorm, int branchtype, AlphaStrong& alphas);
    // wrapper to sample splitting_integrand along z following a 1/(1-x) distribution, for the numerical integration
    double sample_integrand_z(double y,double q, double aijnorm, int branchtype, AlphaStrong& alphas);
    // to integrate a splitting function along z, at set q, using sample_integrand_z
    double integrate_splitting(double q2max, double q2min, int branchtype, double eps, AlphaStrong& alphas);

    public:
};

#endif