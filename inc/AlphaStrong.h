#ifndef ALPHASTRONG_H
#define ALPHASTRONG_H

class AlphaStrong
{
    private:
        /* data */
        // order of computations
        int m_order;
        // tolerance for comparison
        double m_tol;
        // quark masses squared
        double m_quarkmasses2[3];
        // a_s(mu20) related
        double m_mu20;      // reference scale
        double m_as0;       // reference value
        // changing during operation
        int m_regions;
        double m_starts[4];
        double m_stops[4];
        int m_flavs[4];
        int m_direct;
        double m_temp_asplus; // for solve_jump
        double m_temp_mu2;    // for solve_jump
        // get_jump coefficients
        double m_coeffs[2][3];
        // did an error occur ?
        bool m_errorflag;

        /* methods */
        // to indicate which flavour thresholds are in the way
        void get_limits(double mu2ini, double mu2fin);
        // to integrate the beta function via Runge Kutta 4 method
        double integrate_beta(double mu2ini, double mu2fin, double asini, int nflav);
        // to get the number of flavours active at the current scale
        int get_nflav(double mu2, bool& ihit);
        // to recalculate a_s at a flavour treshold
        double do_jump(double asold, double mu2, int nfjump);
        // to compute the jump in a_s (at a flavour threshold)
        double get_jump(double asold, double mu2);
        // function to be solved to get downward jump in flavour
        double solve_jump(double asmin);

    public:
        /* methods */
        // constructor
        AlphaStrong(int order, double mu20 =8315.25, double as0 =0.118, double cmass =1.47, double bmass =4.5, double tmass =173.);
        // to get a_s / 2pi at the desired scale
        double running(double mu2);
};

#endif