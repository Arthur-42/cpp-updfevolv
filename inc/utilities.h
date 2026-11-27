#ifndef UTILITIES_H
#define UTILITIES_H
#include <array>
#include <algorithm>

using std::min;

namespace utilities
{
    // to find the range [x1, x2] in which func(x) = 0
    // taken from QCDNUM 18-00-00
    void brackit(double (*func)(double), double& x1, double& x2);
    // to find a root of func(x) within range [x1, x2]
    // taken from QCDNUM 18-00-00
    double bisect(double (*func)(double), double x1, double x2, double epsi);
    // to determine if a is approximately equal to b, within tolerance
    // taken from QCDNUM 18-00-00
    bool lmb_eq(double a, double b, double tol);
    // to determine if a is approximately equal to b, within tolerance, and with a < b
    // taken from QCDNUM 18-00-00
    bool lmb_lt(double a, double b, double tol);
    // to perform an adaptative Gaussian integration 
    // author T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
    // reference THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
    double dgadap(double a0, double b0, double (*f)(double), double eps);

    /* defined in this header file because of the template */

    // to perform a tabular interpolation using symmetrically placed argument points
    // original from Kernlib Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
    template <typename T, std::size_t N>
    double divdifd(const std::array<T, N>& F, const std::array<T, N>& A, int NN, double X, int MM)
    /*

    */
    {
        int N, M, MPLUS, IX, IY, MID, NPTS, IP, L, ISUB, I, ind, indd, J;
        double SUM;
        bool EXTRA;
        std::array<double,20> TT, DD;

        // START.  FIND SUBSCRIPT IX OF X IN ARRAY A.
        N = NN;
        M = min(MM, min(N-1, 10));
        MPLUS = M+1;
        IX = 0;
        IY = N+1;

        if(A(1) > A(N)){
            // SEARCH DECREASING ARGUMENTS
            do
            {
                MID = (IX+IY)/2;
                if(X <= A(MID)) IX = MID;
                else IY = MID;
            } 
            while (IY-IX > 1);
        }
        else{
            // SEARCH INCREASING ARGUMENTS
            do
            {
                MID = (IX+IY)/2;
                if(X >= A(MID)) IX = MID;
                else IY = MID;
            } 
            while (IY-IX > 1);
        }

        // COPY REORDERED INTERPOLATION POINTS INTO (TT(I),DD(I)), SETTING
        // *EXTRA* TO TRUE IF M+2 POINTS TO BE USED.
        NPTS = M + 2 - M % 2;
        IP = 0;
        L = 0;

        do
        {
            ISUB = IX+L;
            if((1 <= ISUB) && (ISUB <= N)){
                IP = IP+1;
                TT(IP) = A(ISUB);
                DD(IP) = F(ISUB);
            }
            else NPTS = MPLUS;
            L = -L;
            if(L >= 0) L = L+1;
        } 
        while (IP < NPTS);

        EXTRA = (NPTS != MPLUS);

        // REPLACE DD BY THE LEADING DIAGONAL OF A DIVIDED-DIFFERENCE TABLE, 
        // SUPPLEMENTED BY AN EXTRA LINE IF *EXTRA* IS TRUE.
        for(ind = 1; ind <= M; ind++){
            if(EXTRA){
                ISUB = MPLUS - L;
                DD(M+2) = (DD(M+2) - DD(M)) / (TT(M+2) - TT(ISUB));
            }
            I = MPLUS;
            for(J = ind; J <= M; J++){
                ISUB = I - L;
                DD(I) = (DD(I) - DD(I-1)) / (TT(I) - TT(ISUB));
                I = I-1;
            }
        }

        // EVALUATE THE NEWTON INTERPOLATION FORMULA AT X, AVERAGING TWO VALUES
        // OF LAST DIFFERENCE IF *EXTRA* IS TRUE.
        SUM = DD(MPLUS);
        if(EXTRA) SUM = 0.5*(SUM + DD(M+2));
        J = M;
        for(indd = 1; indd <= M; indd++){
            SUM = DD(J) + (X - TT(J))*SUM;
            J = J-1;
        }

        return SUM;
    }
}

#endif 