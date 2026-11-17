#include "specfunc.h"
#include <cmath>

using std::abs, std::floor, std::exp;

// Adapted from GSL. Here is the original header : 
/* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Original author:  G. Jungman */

#define AGSL_LOG_DBL_MIN (-7.0839641853226408e+02)
#define AM_SQRT2 1.41421356237309504880168872421
#define AGSL_SQRT_DBL_EPSILON 1.4901161193847656e-08
#define AGSL_DBL_EPSILON 2.2204460492503131e-16
#define AGSL_LOG_DBL_EPSILON (-3.6043653389117154e+01)
#define AM_LN2      0.69314718055994530941723212146

static double adeb2_data[18] = {
   2.5943810232570770282,
   0.2863357204530719834,
  -0.102062656158046713e-01,
   0.6049109775346844e-03,
  -0.405257658950210e-04,
   0.28633826328811e-05,
  -0.2086394303065e-06,
   0.155237875826e-07,
  -0.11731280087e-08,
   0.897358589e-10,
  -0.69317614e-11,
   0.5398057e-12,
  -0.423241e-13,
   0.33378e-14,
  -0.2645e-15,
   0.211e-16,
  -0.17e-17,
   0.1e-18
};
static cheb_series adeb2_cs = {
  adeb2_data,
  17,
  -1.0, 1.0,
  10
};

static inline double cheb_eval(cheb_series * cs, double x)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += std::abs(y2*temp) + std::abs(dd) + std::abs(cs->c[j]);
    dd = temp;
  }

  { 
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += std::abs(y*temp) + std::abs(dd) + 0.5 * std::abs(cs->c[0]);
  }

  return d;
}

//
double specfunc::debye_2(double x)
/*

*/
{
    const double val_infinity = 4.80822761263837714;
    const double xcut = -AGSL_LOG_DBL_MIN;

    if(x < 2.*AM_SQRT2*AGSL_SQRT_DBL_EPSILON){
        return 1.0 - x/3.0 + x*x/24.0;
    }
    else if(x <= 4.){
        const double t = x*x/8.0 - 1.0;
        double value;
        value = cheb_eval(&adeb2_cs, t);
        return value - x/3.0;
    }

    else if(x < -(AM_LN2 + AGSL_LOG_DBL_EPSILON)){
        const int nexp = std::floor(xcut/x);
        const double ex  = std::exp(-x);
        double xk  = nexp * x;
        double rk  = nexp;
        double sum = 0.0;
        int i;
        for(i=nexp; i>=1; i--) {
            sum *= ex;
            sum += (1.0 + 2.0/xk + 2.0/(xk*xk)) / rk;
            rk -= 1.0;
            xk -= x;
        }
        return val_infinity/(x*x) - 2.0 * sum * ex;
    }
    else if(x < xcut) {
        const double x2  = x*x;
        const double sum = 2.0 + 2.0*x + x2;
        return (val_infinity - 2.0 * sum * exp(-x)) / x2;
    }
    else {
        return (val_infinity/x)/x;
    }
}