/* Parameters.hpp - A class to contain the model parameters and precision tuning
 *   (originally in parameters.c and precision.c)
 *
 * Copyright (C) 2006  Jochen Voss, Andreas Voss.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301 USA.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

// Note: Parameters class now includes precision constants

// Indices for packed parameters array
#define PARAM_a   0
#define PARAM_v   1
#define PARAM_t0  2
#define PARAM_d   3
#define PARAM_szr 4
#define PARAM_sv  5
#define PARAM_st0 6
#define PARAM_zr  7

class Parameters
{
public:
    double a;     // Boundary separation
    double v;     // Mean of the drift
    double t0;    // Non-decision time
    double d;     // Difference between boundaries of non-decision time
    double szr;   // width of zr distribution
    double sv;    // standard deviation of v distribution
    double st0;   // width of t0 distribution
    double zr;    // Mean of diffusion starting point relative to boundaries


    // Precision constants set by SetPrecision()
    double  TUNE_DZ;
    double  TUNE_DV;
    double  TUNE_DT0;

    double  TUNE_PDE_DT_MIN;   // If std=c++11 we can use C++ defaults to set as = 1e-6;
    double  TUNE_PDE_DT_MAX;   // ... we can default to = 1e-6;
    double  TUNE_PDE_DT_SCALE; // ... we can default to = 0.0;

    double  TUNE_INT_T0;
    double  TUNE_INT_Z;

    double  TUNE_SV_EPSILON; // CONVERSION NOTE: See below in SetPrecision()
    double  TUNE_SZ_EPSILON; // CONVERSION NOTE: See below in SetPrecision()

public:
  // Construct the object from the passed in params
    Parameters (NumericVector params, double precision)
    {
        a   = params[PARAM_a];
        v   = params[PARAM_v];
        t0  = params[PARAM_t0];
        d   = params[PARAM_d];
        szr = params[PARAM_szr];
        sv  = params[PARAM_sv];
        st0 = params[PARAM_st0];
        zr  = params[PARAM_zr];

        SetPrecision (precision);
    }

    bool ValidateParams (bool print)
    {
        bool valid = true;
        
        if (a <= 0)                         { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter a = " << a << std::endl;  }
        if (szr < 0 || szr > 1)             { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter szr = " << szr << std::endl; }
        if (st0 < 0)                        { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
        if (sv < 0)                         { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter sv = " << sv << std::endl; }
        if (t0 - fabs(0.5*d) - 0.5*st0 < 0) { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination t0 = " << t0 << ", d = " << d << ", st0 =" << st0 << std::endl; }
        if (zr - 0.5*szr <= 0)              { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
        if (zr + 0.5*szr >= 1)              { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}

        return valid;
    }

private:
    void  SetPrecision (double p)
    {
        // Try to achieve an accuracy of approximately 10^{-p} for the CDF.
        TUNE_PDE_DT_MIN   = pow(10, -0.400825*p-1.422813);
        TUNE_PDE_DT_MAX   = pow(10, -0.627224*p+0.492689);
        TUNE_PDE_DT_SCALE = pow(10, -1.012677*p+2.261668);
        TUNE_DZ  = pow(10, -0.5*p-0.033403);
        TUNE_DV  = pow(10, -1.0*p+1.4);
        TUNE_DT0 = pow(10, -0.5*p-0.323859);

        TUNE_INT_T0 = 0.089045 * exp(-1.037580*p);
        TUNE_INT_Z  = 0.508061 * exp(-1.022373*p);

        // CONVERSION NOTE:
        //     These have been added to optimise code paths by treating very small variances as 0
        //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
        TUNE_SV_EPSILON = pow (10, -(p+2.0));
        TUNE_SZ_EPSILON = pow (10, -(p+2.0));
    }
};

#endif // PARAMETERS_H
