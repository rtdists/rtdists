/* Distribution.hpp - Contains main CDF calls for fast and precise versions 
 * (precise version originally in plot-cdf.c)
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

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <Rcpp.h>
using namespace Rcpp;

//#define _CDF_DEBUG_

// An R-like version which finds the left-hand area (from 0 to RT) - uses Scott Brown's code, pass in boundary to retrieve
NumericVector distribution (NumericVector rts, int boundary)
{
    struct  F_calculator *fc;
    double  mint;
  
    const double *F;
    int     i,j,nz;
  
    double p[1]; // This is CDF(t=Inf)?

#ifdef _CDF_DEBUG_
    Rcpp::Rcout << "distribution (): Entered distribution ()\n";
    
    Rcpp::Rcout << "Using Parameters:\n";
    Rcpp::Rcout << "rts[0] = " << rts[0] << "\n";
    Rcpp::Rcout << "boundary = " << boundary << "\n";
#endif
    
    fc = F_new ();
  
    double scaled = g_Params->zr;// * g_params[p_a];   // SCALING zr BY a (in line with pfastdm and rfastdm)
  
    // Get the value of CDF(t=Inf) (for the upper boundary??) (store in p[0])
    F_start(fc, BOUNDARY_UPPER);
    mint = g_Params->t0 - 0.5*g_Params->st0 ; // Min RT
    F    = F_get_F(fc,mint);
    nz   = F_get_N(fc);
    j    = (int) nz*scaled;
    p[0] = F[j];
  
    int length = rts.length();
    NumericVector out (length);
  
    if (boundary == BOUNDARY_UPPER) // Calc upper boundary
    {
#ifdef _CDF_DEBUG_
        Rcpp::Rcout << "pfastdm_b: Calculating upper boundary\n";
        Rcpp::Rcout << "Using length = "<< length  <<"\n";
        Rcpp::Rcout << "Using mint = "<< mint  <<"\n";
//      Rcpp::Rcout << "Using F = "<< XXXXX  <<"\n";
        Rcpp::Rcout << "Using nz = "<< nz  <<"\n";
        Rcpp::Rcout << "Using j = "<< j  <<"\n";
        Rcpp::Rcout << "Using p[0] = "<< p[0]  <<"\n";
      
#endif
      
        // ASSUMING F_start (,upper) has already been called (to get CDF(t=Inf))
        for (i=0; i < length; i++)
        {
            if (rts[i] <= mint) { out[i] = 0.0; }
            else
            {
                // if (i > 0)
                // {
                //     if (rts[i] <= rts[i-1]) Rcpp::stop("Must call Rfastdm.c with increasing t values (in_RTs[%d]=%g, in_RTs[i-1]=%g.", i, rts[i], rts[i-1]);
                // }
            
                F  = F_get_F(fc, rts[i]);
                nz = F_get_N(fc);
                j  = (int) nz*scaled;
            
                out[i] = F[j]-p[0];
            }
        }
    }
    else     // Calc lower boundary
    {
#ifdef _CDF_DEBUG_
      Rcpp::Rcout << "pfastdm_b: Calculating lower boundary\n";
#endif    
        F_start(fc, BOUNDARY_LOWER); 
        for (i=0; i < length; i++)
        {
            if (rts[i] <= mint) { out[i]=0.0; }
            else
            {
                // if (i > 0)
                // {
                //     if (rts[i] <= rts[i-1]) Rcpp::stop("Must call Rfastdm.c with increasing t values (in_RTs[%d]=%g, in_RTs[i-1]=%g.", i, rts[i], rts[i-1]);
                // }
                
                F  = F_get_F(fc, rts[i]);
                nz = F_get_N(fc);
                j  = (int) nz*scaled;
                out[i] = p[0]-F[j];
            }
          }
    }
    
#ifdef _CDF_DEBUG_
    Rcpp::Rcout << "distribution(): About to delete FC and return.\n";
#endif
    
    F_delete(fc);
    
    return out;
}


// A version of distribution using (slower) F_get_val instead
NumericVector precise_distribution (NumericVector rts, int boundary)
{
    struct  F_calculator *fc;
    int     i;

#ifdef _CDF_DEBUG_
  Rcpp::Rcout << "precise_distribution (): Entered precise_distribution ()\n";
  
  Rcpp::Rcout << "Using Parameters:\n";
  Rcpp::Rcout << "rts[0] = " << rts[0] << "\n";
  Rcpp::Rcout << "length of rts = " << rts.length() << "\n";
  Rcpp::Rcout << "boundary = " << boundary << "\n";
#endif

    double scaled = g_Params->zr * g_Params->a;   // SCALING zr BY a (in line with pfastdm and rfastdm)

    fc = F_new ();

    int length = rts.length();
    NumericVector out (length);
  
    F_start(fc, BOUNDARY_UPPER);
    double offset = F_get_val (fc, 0, scaled); // Subtract the value for t=0
    
    if (boundary == BOUNDARY_UPPER) // Calc upper boundary
    {
#ifdef _CDF_DEBUG_
        Rcpp::Rcout << "Upper bound offset = " << offset << "\n";
#endif
            
        for (i=0; i < length; i++)
        {
            out[i] = F_get_val (fc, rts[i], scaled) - offset;
        }
    
    }
    else
    {
        F_start(fc, BOUNDARY_LOWER);

#ifdef _CDF_DEBUG_
        Rcpp::Rcout << "Lower bound offset = " << offset << "\n";
#endif

        for (i=0; i < length; i++)
        {
          out[i] = offset - F_get_val (fc, rts[i], scaled);
        }
    }

#ifdef _CDF_DEBUG_
    Rcpp::Rcout << "precise_distribution(): About to delete FC and return.\n";
#endif
  
  F_delete(fc);
  
  return out;
}


#endif // DISTRIBUTION_H
