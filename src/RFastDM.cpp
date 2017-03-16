/* RFastDM.cpp - Main source file for the RCpp implementation of fast-dm
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

#include <Rcpp.h>
#include <iostream>
#include <sstream>

#include "RFastDM.hpp"

using namespace Rcpp;

// R-callable PDF for fastdm - pass boundary to retrieve (1 = lower, 2 = upper)
// [[Rcpp::export]]
NumericVector d_fastdm (NumericVector rts, NumericVector params, double precision=3, int boundary=2, bool stop_on_error=true)
{
    int length = rts.length();
    if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }  
  
    if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }
  
    g_Params = new Parameters (params, precision);

    NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..
    
    if (!g_Params->ValidateParams(stop_on_error)) 
    {
        if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
                      else { return out; }
    }
  
    out = density (rts, boundary-1);

    delete g_Params;
    return out;
}


// R-callable CDF which finds the left-hand area (from 0 to RT)  - pass boundary to retrieve (1 = lower, 2 = upper)
// [[Rcpp::export]]
NumericVector p_fastdm (NumericVector rts, NumericVector params, double precision=3, int boundary=2, bool stop_on_error=true)
{
    int length = rts.length();
    if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }  
    
    if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }
  
    g_Params = new Parameters (params, precision);
    
    NumericVector out(length, 0.0);
    
    if (!g_Params->ValidateParams(stop_on_error)) 
    {
        if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
        else return out;
    }

    // Pass through to distribution in Distribution.hpp
    out = distribution (rts, boundary-1);
    
    delete g_Params;
    return out;
}



// More precise R-callable CDF which finds the left-hand area (from 0 to RT) - pass boundary to retrieve (1 = lower, 2 = upper)
// [[Rcpp::export]]
NumericVector p_precise_fastdm (NumericVector rts, NumericVector params, double precision=3, int boundary=2, bool stop_on_error=true)
{
    int length = rts.length();
    if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }  
  
    if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }
  
  
    g_Params = new Parameters (params, precision);
    
    NumericVector out(length, 0.0);
    
    if (!g_Params->ValidateParams(stop_on_error)) 
    {
        if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
                      else return out;
    }
    
    // Pass through to precise_distribution.hpp
    out = precise_distribution (rts, boundary-1);
    delete g_Params;
  
    return out;
}


// R-style sampling from the DM - returns a List consisting of RTs and boundaries
// [[Rcpp::export]]
List r_fastdm (int num_values, NumericVector params, double precision=3, bool stop_on_error=true)
{
    if ((num_values < 1) || (num_values > MAX_INPUT_VALUES))
    {
        Rcpp::stop("Number of samples requested exceeds maximum of %d.\n", MAX_INPUT_VALUES);
    }  
    
    g_Params = new Parameters (params, precision);
    
    if (!g_Params->ValidateParams(stop_on_error)) 
    {
      if (stop_on_error)
      {
          Rcpp::stop("Error validating parameters.\n"); 
      }
      else
      {
          NumericVector out_RTs(num_values, 0.0);
          NumericVector out_bounds(num_values, 0.0);
          return List::create(Named("rt") = out_RTs, Named("boundary") = out_bounds);
      }
    }
    
    // Pass through to Sampling.hpp
    List out = sampling (num_values);

    delete g_Params;
    
    return out;
}
