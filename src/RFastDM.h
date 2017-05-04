/* RFastDM.Hpp - Main header file for the RCpp implementation of fast-dm
 * 
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

#ifndef RFASTDM_H
#define RFASTDM_H

#include <Rcpp.h>
using namespace Rcpp;
#include <assert.h>

#define MAX_INPUT_VALUES 1e+6

// Include all .hpp files 
// Note: This is bad organisation, but Rcpp (and especially RStudio's "sourceCpp" don't seem to handle
//       projects with multiple .cpp files)


// Used by both PDF and CDF (also includes tuning params)
#include "Parameters.h"

// While not enforced, this is the global parameters Singleton 
//   To be created and freed in the d_, p_ and r_ calls in RFastDM.cpp
Parameters *g_Params;

#define BOUNDARY_LOWER 0
#define BOUNDARY_UPPER 1

// PDF
#include "Density.h"



// Placing memory routines used by CDF code here
//   TODO: (relying on C-style realloc until we have a better C++ solution - just vector<double>?)
void *xmalloc (size_t size)
{
    void *ptr;
  
    if (size == 0)  return NULL;
    ptr = malloc (size);
    if (ptr == NULL)  Rcpp::stop("memory exhausted");
    return  ptr;
}

// Originally this was used once as a mem optimization in solve_tridiag () (pde.c, line 43)
//  (relevant code is now in CDF_no_variability.hpp)
void *xrealloc (void *ptr, size_t newsize)
{
    if (newsize == 0) 
    {
        if (ptr)  free (ptr);
        return NULL;
    }
  
    ptr = ptr ? realloc (ptr, newsize) : malloc(newsize);
    if (ptr == NULL)  Rcpp::stop("memory exhausted");
    return  ptr;
}

#define xnew(T,N) ((T *)xmalloc((N)*sizeof(T)))
#define xrenew(T,OLD,N) ((T *)xrealloc(OLD,(N)*sizeof(T)))  

void xfree (void *ptr)
{
    if (ptr)  free (ptr);
}



// CDF code
#include "FCalculator.h"

// Forward declare all FController functions for access from all CDF_*
F_calculator *F_new ();
void          F_delete  (F_calculator *fc);
void          F_start   (F_calculator *fc, int boundary);
int           F_get_N   (const F_calculator *fc);
double        F_get_z   (const F_calculator *fc, int i);
const double *F_get_F   (F_calculator *fc, double t);
double        F_get_val (F_calculator *fc, double t, double z);

#include "FController.h"

#include "Distribution.h"
#include "Sampling.h"


#endif // RFASTDM_H
