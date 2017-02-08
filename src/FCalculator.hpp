/* FCalculator.hpp - A structure containing the data and function pointers 
 *   required to calculate a 'dimension' of the CDF integration.
 *   (originally in cdf.c)
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

#ifndef FCALCULATOR_H
#define FCALCULATOR_H

/**********************************************************************
* struct F_calculator:
* Store information to calculate the cumulative distribution function F.
*
* Usage:
* 1) Allocate a F_calculator structure with the 'F_new' function below.
*    This initialises the appropriate method for the given variabilities.
* 2) Set the initial condition for the PDE with 'F_start'
* 3) Get an array of computed values at time t with 'F_get_F'.
*    The field 'F_calculator.N' gives the length of the array.
* 4) Get the z-value associated with array element 'i' using
*    the function 'F_get_z'.
*
* The values returned by the functions F_get_F and F_get_val are not
* directly the values of the CDF, but they are transformed to ease
* the use of the results by the higher levels of fast-dm.  To get the
* actual values of the CDF, the following transform has to be applied
* to the function results:
*
*    for b_upper: CDF = return value - F^-(\infty)
*    for b_lower: CDF = F^-(\infty) - return value
*
* When all variabilities are zero, F^-(\infty) can be computed using
* the function 'F_limit'.
*/


static double F_limit(double a, double z, double v)
{
    if (fabs(v) < 1e-8) { return 1 - z/a; } 
                   else { return (exp(-2*v*z)-exp(-2*v*a)) / (1-exp(-2*v*a)); }
}

struct F_calculator 
{
    int  N;      
    int  plus;   // boundary
    void *data;

    // To be filled in  
    void          (*start)  (F_calculator *, int plus);   // CONVERSION NOTE: 'plus' was enum boundary (0=lower,1=upper)
    void          (*free)   (F_calculator *);             // CONVERSION NOTE: was "delete" in C version
    const double *(*get_F)  (F_calculator *, double t);
    double        (*get_z)  (const F_calculator *, int i);
};

#endif // FCALCULATOR_H
