/* FController.hpp - Main controller class to manage an FCalculator structure 
 *   and call the function pointers stored within (originally in cdf.c)
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

#ifndef FCONTROLLER_H
#define FCONTROLLER_H

#include "CDF_st0_variability.hpp"

//
// All functions are externally visible
//

F_calculator *F_new ()
// Allocate data required to compute a CDF.
//
// CONVERSION NOTE: while only F_sz_new and F_plain_new rely on copying the global parameters to integrate over
//                  we'll pass them for all *_new calls for orthogonality
//
// The returned structure must be initialised with 'F_start' and must be freed with 'F_delete' after use.
{
    return  F_st0_new (g_Params);  
}
  
void F_delete (F_calculator *fc)
// Free a 'struct F_calculator' and all associated resources. 'fc' must have been allocated by 'F_new'.  After 'F_delete' is
// called, 'fc' cannot be used any longer.
{
    fc->free (fc);
}
  

void F_start (F_calculator *fc, int boundary)
// Set the initial condition for the PDE.
// If upper  boundary prepare to calculate the CDF for hitting a before 0, 
//   otherwise prepare to calculate the CDF for hitting 0 before a.
//
// CONVERSION NOTE: Changed from enum boundary b to int boundary, where 0 = lower and 1 = upper
{
    fc->start (fc, boundary);
}
  
int F_get_N (const F_calculator *fc)
{
    return  fc->N;
}
  
double F_get_z (const F_calculator *fc, int i)
// Get the z-value corresponding to index i.
{
    return  fc->get_z (fc, i);
}
  

const double *F_get_F (F_calculator *fc, double t)
// Get the array of CDF values at time t for all grid points z.
//
// 'F_start' must be used for initialisation before calling 'F_get_F'. Between calls of 'F_start' the calls to 'F_get_F' must have
// increasing values of 't'. 
// The returned array is owned by the 'struct F_calculator' and must not be changed or freed.
//
{
      return  fc->get_F (fc, t);
}

    
double F_get_val (F_calculator *fc, double t, double z)
// Get the value of the CDF for the parameters given when creating fc. 
// The function Uses linear interpolation for z-values between the grid points. 
// Don't use this function for parameter fitting, since it is not very fast (use 'F_get_F' instead).
//
// ? CONVERSION NOTE: Speed difference seems to be within +/- 10% for rt/parameter combinations so far tested
//
{
    const double *F;
    double  z0, z1;
    double  p, x;
    int  N = fc->N;
    int  i;
      
    F = F_get_F (fc, t);
    if (N == 0) 
    {
        x = F[0];
    } 
    else
    {
        z0 = F_get_z (fc, 0);
        z1 = F_get_z (fc, N);
        i = (int)(N*(z-z0)/(z1-z0));
        if (i < N) 
        {
            z0 = F_get_z (fc, i);
            z1 = F_get_z (fc, i+1);
            p = (z1-z) / (z1-z0);
            x = p*F[i] + (1-p)*F[i+1];
        } 
        else
        {
            x = F[N];
        }
    }
    return  x;
}

#endif // FCONTROLLER_H
