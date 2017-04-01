/* Sampling.hpp - Contains main call for random sampling 
 *   (originally in construct-samples.c)
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

#ifndef SAMPLING_H
#define SAMPLING_H

#include <Rcpp.h>
using namespace Rcpp;

// Forward declarations for support functions
static int compare_doubles (const void *a, const void *b);
static int find_slot(double target, const double *value, int l, int r);
  

// RCpp Conversion Note: we've removed the idea of multiple samples (s_count)
List sampling (int s_size)
{
    F_calculator *fc;
    double  *F;
    double  *Fs,   // TODO: Drop this down to a 1D array
            Fs_min, Fs_max;
    double  t_min, t_max, dt;
    int  i, N;

    // get the F-values for the samples
    Fs = xnew(double, s_size);  // Conversion Note: not ready to rewrite entirely using Rcpp NumericVector etc. 
    
    Fs_min = 1;
    Fs_max = 0;

    // Create Fs[], an s_size matrix of random or sequential numbers (depending on random_flag), 
    //   and set Fs_min and Fs_max, the min and max produced
    bool random_flag = true;  // false = allow non-random sampling for testing 
    if (random_flag)
    {

        // ? TODO: vectorise this (though first attempt when quite badly performance-wise...)
        for (i=0; i<s_size; ++i)
        {
             Fs[i] = runif (1,0,1)[0];  // CHECK: Omitting n from runif uses scalar version but still requires subsetting?
                                        //        Compare to Rf_runif performance (but check if that needs get/set state)

            // Maintain the min and max of the sample
            if (Fs[i] > Fs_max)  Fs_max = Fs[i];
            if (Fs[i] < Fs_min)  Fs_min = Fs[i];
        }
    }
    else
    {
        // Generate equally-spaced samples
        for (i=0; i<s_size; ++i) { Fs[i] = (i+0.5)/s_size; }
        Fs_min = Fs[0];
        Fs_max = Fs[s_size-1];
    }

    // Create new F calculator
    fc = F_new ();
  
    // get the required t-range
    t_max = 0.5;
    F_start (fc, BOUNDARY_UPPER);
    double scaled_z = g_Params->zr * g_Params->a;
    while (F_get_val (fc, t_max, scaled_z) < Fs_max)	t_max += 0.1;

    t_min = -0.5;
    F_start (fc, BOUNDARY_LOWER);
    while (F_get_val (fc, -t_min, scaled_z) > Fs_min)	t_min -= 0.1;
  
    // get a table of F-values
    N = (int)((t_max-t_min)/0.001 + 0.5);
    dt = (t_max-t_min)/N;
    
    F = xnew(double, N+1);
  
    F_start (fc, BOUNDARY_UPPER);
    for (i=0; i<=N; ++i)
    {
        double  t = t_min+i*dt;
        if (t < 0)  continue;
        F[i] = F_get_val (fc, t, scaled_z);
    }
    F_start (fc, BOUNDARY_LOWER);
    for (i=N; i>=0; --i)
    {
        double  t = -(t_min+i*dt);
        if (t < 0)  continue;
        F[i] = F_get_val (fc, t, scaled_z);
    }

    F_delete (fc);
  
    // protect against rounding errors: make F increasing and restrict to valid range
    for (i=0; i<=N; ++i)
    {
        if (F[i] < 0)	F[i] = 0;
        if (F[i] > 1)	F[i] = 1;
    }

    qsort(F, N+1, sizeof(double), compare_doubles);
    if (F[0] > Fs_min)		F[0] = Fs_min;
    if (F[N] < Fs_max)		F[N] = Fs_max;

    NumericVector out_RTs(s_size);
    NumericVector out_bounds(s_size);
    
    for (i=0; i<s_size; ++i)
    {
        double  y = Fs[i];
        int  k = find_slot(y, F, 0, N);
        double  t = t_min + (k + (y-F[k])/(F[k+1]-F[k]))*dt;
      
        assert (F[k]<=y && y<=F[k+1]);
      
        out_bounds[i] = t >= 0;
        out_RTs[i] = fabs(t);
    }
  
    // Free up allocated objects
    xfree(F);
    xfree(Fs);

    return List::create(Named("rt") = out_RTs, Named("boundary") = out_bounds);
}



// Support functions for rfastdm ()
static int compare_doubles (const void *a, const void *b)
{
  const double *da = (double *)a;
  const double *db = (double *)b;
  
  if (*da < *db)  return -1;
  if (*da > *db)  return 1;
  return  0;
}

static int find_slot(double target, const double *value, int l, int r)
{
  int m = (l+r)/2;
  if (m == l) { return l; }
  else if ( value[m] > target) { return find_slot(target, value, l, m); }
  else    { return find_slot(target, value, m, r); }
}

#endif // SAMPLING_H



// OLDER VERSION... TODO: Delete this
//
// List sampling (int s_size)
// {
//   F_calculator *fc;
//   double  *F;
//   double  **Fs,   // TODO: Drop this down to a 1D array
//   Fs_min, Fs_max;
//   double  t_min, t_max, dt;
//   int  i, j, N;
//   
//   bool random_flag = true;  // false = allow non-random sampling for testing 
//   
//   // get the F-values for the samples
//   Fs = xnew(double *, 1);
//   Fs_min = 1;
//   Fs_max = 0;
//   
//   // Create Fs[][], an s_count X s_size matrix of random or sequential numbers,
//   //   and set Fs_min and Fs_max, the min and max produced
//   
//   j = 0;  // !! Interim HACK for s_count = 1
//   
//   Fs[j] = xnew(double, s_size);
//   
//   if (random_flag)
//   {
//     // ! TODO: vectorise this...
//     for (i=0; i<s_size; ++i)
//     {
//       Fs[j][i] = runif (1,0,1)[0];  // CHECK: Omitting n from runif uses scalar version but still requires subsetting?
//       //        Compare to Rf_runif performance (but check if that needs get/set state)
//       
//       //            Rcout << "Generated random uniform of " << Fs[j][i] << "\n";
//       // Get the min and max of the sample
//       if (Fs[j][i] > Fs_max)  Fs_max = Fs[j][i];
//       if (Fs[j][i] < Fs_min)  Fs_min = Fs[j][i];
//     }
//   }
//   else
//   {
//     // Generate equal-spaced samples
//     for (i=0; i<s_size; ++i)
//     {
//       Fs[j][i] = (i+0.5)/s_size;
//     }
//     Fs_min = Fs[j][0];
//     Fs_max = Fs[j][s_size-1];
//   }
//   
//   // Create new F calculator
//   fc = F_new ();
//   
//   // get the required t-range
//   t_max = 0.5;
//   F_start (fc, BOUNDARY_UPPER);
//   double scaled_z = g_Params->zr * g_Params->a;
//   while (F_get_val (fc, t_max, scaled_z) < Fs_max)	t_max += 0.1;
//   
//   t_min = -0.5;
//   F_start (fc, BOUNDARY_LOWER);
//   while (F_get_val (fc, -t_min, scaled_z) > Fs_min)	t_min -= 0.1;
//   
//   // get a table of F-values
//   N = (int)((t_max-t_min)/0.001 + 0.5);
//   dt = (t_max-t_min)/N;
//   F = xnew(double, N+1);
//   
//   F_start (fc, BOUNDARY_UPPER);
//   for (i=0; i<=N; ++i)
//   {
//     double  t = t_min+i*dt;
//     if (t < 0)  continue;
//     F[i] = F_get_val (fc, t, scaled_z);
//   }
//   F_start (fc, BOUNDARY_LOWER);
//   for (i=N; i>=0; --i)
//   {
//     double  t = -(t_min+i*dt);
//     if (t < 0)  continue;
//     F[i] = F_get_val (fc, t, scaled_z);
//   }
//   
//   F_delete (fc);
//   
//   // protect against rounding errors: make F increasing and restrict to valid range
//   for (i=0; i<=N; ++i)
//   {
//     if (F[i] < 0)	F[i] = 0;
//     if (F[i] > 1)	F[i] = 1;
//   }
//   
//   qsort(F, N+1, sizeof(double), compare_doubles);
//   if (F[0] > Fs_min)		F[0] = Fs_min;
//   if (F[N] < Fs_max)		F[N] = Fs_max;
//   
//   NumericVector out_RTs(s_size);
//   NumericVector out_bounds(s_size);
//   
//   for (i=0; i<s_size; ++i)
//   {
//     double  y = Fs[j][i];
//     int  k = find_slot(y, F, 0, N);
//     double  t = t_min + (k + (y-F[k])/(F[k+1]-F[k]))*dt;
//     
//     assert (F[k]<=y && y<=F[k+1]);
//     
//     out_bounds[(j*s_size)+i] = t >= 0;
//     out_RTs[(j*s_size)+i] = fabs(t);
//   }
//   
//   // Free up allocated objects
//   xfree(F);
//   xfree(Fs[j]);
//   xfree(Fs);
//   
//   return List::create(Named("rt") = out_RTs, Named("boundary") = out_bounds);
// }
