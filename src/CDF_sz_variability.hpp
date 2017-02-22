/* CDF_sz_variability.hpp - Functions to calculate CDF when there is variance 
 *   in the start point parameter (szr != 0) (originally in cdf.c and phi.c)
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

#ifndef CDF_SZ_VARIABILITY_H
#define CDF_SZ_VARIABILITY_H
//
// sz: variability in z
//
#include "CDF_no_variability.hpp"


static void          F_sz_delete (F_calculator *fc);
static void          F_sz_start  (F_calculator *fc, int plus);
static const double *F_sz_get_F  (F_calculator *fc, double t);
static double        F_sz_get_z  (const F_calculator *fc, int i);








struct F_sz_data 
{
    F_calculator *base_fc;    // gives the values we average over 
    double *avg;              // the computed averages 
    int  k;                   // the average involves 2*k+1 cells 
    double  q;                // unused part of the outermost cells 
    double  f;                // scale factor for the integration 
};






static F_calculator *F_sz_new (Parameters *params)
// Allocate a new 'struct F_calculator' (with sv == 0 and st == 0).
//
// This function can deal with variabilities in z.
// If 'sz == 0', it just returns the result of 'F_plain_new'.
//
//
{
    struct F_calculator *base_fc;
    struct F_calculator *fc;
    struct F_sz_data *data;
    double  sz, tmp, dz;
    int  N, k;
    
    base_fc = F_plain_new (params);
    
    sz = params->szr*params->a;
    if (sz < params->TUNE_SZ_EPSILON)  return base_fc;
    
    N = base_fc->N;
    dz = F_get_z (base_fc, 1) - F_get_z (base_fc, 0);
    tmp = sz/(2*dz);
    k = (int)(ceil(tmp) + 0.5);
    assert (2*k <= N);
    
    fc = xnew (struct F_calculator, 1);         // NOTE: MEMORY ALLOCATION
    fc->N = N-2*k;
    fc->plus = -1;
    data = xnew (struct F_sz_data, 1);          // NOTE: MEMORY ALLOCATION
    data->base_fc = base_fc;
    data->avg = xnew (double, fc->N+1);         // NOTE: MEMORY ALLOCATION
    data->k = k;
    data->q = k - tmp;
    data->f = dz/sz;
    fc->data = data;
    
    fc->start = F_sz_start;
    fc->free = F_sz_delete;
    fc->get_F = F_sz_get_F;
    fc->get_z = F_sz_get_z;
    
    return  fc;
}




static void F_sz_delete (F_calculator *fc)
{
    F_sz_data *data = (F_sz_data *)fc->data;
  
    F_delete (data->base_fc);
    xfree (data->avg);
    xfree (data);
    xfree (fc);
}






static void F_sz_start (F_calculator *fc, int plus)
{
  F_sz_data *data = (F_sz_data *)fc->data;
  fc->plus = plus;
  F_start (data->base_fc, plus);
}














static const double *F_sz_get_F (F_calculator *fc, double t)
{
    F_sz_data *data = (F_sz_data *)fc->data;
  
    const double *F;
    double  tmp, q, f;
    int  i, j, m;
    
    F = F_get_F (data->base_fc, t);
    m = 2*data->k;
    q = data->q;
    f = data->f;
    if (m >= 3)
    {
        for (i=0; i<=fc->N; ++i) 
        {
            tmp  = F[i] * 0.5*(1-q)*(1-q);
            tmp += F[i+1] * (1-0.5*q*q);
            for (j=i+2; j<i+m-1; ++j)  tmp += F[j];
            tmp += F[i+m-1] * (1-0.5*q*q);
            tmp += F[i+m] * 0.5*(1-q)*(1-q);
            data->avg[i] = tmp * f;
        }
    } 
    else
    {
        /* m == 2 */
        for (i=0; i<=fc->N; ++i) 
        {
            tmp = F[i] * 0.5*(1-q)*(1-q);
            tmp += F[i+1] * (1-q*q);
            tmp += F[i+2] * 0.5*(1-q)*(1-q);
            data->avg[i] = tmp * f;
        }
    }
    /* m == 1 is impossible here */
    
    return  data->avg;
}


static double F_sz_get_z (const F_calculator *fc, int i)
{
    F_sz_data *data = (F_sz_data *)fc->data;
    return  F_get_z (data->base_fc, i+data->k);
}

#endif // CDF_SZ_VARIABILITY_H

