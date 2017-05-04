/* CDF_sv_variability.hpp - Functions to calculate CDF when there is variance
 *   in the drift rate parameter (sv != 0) (originally in cdf.c and phi.c)
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

#ifndef CDF_SV_VARIABILITY_H
#define CDF_SV_VARIABILITY_H
//
// sv: variability in v
//
#include "CDF_sz_variability.h"

// Includes normal and inverse normal CDF code required only here (originally in phi.c)
double Phi (double x);
double Phi_inverse (double y);

static void          F_sv_delete (F_calculator *fc);
static void          F_sv_start  (F_calculator *fc, int plus);
static const double *F_sv_get_F  (F_calculator *fc, double t);
static double        F_sv_get_z  (const F_calculator *fc, int i);





struct F_sv_data
{
    int  nv;                // number of points in integration
    F_calculator **base_fc; // F_calculators for different v
    double *avg;
};








static F_calculator *F_sv_new (Parameters *params)
// Allocate a new 'struct F_calculator'.
//
// This initialises the PDE and prepares most things for the calculation.  The initial condition for the returned
// 'struct F_calculator' has to be set using 'F_start'.
//
// This function can deal with variabilities in all parameters. If 'sv == 0', it just return the result of 'F_st_new'.
{
    struct F_calculator **base_fc;
    struct F_calculator *fc;
    struct F_sv_data *data;

    int  nv, j;

    double sv = params->sv; // convenience only: cache sv locally

    if (sv < params->TUNE_SV_EPSILON)  return  F_sz_new (params);  // No need to integrate

    nv = (int)(sv/params->TUNE_DV + 0.5);
    if (nv < 3)  nv = 3;

    // Create a temp copy of the parameters
    Parameters temp_params = *params;    // SHOULD WORK, BUT CHECK THIS

    // Integrate across svs
    temp_params.sv = 0;
    base_fc = xnew (struct F_calculator *, nv);         // NOTE: MEMORY ALLOCATION
    for (j=0; j<nv; ++j)
    {
        double  x = Phi_inverse ((0.5+j)/nv);
        temp_params.v = sv*x + params->v;
        base_fc[j] = F_sz_new (&temp_params);
    }

    fc = xnew (struct F_calculator, 1);                   // NOTE: MEMORY ALLOCATION
    fc->N = base_fc[0]->N;
    fc->plus = -1;
    data = xnew (struct F_sv_data, 1);                    // NOTE: MEMORY ALLOCATION
    data->nv = nv;
    data->base_fc = base_fc;
    data->avg = xnew (double, fc->N+1);                   // NOTE: MEMORY ALLOCATION
    fc->data = data;

    fc->start = F_sv_start;
    fc->free = F_sv_delete;
    fc->get_F = F_sv_get_F;
    fc->get_z = F_sv_get_z;

    return  fc;
  }



static void F_sv_delete (F_calculator *fc)
{
    F_sv_data *data = (F_sv_data *)fc->data;
    int  j;

    for (j=0; j<data->nv; ++j) F_delete (data->base_fc[j]);
    xfree (data->base_fc);
    xfree (data->avg);
    xfree (data);
    xfree (fc);
}



static void F_sv_start (F_calculator *fc, int plus)
{
    F_sv_data *data = (F_sv_data *)fc->data;
    int  j;

    fc->plus = plus;
    for (j=0; j<data->nv; ++j) F_start (data->base_fc[j], plus);
}











static const double *F_sv_get_F (F_calculator *fc, double t)
{
    F_sv_data *data = (F_sv_data *)fc->data;
    const double *F;
    double *avg = data->avg;
    int  i, j;

    F = F_get_F(data->base_fc[0], t);
    for (i=0; i<=fc->N; ++i)  avg[i] = F[i];
    for (j=1; j<data->nv; ++j)
    {
        F = F_get_F(data->base_fc[j], t);
        for (i=0; i<=fc->N; ++i)  avg[i] += F[i];
    }
    for (i=0; i<=fc->N; ++i)  avg[i] /= data->nv;

    return  avg;
}


static double F_sv_get_z (const F_calculator *fc, int i)
{
    F_sv_data *data = (F_sv_data *)fc->data;
    return  F_get_z (data->base_fc[0], i);
}



double Phi (double x)
/* The distribution function of the standard normal distribution.  */
{
    return  0.5*(1+erf (x/M_SQRT2));
}

double Phi_inverse (double y)
/* The inverse of Phi, calculated using the bisection method */
{
    double  l, r;

    if (y<=0.5)
    {
        l = -1;
        while (Phi(l)>=y)  l -= 1;
        r = l+1;
    }
    else
    {
        r = 0;
        while (Phi(r)<y)  r += 1;
        l = r-1;
    }

    do
    {
        double m = 0.5*(l+r);
        if (Phi(m) < y) { l = m; }
                   else { r = m; }
    } while (r-l > 1e-8);

    return  0.5*(l+r);
}

#endif // CDF_SV_VARIABILITY_H
