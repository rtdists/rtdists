/* CDF_st0_variability.hpp - Functions to calculate CDF when there is variance
 *   in the non-decision-time parameter (st0 != 0) (originally in cdf.c)
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

#ifndef CDF_ST0_VARIABILITY_H
#define CDF_ST0_VARIABILITY_H
//
// st0: variability in t0
//
// This implements numerical integration in t-direction, using the trapezoidal rule.  Since computing the CDF is slow and since we can
// solve the PDE only forward in time, we cache old values of the CDF.
//
// The cached values form a grid of M different t-values such that the smallest cached t-value is smaller or equal than t-0.5*st0, the
// biggest cached t-value is bigger or equal than t+0.5*st0.  The total length of the cached time interval is (M-1)*dt where M and dt
// are chosen such that st0 = (M-2)*dt.

#include "CDF_sv_variability.hpp"

// Forward declarations
static void F_st0_start (F_calculator *fc, int plus);
static const double *F_st0_get_F (F_calculator *fc, double t);
static double F_st0_get_z (const F_calculator *fc, int i);
static void F_st0_delete (F_calculator *fc);

struct F_st0_data
{
    struct F_calculator *base_fc;
    double  st0;		// variability of t0
    int     M;			// number of stored grid lines
    double  start;	// t-value of first stored grid line
    double  dt;		  // t-spacing of stored grid lines
    double *values;	// array: stored grid lines (length M*(N+1))
    char   *valid;	// which lines in 'values' are valid
    int     base;		// first grid line starts at pos. base*(N+1)
    double *avg;		// the computed average (size N+1)
};


static struct F_calculator *F_st0_new (Parameters *params)
// Allocate a new 'struct F_calculator' (with sv == 0).
//
// This function can deal with variabilities in z and t.
// If 'st0 == 0', it just returns the result of 'F_sz_new'.
//
//
{
    F_calculator *base_fc;
    F_calculator *fc;
    F_st0_data *data;
    double  dt;
    int  M, N;

    double st0 = params->st0; // convenience only: cache st0 locally

    base_fc = F_sv_new (params);
    if (st0 <= params->TUNE_DT0*1e-6)  return base_fc;

    M = (int)(st0/params->TUNE_DT0 + 1.5);
    if (M < 3)  M = 3;
    dt = st0/(M-2);
    N = base_fc->N;

    fc = xnew (F_calculator, 1);     // NOTE: MEMORY ALLOCATION: struct F_calculator
    fc->N = N;
    fc->plus = -1;
    data = xnew (F_st0_data, 1);     // NOTE: MEMORY ALLOCATION: struct F_st0_data
    data->st0 = st0;
    data->base_fc = base_fc;
    data->M = M;
    /* data->start is set in F_st0_start */
    data->dt = dt;
    data->values = xnew (double, M*(N+1));  // NOTE: MEMORY ALLOCATION: M*(N+1) doubles
    data->valid = xnew (char, M);           // NOTE: MEMORY ALLOCATION: M chars
    data->base = 0;
    data->avg = xnew (double, N+1);         // NOTE: MEMORY ALLOCATION: N+1 doubles
    fc->data = data;

    fc->start = F_st0_start;
    fc->free = F_st0_delete;
    fc->get_F = F_st0_get_F;
    fc->get_z = F_st0_get_z;

    return  fc;
  }


static void F_st0_delete (F_calculator *fc)
{
    F_st0_data *data = (F_st0_data *)fc->data;

    F_delete (data->base_fc);
    xfree (data->valid);
    xfree (data->values);
    xfree (data->avg);
    xfree (data);
    xfree (fc);
}




static void F_st0_start (F_calculator *fc, int plus)
{
    F_st0_data *data = (F_st0_data *)fc->data;
    int j;

    fc->plus = plus;
    F_start (data->base_fc, plus);
    data->start = -DBL_MAX;

    // initially mark all of the cache as invalid
    for (j = 0; j < data->M; ++j) data->valid[j] = 0;
}







static const double *F_st0_get_row(const F_calculator *fc, int j)
// Get a pointer to one of the stored grid lines.
// The value j is the grid line index (range 0, ..., M-1).
// The returned value is a pointer to an array of length N+1.
{
    const F_st0_data *data = (F_st0_data *)fc->data;

    int  M, N, idx;
    double *row;

    M = data->M;
    N = fc->N;
    assert(0 <= j && j < M);
    idx = (data->base + j)%M;
    row = data->values + idx*(N+1);

    if (! data->valid[idx])
    {
        double t;
        const double *F;

        t = data->start + j*data->dt;
        F = F_get_F(data->base_fc, t);
        memcpy(row, F, (N+1)*sizeof(double));
        data->valid[idx] = 1;
    }

    return row;
}

static void add_vec(long n, double a, const double *x, double *y)
// the vector operation y += a*x
{
#ifdef HAVE_LIBBLAS
    extern void daxpy_(long *Np, double *DAp, const double *X, long *INCXp, double *Y, long *INCYp);
    long inc = 1;
    daxpy_(&n, &a, x, &inc, y, &inc);
#else /* ! HAVE_LIBBLAS */
    int i;
    if (a == 1) { for (i=0; i<n; ++i) y[i] += x[i]; }
           else { for (i=0; i<n; ++i) y[i] += a*x[i]; }
#endif /* ! HAVE_LIBBLAS */
}

static const double *F_st0_get_F (struct F_calculator *fc, double t)
{
    F_st0_data *data = (F_st0_data*)fc->data;
    double  a, b, dt;
    const double *row;
    double  q, r, *avg;
    int  M, N, shift;
    int  i, j, m;

    a = t - 0.5*data->st0;
    b = t + 0.5*data->st0;
    dt = data->dt;
    M = data->M;
    N = fc->N;

    // how many of the precalculated rows can we keep?
    if (a - data->start >= M*dt)
    {
        // beware of integer overflows for small dt
        shift = M;
    }
    else
    {
        shift = (int)((a - data->start)/dt);
        assert (shift >= 0);
    }

    for (j=0; j<shift; ++j) data->valid[(data->base+j)%M] = 0;

    if (shift < M)
    {
        data->start += shift*dt;
        data->base = (data->base+shift)%M;
    }
    else
    {
        data->start = a;
    }

    /* compute the average over the rows from a to b */
    avg = data->avg;
    for (i=0; i<=N; ++i)  avg[i] = 0;
    {
        double  tmp = (b - data->start)/dt;
        m = (int)(ceil (tmp) + 0.5);
        if (m >= M)  m = M-1; /* protect against rounding errors */
        q = (a - data->start)/dt;
        r = m - tmp;
    }

    if (m >= 3)
    {
        row = F_st0_get_row(fc, 0);
        add_vec(N+1, 0.5*(1-q)*(1-q), row, avg);

        row = F_st0_get_row(fc, 1);
        add_vec(N+1, 1-0.5*q*q, row, avg);

        for (j=2; j<m-1; ++j)
        {
          row = F_st0_get_row(fc, j);
          add_vec(N+1, 1, row, avg);
        }

        row = F_st0_get_row(fc, m-1);
        add_vec(N+1, 1-0.5*r*r, row, avg);

        row = F_st0_get_row(fc, m);
        add_vec(N+1, 0.5*(1-r)*(1-r), row, avg);
    }
    else if (m == 2)
    {
        row = F_st0_get_row(fc, 0);
        add_vec(N+1, 0.5*(1-q)*(1-q), row, avg);

        row = F_st0_get_row(fc, 1);
        add_vec(N+1, 1-0.5*(q*q+r*r), row, avg);

        row = F_st0_get_row(fc, 2);
        add_vec(N+1, 0.5*(1-r)*(1-r), row, avg);
    }
    else if (m == 1)
    {
        row = F_st0_get_row(fc, 0);
        add_vec(N+1, 0.5*((1-q)*(1-q)-r*r), row, avg);

        row = F_st0_get_row(fc, 1);
        add_vec(N+1, 0.5*((1-r)*(1-r)-q*q), row, avg);
    }

    for (i=0; i<=N; ++i)  avg[i] *= dt/(b-a);
    return  avg;
}

static double F_st0_get_z (const F_calculator *fc, int i)
{
    F_st0_data *data = (F_st0_data*)fc->data;
    return  F_get_z (data->base_fc, i);
}

#endif // CDF_ST0_VARIABILITY
