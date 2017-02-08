/* CDF_no_variability.hpp - Functions to calculate CDF when all variance 
 *   parameters are 0 (originally in cdf.c and pde.c)
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

#ifndef CDF_NO_VARIABILITY_H
#define CDF_NO_VARIABILITY_H

#include "Parameters.hpp"
#include "FController.hpp"

// Forward declarations
static void          F_plain_delete (F_calculator *fc);  
static void          F_plain_start  (F_calculator *fc, int plus);
static const double *F_plain_get_F  (F_calculator *fc, double t);
static double        F_plain_get_z  (const F_calculator *fc, int i);

void advance_to (int N, double *vector, double t0, double t1, double dz, double v);







struct F_plain_data 
{
    double  a, v, t0, d; // parameters (except z)
    double  dz;				   // z step-size
    double  t_offset;		 // time adjustment, resulting from t0 and d
    double  t;				   // adjusted time, corresponding to the vector F
    double *F;				   // state at time t + t_offset 
};






static F_calculator *F_plain_new (Parameters *params)
// Allocate a new 'struct F_calculator' (without variabilities).
{
    F_calculator *fc;
    F_plain_data *data;
    int  N;
    
    // N must be even, otherwise the case szr == 1 fails 
    N = 2*(int)(params->a*0.5/params->TUNE_DZ+0.5);
    if (N<4)  N = 4; 
    
    fc = xnew (struct F_calculator, 1);         // NOTE: MEMORY ALLOCATION
    fc->N = N;
    fc->plus = -1;
    
    data = xnew (struct F_plain_data, 1);       // NOTE: MEMORY ALLOCATION
    data->a  = params->a;
    data->v  = params->v;
    data->t0 = params->t0;
    data->d  = params->d;
    data->dz = params->a/N;
    data->F  = xnew (double, N+1);              // NOTE: MEMORY ALLOCATION
    fc->data = data;
    
    fc->start = F_plain_start;
    fc->free  = F_plain_delete;
    fc->get_F = F_plain_get_F;
    fc->get_z = F_plain_get_z;
    
    return  fc;
}










static void F_plain_delete (F_calculator *fc)
{
    F_plain_data *data = (F_plain_data*)fc->data;
    
    xfree (data->F);
    xfree (data);
    xfree (fc);
}







static void F_plain_start (F_calculator *fc, int plus)
{
  F_plain_data *data = (F_plain_data*)fc->data;
  double  a = data->a;
  double  v = data->v;
  int  N = fc->N;
  int  i;
  
  fc->plus = plus;
  data->t_offset = data->t0 - data->d * (plus == 1? 0.5 : -0.5);
  data->t = 0;
  
  data->F[0] = (plus == 1) ? 1 : 0;
  for (i=1; i<N; i++)
  {
      double  z = F_get_z (fc, i); 
      data->F[i] = F_limit(a, z, v);
  }
  data->F[N] = (plus == 1) ? 1 : 0;
}

static const double *F_plain_get_F (F_calculator *fc, double t)
{
    F_plain_data *data = (F_plain_data*)fc->data;
    
    t -= data->t_offset;
    
    if (t > data->t) 
    {
      advance_to (fc->N, data->F, data->t, t, data->dz, data->v);
      data->t = t;
    }
    return  data->F;
}


static double F_plain_get_z (const F_calculator *fc, int i)
{
    F_plain_data *data = (F_plain_data*)fc->data;
    return  i * data->dz;
}





// NOTE: Functions from pde.c 
//       (note the use of a static array of doubles and a realloc for optimisation)
//       TODO: Investigate if there's a better/faster/more robust/more C++ way of doing this
// NOTE: *res is in/out
//
static void solve_tridiag(int n, const double *rhs, double *res, double left, double mid, double right)
// Solve an n by n tridiagonal system of linear equations.
// The matrix has 'mid' on the diagonal, 'left' on the subdiagonal and 'right' on the superdiagonal.
//
{
    static double *tmp = NULL;
    static int tmp_len = 0;
    double p, old_res, old_tmp;
    int  i;
    
    if (n-1 > tmp_len) 
    {
        // Reallocating/freeing the scratch buffer for every call to 'solve_tridiag' caused about 10% of the total CPU load during
        // some fast-dm runs.  To avoid this problem, re-use the same buffer between runs if possible.
        tmp = xrenew(double, tmp, n-1);
        tmp_len = n-1;
    }
    
    /* step 1: solving forward */
    tmp[0] = old_tmp = right / mid;
    res[0] = old_res = rhs[0] / mid;
    for (i=1; i<n-1; ++i) 
    {
        p = 1.0/(mid-left*old_tmp);
        res[i] = old_res = (rhs[i] - left*old_res)*p;
        tmp[i] = old_tmp = right*p;
    }
    p = 1.0/(mid-left*old_tmp);
    res[n-1] = (rhs[n-1] - left*old_res)*p;
    
    // step 2: solving backward
    for (i=n-1; i>0; --i)  res[i-1] -= tmp[i-1]*res[i];
}

// NOTE: *vector is in/out
static void make_step (int N, double *vector, double dt, double dz, double v)
// Advance the numerical solution of the PDE by one step in time, using the Crank-Nicolson scheme.
// The time step size is 'dt', the space grid size is 'dz'.  */
{
    double	*tmp_vector;
    double  left, mid, right;
    int	i;
    
    tmp_vector = xnew (double, N+1);
    
    left  =  (1-dz*v) / (2*dz*dz);
    mid   =  -1 / (dz*dz);
    right =  (1+dz*v) / (2*dz*dz);
    
    tmp_vector[1] = (dt*left       * vector[0] + 
                    (1+0.5*dt*mid) * vector[1] + 
                    0.5*dt*right   * vector[2]);
    
    for (i=2; i<N-1; i++) 
    {
        tmp_vector[i] = (0.5*dt*left   * vector[i-1] +
                        (1+0.5*dt*mid) * vector[i] +
                        0.5*dt*right   * vector[i+1]);
    }
    tmp_vector[N-1] = (0.5*dt*left   * vector[N-2] +
                      (1+0.5*dt*mid) * vector[N-1] +
                      dt*right       * vector[N]);
    
    solve_tridiag(N-1, tmp_vector+1, vector+1, -0.5*dt*left, 1-0.5*dt*mid, -0.5*dt*right);
    
    xfree (tmp_vector);
}

// NOTE: *vector is in/out. Passed in as data->F from F_plain_get_F
void advance_to (int N, double *vector, double t0, double t1, double dz, double v)
// Advance the state 'vector' of the PDE from time 't0' to time 't1' 
{
    int  done = 0;
    
    do 
    {
        double  dt = g_Params->TUNE_PDE_DT_MIN + g_Params->TUNE_PDE_DT_SCALE*t0;
        if (dt > g_Params->TUNE_PDE_DT_MAX)  dt = g_Params->TUNE_PDE_DT_MAX;
        if (t0 + dt >= t1)
        {
            dt = t1 - t0;
            done = 1;
        }
        make_step (N, vector, dt, dz, v);
        t0 += dt;
    } while (!done);
}

#endif // CDF_NO_VARIABILITY_H
