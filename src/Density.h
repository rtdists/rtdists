/* Density.hpp - Functions for PDF calculation (originally in density.c)
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

#ifndef DENSITY_H
#define DENSITY_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp;

#define EPSILON 1e-6

// Forward declarations
double g_minus (double t);
double g_plus  (double t);

static double integral_t0_g_minus (double t, Parameters *params);
static double integral_z_g_minus  (double t, Parameters *params);
static double integral_v_g_minus  (double t, double zr, Parameters *params);

static double g_minus_no_var     (double t, double a, double zr, double v);
static double g_minus_small_time (double t, double zr, int N);
static double g_minus_large_time (double t, double zr, int N);

// TODO: Make sure these function names are accurate
static double integrate_z_over_t  (Parameters *params, double a, double b, double step_width);
static double integrate_v_over_zr (Parameters *params, double a, double b, double t, double step_width);


// Main calls
NumericVector density (NumericVector rts, int boundary)
{
    int length = rts.length();
    NumericVector out(length);

    if (boundary == 1) { for (int i = 0; i < length; i++) { out[i] =  g_plus(rts[i]);  } } // Calc upper
                  else { for (int i = 0; i < length; i++) { out[i] = -g_minus(rts[i]); } } // Calc lower
  
    return out;
}

double g_minus(double t)
{
    return integral_t0_g_minus (t - g_Params->t0 - 0.5*g_Params->d, g_Params);
}

double g_plus(double t)
{
    // Make a copy so we don't disturb our params 
    // (?TODO: we could optimise the object creation out and just set them back after the call)     
    Parameters new_params(*g_Params);
    new_params.zr = 1 - g_Params->zr;
    new_params.v = -g_Params->v;
  
    return integral_t0_g_minus (t - new_params.t0 + 0.5*new_params.d, &new_params);
}



static double integral_t0_g_minus (double t, Parameters *params)
{
    double res;
    
    if (params->st0 < params->TUNE_ST0_EPSILON) // 170501   was == 0) 
    {
        res = integral_z_g_minus (t, params);
    } 
    else 
    {
        res = integrate_z_over_t(params, 
                        t - .5*params->st0, 
                        t + .5*params->st0, params->TUNE_INT_T0) / params->st0;
    }

    return res;
}


static double integral_z_g_minus (double t, Parameters *params)
{
    double res;
    
    if (t <= 0) return 0;
    
    if (params->szr < params->TUNE_SZ_EPSILON) // 170501   was == 0) 
    {
        res = integral_v_g_minus (t, params->zr, params);
    } 
    else 
    {
        res = integrate_v_over_zr(params, params->zr - .5*params->szr, params->zr + .5*params->szr, 
                                  t, params->TUNE_INT_Z) / params->szr;
    }
    return res;
}


static double integral_v_g_minus (double t, double zr, Parameters *params)
{
    double a = params->a;
    double v = params->v;
    double sv = params->sv;

    int N_small, N_large;
    double simple, factor, eps;
    
    double ta = t/(a*a);
    
    factor = 1 / (a*a * sqrt(t * sv*sv + 1)) * exp(-0.5 * (v*v*t + 2*v*a*zr - a*a * zr*zr * sv*sv) / (t*sv*sv+1));
    
    if (std::isinf(factor)) 
    {
        return 0;  
    }
    
    eps = EPSILON / factor;
    
    if (params->sv == 0) 
    {
        return g_minus_no_var(t, a, zr, v);
    }
    
    N_large = (int)ceil(1 / (M_PI*sqrt(t)));
    if (M_PI*ta*eps < 1) 
    {
        N_large = std::max(N_large, (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
    }
    
    if (2*sqrt(2*M_PI*ta)*eps < 1) 
    {
        N_small = (int)ceil(fmax(sqrt(ta) + 1, 2 + sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
    } 
    else
    {
        N_small = 2;
    }
    
    if (N_small < N_large) 
    {
        simple = g_minus_small_time(t/(a*a), zr, N_small);
    } 
    else 
    {
        simple = g_minus_large_time(t/(a*a), zr, N_large);
    }
    return factor * simple;
}


static double g_minus_no_var(double t, double a, double zr, double v)
{
    int N_small, N_large;
    double simple, factor, eps;
    double ta = t/(a*a);
    
    factor = exp(-a*zr*v - 0.5*v*v*t) / (a*a);
    if (std::isinf(factor)) { return 0; }
    
    eps = EPSILON / factor;
    
    N_large = (int)ceil (1/ (M_PI*sqrt(t)));
    if (M_PI*ta*eps < 1) 
    {
        N_large = std::max(N_large, (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
    }
    
    if (2*sqrt(2*M_PI*ta)*eps < 1) 
    {
        N_small = (int)ceil(fmax(sqrt(ta) + 1, 2 + sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
    } 
    else
    {
        N_small = 2;
    }
    
    if (N_small < N_large) 
    {
        simple = g_minus_small_time(t/(a*a), zr, N_small);
    } 
    else 
    {
        simple = g_minus_large_time(t/(a*a), zr, N_large);
    }
    return factor * simple;
}


static double g_minus_small_time(double t, double zr, int N)
{
    int i;
    double sum = 0;
    double d;
    
    for(i = -N/2; i <= N/2; i++) 
    {
        d = 2*i + zr;
        sum += exp(-d*d / (2*t)) * d;
    }
    
    return sum / sqrt(2*M_PI*t*t*t);
}

static double g_minus_large_time(double t, double zr, int N)
{
    int i;
    double sum = 0;
    double d;
    
    for(i = 1; i <= N; i++) 
    {
        d = i * M_PI;
        sum += exp(-0.5 * d*d * t) * sin(d*zr) * i;
    }
    
    return sum * M_PI;
}

// --- Rectangular (midpoint) rule fallback ---
// Kept as a fallback when GSL adaptive quadrature fails to converge.
static double midpoint_integrate(double (*func)(double, void*), void *params, double a, double b, double step_width)
{
    double width = b - a;
    if (width <= 0) return 0;
    double tmp_N = width / step_width;
    if (std::isnan(tmp_N)) tmp_N = 20;
    int N = std::max(4, static_cast<int>(tmp_N));
    double step = std::max(width / N, EPSILON);
    double x, result = 0;
    for (x = a + 0.5*step; x < b; x += step)
        result += step * func(x, params);
    return result;
}

// GSL callback for integrate_z_over_t: integrates integral_z_g_minus over st0 variability
static double gsl_integrand_z_over_t(double x, void *p)
{
    Parameters *params = static_cast<Parameters*>(p);
    return integral_z_g_minus(x, params);
}

// GSL callback for integrate_v_over_zr: integrates integral_v_g_minus over sz variability
struct IntegrateVContext { double t; Parameters *params; };
static double gsl_integrand_v_over_zr(double x, void *ctx)
{
    IntegrateVContext *c = static_cast<IntegrateVContext*>(ctx);
    return integral_v_g_minus(c->t, x, c->params);
}

// Wraps GSL adaptive quadrature (QAGS for singularity handling) with fallback to midpoint rule
static double gsl_integrate_safe(double (*func)(double, void*), void *params,
                                 double a, double b, double epsabs, double epsrel,
                                 double fallback_step)
{
    if (a >= b) return 0;

    gsl_function f;
    f.function = func;
    f.params = params;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    // QAGS handles singularities/discontinuities (e.g. integrand is 0 for t<=0)
    int status = gsl_integration_qags(&f, a, b, epsabs, epsrel, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS)
    {
        // Fall back to the original midpoint rule
        result = midpoint_integrate(func, params, a, b, fallback_step);
    }
    return result;
}

static double integrate_z_over_t (Parameters *params, double a, double b, double step_width)
{
    return gsl_integrate_safe(&gsl_integrand_z_over_t, params,
                              a, b, params->TUNE_INT_EPSABS, params->TUNE_INT_EPSREL,
                              step_width);
}

static double integrate_v_over_zr (Parameters *params, double a, double b, double t, double step_width)
{
    IntegrateVContext ctx = { t, params };
    return gsl_integrate_safe(&gsl_integrand_v_over_zr, &ctx,
                              a, b, params->TUNE_INT_EPSABS, params->TUNE_INT_EPSREL,
                              step_width);
}

#endif // DENSITY_H
