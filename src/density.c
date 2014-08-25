/* density.c - compute the densities g- and g+ of the first exit time
 *
 * Copyright (C) 2012  Andreas Voss, Jochen Voss.
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

#define _USE_MATH_DEFINES
#include <math.h>

#include "fast-dm.h"


#define EPSILON 1e-6


#ifndef HAVE_FABS
double
fmax(double a, double b)
{
	return a > b ? a : b;
}
#endif

// [MG] This was causing a weird compiler fail under current gcc versions, which (now??) has its own isinf()
// If we see any more issues in the wild, we should try using R's ISINF function
/* 
#ifdef _WIN32
static int
isinf(double x)
{
	return (x-x != 0);
}
#endif
*/


static int
imax(int a, int b)
{
	return a > b ? a : b;
}


struct para {
	double t, a, zr, v, st0, szr, sv;
};

struct function {
	double (*f) (double x, void *data);
	void *data;
};

static double
integrate(struct function *F, double a, double b, double step_width)
{
	double width = b-a;
	int N = imax(4, (int) (width / step_width));
	double step = width / N;
	double x;
	double result = 0;

	for(x = a+0.5*step; x < b; x += step) {
		result += step * F->f(x, F->data);
	}
	return result;
}


static double
g_minus_small_time(double t, double zr, int N)
{
	int i;
	double sum = 0;

	for(i = -N/2; i <= N/2; i++) {
		double d = 2*i + zr;
		sum += exp(-d*d / (2*t)) * d;
	}

	return sum / sqrt(2*M_PI*t*t*t);
}

static double
g_minus_large_time(double t, double zr, int N)
{
	int i;
	double sum = 0;

	for(i = 1; i <= N; i++) {
		double d = i * M_PI;
		sum += exp(-0.5 * d*d * t) * sin(d*zr) * i;
	}

	return sum * M_PI;
}

static double
g_minus_no_var(double t, double a, double zr, double v)
{
	int N_small, N_large;
	double simple, factor, eps;
	double ta = t/(a*a);

	factor = exp(-a*zr*v - 0.5*v*v*t) / (a*a);
	if (isinf(factor)) {
		return 0;
	}
	eps = EPSILON / factor;

	N_large = (int)ceil(1/(M_PI*sqrt(t)));
	if (M_PI*ta*eps < 1) {
		N_large = imax(N_large,
					   (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
	}

	if (2*sqrt(2*M_PI*ta)*eps < 1) {
		N_small = (int)ceil(fmax(sqrt(ta) + 1,
								 2 + sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
	} else {
		N_small = 2;
	}

	if (N_small < N_large) {
		simple = g_minus_small_time(t/(a*a), zr, N_small);
	} else {
		simple = g_minus_large_time(t/(a*a), zr, N_large);
	}
	return factor * simple;
}

static double
integral_v_g_minus(double zr, void *data)
{
	struct para *P = data;
	double t = P->t;
	double a = P->a;
	double v = P->v;
	double sv = P->sv;
	int N_small, N_large;
	double simple, factor, eps;
	double ta = t/(a*a);

	factor = 1 / (a*a * sqrt(t * sv*sv + 1)) * exp(-0.5 * (v*v*t + 2*v*a*zr - a*zr*a*zr*sv*sv) / (t*sv*sv+1));
	if (isinf(factor)) {
		return 0;
	}
	eps = EPSILON / factor;

	if (P->sv == 0) {
		return g_minus_no_var(P->t, P->a, zr, P->v);
	}

	N_large = (int)ceil(1 / (M_PI*sqrt(t)));
	if (M_PI*ta*eps < 1) {
		N_large = imax(N_large,
					   (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
	}

	if (2*sqrt(2*M_PI*ta)*eps < 1) {
		N_small = (int)ceil(fmax(sqrt(ta)+1,
								 2+sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
	} else {
		N_small = 2;
	}

	if (N_small < N_large) {
		simple = g_minus_small_time(t/(a*a), zr, N_small);
	} else {
		simple = g_minus_large_time(t/(a*a), zr, N_large);
	}
	return factor * simple;
}

static double
integral_z_g_minus(double t, void *data)
{
	struct para *P = data;
	double res;

	if (t <= 0) return 0;

	P->t = t;
	if (P->szr == 0) {
		res = integral_v_g_minus(P->zr, P);
	} else {
		struct function F = {
			integral_v_g_minus,
			P
		};
		res = integrate(&F, P->zr - .5*P->szr, P->zr + .5*P->szr,
						TUNE_INT_Z) / P->szr;
	}
	return res;
}

static double
integral_t0_g_minus(double t, void *data)
{
	struct para *P = data;
	double res;

	if (P->st0 == 0) {
		res = integral_z_g_minus(t, P);
	} else {
		struct function F = {
			integral_z_g_minus,
			P
		};
		res = integrate(&F, t - .5*P->st0, t + .5*P->st0, TUNE_INT_T0) / P->st0;
	}
	return res;
}

double
g_minus(double t, const double *para)
{
	struct para P;

	P.a = para[p_a];
	P.zr = para[p_count];
	P.v = para[p_v];
	P.szr = para[p_szr];
	P.sv = para[p_sv];
	P.st0 = para[p_st0];

	return integral_t0_g_minus(t - para[p_t0] - 0.5*para[p_d], &P);
}

double
g_plus(double t, const double *para)
{
	struct para P;

	P.a = para[p_a];
	P.zr = 1 - para[p_count];
	P.v = -para[p_v];
	P.szr = para[p_szr];
	P.sv = para[p_sv];
	P.st0 = para[p_st0];

	return integral_t0_g_minus(t - para[p_t0] + 0.5*para[p_d], &P);
}
