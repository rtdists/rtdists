/*****************************************************************************/
//
// R Wrapper for fastdm.
// Contains functions for calculating pdf, cdf, and random sampling
//
// fast-dm - a fast method for diffusion model analysis by:
//    Jochen Voss <voss@seehuhn.de>
//    Andreas Voss <andreas.voss@psychologie.uni-heidelberg.de>
//
// GNU General Public License v2:
/*
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
//
/*****************************************************************************/
// Any questions or problems with this wrapper?
//   Contact: Matthew Gretton (matthew.gretton@uon.edu.au)
//
// Compile via R:
/*
     R CMD SHLIB Rfastdm2.c density.c cdf.c pde.c phi.c precision.c xmalloc.c
*/
//   NOTE: To enable compilation options I added a home/.R/Makevars file containing "CFLAGS=CFLAGS=-g -O2 -Wall"
//   Do -c for clean and --preclean to clean previous versions
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "fast-dm.h"

#define FASTDM_DEBUG  1  // Undef below for release
#define SHOW_INFO     1  // Show warnings / extra info, etc. (undef below)
#undef  FASTDM_DEBUG
#undef  SHOW_INFO

#define USE_R_PRNG  1 // Turn off to use internal PRNG

#ifndef USE_R_PRNG    // Don't bother including if we're using R's PRNG
  #include "jvrandom.h"
#endif

// Andrew Heathcote's additional bounds (no longer needed?)
//#define ADDITIONAL_BOUNDS 1  
#undef ADDITIONAL_BOUNDS    

// Maximum number of values that each function can return (as per Voss&Voss implementation)
#define MAX_VALUES 1000000

// Forward declarations for R-callable functions
void dfastdm   (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, double *out_densities_a, double *out_densities_0);
void dfastdm_b (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, int *boundary, double *out_densities);
void pfastdm   (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, double *out_pvalues_upper, double *out_pvalues_lower);
void pfastdm_b (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, int *boundary, double *out_pvalues);
void rfastdm   (int *in_numvalues, double *in_params, double *in_precision, double *out_RTs, double *out_bounds);


// Support function forward declarations for rfastdm ()
static int compare_doubles (const void *a, const void *b);
static int find_slot(double target, const double *value, int l, int r);

// Reimplementing these for R compatibility
extern void params_write(double para[p_count], double zr, double precision, int n);
extern void params_check(double para[p_count], double zr);

// Order of parameters (note: doesn't include zr..)
// enum parameter_index {
//	p_a,       // Boundary separation
//	p_v,       // Mean of the drift
//	p_t0,      // Non-decision time
//	p_d,       // Difference between boundaries of non-decision time IS THIS NEW AS WELL???
//	p_szr,     // width of zr distribution
//	p_sv,      // standard deviation of v distribution
//	p_st0,     // width of t0 distribution
//	p_count    // Num parameters but often used internally to refer to zr (...!?)
//};
//  zr         // Mean of diffusion starting point relative to boundaries

// Global variables for common stuff (nasty shortcut, used more to ensure correctness between functions)
double g_precision;
double g_params[7];  // An array of 8 parameters is passed in from R, but internally the code uses 7 + zr
double g_zr;
int    g_num_values;

// Internal helper function to set up common elements for each function
void _setup (double num_values, double *params, double precision)
{
#ifdef SHOW_INFO
  #ifdef HAVE_LIBBLAS
    Rprintf ("Note: Using LIBBLAS.\n");
  #else
    Rprintf ("Note: Not using LIBBLAS.\n");
  #endif
#endif

#ifdef FASTDM_DEBUG
    Rprintf ("In _setup.\n");
#endif
    // Set precision
    g_precision = precision;
    set_precision (g_precision);

    // Set up the parameters
    for (int i=0; i<7;i++)
    {
        g_params[i] = params[i];
    }
    g_zr = params[7];

#ifdef SHOW_INFO
    Rprintf ("Writing out parameters.\n");
    params_write(g_params, g_zr, g_precision, 0);
#endif

	params_check(g_params, g_zr);

    // Get number of values we're calculating
    g_num_values = (int)num_values;
    if ((g_num_values <= 0) || (g_num_values > MAX_VALUES))
    {
        Rf_error ("Number of values requested is either <= 0 or exceeds maximum of %d\n", MAX_VALUES);
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// R-callable probability density function (PDF) for fastdm
void dfastdm   (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, double *out_densities_a, double *out_densities_0)
{
#ifdef FASTDM_DEBUG
    Rprintf ("In dfastdm.\n");
#endif

    _setup (*in_numvalues, in_params, *in_precision);

    double *RTs = in_RTs;

    // Coercing params into format g_plus/g_minus wants
    double params[8];
    for (int i=0; i<7; i++) { params[i] = g_params[i]; }
    params[7] = g_zr;

    // Loop through each RT and add to the output vector
    for (int i = 0; i < g_num_values; i++)
    {
        out_densities_a[i] = g_plus(RTs[i], params);
        out_densities_0[i] = -g_minus(RTs[i], params);
#ifdef SHOW_INFO
        Rprintf ("Calculating density for RT[%2d] = %3.4g: upper: %.4g, lower: %.4g.\n", i, RTs[i], out_densities_a[i], out_densities_0[i]);
#endif
	}
}

// R-callable probability density function (PDF) for fastdm - pass boundary to retrieve (1 = lower, 2 = upper)
void dfastdm_b (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, int *boundary, double *out_densities)
{
#ifdef FASTDM_DEBUG
    Rprintf ("In dfastdm.\n");
#endif

    _setup (*in_numvalues, in_params, *in_precision);

    int bound = *boundary;
    if ((bound < 1) || (bound > 2)) Rf_error ("Error: invalid boundary!\n");

    // Coercing params into format g_plus/g_minus wants
    double params[8];
    for (int i=0; i<7; i++) { params[i] = g_params[i]; }
    params[7] = g_zr;

    if (bound == 2)  // Calc upper
    {
        for (int i = 0; i < g_num_values; i++)
        {
            out_densities[i] = g_plus(in_RTs[i], params);
    	}
    }
    else // Calc lower
    {
        for (int i = 0; i < g_num_values; i++)
        {
            out_densities[i] = -g_minus(in_RTs[i], params);
    	}
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// An R-like version which finds the left-hand area (from 0 to RT)
void pfastdm   (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, double *out_pvalues_upper, double *out_pvalues_lower)
{
#ifdef FASTDM_DEBUG
    Rprintf ("In pfastdm.\n");
#endif

    _setup (*in_numvalues, in_params, *in_precision);

	struct F_calculator *fc;

    double z;
	z = g_zr * g_params[p_a];
    double *RTs = in_RTs;
	fc = F_new (g_params);
    F_start (fc, b_upper);

    double val;
    double Fv;
#ifdef SHOW_INFO
    double Fz;
    const double *FF;
    double n;
#endif

    for (int i = 0; i < g_num_values; i++)
    {
        val = RTs[i];
        Fv = F_get_val (fc, val, z);

#ifdef SHOW_INFO
        Fz = F_get_z (fc, val);
        FF = F_get_F (fc, val);
        n = F_get_N (fc);

        Rprintf ("Calculating upper p-value for RT %d = %g.\n", i, val);
        Rprintf ("  F_get_z returns %g.\n", Fz);
        Rprintf ("  F_get_val returns %g.\n", Fv);
        Rprintf ("  F_get_F returns ");
for (int j = 0; j < n; j++)
{
        Rprintf (" %g,", FF[j]);
}
        Rprintf ("\n");
#endif
        out_pvalues_upper[i] = Fv;
	}

    // Do lower boundary
    F_start (fc, b_lower);

#ifdef FASTDM_DEBUG
    Rprintf ("F_get_N returns %d.\n", F_get_N (fc));
#endif

    for (int i = 0; i < g_num_values; i++)
    {
        val = RTs[i];
#ifdef SHOW_INFO
        Fv = F_get_val (fc, val, z);
        Fz = F_get_z (fc, val);
        FF = F_get_F (fc, val);
        double n = F_get_N (fc);

	     	Rprintf("Calculating lower p-value for RT %d = %g.\n", i, val);
        Rprintf ("  F_get_z returns %g.\n", Fz);
        Rprintf ("  F_get_val returns %g.\n", Fv);
        Rprintf ("  F_get_F returns ");
for (int j = 0; j < n; j++)
{
        Rprintf (" %g,", FF[j]);
}
        Rprintf ("\n");

#endif
	    out_pvalues_lower[g_num_values-i] = Fv;
	}

	F_delete (fc);
}

// An R-like version which finds the left-hand area (from 0 to RT) - uses Scott Brown's code, pass in boundary to retrieve
void pfastdm_b (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, int *boundary, double *out_pvalues)
{
	struct  F_calculator *fc;
	double  mint;

	const double *F;
	int     i,j,nz;

#ifdef FASTDM_DEBUG
    Rprintf ("In pfastdm_b. Calling _setup ()\n");
#endif
    double p[1]; // This is CDF(t=Inf)?

#ifdef SHOW_INFO
    Rprintf ("p = %g\n", p);
#endif
    _setup (*in_numvalues, in_params, *in_precision);

    int bound = *boundary;
    if ((bound < 1) || (bound > 2)) Rf_error ("Error: invalid boundary!\n");

#ifdef FASTDM_DEBUG
    Rprintf ("Creating F object.\n");
#endif

	fc = F_new(g_params);

  double scaled = g_zr;// * g_params[p_a];   // SCALING zr BY a (in line with pfastdm and rfastdm)

	// Get the value of CDF(t=Inf) (for the upper boundary??) (store in p[0])
	F_start(fc, b_upper);
	mint = g_params[p_t0] - 0.5*g_params[p_st0] ; // Min RT
	F    = F_get_F(fc,mint);
	nz   = F_get_N(fc);
	j    = (int) nz*scaled;
	p[0] = F[j];

    if (bound == 2) // Calc upper boundary
    {
        #ifdef FASTDM_DEBUG
            Rprintf ("Calculating upper boundary\n");
        #endif

        // ASSUMING F_start (,upper) has already been called (to get CDF(t=Inf))
		for (i=0; i < g_num_values; i++)
        {
			if (in_RTs[i] <= mint) { out_pvalues[i] = 0.0; }
            else
            {
				if (i > 0)
                {
                    if (in_RTs[i] <= in_RTs[i-1]) Rf_error("Must call Rfastdm.c with increasing t values (in_RTs[%d]=%g, in_RTs[i-1]=%g.", i, in_RTs[i], in_RTs[i-1]);
                }

				F  = F_get_F(fc, in_RTs[i]);
				nz = F_get_N(fc);
				j  = (int) nz*scaled;

				out_pvalues[i] = F[j]-p[0];
			}
		}
	}
    else     // Calc lower boundary
    {
        #ifdef FASTDM_DEBUG
            Rprintf ("Calculating lower boundary\n");
        #endif

		F_start(fc, b_lower);
		for (i=0; i < g_num_values; i++)
        {
			if (in_RTs[i] <= mint) { out_pvalues[i]=0.0; }
            else
            {
				if (i > 0)
                {
                    if (in_RTs[i] <= in_RTs[i-1]) Rf_error("Must call Rfastdm.c with increasing t values (in_RTs[%d]=%g, in_RTs[i-1]=%g.", i, in_RTs[i], in_RTs[i-1]);
                }

				F  = F_get_F(fc, in_RTs[i]);
				nz = F_get_N(fc);
				j  = (int) nz*scaled;
				out_pvalues[i] = p[0]-F[j];
			}
		}
	}

	F_delete(fc);
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// rfastdm - uses construct-samples method
//
void rfastdm (int *in_numvalues, double *in_params, double *in_precision, double *out_RTs, double *out_bounds)
{
	struct F_calculator *fc;
	double  *F;
	double  **Fs,
            Fs_min,
            Fs_max;
	double  t_min, t_max, dt;

	int  random_flag = 0;
	int  i, j, N, s_size, s_count;

    // Will be deterministic if this isn't set
    random_flag = 1;

#ifdef FASTDM_DEBUG
    Rprintf ("In rfastdm.\n");
#endif

    _setup (*in_numvalues, in_params, *in_precision);

    s_size = g_num_values;
    if ((s_size < 1) || (s_size > MAX_VALUES))
    {
        Rf_error ("Number of samples requested exceeds maximum of %d\n", MAX_VALUES);
    }

	s_count = 1; // TODO: NOT SURE WHAT THE DIFFERENCE IS HERE
    if ((s_count < 1) || (s_count > MAX_VALUES))
    {
        Rf_error ("Number of samples requested exceeds maximum of %d\n", MAX_VALUES);
    }

	if ((s_count > 1) && (!random_flag))
    {
        Rf_error ("Error: generating multiple samples but -r not given\n");
	}

#ifdef SHOW_INFO
	Rprintf("Generating %d %s sample%s of size %d\n", s_count, random_flag?"random":"deterministic", s_count!=1?"s":"", s_size);
#endif

	// get the F-values for the samples
	Fs = xnew(double *, s_count);
	Fs_min = 1;
	Fs_max = 0;

	if (random_flag)
	{
#ifdef USE_R_PRNG
        GetRNGstate();  // Retrieve R's PRNG state
#else
        unsigned long c = clock();
  	    init_noise(c); 
#endif
	}

    // Create Fs[][], an s_count X s_size matrix of random or sequential numbers,
    //   and set Fs_min and Fs_max, the min and max produced
	for (j=0; j<s_count; ++j)
    {
		Fs[j] = xnew(double, s_size);

		if (random_flag)
        {
			for (i=0; i<s_size; ++i)
            {
#ifdef USE_R_PRNG
				Fs[j][i] = unif_rand ();
#else
				Fs[j][i] = jvrand_real2();
#endif

				
				
				if (Fs[j][i] > Fs_max)  Fs_max = Fs[j][i];
				if (Fs[j][i] < Fs_min)  Fs_min = Fs[j][i];
			}
		}
        else
        {
			for (i=0; i<s_size; ++i)
            {
				Fs[j][i] = (i+0.5)/s_size;
			}
			Fs_min = Fs[j][0];
			Fs_max = Fs[j][s_size-1];
		}
	}

    // Create new F calculator
	fc = F_new (g_params);

	// get the required t-range
	t_max = 0.5;
	F_start (fc, b_upper);
	while (F_get_val (fc, t_max, g_zr*g_params[p_a]) < Fs_max)	t_max += 0.1;

	t_min = -0.5;
	F_start (fc, b_lower);
	while (F_get_val (fc, -t_min, g_zr*g_params[p_a]) > Fs_min)	t_min -= 0.1;

	// get a table of F-values
	N = (int)((t_max-t_min)/0.001 + 0.5);
	dt = (t_max-t_min)/N;
	F = xnew(double, N+1);

	F_start (fc, b_upper);
	for (i=0; i<=N; ++i)
    {
		double  t = t_min+i*dt;
		if (t < 0)  continue;
		F[i] = F_get_val (fc, t, g_zr*g_params[p_a]);
	}
	F_start (fc, b_lower);
	for (i=N; i>=0; --i)
    {
		double  t = -(t_min+i*dt);
		if (t < 0)  continue;
		F[i] = F_get_val (fc, t, g_zr*g_params[p_a]);
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

	for (j=0; j<s_count; ++j)
    {
		for (i=0; i<s_size; ++i)
        {
			double  y = Fs[j][i];
			int  k = find_slot(y, F, 0, N);
			double  t = t_min + (k + (y-F[k])/(F[k+1]-F[k]))*dt;

			assert (F[k]<=y && y<=F[k+1]);

            // MG: GUTS
			//fprintf (out, "%3d %6.3f\n", t >= 0, fabs(t));
            out_bounds[(j*s_size)+i] = t >= 0;
            out_RTs[(j*s_size)+i] = fabs(t);
		}
	}

    // Free up allocated objects
	xfree(F);
	for (j=0; j<s_count; ++j) xfree(Fs[j]);
	xfree(Fs);
	
#ifdef USE_R_PRNG
	if (random_flag)
	{
        PutRNGstate();  // Must put back R's PRNG state
    }
#endif

}

// Support functions for rfastdm ()
static int compare_doubles (const void *a, const void *b)
{
	const double *da = a;
	const double *db = b;

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


// Alternative CDF method adapted from code by Scott Brown 
// int fastdmcdf(double *parain,            -- Parameter set
//              double *nhi, double *nlo,   -- Number of upper and lower boundary values...?
//              double *thi, double *tlo,   -- RT arrays for upper and lower boundaries
//              double *phi, double *plo,   -- Upper and lower p values
//              double *p,                  -- CDF(t=Inf) ... ?
//              double *precision)
int Scott_fastdmcdf(double *in_params,
              double *nhi, double *nlo,
              double *thi, double *tlo,
              double *phi, double *plo,
              double *p, double *precision)
{
	struct  F_calculator *fc;
	double  mint;

	const double *F;
	int     inhi,
            inlo,

            i,j,nz;
#ifdef FASTDM_DEBUG
    Rprintf ("In Scott_fastdmcdf. Calling _setup ()\n");
#endif

    // PREVIOUSLY: Assume that the "parain" vector has parameters in this order: a,v,t0,sz,sv,st,z
    // BUT NOTE: ADDITIONAL PARAMETERS HAVE BEEN ADDED
    _setup (*nhi, in_params, *precision);

    // Set number of upper and lower boundary values
	inhi = (int) *nhi ;
	inlo = (int) *nlo ;

    // start doing maths
#ifdef FASTDM_DEBUG
    Rprintf ("Creating F object.\n");
#endif

	fc = F_new(g_params);

    double scaled = g_zr /g_params[p_a];

	// Get the value of CDF(t=Inf)  (store in p[0])
	F_start(fc, b_upper);
	mint = g_params[p_t0] - 0.5*g_params[p_st0] ; // Min RT
	F    = F_get_F(fc,mint);

	nz   = F_get_N(fc);
	j    = (int) nz*scaled;
	p[0] = F[j];

#ifdef FASTDM_DEBUG
    Rprintf ("Calculating upper boundary\n");
#endif
    // Calc upper boundary
	if (inhi > 0)
    {
		for (i=0; i<inhi; i++)
        {
			if (thi[i]<=mint) { phi[i] = 0.0; }
            else
            {
				if (i > 0)
                {
                    if (thi[i] <= thi[i-1]) Rf_error("Must call Rfastdm.c with increasing t values (thi[%d]=%g, thi[i-1]=%g.", i, thi[i], thi[i-1]);
                }

				F  = F_get_F(fc, thi[i]);
				nz = F_get_N(fc);
				j  = (int) nz*scaled;

				phi[i] = F[j]-p[0];
			}
		}
	}

#ifdef FASTDM_DEBUG
    Rprintf ("Calculating lower boundary\n");
#endif
    // Calc lower boundary
	if (inlo > 0)
    {
		F_start(fc, b_lower);
		for (i=0; i<inlo; i++)
        {
			if (tlo[i]<=mint) { plo[i]=0.0; }
            else
            {
				if (i > 0)
                {
                    if (tlo[i] <= tlo[i-1]) Rf_error("Must call Rfastdm.c with increasing t values (tlo[%d]=%g, tlo[i-1]=%g.", i, tlo[i], tlo[i-1]);
                }

				F  = F_get_F(fc, tlo[i]);
				nz = F_get_N(fc);
				j  = (int) nz*scaled;
				plo[i] = p[0]-F[j];
			}
		}
	}

	F_delete(fc);
	return  0;
}



// Alternative method: pass through to Scott's code
// An R-like version which finds the left-hand area (from 0 to RT)
void pfastdm_alt (double *in_numvalues, double *in_params, double *in_RTs, double *in_precision, double *out_pvalues_upper, double *out_pvalues_lower)
{
    double p[1]; // This is CDF(t=Inf)?

    Scott_fastdmcdf(in_params,
                    in_numvalues, in_numvalues,
                    in_RTs, in_RTs,

                    out_pvalues_upper, out_pvalues_lower,
                    p, in_precision);
#ifdef SHOW_INFO
    Rprintf ("p = %g\n", p);
#endif
}


void params_write(double para[p_count], double zr, double precision, int n)
{
	Rprintf("a = %g\n", para[p_a]);
	Rprintf("zr = %g\n", zr);
	Rprintf("v = %g\n", para[p_v]);
	Rprintf("t0 = %g\n", para[p_t0]);
	Rprintf("d = %g\n", para[p_d]);
	Rprintf("szr = %g\n", para[p_szr]);
	Rprintf("sv = %g\n", para[p_sv]);
	Rprintf("st0 = %g\n", para[p_st0]);
	Rprintf("precision = %g\n", precision);
	if (n > 0) {  Rprintf("n = %d\n", n);}
}

void params_check(double para[p_count], double zr)
{
	if (para[p_a] <= 0)                     { Rf_error ("error: invalid parameter a=%g\n", para[p_a]); }
	if (para[p_szr] < 0 || para[p_szr] > 1) { Rf_error ("error: invalid parameter szr=%g\n", para[p_szr]); }
	if (para[p_st0] < 0)                    { Rf_error ("error: invalid parameter st0=%g\n", para[p_st0]); }
	if (para[p_sv] < 0)                     { Rf_error ("error: invalid parameter sv=%g\n",	para[p_sv]); }
	if (para[p_t0] - fabs(0.5*para[p_d]) - 0.5*para[p_st0] < 0) { Rf_error ("error: invalid parameter combination t0=%g, d=%g, st0=%g\n", para[p_t0], para[p_d], para[p_st0]); }
	if (zr - 0.5*para[p_szr] <= 0)          { Rf_error ("error: invalid parameter combination zr=%g, szr=%g\n", zr, para[p_szr]); }
	if (zr + 0.5*para[p_szr] >= 1)          { Rf_error ("error: invalid parameter combination zr=%g, szr=%g\n",	zr, para[p_szr]); }

#ifdef ADDITIONAL_BOUNDS
    if (para[p_a] > 2)                      { Rf_error ("Additional Bounds error: invalid parameter a=%g\n", para[p_a]); }
    if (para[p_sv] > 2)                     { Rf_error ("Additional Bounds error: invalid parameter sv=%g\n", para[p_sv]); }
	if (para[p_v] < -5 || para[p_v] > 5)    { Rf_error ("Additional Bounds error: invalid parameter v=%g\n", para[p_v]); }
#endif
}

