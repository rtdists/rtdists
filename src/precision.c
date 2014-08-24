/* precision.c - tunable parameter for the precision of CDF and density
 *
 * Copyright (C) 2013  Jochen Voss.
 */

#include <math.h>

#include "fast-dm.h"


double  TUNE_DZ;
double  TUNE_DV;
double  TUNE_DT0;

double  TUNE_PDE_DT_MIN = 1e-6;
double  TUNE_PDE_DT_MAX = 1e-6;
double  TUNE_PDE_DT_SCALE = 0.0;

double  TUNE_INT_T0;
double  TUNE_INT_Z;

int  precision_set = 0;


void
set_precision (double p)
/* Try to achieve an accuracy of approximately 10^{-p} for the CDF.  */
{
	TUNE_PDE_DT_MIN = pow(10, -0.400825*p-1.422813);
	TUNE_PDE_DT_MAX = pow(10, -0.627224*p+0.492689);
	TUNE_PDE_DT_SCALE = pow(10, -1.012677*p+2.261668);
	TUNE_DZ = pow(10, -0.5*p-0.033403);
	TUNE_DV = pow(10, -1.0*p+1.4);
	TUNE_DT0 = pow(10, -0.5*p-0.323859);

	TUNE_INT_T0 = 0.089045 * exp(-1.037580*p);
	TUNE_INT_Z = 0.508061 * exp(-1.022373*p);

	precision_set = 1;
}
