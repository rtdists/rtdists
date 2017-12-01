#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rtdists_d_fastdm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rtdists_p_fastdm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rtdists_p_precise_fastdm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rtdists_r_fastdm(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_rtdists_d_fastdm",         (DL_FUNC) &_rtdists_d_fastdm,         5},
    {"_rtdists_p_fastdm",         (DL_FUNC) &_rtdists_p_fastdm,         5},
    {"_rtdists_p_precise_fastdm", (DL_FUNC) &_rtdists_p_precise_fastdm, 5},
    {"_rtdists_r_fastdm",         (DL_FUNC) &_rtdists_r_fastdm,         4},
    {NULL, NULL, 0}
};

void R_init_rtdists(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

