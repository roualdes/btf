#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP btf_dexp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP btf_gdp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP btf_tf_approx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"btf_dexp",      (DL_FUNC) &btf_dexp,      7},
  {"btf_gdp",       (DL_FUNC) &btf_gdp,       7},
  {"btf_tf_approx", (DL_FUNC) &btf_tf_approx, 7},
  {NULL, NULL, 0}
};

void R_init_btf(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
