#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void SBDE(double *, double *, int *, double *, double *, int *, double *, double *, double *, double *, int *, int *, double *, int *, double *, double *, double *, int *);
extern void DEV(double *, double *, int *, double *, double *, int *, double *, double *, double *, double *, double *);
extern void PRED(double *, double *, double *, int *, double *, double *, double *);

static const R_CMethodDef CEntries[] = {
    {"SBDE", (DL_FUNC) &SBDE, 18},
    {"DEV",  (DL_FUNC) &DEV,  11},
    {"PRED", (DL_FUNC) &PRED,  7},
    {NULL, NULL, 0}
};

void R_init_sbde(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
