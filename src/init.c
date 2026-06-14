#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP wrap__scclr_overdispersion(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP wrap__scclr_normalize(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP wrap__scclr_pca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP wrap__scclr_normalize_pca(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"wrap__scclr_overdispersion", (DL_FUNC) &wrap__scclr_overdispersion, 5},
    {"wrap__scclr_normalize", (DL_FUNC) &wrap__scclr_normalize, 10},
    {"wrap__scclr_pca", (DL_FUNC) &wrap__scclr_pca, 11},
    {"wrap__scclr_normalize_pca", (DL_FUNC) &wrap__scclr_normalize_pca, 13},
    {NULL, NULL, 0}
};

void R_init_scclrR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
