#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */

extern "C" {
extern SEXP tgs_cor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_cor_blas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_cor_graph(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_cor_knn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_cross_cor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_cross_cor_blas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_cross_cor_knn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_dist(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_dist_blas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_finite(SEXP, SEXP);
extern SEXP tgs_graph2cluster(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_graph2cluster_multi_edges(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_graph2cluster_multi_full(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_graph2cluster_multi_hash(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_knn(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tgs_matrix_tapply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

static const R_CallMethodDef CallEntries[] = {
    {"tgs_cor",                       (DL_FUNC) &tgs_cor,                       6},
    {"tgs_cor_blas",                  (DL_FUNC) &tgs_cor_blas,                  6},
    {"tgs_cor_graph",                 (DL_FUNC) &tgs_cor_graph,                 5},
    {"tgs_cor_knn",                   (DL_FUNC) &tgs_cor_knn,                   6},
    {"tgs_cross_cor",                 (DL_FUNC) &tgs_cross_cor,                 7},
    {"tgs_cross_cor_blas",            (DL_FUNC) &tgs_cross_cor_blas,            7},
    {"tgs_cross_cor_knn",             (DL_FUNC) &tgs_cross_cor_knn,             7},
    {"tgs_dist",                      (DL_FUNC) &tgs_dist,                      6},
    {"tgs_dist_blas",                 (DL_FUNC) &tgs_dist_blas,                 6},
    {"tgs_finite",                    (DL_FUNC) &tgs_finite,                    2},
    {"tgs_graph2cluster",             (DL_FUNC) &tgs_graph2cluster,             5},
    {"tgs_graph2cluster_multi_edges", (DL_FUNC) &tgs_graph2cluster_multi_edges, 8},
    {"tgs_graph2cluster_multi_full",  (DL_FUNC) &tgs_graph2cluster_multi_full,  8},
    {"tgs_graph2cluster_multi_hash",  (DL_FUNC) &tgs_graph2cluster_multi_hash,  8},
    {"tgs_knn",                       (DL_FUNC) &tgs_knn,                       5},
    {"tgs_matrix_tapply",             (DL_FUNC) &tgs_matrix_tapply,             6},
    {NULL, NULL, 0}
};

void R_init_tgstat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
