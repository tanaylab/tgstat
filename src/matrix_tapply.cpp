#include <iterator>

#include <R.h>
#include <Rinternals.h>

#include "tgstat.h"

extern "C" {

SEXP matrix_tapply(SEXP _x, SEXP _index, SEXP _fn, SEXP _envir)
{
    SEXP answer = R_NilValue;
//  double *res = (double *)MAP_FAILED;
//  size_t res_sizeof = 0;

	try {
        TGStat tgstat(_envir);

        SEXP _rdims;

        if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1 || !isInteger(_rdims = getAttrib(_x, R_DimSymbol)) || xlength(_rdims) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        int num_rows = nrows(_x);
        int num_cols = ncols(_x);

        if (!isInteger(_index) || xlength(_index) < 1)
            verror("\"index\" argument must be a factor");

        if (xlength(_index) != num_cols)
            verror("Arguments \"x\" and \"index\" must have same length");

        if (!isFunction(_fn) || xlength(_fn))
            verror("\"fn\" argument must be a function");

        SEXP rcall;
        SEXP rgroup;
        int cols[] = { 2, 3, 5, 6, 8 };
        int group_size = sizeof(cols) / sizeof(cols[0]);

        if (isReal(_x)) {
            rprotect(rgroup = allocVector(REALSXP, group_size));
            for (auto col : cols)
                REAL(rgroup)[0] = REAL(_x)[num_rows * col];
        } else {
            rprotect(rgroup = allocVector(INTSXP, group_size));
            for (auto col : cols)
                INTEGER(rgroup)[0] = INTEGER(_x)[num_rows * col];
        }

        rprotect(rcall = lang2(_fn, R_NilValue));
        SETCADR(rcall, rgroup);

        return eval_in_R(rcall, _envir);
        
    } catch (TGLException &e) {
//      if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
//          munmap(res, res_sizeof);
//          res = (double *)MAP_FAILED;
//      }
        rerror("%s", e.msg());
    }

    rreturn(answer);
}

}
