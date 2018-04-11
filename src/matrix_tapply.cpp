#include <R.h>
#include <Rinternals.h>

#include "tgstat.h"

extern "C" {

SEXP matrix_tapply(SEXP _x, SEXP _index, SEXP _fun, SEXP _envir)
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

        if (!isInteger(_index) || xlength(_index) <= 1)
            verror("\"index\" argument must be a factor");

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
