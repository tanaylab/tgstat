#include <R.h>
#include <Rinternals.h>

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "tgstat.h"

extern "C" {

SEXP tgs_finite(SEXP _x, SEXP _envir)
{
	try {
        TGStat tgstat(_envir);
        size_t len = xlength(_x);

//char buf[1000];
//sprintf(buf, "set.seed(%d)", asInteger(_x));
//sprintf(buf, "set.seed(%d)", (int)(unif_rand() * 10000));
//run_in_R(buf, _envir);
//GetRNGstate();
//for (int i = 0; i < 10; ++i)
//printf("%g ", unif_rand());
//PutRNGstate();
//printf("\n");
//return R_NilValue;

        if (!isReal(_x) && !isInteger(_x))
            verror("\"x\" argument must be numeric or integer");

        for (size_t i = 0; i < len; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
                rreturn(ScalarLogical(false));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    }

    rreturn(ScalarLogical(true));
}

}
