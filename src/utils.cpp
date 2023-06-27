#include <cstdint>
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
        uint64_t len = xlength(_x);

        if (!isReal(_x) && !isInteger(_x))
            verror("\"x\" argument must be numeric or integer");

        for (uint64_t i = 0; i < len; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
                rreturn(ScalarLogical(false));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(ScalarLogical(true));
}

}
