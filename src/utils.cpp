#include <cstdint>
#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
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
        uint64_t len = Rf_xlength(_x);

        if (!Rf_isReal(_x) && !Rf_isInteger(_x))
            verror("\"x\" argument must be numeric or integer");

        for (uint64_t i = 0; i < len; ++i) {
            if ((Rf_isReal(_x) && !R_FINITE(REAL(_x)[i])) || (Rf_isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
                rreturn(Rf_ScalarLogical(false));
        }
    } catch (TGLException &e) {
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(Rf_ScalarLogical(true));
}

}
