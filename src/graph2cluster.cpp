#include <algorithm>
#include <limits>

#include <R.h>
#include <Rinternals.h>

#ifdef length
#undef length
#endif
#ifdef error
#undef error
#endif

#include "ProgressReporter.h"
#include "tgstat.h"

extern "C" {

SEXP tgs_graph2cluster(SEXP _graph, SEXP _envir)
{
	try {
        TGStat tgstat(_envir);

        int *pcol1;
        int *pcol2;
        double *pweight;
        size_t num_edges;

        {
            enum { COL1, COL2, WEIGHT, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "weight" };

            SEXP rnames = getAttrib(_graph, R_NamesSymbol);

            if (!isVector(_graph) || xlength(_graph) != NUM_COLS || xlength(rnames) != NUM_COLS ||
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || !isInteger(VECTOR_ELT(_graph, COL1)) && !isFactor(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || !isInteger(VECTOR_ELT(_graph, COL2)) && !isFactor(VECTOR_ELT(_graph, COL2)) ||
                xlength(VECTOR_ELT(_graph, COL2)) != xlength(VECTOR_ELT(_graph, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, WEIGHT)), COL_NAMES[WEIGHT]) || !isReal(VECTOR_ELT(_graph, WEIGHT)) || xlength(VECTOR_ELT(_graph, WEIGHT)) != xlength(VECTOR_ELT(_graph, COL1)))
                verror("\"graph\" argument must be in the format that is returned by tgs_cor_graph function");

            pcol1 = INTEGER(VECTOR_ELT(_graph, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_graph, COL2));
            pweight = REAL(VECTOR_ELT(_graph, WEIGHT));
            num_edges = xlength(VECTOR_ELT(_graph, COL1));
        }

        SEXP ranswer;

        rreturn(ranswer);
    } catch (TGLException &e) {
		rerror("%s", e.msg());
	}

    rreturn(R_NilValue);
}

}


