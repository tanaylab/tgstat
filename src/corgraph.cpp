#include <algorithm>
#include <limits>
#include <unordered_map>

#include "HashFunc.h"

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

SEXP tgs_cor_graph(SEXP _ranks, SEXP _knn, SEXP _k_expand, SEXP _k_beta, SEXP _envir)
{
	try {
        TGStat tgstat(_envir);

        int *pcol1;
        int *pcol2;
        int *prank;
        uint64_t num_ranks;

        if ((!isReal(_k_beta) && !isInteger(_k_beta)) || xlength(_k_beta) != 1)
            verror("\"k_beta\" argument must be a numeric value");

        {
            enum { COL1, COL2, COR, RANK, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "cor", "rank" };

            SEXP rnames = getAttrib(_ranks, R_NamesSymbol);

    		if (!isVector(_ranks) || xlength(_ranks) != NUM_COLS || xlength(rnames) != NUM_COLS ||
                strcmp(CHAR(STRING_ELT(rnames, COL1)), COL_NAMES[COL1]) || (!isInteger(VECTOR_ELT(_ranks, COL1)) && !isFactor(VECTOR_ELT(_ranks, COL1))) ||
                strcmp(CHAR(STRING_ELT(rnames, COL2)), COL_NAMES[COL2]) || (!isInteger(VECTOR_ELT(_ranks, COL2)) && !isFactor(VECTOR_ELT(_ranks, COL2))) ||
                xlength(VECTOR_ELT(_ranks, COL2)) != xlength(VECTOR_ELT(_ranks, COL1)) ||
                !isReal(VECTOR_ELT(_ranks, COR)) || xlength(VECTOR_ELT(_ranks, COR)) != xlength(VECTOR_ELT(_ranks, COL1)) ||
                strcmp(CHAR(STRING_ELT(rnames, RANK)), COL_NAMES[RANK]) || !isInteger(VECTOR_ELT(_ranks, RANK)) || xlength(VECTOR_ELT(_ranks, RANK)) != xlength(VECTOR_ELT(_ranks, COL1)))
    			verror("\"ranks\" argument must be in the format that is returned by tgs_knn function");

            pcol1 = INTEGER(VECTOR_ELT(_ranks, COL1));
            pcol2 = INTEGER(VECTOR_ELT(_ranks, COL2));
            prank = INTEGER(VECTOR_ELT(_ranks, RANK));
            num_ranks = xlength(VECTOR_ELT(_ranks, RANK));
        }

        if (!isNull(_knn) && ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1))
            verror("\"knn\" argument must be a numeric value");

        if (!isNull(_k_expand) && ((!isReal(_k_expand) && !isInteger(_k_expand)) || xlength(_k_expand) != 1))
            verror("\"k_expand\" argument must be a numeric value");

        double knn_d = isNull(_knn) ? 0 : asReal(_knn);
        double k_expand = asReal(_k_expand);
        double k_beta = asReal(_k_beta);

        if (knn_d < 1)
            verror("\"knn\" argument must be a positive integer");

        if (k_expand <= 0)
            verror("\"k_expand\" argument must be a positive number");

        if (k_beta <= 0)
            verror("\"k_beta\" argument must be a positive number");

        vdebug("Building the graph\n");

        uint64_t knn = (uint64_t)knn_d;
        unsigned num_points = 0;
        unordered_map<pair<unsigned, unsigned>, uint64_t> ij2weight;
        unordered_map<pair<unsigned, unsigned>, uint64_t> ij2rank;
        uint64_t max_weight = knn * knn * k_expand;

        vdebug("Reading ranks\n");
        ij2rank.reserve(num_ranks);
        for (uint64_t i = 0; i < num_ranks; ++i)
            ij2rank[{pcol1[i], pcol2[i]}] = prank[i];

        vdebug("Building edges weights\n");
        ij2weight.reserve(num_ranks * 2);
        for (const auto &r : ij2rank) {
            const auto &ij = r.first;
            num_points = max(num_points, ij.first + 1);
            num_points = max(num_points, ij.second + 1);
            if (ij.first < ij.second) {
                auto itr = ij2rank.find({ij.second, ij.first});
                if (itr != ij2rank.end()) {
                    uint64_t weight = r.second * itr->second;    // weight = rank[i,j] * rank[j,i]
                    if (weight <= max_weight)
                        ij2weight[ij] = ij2weight[{ij.second, ij.first}] = weight;
                }
            }
        }

        {
            decltype(ij2rank) cleaner;
            ij2rank.swap(cleaner);     // force ij2rank to release memory
        }

        struct Edge {
            unsigned node;
            uint64_t weight;
            Edge(unsigned _node, unsigned _weight) : node(_node), weight(_weight) {}
            bool operator<(const Edge &o) const { return weight < o.weight || (weight == o.weight && node < o.node); }
        };

        vdebug("Filter out by incoming edges\n");

        // leave max k_beta * knn incoming edges
        vector<vector<Edge>> incoming(num_points);
        for (const auto &w : ij2weight) {
            unsigned i = w.first.first;
            unsigned j = w.first.second;
            incoming[j].push_back(Edge(i, w.second));
        }

        for (auto iedges = incoming.begin(); iedges < incoming.end(); ++iedges) {
            if (iedges->size() > k_beta * knn) {
                partial_sort(iedges->begin(), iedges->begin() + k_beta * knn, iedges->end());
                for (auto iedge = iedges->begin() + k_beta * knn; iedge < iedges->end(); ++iedge)
                    ij2weight.erase({iedge->node, iedges - incoming.begin()});
            }
        }

        {
            decltype(incoming) cleaner;
            incoming.swap(cleaner);     // force incoming to release memory
        }

        vdebug("Filter out by outgoing edges\n");

        // leave max knn outgoing edges and pack the answer
        vector<vector<Edge>> outgoing(num_points);

        for (const auto &w : ij2weight) {
            unsigned i = w.first.first;
            unsigned j = w.first.second;
            outgoing[i].push_back(Edge(j, w.second));
        }

        {
            decltype(ij2weight) cleaner;
            ij2weight.swap(cleaner);     // force ij2weight to release memory
        }

        vdebug("PACKING\n");

        enum { COL1, COL2, WEIGHT, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "weight" };

        uint64_t answer_size = 0;
        SEXP ranswer, rcol1, rcol2, rweight, rrownames, rcolnames, rlevels;

        for (const auto &edges : outgoing)
            answer_size += min((uint64_t)edges.size(), knn);

        rprotect(ranswer = RSaneAllocVector(VECSXP, NUM_COLS));
        rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rweight = RSaneAllocVector(REALSXP, answer_size));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

        rlevels = getAttrib(VECTOR_ELT(_ranks, 0), R_LevelsSymbol);
        if (rlevels != R_NilValue) {
            setAttrib(rcol1, R_LevelsSymbol, rlevels);
            setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
        }

        rlevels = getAttrib(VECTOR_ELT(_ranks, 1), R_LevelsSymbol);
        if (rlevels != R_NilValue) {
            setAttrib(rcol2, R_LevelsSymbol, rlevels);
            setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
        }

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        uint64_t idx = 0;
        for (auto iedges = outgoing.begin(); iedges < outgoing.end(); ++iedges) {
            double rank = 0;

            if (iedges->size() <= knn)
                sort(iedges->begin(), iedges->end());
            else
                partial_sort(iedges->begin(), iedges->begin() + knn, iedges->end());

            auto iedge_end = iedges->size() <= knn ? iedges->end() : iedges->begin() + knn;

            for (auto iedge = iedges->begin(); iedge < iedge_end; ++iedge) {
                int i = iedges - outgoing.begin();
                int j = iedge->node;
                INTEGER(rcol1)[idx] = i;
                INTEGER(rcol2)[idx] = j;
                REAL(rweight)[idx] = 1. - rank / knn;
                INTEGER(rrownames)[idx] = idx + 1;
                ++idx;
                ++rank;
            }
        }

        vdebug("END\n");

        SET_VECTOR_ELT(ranswer, COL1, rcol1);
        SET_VECTOR_ELT(ranswer, COL2, rcol2);
        SET_VECTOR_ELT(ranswer, WEIGHT, rweight);

        setAttrib(ranswer, R_NamesSymbol, rcolnames);
        setAttrib(ranswer, R_ClassSymbol, mkString("data.frame"));
        setAttrib(ranswer, R_RowNamesSymbol, rrownames);

        rreturn(ranswer);
    } catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(R_NilValue);
}

}

