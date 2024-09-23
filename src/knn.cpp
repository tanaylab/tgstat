#include <algorithm>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>

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

#include "ProgressReporter.h"
#include "tgstat.h"

extern "C" {

SEXP tgs_knn(SEXP _x, SEXP _knn, SEXP _diag, SEXP _threshold, SEXP _envir)
{
    SEXP answer = R_NilValue;
    uint64_t *res = (uint64_t *)MAP_FAILED;
    uint64_t res_sizeof = 0;

	try {
        TGStat tgstat(_envir);

        if ((!Rf_isReal(_threshold) && !Rf_isInteger(_threshold)) || Rf_xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        if (!Rf_isLogical(_diag) || Rf_xlength(_diag) != 1)
            verror("\"diag\" argument must be a logical value");

        if ((!Rf_isReal(_knn) && !Rf_isInteger(_knn)) || Rf_xlength(_knn) != 1 || Rf_asReal(_knn) < 1 || Rf_asReal(_knn) != (double)Rf_asInteger(_knn))
            verror("\"knn\" argument must be a positive integer");

        double threshold = Rf_asReal(_threshold);
        uint64_t knn = Rf_asInteger(_knn);
        bool diag = Rf_asLogical(_diag);
        bool tidy = false;
        uint64_t num_pairs = 0;
        vector<uint64_t> point2size;
        int *pcol1 = NULL;
        int *pcol2 = NULL;
        double *data;
        auto cmp = [&data](uint64_t idx1, uint64_t idx2) { return data[idx1] > data[idx2] || (data[idx1] == data[idx2] && idx1 < idx2); };

        vector<uint64_t> sorted_rows;     // contains number of rows of _x sorted by col1
        uint64_t num_rows = 0;
        uint64_t num_cols = 0;

        vdebug("Preparing for multitasking...\n");
        if (Rf_isReal(_x) && Rf_xlength(_x) >= 1) {    // x is a matrix
            SEXP rdim = Rf_getAttrib(_x, R_DimSymbol);

            if (!Rf_isInteger(rdim) || Rf_xlength(rdim) != 2)
                verror("Invalid format of \"x\" argument");

            num_rows = Rf_nrows(_x);
            num_cols = Rf_ncols(_x);

            tidy = false;
            data = REAL(_x);

            point2size.resize(num_cols, 0);

            {
                uint64_t idx = 0;
                for (uint64_t icol = 0; icol < num_cols; ++icol) {
                    for (uint64_t irow = 0; irow < num_rows; ++irow) {
                        if ((irow != icol || diag) && R_FINITE(data[idx]) && data[idx] >= threshold) {
                            ++point2size[icol];
                            ++num_pairs;
                        }
                        ++idx;
                    }
                }
            }

            if (num_pairs) {
                res_sizeof = sizeof(uint64_t) * num_pairs;
                res = (uint64_t *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

                if (res == (uint64_t *)MAP_FAILED)
                    verror("Failed to allocate shared memory: %s", strerror(errno));

                {
                    uint64_t idx1 = 0;
                    uint64_t idx2 = 0;

                    for (uint64_t icol = 0; icol < num_cols; ++icol) {
                        for (uint64_t irow = 0; irow < num_rows; ++irow) {
                            if ((irow != icol || diag) && R_FINITE(data[idx1]) && data[idx1] >= threshold)
                                res[idx2++] = idx1;
                            ++idx1;
                        }
                    }

                }
            }
        } else if (Rf_isVector(_x) && Rf_xlength(_x) >= 3) {
            tidy = true;

            SEXP rcol1 = VECTOR_ELT(_x, 0);
            SEXP rcol2 = VECTOR_ELT(_x, 1);
            SEXP rval = VECTOR_ELT(_x, 2);

            num_rows = Rf_xlength(rcol1);

            if ((!Rf_isInteger(rcol1) && !Rf_isFactor(rcol1)) || (!Rf_isInteger(rcol2) && !Rf_isFactor(rcol1)) || (unsigned)Rf_xlength(rcol2) != num_rows ||
                (!Rf_isInteger(rval) && !Rf_isReal(rval)) || (unsigned)Rf_xlength(rval) != num_rows)
                verror("Invalid format of \"x\" argument");

            bool sort_needed = false;

            pcol1 = INTEGER(rcol1);
            pcol2 = INTEGER(rcol2);
            data = REAL(rval);

            auto cmp_rows = [&pcol1, &pcol2](uint64_t idx1, uint64_t idx2) {
                                return pcol1[idx1] < pcol1[idx2] || (pcol1[idx1] == pcol1[idx2] && pcol2[idx1] < pcol2[idx2]);
                            };

            sorted_rows.resize(num_rows);
            for (uint64_t i = 0; i < num_rows; ++i) {
                if (i && pcol1[i - 1] > pcol1[i])
                    sort_needed = true;

                sorted_rows[i] = i;
            }

            if (sort_needed)
                sort(sorted_rows.begin(), sorted_rows.end(), cmp_rows);

            for (auto irow = sorted_rows.begin(); irow != sorted_rows.end(); ++irow) {
                if (irow != sorted_rows.begin() && pcol1[*irow] == pcol1[*(irow - 1)] && pcol2[*irow] == pcol2[*(irow - 1)])
                    verror("\"x\" contains a duplicated column pair [%d, %d]", pcol1[*irow], pcol2[*irow]);

                if (irow == sorted_rows.begin() || pcol1[*irow] != pcol1[*(irow - 1)])
                    point2size.push_back(0);

                if (R_FINITE(data[*irow]) && data[*irow] >= threshold) {
                    ++point2size.back();
                    ++num_pairs;
                }
            }

            if (num_pairs) {
                res_sizeof = sizeof(uint64_t) * num_pairs;
                res = (uint64_t *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

                if (res == (uint64_t *)MAP_FAILED)
                    verror("Failed to allocate shared memory: %s", strerror(errno));

                uint64_t idx = 0;
                for (auto irow = sorted_rows.begin(); irow != sorted_rows.end(); ++irow) {
                    if (R_FINITE(data[*irow]) && data[*irow] >= threshold)
                        res[idx++] = *irow;
                }
            }
        } else {
            verror("Invalid format of \"x\" argument");
        }

        int num_processes = (int)min(num_pairs / 1000000, (uint64_t)g_tgstat->num_processes());

        if (num_pairs)
            num_processes = max(num_processes, 1);

        ProgressReporter progress;
        progress.init(point2size.size(), 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                uint64_t spoint = point2size.size() * (iprocess / (double)num_processes);
                uint64_t epoint = point2size.size() * ((iprocess + 1) / (double)num_processes);
                uint64_t idx = 0;

                for (uint64_t i = 0; i < spoint; ++i)
                    idx += point2size[i];

                for (uint64_t ipoint = spoint; ipoint < epoint; ++ipoint) {
                    if (point2size[ipoint] > knn)
                        partial_sort(res + idx, res + idx + knn, res + idx + point2size[ipoint], cmp);
                    else
                        sort(res + idx, res + idx + point2size[ipoint], cmp);

                    idx += point2size[ipoint];
                    TGStat::itr_idx(ipoint - spoint + 1);
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        enum { COL1, COL2, VAL, RANK, NUM_COLS };
        const char *COL_NAMES[NUM_COLS] = { "col1", "col2", "val", "rank" };

        SEXP rcol1, rcol2, rval, rrank, rrownames, rcolnames, rold_colnames;
        uint64_t answer_size = 0;

        vdebug("PACKING\n");
        for (auto point_size : point2size)
            answer_size += min(knn, point_size);

        if (tidy)
            rold_colnames = Rf_getAttrib(VECTOR_ELT(_x, 0), R_LevelsSymbol);
        else {
            SEXP rold_dimnames = Rf_getAttrib(_x, R_DimNamesSymbol);
            rold_colnames = !Rf_isNull(rold_dimnames) && Rf_xlength(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;
        }

        rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));
        rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rval = RSaneAllocVector(REALSXP, answer_size));
        rprotect(rrank = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, Rf_mkChar(COL_NAMES[i]));

        if (tidy) {
            rold_colnames = Rf_getAttrib(_x, R_NamesSymbol);
            if (Rf_isString(rold_colnames) && Rf_xlength(rold_colnames) >= 2)
                SET_STRING_ELT(rcolnames, VAL, STRING_ELT(rold_colnames, 2));
        }

        uint64_t row = 0;
        uint64_t idx = 0;
        int rank = 1;

        if (tidy) {
            for (uint64_t ipoint = 0; ipoint < point2size.size(); ++ipoint) {
                uint64_t point_knn = min(knn, point2size[ipoint]);

                for (uint64_t i = 0; i < point_knn; ++i) {
                    INTEGER(rcol1)[row] = pcol1[res[idx + i]];
                    INTEGER(rcol2)[row] = pcol2[res[idx + i]];
                    REAL(rval)[row] = data[res[idx + i]];

                    if (row && INTEGER(rcol1)[row - 1] != INTEGER(rcol1)[row])
                        rank = 1;

                    INTEGER(rrank)[row] = rank;
                    INTEGER(rrownames)[row] = row + 1;
                    ++rank;
                    ++row;
                }
                idx += point2size[ipoint];
            }
        } else {
            for (uint64_t icol = 0; icol < num_cols; ++icol) {
                uint64_t point_knn = min(knn, point2size[icol]);

                for (uint64_t i = 0; i < point_knn; ++i) {
                    INTEGER(rcol1)[row] = icol + 1;
                    INTEGER(rcol2)[row] = res[idx + i] % num_rows + 1;
                    REAL(rval)[row] = data[res[idx + i]];

                    if (row && INTEGER(rcol1)[row - 1] != INTEGER(rcol1)[row])
                        rank = 1;

                    INTEGER(rrank)[row] = rank;
                    INTEGER(rrownames)[row] = row + 1;
                    ++rank;
                    ++row;
                }
                idx += point2size[icol];
            }
        }

        if (rold_colnames != R_NilValue) {
            Rf_setAttrib(rcol1, R_LevelsSymbol, rold_colnames);
            Rf_setAttrib(rcol1, R_ClassSymbol, Rf_mkString("factor"));
            Rf_setAttrib(rcol2, R_LevelsSymbol, rold_colnames);
            Rf_setAttrib(rcol2, R_ClassSymbol, Rf_mkString("factor"));
        }

        SET_VECTOR_ELT(answer, COL1, rcol1);
        SET_VECTOR_ELT(answer, COL2, rcol2);
        SET_VECTOR_ELT(answer, VAL, rval);
        SET_VECTOR_ELT(answer, RANK, rrank);

        Rf_setAttrib(answer, R_NamesSymbol, rcolnames);
        Rf_setAttrib(answer, R_ClassSymbol, Rf_mkString("data.frame"));
        Rf_setAttrib(answer, R_RowNamesSymbol, rrownames);

        vdebug("END\n");
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (uint64_t *)MAP_FAILED) {
            munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
            res = (uint64_t *)MAP_FAILED;
        }
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (uint64_t *)MAP_FAILED) {
        munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
        res = (uint64_t *)MAP_FAILED;
    }
	rreturn(answer);
}

}
