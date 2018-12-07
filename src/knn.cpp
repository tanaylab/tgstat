#include <algorithm>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>

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
    size_t *res = (size_t *)MAP_FAILED;
    size_t res_sizeof = 0;

	try {
        TGStat tgstat(_envir);

        if ((!isReal(_threshold) && !isInteger(_threshold)) || xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        if (!isLogical(_diag) || xlength(_diag) != 1)
            verror("\"diag\" argument must be a logical value");

        if ((!isReal(_knn) && !isInteger(_knn)) || xlength(_knn) != 1 || asReal(_knn) < 1 || asReal(_knn) != (double)asInteger(_knn))
            verror("\"knn\" argument must be a positive integer");

        double threshold = asReal(_threshold);
        size_t knn = asInteger(_knn);
        bool diag = asLogical(_diag);
        bool tidy;
        size_t num_pairs = 0;
        vector<size_t> point2size;
        int *pcol1 = NULL;
        int *pcol2 = NULL;
        double *data;
        auto cmp = [&data](size_t idx1, size_t idx2) { return data[idx1] > data[idx2] || (data[idx1] == data[idx2] && idx1 < idx2); };

        vector<size_t> sorted_rows;     // contains number of rows of _x sorted by col1
        size_t num_rows;
        size_t num_cols;

        vdebug("Preparing for multitasking...\n");
        if (isReal(_x) && xlength(_x) >= 1) {    // x is a matrix
            SEXP rdim = getAttrib(_x, R_DimSymbol);

            if (!isInteger(rdim) || xlength(rdim) != 2)
                verror("Invalid format of \"x\" argument");

            num_rows = nrows(_x);
            num_cols = ncols(_x);

            tidy = false;
            data = REAL(_x);

            point2size.resize(num_cols, 0);

            {
                size_t idx = 0;
                for (size_t icol = 0; icol < num_cols; ++icol) {
                    for (size_t irow = 0; irow < num_rows; ++irow) {
                        if ((irow != icol || diag) && R_FINITE(data[idx]) && data[idx] >= threshold) {
                            ++point2size[icol];
                            ++num_pairs;
                        }
                        ++idx;
                    }
                }
            }

            if (num_pairs) {
                res_sizeof = sizeof(size_t) * num_pairs;
                res = (size_t *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

                if (res == (size_t *)MAP_FAILED)
                    verror("Failed to allocate shared memory: %s", strerror(errno));

                {
                    size_t idx1 = 0;
                    size_t idx2 = 0;

                    for (size_t icol = 0; icol < num_cols; ++icol) {
                        for (size_t irow = 0; irow < num_rows; ++irow) {
                            if ((irow != icol || diag) && R_FINITE(data[idx1]) && data[idx1] >= threshold)
                                res[idx2++] = idx1;
                            ++idx1;
                        }
                    }

                }
            }
        } else if (isVector(_x) && xlength(_x) >= 3) {
            tidy = true;

            SEXP rcol1 = VECTOR_ELT(_x, 0);
            SEXP rcol2 = VECTOR_ELT(_x, 1);
            SEXP rval = VECTOR_ELT(_x, 2);

            num_rows = xlength(rcol1);

            if ((!isInteger(rcol1) && !isFactor(rcol1)) || (!isInteger(rcol2) && !isFactor(rcol1)) || xlength(rcol2) != num_rows ||
                (!isInteger(rval) && !isReal(rval)) || xlength(rval) != num_rows)
                verror("Invalid format of \"x\" argument");

            bool sort_needed = false;

            pcol1 = INTEGER(rcol1);
            pcol2 = INTEGER(rcol2);
            data = REAL(rval);

            auto cmp_rows = [&pcol1, &pcol2](size_t idx1, size_t idx2) {
                                return pcol1[idx1] < pcol1[idx2] || (pcol1[idx1] == pcol1[idx2] && pcol2[idx1] < pcol2[idx2]);
                            };

            sorted_rows.resize(num_rows);
            for (size_t i = 0; i < num_rows; ++i) {
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
                res_sizeof = sizeof(size_t) * num_pairs;
                res = (size_t *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

                if (res == (size_t *)MAP_FAILED)
                    verror("Failed to allocate shared memory: %s", strerror(errno));

                size_t idx = 0;
                for (auto irow = sorted_rows.begin(); irow != sorted_rows.end(); ++irow) {
                    if (R_FINITE(data[*irow]) && data[*irow] >= threshold)
                        res[idx++] = *irow;
                }
            }
        } else
            verror("Invalid format of \"x\" argument");

        int num_processes = (int)min(num_pairs / 1000000, (size_t)g_tgstat->num_processes());

        if (num_pairs)
            num_processes = max(num_processes, 1);

        ProgressReporter progress;
        progress.init(point2size.size(), 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                size_t spoint = point2size.size() * (iprocess / (double)num_processes);
                size_t epoint = point2size.size() * ((iprocess + 1) / (double)num_processes);
                size_t idx = 0;

                for (size_t i = 0; i < spoint; ++i)
                    idx += point2size[i];

                for (size_t ipoint = spoint; ipoint < epoint; ++ipoint) {
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
        size_t answer_size = 0;

        vdebug("PACKING\n");
        for (auto point_size : point2size)
            answer_size += min(knn, point_size);

        if (tidy)
            rold_colnames = getAttrib(VECTOR_ELT(_x, 0), R_LevelsSymbol);
        else {
            SEXP rold_dimnames = getAttrib(_x, R_DimNamesSymbol);
            rold_colnames = !isNull(rold_dimnames) && xlength(rold_dimnames) == 2 ? VECTOR_ELT(rold_dimnames, 1) : R_NilValue;
        }

        rprotect(answer = RSaneAllocVector(VECSXP, NUM_COLS));
        rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rval = RSaneAllocVector(REALSXP, answer_size));
        rprotect(rrank = RSaneAllocVector(INTSXP, answer_size));
        rprotect(rcolnames = RSaneAllocVector(STRSXP, NUM_COLS));
        rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

        for (int i = 0; i < NUM_COLS; i++)
            SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

        if (tidy) {
            rold_colnames = getAttrib(_x, R_NamesSymbol);
            if (isString(rold_colnames) && xlength(rold_colnames) >= 2)
                SET_STRING_ELT(rcolnames, VAL, STRING_ELT(rold_colnames, 2));
        }

        size_t row = 0;
        size_t idx = 0;
        int rank = 1;

        if (tidy) {
            for (size_t ipoint = 0; ipoint < point2size.size(); ++ipoint) {
                size_t point_knn = min(knn, point2size[ipoint]);

                for (size_t i = 0; i < point_knn; ++i) {
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
            for (size_t icol = 0; icol < num_cols; ++icol) {
                size_t point_knn = min(knn, point2size[icol]);

                for (size_t i = 0; i < point_knn; ++i) {
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
            setAttrib(rcol1, R_LevelsSymbol, rold_colnames);
            setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
            setAttrib(rcol2, R_LevelsSymbol, rold_colnames);
            setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
        }

        SET_VECTOR_ELT(answer, COL1, rcol1);
        SET_VECTOR_ELT(answer, COL2, rcol2);
        SET_VECTOR_ELT(answer, VAL, rval);
        SET_VECTOR_ELT(answer, RANK, rrank);

        setAttrib(answer, R_NamesSymbol, rcolnames);
        setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
        setAttrib(answer, R_RowNamesSymbol, rrownames);

        vdebug("END\n");
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (size_t *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (size_t *)MAP_FAILED;
        }
		rerror("%s", e.msg());
	} catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (size_t *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (size_t *)MAP_FAILED;
    }
	rreturn(answer);
}

}
