
#include <algorithm>
#include <limits>
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

SEXP tgs_dist(SEXP _x, SEXP _attrs, SEXP _tidy, SEXP _threshold, SEXP _rrownames, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    size_t res_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if (!isReal(_x) && !isInteger(_x) || Rf_length(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!isLogical(_tidy) || Rf_length(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if (!isReal(_threshold) && !isInteger(_threshold) || Rf_length(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = getAttrib(_x, R_DimSymbol);

        if (!isInteger(rdim) || Rf_length(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));
        int num_rows = nrows(_x);
        int num_cols = ncols(_x);

        if (num_rows <= 1 || num_cols <= 1)
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_vals = num_rows * num_cols;
        vector<double> vals;

        vals.reserve(num_vals);

        for (size_t i = 0; i < num_vals; ++i) {
            if (isReal(_x) && !R_FINITE(REAL(_x)[i]) || isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)
                vals.push_back(numeric_limits<double>::quiet_NaN());
            else
                vals.push_back(isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
        }

        size_t res_size = num_rows * num_rows;
        res_sizeof = sizeof(double) * res_size;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        int num_cores = max(1, (int)sysconf(_SC_NPROCESSORS_ONLN));
        int num_processes = min(num_rows / 2, num_cores);
        double num_row4process = num_rows / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_rows * num_rows / 2 - num_rows, 1);

        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                int srow[2] = { (int)(iprocess * num_row4process / 2.), (int)(num_rows - (iprocess + 1) * num_row4process / 2.) };
                int erow[2] = { (int)((iprocess + 1) * num_row4process / 2.), (int)(num_rows - iprocess * num_row4process / 2.) };
                size_t itr_idx = 0;

                for (int ipart = 0; ipart < 2; ipart++) {
                    for (int irow1 = 0; irow1 < num_rows; ++irow1) {
                        int start = max(srow[ipart], irow1 + 1);
                        if (start >= erow[ipart])
                            break;

                        size_t idx = start * num_rows + irow1;

                        for (int irow2 = start; irow2 < erow[ipart]; ++irow2) {
                            double dist = 0;
                            double dev;
                            int count = 0;
                            size_t idx1 = irow1;
                            size_t idx2 = irow2;

                            for (int icol = 0; icol < num_cols; ++icol) {
                                if (!isnan(vals[idx1]) && !isnan(vals[idx2])) {
                                    dev = vals[idx1] - vals[idx2];
                                    dist += dev * dev;
                                    ++count;
                                }
                                idx1 += num_rows;
                                idx2 += num_rows;
                            }

                            if (count) {
                                if (count == num_cols)
                                    res[idx] = sqrt(dist);
                                else
                                    res[idx] = sqrt(dist * num_cols / (double)count);
                            } else
                                res[idx] = NA_REAL;

                            idx += num_rows;
                            ++itr_idx;
                            TGStat::itr_idx(itr_idx);
                        }
                    }
                }
                exit(0);
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        if (tidy) {
            enum { ROW1, ROW2, DIST, NUM_COLS };
            const char *COL_NAMES[NUM_COLS] = { "row1", "row2", "dist" };

            rprotect(answer = allocVector(VECSXP, NUM_COLS));

            size_t answer_size = 0;

            for (int irow1 = 0; irow1 < num_rows; ++irow1) {
                size_t idx2 = (irow1 + 1) * num_rows + irow1;
                for (int irow2 = irow1 + 1; irow2 < num_rows; ++irow2) {
                    if (res[idx2] <= threshold)
                        ++answer_size;
                    idx2 += num_rows;
                }
            }

            SEXP rcol1, rcol2, rdist, rrownames, rcolnames;

            SET_VECTOR_ELT(answer, ROW1, (rcol1 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, ROW2, (rcol2 = allocVector(INTSXP, answer_size)));
            SET_VECTOR_ELT(answer, DIST, (rdist = allocVector(REALSXP, answer_size)));

            if (_rrownames != R_NilValue) {
                setAttrib(rcol1, R_LevelsSymbol, _rrownames);
                setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
                setAttrib(rcol2, R_LevelsSymbol, _rrownames);
                setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
            }

            setAttrib(answer, R_NamesSymbol, (rcolnames = allocVector(STRSXP, NUM_COLS)));
            setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
            setAttrib(answer, R_RowNamesSymbol, (rrownames = allocVector(INTSXP, answer_size)));

            for (int i = 0; i < NUM_COLS; i++)
                SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

            size_t idx1 = 0;
            for (int irow1 = 0; irow1 < num_rows; ++irow1) {
                size_t idx2 = (irow1 + 1) * num_rows + irow1;
                for (int irow2 = irow1 + 1; irow2 < num_rows; ++irow2) {
                    if (res[idx2] <= threshold) {
                        INTEGER(rcol1)[idx1] = irow1 + 1;
                        INTEGER(rcol2)[idx1] = irow2 + 1;
                        REAL(rdist)[idx1] = res[idx2];
                        INTEGER(rrownames)[idx1] = idx1 + 1;
                        ++idx1;
                    }
                    idx2 += num_rows;
                }
            }
        } else {
            rprotect(answer = allocVector(REALSXP, (size_t)num_rows * (num_rows - 1) / 2));

            size_t idx1 = 0;
            double *d = REAL(answer);

            for (int irow1 = 0; irow1 < num_rows; ++irow1) {
                size_t idx2 = (irow1 + 1) * num_rows + irow1;
                for (int irow2 = irow1 + 1; irow2 < num_rows; ++irow2) {
                    d[idx1++] = res[idx2];
                    idx2 += num_rows;
                }
            }

            SEXP names = getAttrib(_attrs, R_NamesSymbol);
            for (int i = 0; i < LENGTH(_attrs); i++)
                setAttrib(answer, install(translateChar(STRING_ELT(names, i))), VECTOR_ELT(_attrs, i));

        }

    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (double *)MAP_FAILED;
        }
		rerror("%s", e.msg());
	}

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (double *)MAP_FAILED;
    }
	rreturn(answer);
}

}

