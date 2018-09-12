#include <iterator>
#include <limits>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>

#include <R.h>
#include <Rinternals.h>

#include "ProgressReporter.h"
#include "tgstat.h"

extern "C" {

SEXP tgs_matrix_tapply(SEXP _x, SEXP _index, SEXP _fn, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    size_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        SEXP _rdims = R_NilValue;
        SEXP _xclass = getAttrib(_x, R_ClassSymbol);
        SEXP _xx = getAttrib(_x, install("x"));   // non zero values of sparse matrix
        SEXP _xi = getAttrib(_x, install("i"));   // row number within sparse matrix
        SEXP _xp = getAttrib(_x, install("p"));   // index offset of the column within sparse matrix

        if ((!isReal(_x) && !isInteger(_x) && !isObject(_x)) ||
            ((isReal(_x) || isInteger(_x)) && (xlength(_x) < 1 || !isInteger(_rdims = getAttrib(_x, R_DimSymbol)) || xlength(_rdims) != 2)) ||
            (isObject(_x) && (!isString(_xclass) || xlength(_xclass) < 1 || strcmp(CHAR(STRING_ELT(_xclass, 0)), "dgCMatrix") ||
                              !isInteger(_rdims = getAttrib(_x, install("Dim"))) || xlength(_rdims) != 2 ||
                              (!isInteger(_xx) && !isReal(_xx)) || xlength(_xx) < 1 ||
                              !isInteger(_xi) || xlength(_xi) != xlength(_xx) ||
                              !isInteger(_xp) || xlength(_xp) != INTEGER(_rdims)[1] + 1)))
            verror("\"x\" argument must be a matrix of numeric values");

        size_t num_rows = INTEGER(_rdims)[0];   // do not use nrows(): it doesn't work for sparse (i.e. dgCMatrix) matrix
        size_t num_cols = INTEGER(_rdims)[1];

        if (!isFactor(_index) || xlength(_index) < 1)
            verror("\"index\" argument must be a factor");

        if (xlength(_index) != num_cols)
            verror("Arguments \"x\" and \"index\" must have same length");

        if (!isFunction(_fn) || xlength(_fn) != 1)
            verror("\"fn\" argument must be a function");

        SEXP rcall;
        SEXP rindex_levels = getAttrib(_index, R_LevelsSymbol);
        size_t num_groups = xlength(rindex_levels);

        res_sizeof = sizeof(double) * num_rows * num_groups;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        fill_n(res, num_rows * num_groups, NA_REAL);

        vector<vector<int>> group2cols(num_groups);

        for (int i = 0; i < num_cols; ++i) {
            int group = INTEGER(_index)[i] - 1;
            if (group < 0 || group >= num_groups)
                verror("Invalid group index %d, must be in [0, %d] range", group, (int)num_groups);
            group2cols[group].push_back(i);
        }

        vdebug("Preparing for multitasking...\n");
        int num_processes = (int)min(num_rows / 10, (size_t)g_tgstat->num_processes());

        if (num_rows)
            num_processes = max(num_processes, 1);

        ProgressReporter progress;
        progress.init(num_rows * num_groups, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                size_t srow = num_rows * (iprocess / (double)num_processes);
                size_t erow = num_rows * ((iprocess + 1) / (double)num_processes);
                int *int_vals = _xclass == R_NilValue ? (isInteger(_x) ? INTEGER(_x) : NULL) : (isInteger(_xx) ? INTEGER(_xx) : NULL);
                double *dbl_vals = _xclass == R_NilValue ? (isReal(_x) ? REAL(_x) : NULL) : (isReal(_xx) ? REAL(_xx) : NULL);
                int *rows_sparse = _xclass == R_NilValue ? NULL : INTEGER(_xi);
                int *col_offsets_sparse = _xclass == R_NilValue ? NULL : INTEGER(_xp);
                vector<int> gcol2offset_sparse;   // translates between idx within the group and index within "x"/"i" attributes of dgCMatrix

                for (size_t igroup = 0; igroup < num_groups; ++igroup) {
                    auto &group_cols = group2cols[igroup];

                    if (group_cols.empty())
                        continue;

                    SEXP rgroup;
                    int *int_group = NULL;
                    double *dbl_group = NULL;

                    if ((_xclass == R_NilValue && isReal(_x)) || (_xclass != R_NilValue && isReal(_xx))) {
                        rprotect(rgroup = RSaneAllocVector(REALSXP, group_cols.size()));
                        dbl_group = REAL(rgroup);
                    } else {
                        rprotect(rgroup = RSaneAllocVector(INTSXP, group_cols.size()));
                        int_group = INTEGER(rgroup);
                    }

                    if (_xclass != R_NilValue) {
                        gcol2offset_sparse.resize(group_cols.size());
                        for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol)
                            gcol2offset_sparse[icol - group_cols.begin()] =
                                lower_bound(rows_sparse + col_offsets_sparse[*icol], rows_sparse + col_offsets_sparse[*icol + 1], srow) - rows_sparse;
                    }

                    rprotect(rcall = lang2(_fn, R_NilValue));
                    SETCADR(rcall, rgroup);

                    for (int irow = srow; irow < erow; ++irow) {
                        if (_xclass == R_NilValue) {     // _x is regular matrix
                            if (isReal(_x)) {
                                for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol)
                                    dbl_group[icol - group_cols.begin()] = dbl_vals[num_rows * *icol + irow];
                            } else {
                                for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol)
                                    int_group[icol - group_cols.begin()] = int_vals[num_rows * *icol + irow];
                            }
                        } else {    // _x is sparse matrix of dgCMatrix type
                            if (isReal(_xx)) {
                                for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                    size_t idx = icol - group_cols.begin();
                                    if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == irow) {
                                        dbl_group[idx] = dbl_vals[gcol2offset_sparse[idx]];
                                        ++gcol2offset_sparse[idx];
                                    } else
                                        dbl_group[idx] = 0;
                                }
                            } else {
                                for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                    size_t idx = icol - group_cols.begin();
                                    if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == irow) {
                                        int_group[idx] = int_vals[gcol2offset_sparse[idx]];
                                        ++gcol2offset_sparse[idx];
                                    } else
                                        int_group[idx] = 0;
                                }
                            }
                        }

                        SEXP retv = eval_in_R(rcall, _envir);

                        if (xlength(retv) == 1) {
                            if (isReal(retv))
                                res[num_groups * irow + igroup] = REAL(retv)[0];
                            else if (isInteger(retv))
                                res[num_groups * irow + igroup] = INTEGER(retv)[0];
                            else
                                verror("Evaluation returned neither numeric nor integer");
                        } else
                            verror("Evaluation returned not a scalar");

                        runprotect(1);
                    }
                    TGStat::itr_idx((erow - srow) * (igroup + 1));

                    runprotect(2);
                }

                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        SEXP dim, dimnames;

        rprotect(answer = RSaneAllocVector(REALSXP, num_rows * num_groups));
        memcpy(REAL(answer), res, res_sizeof);

        rprotect(dim = RSaneAllocVector(INTSXP, 2));
        INTEGER(dim)[0] = num_groups;
        INTEGER(dim)[1] = num_rows;
        setAttrib(answer, R_DimSymbol, dim);

        rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, rindex_levels);
        SET_VECTOR_ELT(dimnames, 1, R_NilValue);
        setAttrib(answer, R_DimNamesSymbol, dimnames);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap(res, res_sizeof);
            res = (double *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap(res, res_sizeof);
        res = (double *)MAP_FAILED;
    }

    rreturn(answer);
}

}
