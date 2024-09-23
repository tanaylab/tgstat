#include <iterator>
#include <limits>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>
#include <stdint.h>
#include <string>
#include <stdlib.h>
#include <set>
#include <pthread.h>
#include <semaphore.h>
#include <signal.h>
#include <time.h>
#include <sys/types.h>
#include <cstdint>

#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>

#include "ProgressReporter.h"
#include "tgstat.h"

struct SumData {
    double pre_eval_sum{0};
    bool na_rm{false};
};

struct MeanData {
    double trim{0};
    bool na_rm{false};
};

void init_sum_data(SEXP rargs, SEXP rarg_names, SEXP renvir, SumData *sum_data)
{
    vdebug("overriding R's \"sum\" function");
    int num_unnamed_args = 0;

    if (Rf_isNull(rarg_names))
        num_unnamed_args = Rf_length(rargs);
    else
    {
        for (int i = 0; i < Rf_length(rargs); ++i)
        {
            const char *arg_name = CHAR(STRING_ELT(rarg_names, i));
            if (!*arg_name)
                ++num_unnamed_args;
            else if (!strcmp(arg_name, "na.rm"))
            {
                SEXP retv = eval_in_R(VECTOR_ELT(rargs, i), renvir);
                sum_data->na_rm = Rf_asLogical(retv);
                runprotect(1);
            }
        }
    }

    // Sum may accept more than one unnamed argument (vectors).
    // Calculate the sum now and add the value later to our matrix sum.
    if (num_unnamed_args)
    {
        SEXP parsed_expr;
        SEXP cmd;
        ParseStatus status;

        rprotect(cmd = Rf_ScalarString(Rf_mkChar("sum")));
        rprotect(parsed_expr = R_ParseVector(cmd, -1, &status, R_NilValue));

        if (status != PARSE_OK)
            verror("R parsing of expression \"sum\" failed");

        SEXP eval_expr = VECTOR_ELT(parsed_expr, 0);
        SEXP rcall = R_NilValue;

        // Build the call using Rf_lcons
        rprotect(rcall = Rf_lcons(eval_expr, R_NilValue));
        SEXP tail = rcall;

        for (int i = 0; i < Rf_length(rargs); ++i)
        {
            SEXP next = Rf_lcons(VECTOR_ELT(rargs, i), R_NilValue);
            SETCDR(tail, next);
            tail = next;
            if (!Rf_isNull(rarg_names))
            {
                const char *arg_name = CHAR(STRING_ELT(rarg_names, i));
                if (*arg_name)
                    SET_TAG(next, Rf_install(arg_name));
            }
        }

        SEXP retv = eval_in_R(rcall, renvir);
        if (Rf_xlength(retv) != 1)
            verror("Evaluation of \"sum\" did not return a scalar");

        sum_data->pre_eval_sum = Rf_asReal(retv);
        runprotect(3); // Adjust protect count if necessary
    }
}

void init_mean_data(SEXP rargs, SEXP rarg_names, SEXP renvir, MeanData *mean_data)
{
    vdebug("overriding R's \"mean\" function");
    bool trim_set = false;
    bool na_rm_set = false;

    // first set the named arguments
    if (!Rf_isNull(rarg_names)) {
        for (int i = 0; i < Rf_length(rargs); ++i) {
            const char *arg_name = CHAR(STRING_ELT(rarg_names, i));
            if (*arg_name) {
                if (!strcmp(arg_name, "trim")) {
                    SEXP retv = eval_in_R(VECTOR_ELT(rargs, i), renvir);
                    mean_data->trim = Rf_asReal(retv);
                    trim_set = true;
                    runprotect(1);
                } else if (!strcmp(arg_name, "na.rm")) {
                    SEXP retv = eval_in_R(VECTOR_ELT(rargs, i), renvir);
                    mean_data->na_rm = Rf_asLogical(retv);
                    na_rm_set = true;
                    runprotect(1);
                }
            }
        }
    }

    // now set the unnamed if the named had not been set yet
    for (int i = 0; i < Rf_length(rargs); ++i) {
        if (Rf_isNull(rarg_names) || !*CHAR(STRING_ELT(rarg_names, i))) {
            if (!trim_set) {
                SEXP retv = eval_in_R(VECTOR_ELT(rargs, i), renvir);
                mean_data->trim = Rf_asReal(retv);
                trim_set = true;
                runprotect(1);
            } else if (!na_rm_set) {
                SEXP retv = eval_in_R(VECTOR_ELT(rargs, i), renvir);
                mean_data->na_rm = Rf_asLogical(retv);
                na_rm_set = true;
                runprotect(1);
            }
        }
    }
}

extern "C" {

SEXP tgs_matrix_tapply(SEXP _x, SEXP _index, SEXP _fn, SEXP _fn_name, SEXP _args, SEXP _envir)
{
    SEXP answer = R_NilValue;
    double *res = (double *)MAP_FAILED;
    uint64_t res_sizeof = 0;

    try {
        TGStat tgstat(_envir);

        SEXP _rdims = R_NilValue;
        SEXP _xclass = Rf_getAttrib(_x, R_ClassSymbol);
        SEXP _xx = Rf_getAttrib(_x, Rf_install("x"));   // non zero values of sparse matrix
        SEXP _xi = Rf_getAttrib(_x, Rf_install("i"));   // row number within sparse matrix
        SEXP _xp = Rf_getAttrib(_x, Rf_install("p"));   // index offset of the column within sparse matrix

        if ((!Rf_isReal(_x) && !Rf_isInteger(_x) && !Rf_isObject(_x)) ||
            ((Rf_isReal(_x) || Rf_isInteger(_x)) && (Rf_xlength(_x) < 1 || !Rf_isInteger(_rdims = Rf_getAttrib(_x, R_DimSymbol)) || Rf_xlength(_rdims) != 2)) ||
            (Rf_isObject(_x) && (!Rf_isString(_xclass) || Rf_xlength(_xclass) < 1 || strcmp(CHAR(STRING_ELT(_xclass, 0)), "dgCMatrix") ||
                              !Rf_isInteger(_rdims = Rf_getAttrib(_x, Rf_install("Dim"))) || Rf_xlength(_rdims) != 2 ||
                              (!Rf_isInteger(_xx) && !Rf_isReal(_xx)) || Rf_xlength(_xx) < 1 ||
                              !Rf_isInteger(_xi) || Rf_xlength(_xi) != Rf_xlength(_xx) ||
                              !Rf_isInteger(_xp) || Rf_xlength(_xp) != INTEGER(_rdims)[1] + 1)))
            verror("\"x\" argument must be a matrix of numeric values");

        SEXP _xdimnames = _xclass == R_NilValue ? Rf_getAttrib(_x, R_DimNamesSymbol) : Rf_getAttrib(_x, Rf_install("Dimnames"));
        uint64_t num_rows = INTEGER(_rdims)[0];   // do not use Rf_nrows(): it doesn't work for sparse (i.e. dgCMatrix) matrix
        uint64_t num_cols = INTEGER(_rdims)[1];

        if (!Rf_isFactor(_index) || Rf_xlength(_index) < 1)
            verror("\"index\" argument must be a factor");

        if ((uint64_t)Rf_xlength(_index) != num_cols)
            verror("Arguments \"x\" and \"index\" must have same length");

        if (!Rf_isFunction(_fn) || Rf_xlength(_fn) != 1)
            verror("\"fn\" argument must be a function");

        if (!Rf_isString(_fn_name))
            verror("\"fn_name\" argument must be a string");

        // in case "fun=function(x) mean(x)", _fn_name will be a vector of strings
        string fn_name(Rf_xlength(_fn_name) == 1 ? CHAR(Rf_asChar(_fn_name)) : "");
        SEXP rarg_names = Rf_getAttrib(_args, R_NamesSymbol);
        SEXP rindex_levels = Rf_getAttrib(_index, R_LevelsSymbol);
        uint64_t num_groups = (uint64_t)Rf_xlength(rindex_levels);
        bool is_sum = fn_name == "sum";
        SumData sum_data;
        bool is_mean = fn_name == "mean";
        MeanData mean_data;

        if (is_sum)
            init_sum_data(_args, rarg_names, _envir, &sum_data);
        else if (is_mean) {
            init_mean_data(_args, rarg_names, _envir, &mean_data);
            if (mean_data.trim) {
                vdebug("mean's \"trim\" argument prevents override of \"mean\"\n");
                is_mean = false;
            }
        }

        res_sizeof = sizeof(double) * num_rows * num_groups;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        fill_n(res, num_rows * num_groups, NA_REAL);

        vector<vector<int>> group2cols(num_groups);

        for (uint64_t i = 0; i < num_cols; ++i) {
            int group = INTEGER(_index)[i] - 1;
            if (group < 0 || (unsigned)group >= num_groups)
                verror("Invalid group index %d, must be in [0, %d] range", group, (int)num_groups);
            group2cols[group].push_back(i);
        }

        vdebug("Preparing for multitasking...\n");
        int num_processes = (int)min(num_rows / 10, (uint64_t)(g_tgstat->num_processes() / 2));

        // sum and mean are simple functions with complexity O(N). The overhead of extreme parallelization might obliterate the gain.
        if (is_sum || is_mean)
            num_processes = min(5, num_processes);

        if (num_rows)
            num_processes = max(num_processes, 1);

        ProgressReporter progress;
        progress.init(num_rows * num_groups, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                uint64_t srow = num_rows * (iprocess / (double)num_processes);
                uint64_t erow = num_rows * ((iprocess + 1) / (double)num_processes);
                int *int_vals = _xclass == R_NilValue ? (Rf_isInteger(_x) ? INTEGER(_x) : NULL) : (Rf_isInteger(_xx) ? INTEGER(_xx) : NULL);
                double *dbl_vals = _xclass == R_NilValue ? (Rf_isReal(_x) ? REAL(_x) : NULL) : (Rf_isReal(_xx) ? REAL(_xx) : NULL);
                int *rows_sparse = _xclass == R_NilValue ? NULL : INTEGER(_xi);
                int *col_offsets_sparse = _xclass == R_NilValue ? NULL : INTEGER(_xp);
                vector<int> gcol2offset_sparse;   // translates between idx within the group and index within "x"/"i" attributes of dgCMatrix

                for (uint64_t igroup = 0; igroup < num_groups; ++igroup) {
                    auto &group_cols = group2cols[igroup];

                    if (group_cols.empty())
                        continue;

                    if (_xclass != R_NilValue) {
                        gcol2offset_sparse.resize(group_cols.size());
                        for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol)
                            gcol2offset_sparse[icol - group_cols.begin()] =
                                lower_bound(rows_sparse + col_offsets_sparse[*icol], rows_sparse + col_offsets_sparse[*icol + 1], srow) - rows_sparse;
                    }

                    if (is_sum) {
                        for (auto irow = srow; irow < erow; ++irow) {
                            double sum = sum_data.pre_eval_sum;

                            if (_xclass == R_NilValue) {     // _x is regular matrix
                                if (Rf_isReal(_x)) {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        double val = dbl_vals[num_rows * *icol + irow];
                                        if (!sum_data.na_rm || !std::isnan(val))
                                            sum += val;
                                    }
                                } else {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol)
                                        sum += int_vals[num_rows * *icol + irow];
                                }
                            } else {    // _x is sparse matrix of dgCMatrix type
                                if (Rf_isReal(_xx)) {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        uint64_t idx = icol - group_cols.begin();
                                        if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == (int)irow) {
                                            double val = dbl_vals[gcol2offset_sparse[idx]];
                                            if (!std::isnan(val) || !sum_data.na_rm)
                                                sum += val;
                                            ++gcol2offset_sparse[idx];
                                        }
                                    }
                                } else {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        uint64_t idx = icol - group_cols.begin();
                                        if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == (int)irow) {
                                            sum += int_vals[gcol2offset_sparse[idx]];
                                            ++gcol2offset_sparse[idx];
                                        }
                                    }
                                }
                            }

                            res[num_groups * irow + igroup] = std::isnan(sum) ? NA_REAL : sum;
                        }
                    } else if (is_mean) {
                        for (auto irow = srow; irow < erow; ++irow) {
                            double sum = 0;
                            uint64_t num_vals = 0;

                            if (_xclass == R_NilValue) {     // _x is regular matrix
                                if (Rf_isReal(_x)) {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        double val = dbl_vals[num_rows * *icol + irow];
                                        if (!mean_data.na_rm || !std::isnan(val)) {
                                            sum += val;
                                            ++num_vals;
                                        }
                                    }
                                } else {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        sum += int_vals[num_rows * *icol + irow];
                                        ++num_vals;
                                    }
                                }
                            } else {    // _x is sparse matrix of dgCMatrix type
                                if (Rf_isReal(_xx)) {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        uint64_t idx = icol - group_cols.begin();
                                        if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == (int)irow) {
                                            double val = dbl_vals[gcol2offset_sparse[idx]];
                                            if (!std::isnan(val) || !mean_data.na_rm) {
                                                sum += val;
                                                ++num_vals;
                                            }
                                            ++gcol2offset_sparse[idx];
                                        } else
                                            ++num_vals;
                                    }
                                } else {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        uint64_t idx = icol - group_cols.begin();
                                        if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == (int)irow) {
                                            sum += int_vals[gcol2offset_sparse[idx]];
                                            ++gcol2offset_sparse[idx];
                                            ++num_vals;
                                        }
                                    }
                                }
                            }

                            res[num_groups * irow + igroup] = std::isnan(sum) || !num_vals ? NA_REAL : sum / num_vals;
                        }
                    } else {
                        SEXP rgroup;
                        int *int_group = NULL;
                        double *dbl_group = NULL;

                        if ((_xclass == R_NilValue && Rf_isReal(_x)) || (_xclass != R_NilValue && Rf_isReal(_xx))) {
                            rprotect(rgroup = RSaneAllocVector(REALSXP, group_cols.size()));
                            dbl_group = REAL(rgroup);
                        } else {
                            rprotect(rgroup = RSaneAllocVector(INTSXP, group_cols.size()));
                            int_group = INTEGER(rgroup);
                        }

                        // construct language expression for eval:
                        // function, vector of values, additional arguments in _args
                        SEXP rcall = R_NilValue;
                        rprotect(rcall = Rf_lcons(_fn, Rf_lcons(rgroup, R_NilValue)));

                        SEXP tail = CDR(CDR(rcall));
                        for (int i = 0; i < Rf_length(_args); ++i)
                        {
                            SEXP next = Rf_lcons(VECTOR_ELT(_args, i), R_NilValue);
                            SETCDR(tail, next);
                            tail = next;
                            if (!Rf_isNull(rarg_names))
                            {
                                const char *arg_name = CHAR(STRING_ELT(rarg_names, i));
                                if (*arg_name)
                                    SET_TAG(next, Rf_install(arg_name));
                            }
                        }

                        for (auto irow = srow; irow < erow; ++irow) {
                            if (_xclass == R_NilValue) {     // _x is regular matrix
                                if (Rf_isReal(_x)) {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol)
                                        dbl_group[icol - group_cols.begin()] = dbl_vals[num_rows * *icol + irow];
                                } else {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol)
                                        int_group[icol - group_cols.begin()] = int_vals[num_rows * *icol + irow];
                                }
                            } else {    // _x is sparse matrix of dgCMatrix type
                                if (Rf_isReal(_xx)) {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        uint64_t idx = icol - group_cols.begin();
                                        if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == (int)irow) {
                                            dbl_group[idx] = dbl_vals[gcol2offset_sparse[idx]];
                                            ++gcol2offset_sparse[idx];
                                        } else
                                            dbl_group[idx] = 0;
                                    }
                                } else {
                                    for (auto icol = group_cols.begin(); icol != group_cols.end(); ++icol) {
                                        uint64_t idx = icol - group_cols.begin();
                                        if (gcol2offset_sparse[idx] < col_offsets_sparse[*icol + 1] && rows_sparse[gcol2offset_sparse[idx]] == (int)irow) {
                                            int_group[idx] = int_vals[gcol2offset_sparse[idx]];
                                            ++gcol2offset_sparse[idx];
                                        } else
                                            int_group[idx] = 0;
                                    }
                                }
                            }

                            SEXP retv = eval_in_R(rcall, _envir);

                            if (Rf_xlength(retv) == 1) {
                                if (Rf_isReal(retv))
                                    res[num_groups * irow + igroup] = REAL(retv)[0];
                                else if (Rf_isInteger(retv))
                                    res[num_groups * irow + igroup] = INTEGER(retv)[0];
                                else
                                    verror("Evaluation returned neither numeric nor integer");
                            } else
                                verror("Evaluation did not return a scalar");

                            runprotect(1);
                        }
                        runprotect(2);
                    }
                    TGStat::itr_idx((erow - srow) * (igroup + 1));
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
        Rf_setAttrib(answer, R_DimSymbol, dim);

        rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, rindex_levels);
        SET_VECTOR_ELT(dimnames, 1, !Rf_isNull(_xdimnames) && Rf_xlength(_xdimnames) > 0 ? VECTOR_ELT(_xdimnames, 0) : R_NilValue);
        Rf_setAttrib(answer, R_DimNamesSymbol, dimnames);
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
            munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
            res = (double *)MAP_FAILED;
        }
        rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && res != (double *)MAP_FAILED) {
        munmap((char *)res, res_sizeof);  // needs to be char * for some versions of Solaris
        res = (double *)MAP_FAILED;
    }

    rreturn(answer);
}

}
