#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <cstring>
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

// Sets na_rm from a named "na.rm" argument.
// Returns false if there are other arguments or if "na.rm" is neither TRUE nor FALSE.
bool init_na_rm(SEXP rargs, SEXP rarg_names, SEXP renvir, bool *na_rm)
{
    for (int i = 0; i < Rf_length(rargs); ++i) {
        if (Rf_isNull(rarg_names) || strcmp(CHAR(STRING_ELT(rarg_names, i)), "na.rm"))
            return false;

        SEXP retv = eval_in_R(VECTOR_ELT(rargs, i), renvir);
        int val = Rf_asLogical(retv);
        runprotect(1);
        if (val == NA_LOGICAL)
            return false;
        *na_rm = val;
    }
    return true;
}

// Returns the k-th smallest (0-based) of the n values of v together with num_zeros zeros. v must not contain NaN and is reordered.
static double kth_smallest(double *v, uint64_t n, uint64_t num_zeros, uint64_t k)
{
    if (num_zeros) {
        uint64_t num_neg = partition(v, v + n, [](double a) { return a < 0; }) - v;
        if (k < num_neg)
            n = num_neg;
        else {
            uint64_t num_nonpos = partition(v + num_neg, v + n, [](double a) { return a == 0; }) - v;
            if (k < num_nonpos + num_zeros)
                return 0;
            v += num_nonpos;
            n -= num_nonpos;
            k -= num_nonpos + num_zeros;
        }
    }
    nth_element(v, v + k, v + n);
    return v[k];
}

// Returns the median of the n values of v together with num_zeros zeros, as R's median() computes it: NA if a value is
// missing (unless na_rm) or if there are no values, the mean of the two middle values for an even count. v is reordered.
static double median(double *v, uint64_t n, uint64_t num_zeros, bool na_rm)
{
    double *end = v + n;

    if (na_rm)
        end = remove_if(v, end, [](double a) { return std::isnan(a); });
    else if (any_of(v, end, [](double a) { return std::isnan(a); }))
        return NA_REAL;

    n = end - v;
    uint64_t num_vals = n + num_zeros;

    if (!num_vals)
        return NA_REAL;

    if (num_vals % 2)
        return kth_smallest(v, n, num_zeros, num_vals / 2);

    double a = kth_smallest(v, n, num_zeros, num_vals / 2 - 1);
    double b = kth_smallest(v, n, num_zeros, num_vals / 2);

    // R's mean() of a and b (real_mean in R's summary.c)
    long double s = (long double)a + b;
    bool finite_s = R_FINITE((double)s);

    s = finite_s ? s / 2 : (long double)(a / 2) + b / 2;
    if (finite_s && R_FINITE((double)s))
        s += ((a - s) + (b - s)) / 2;
    else if (R_FINITE((double)s))
        s += (a - s) / 2 + (b - s) / 2;
    return (double)s;
}

enum Fn { GENERIC, SUM, MEAN, MAX, MIN, MEDIAN };

extern "C" {

SEXP tgs_matrix_tapply(SEXP _x, SEXP _index, SEXP _fn, SEXP _fn_name, SEXP _args, SEXP _margin, SEXP _envir)
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

        bool is_sparse = _xclass != R_NilValue;
        SEXP _xdimnames = is_sparse ? Rf_getAttrib(_x, Rf_install("Dimnames")) : Rf_getAttrib(_x, R_DimNamesSymbol);
        uint64_t num_rows = INTEGER(_rdims)[0];   // do not use Rf_nrows(): it doesn't work for sparse (i.e. dgCMatrix) matrix
        uint64_t num_cols = INTEGER(_rdims)[1];

        if (!Rf_isInteger(_margin) || Rf_xlength(_margin) != 1 || (INTEGER(_margin)[0] != 1 && INTEGER(_margin)[0] != 2))
            verror("\"margin\" argument must be 1 or 2");

        // margin 1: "index" groups the rows of x and the result has a column for each column of x;
        // margin 2: "index" groups the columns of x and the result has a column for each row of x
        bool by_rows = INTEGER(_margin)[0] == 1;
        uint64_t num_items = by_rows ? num_cols : num_rows;

        if (!Rf_isFactor(_index) || Rf_xlength(_index) < 1)
            verror("\"index\" argument must be a factor");

        if ((uint64_t)Rf_xlength(_index) != (by_rows ? num_rows : num_cols))
            verror("The length of \"index\" must be equal to the number of %s of \"x\"", by_rows ? "rows" : "columns");

        if (!Rf_isFunction(_fn) || Rf_xlength(_fn) != 1)
            verror("\"fn\" argument must be a function");

        if (!Rf_isString(_fn_name))
            verror("\"fn_name\" argument must be a string");

        // in case "fun=function(x) mean(x)", _fn_name will be a vector of strings
        string fn_name(Rf_xlength(_fn_name) == 1 ? CHAR(Rf_asChar(_fn_name)) : "");
        SEXP rarg_names = Rf_getAttrib(_args, R_NamesSymbol);
        SEXP rindex_levels = Rf_getAttrib(_index, R_LevelsSymbol);
        uint64_t num_groups = (uint64_t)Rf_xlength(rindex_levels);
        Fn fn = GENERIC;
        bool na_rm = false;
        double init = 0;   // initial value of sum, mean, max and min

        if (fn_name == "sum") {
            SumData sum_data;
            init_sum_data(_args, rarg_names, _envir, &sum_data);
            fn = SUM;
            na_rm = sum_data.na_rm;
            init = sum_data.pre_eval_sum;
        } else if (fn_name == "mean") {
            MeanData mean_data;
            init_mean_data(_args, rarg_names, _envir, &mean_data);
            if (mean_data.trim)
                vdebug("mean's \"trim\" argument prevents override of \"mean\"\n");
            else {
                fn = MEAN;
                na_rm = mean_data.na_rm;
            }
        } else if (fn_name == "max" || fn_name == "min" || fn_name == "median") {
            vdebug("overriding R's \"%s\" function", fn_name.c_str());
            if (init_na_rm(_args, rarg_names, _envir, &na_rm)) {
                fn = fn_name == "max" ? MAX : fn_name == "min" ? MIN : MEDIAN;
                init = fn == MAX ? -numeric_limits<double>::infinity() : numeric_limits<double>::infinity();
            } else
                vdebug("arguments other than \"na.rm\" prevent override of \"%s\"\n", fn_name.c_str());
        }

        res_sizeof = sizeof(double) * num_items * num_groups;
        res = (double *)mmap(NULL, res_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (res == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        fill_n(res, num_items * num_groups, NA_REAL);

        const int *codes = INTEGER(_index);
        vector<uint64_t> group_sizes(num_groups, 0);
        // the columns (margin 2) or rows (margin 1) of each group, in increasing order; not needed by the native functions with margin 1
        vector<vector<int>> group2idx(by_rows && fn != GENERIC ? 0 : num_groups);

        for (uint64_t i = 0; i < (uint64_t)Rf_xlength(_index); ++i) {
            int group = codes[i] - 1;
            if (group < 0 || (unsigned)group >= num_groups)
                verror("Invalid group index %d, must be in [0, %d] range", group, (int)num_groups);
            ++group_sizes[group];
            if (!group2idx.empty())
                group2idx[group].push_back(i);
        }

        vdebug("Preparing for multitasking...\n");
        int num_processes = (int)min(num_items / 10, (uint64_t)(g_tgstat->num_processes() / 2));

        // sum, mean, max and min are simple functions with complexity O(N). The overhead of extreme parallelization might obliterate the gain.
        if (fn == SUM || fn == MEAN || fn == MAX || fn == MIN)
            num_processes = min(5, num_processes);

        if (num_items)
            num_processes = max(num_processes, 1);

        ProgressReporter progress;
        progress.init(num_items * num_groups, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                // the rows (margin 2) or columns (margin 1) of x handled by this process
                uint64_t sitem = iprocess * num_items / num_processes;
                uint64_t eitem = (iprocess + 1) * num_items / num_processes;
                const int *int_vals = is_sparse ? (Rf_isInteger(_xx) ? INTEGER(_xx) : NULL) : (Rf_isInteger(_x) ? INTEGER(_x) : NULL);
                const double *dbl_vals = is_sparse ? (Rf_isReal(_xx) ? REAL(_xx) : NULL) : (Rf_isReal(_x) ? REAL(_x) : NULL);
                const int *rows_sparse = is_sparse ? INTEGER(_xi) : NULL;
                const int *col_offsets_sparse = is_sparse ? INTEGER(_xp) : NULL;

                if (fn != GENERIC) {
                    // i-th value of x (dense) or of its non zero values (sparse) as double
                    auto val = [&](uint64_t i) { return dbl_vals ? dbl_vals[i] : int_vals[i] == NA_INTEGER ? NA_REAL : (double)int_vals[i]; };

                    // sum, mean, max and min are computed in one pass over the values, in accs and counts;
                    // counts is the number of values skipped due to na.rm for mean and the number of values seen for max and min
                    vector<double> accs;
                    vector<uint64_t> counts;
                    vector<double> vals;        // median: the values of all the slots, slot after slot
                    vector<uint64_t> starts;    // median: the offset of each slot in vals

                    auto update = [&](double &acc, uint64_t &count, double v) {
                        if (fn == SUM) {
                            if (!na_rm || !std::isnan(v))
                                acc += v;
                        } else if (fn == MEAN) {
                            if (!na_rm || !std::isnan(v))
                                acc += v;
                            else
                                ++count;
                        } else {   // max, min: as in R, NA wins over NaN which wins over any number
                            ++count;
                            if (std::isnan(v)) {
                                if (!na_rm && !R_IsNA(acc))
                                    acc = v;
                            } else if (fn == MAX ? v > acc : v < acc)
                                acc = v;
                        }
                    };

                    // size: the number of values of the slot, including the implicit zeros of a sparse matrix
                    auto finish = [&](double acc, uint64_t count, uint64_t size) {
                        if (fn == SUM)
                            return std::isnan(acc) ? NA_REAL : acc;
                        if (fn == MEAN)
                            return std::isnan(acc) || count == size ? NA_REAL : acc / (size - count);
                        if (count < size)  // implicit zeros
                            update(acc, count, 0);
                        return acc;
                    };

                    // Applies the function to num_slots (item, group) slots. each(f) calls f(slot, value) for each value of x in the slots
                    // except the implicit zeros, ordered by column (margin 2) or row (margin 1) within each slot. size(slot) is the
                    // number of values of the slot including the implicit zeros and res_idx(slot) its index in res.
                    auto reduce = [&](uint64_t num_slots, auto each, auto size, auto res_idx) {
                        if (fn == MEDIAN) {
                            starts.assign(num_slots + 1, 0);
                            if (is_sparse)
                                each([&](uint64_t slot, double) { ++starts[slot + 1]; });
                            else {
                                for (uint64_t slot = 0; slot < num_slots; ++slot)
                                    starts[slot + 1] = size(slot);
                            }
                            partial_sum(starts.begin(), starts.end(), starts.begin());
                            vals.resize(starts.back());
                            counts.assign(starts.begin(), starts.end() - 1);
                            each([&](uint64_t slot, double v) { vals[counts[slot]++] = v; });
                            for (uint64_t slot = 0; slot < num_slots; ++slot) {
                                uint64_t n = starts[slot + 1] - starts[slot];
                                if (size(slot))
                                    res[res_idx(slot)] = median(vals.data() + starts[slot], n, size(slot) - n, na_rm);
                            }
                        } else {
                            accs.assign(num_slots, init);
                            counts.assign(num_slots, 0);
                            each([&](uint64_t slot, double v) { update(accs[slot], counts[slot], v); });
                            for (uint64_t slot = 0; slot < num_slots; ++slot) {
                                if (size(slot))
                                    res[res_idx(slot)] = finish(accs[slot], counts[slot], size(slot));
                            }
                        }
                    };

                    if (by_rows) {
                        auto size = [&](uint64_t igroup) { return group_sizes[igroup]; };

                        for (uint64_t icol = sitem; icol < eitem; ++icol) {
                            auto res_idx = [&](uint64_t igroup) { return num_groups * icol + igroup; };

                            if (is_sparse) {
                                reduce(num_groups, [&](auto f) {
                                    for (int i = col_offsets_sparse[icol]; i < col_offsets_sparse[icol + 1]; ++i)
                                        f(codes[rows_sparse[i]] - 1, val(i));
                                }, size, res_idx);
                            } else {
                                reduce(num_groups, [&](auto f) {
                                    for (uint64_t irow = 0; irow < num_rows; ++irow)
                                        f(codes[irow] - 1, val(num_rows * icol + irow));
                                }, size, res_idx);
                            }
                            TGStat::itr_idx((icol - sitem + 1) * num_groups);
                        }
                    } else {
                        for (uint64_t igroup = 0; igroup < num_groups; ++igroup) {
                            const auto &group_cols = group2idx[igroup];
                            auto size = [&](uint64_t) { return group_cols.size(); };

                            if (is_sparse) {
                                // one slot per row of the process
                                reduce(eitem - sitem, [&](auto f) {
                                    for (int icol : group_cols) {
                                        const int *end = rows_sparse + col_offsets_sparse[icol + 1];
                                        for (const int *p = lower_bound(rows_sparse + col_offsets_sparse[icol], end, (int)sitem); p < end && *p < (int)eitem; ++p)
                                            f(*p - sitem, val(p - rows_sparse));
                                    }
                                }, size, [&](uint64_t slot) { return num_groups * (sitem + slot) + igroup; });
                            } else if (!group_cols.empty()) {
                                // one slot per row, in blocks of rows that limit the values held for median to about 1M
                                uint64_t block_size = max((uint64_t)1, ((uint64_t)1 << 20) / group_cols.size());

                                for (uint64_t srow = sitem; srow < eitem; srow += block_size) {
                                    uint64_t erow = min(eitem, srow + block_size);

                                    reduce(erow - srow, [&](auto f) {
                                        for (int icol : group_cols) {
                                            for (uint64_t irow = srow; irow < erow; ++irow)
                                                f(irow - srow, val(num_rows * icol + irow));
                                        }
                                    }, size, [&](uint64_t slot) { return num_groups * (srow + slot) + igroup; });
                                }
                            }
                            TGStat::itr_idx((eitem - sitem) * (igroup + 1));
                        }
                    }
                } else {
                    vector<int> gcol2offset_sparse;   // translates between idx within the group and index within "x"/"i" attributes of dgCMatrix

                    for (uint64_t igroup = 0; igroup < num_groups; ++igroup) {
                        auto &group_idx = group2idx[igroup];

                        if (group_idx.empty()) {
                            TGStat::itr_idx((eitem - sitem) * (igroup + 1));
                            continue;
                        }

                        if (is_sparse && !by_rows) {
                            gcol2offset_sparse.resize(group_idx.size());
                            for (auto icol = group_idx.begin(); icol != group_idx.end(); ++icol)
                                gcol2offset_sparse[icol - group_idx.begin()] =
                                    lower_bound(rows_sparse + col_offsets_sparse[*icol], rows_sparse + col_offsets_sparse[*icol + 1], (int)sitem) - rows_sparse;
                        }

                        SEXP rgroup;
                        int *int_group = NULL;
                        double *dbl_group = NULL;

                        if (dbl_vals) {
                            rprotect(rgroup = RSaneAllocVector(REALSXP, group_idx.size()));
                            dbl_group = REAL(rgroup);
                        } else {
                            rprotect(rgroup = RSaneAllocVector(INTSXP, group_idx.size()));
                            int_group = INTEGER(rgroup);
                        }

                        // construct language expression for eval:
                        // function, vector of values, additional arguments in _args
                        SEXP rcall = R_NilValue;
                        rprotect(rcall = Rf_lcons(_fn, Rf_lcons(rgroup, R_NilValue)));

                        SEXP tail = CDR(rcall);
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

                        for (auto item = sitem; item < eitem; ++item) {
                            const int *col_start = by_rows && is_sparse ? rows_sparse + col_offsets_sparse[item] : NULL;
                            const int *col_end = by_rows && is_sparse ? rows_sparse + col_offsets_sparse[item + 1] : NULL;

                            for (uint64_t k = 0; k < group_idx.size(); ++k) {
                                int64_t i;   // index of the value in x (dense) or in its non zero values (sparse), -1 for an implicit zero

                                if (!is_sparse)
                                    i = by_rows ? num_rows * item + group_idx[k] : num_rows * group_idx[k] + item;
                                else if (!by_rows) {
                                    i = -1;
                                    if (gcol2offset_sparse[k] < col_offsets_sparse[group_idx[k] + 1] && rows_sparse[gcol2offset_sparse[k]] == (int)item)
                                        i = gcol2offset_sparse[k]++;
                                } else {
                                    col_start = lower_bound(col_start, col_end, group_idx[k]);
                                    i = col_start < col_end && *col_start == group_idx[k] ? col_start - rows_sparse : -1;
                                }

                                if (dbl_group)
                                    dbl_group[k] = i < 0 ? 0 : dbl_vals[i];
                                else
                                    int_group[k] = i < 0 ? 0 : int_vals[i];
                            }

                            SEXP retv = eval_in_R(rcall, _envir);

                            if (Rf_xlength(retv) == 1) {
                                if (Rf_isReal(retv))
                                    res[num_groups * item + igroup] = REAL(retv)[0];
                                else if (Rf_isInteger(retv))
                                    res[num_groups * item + igroup] = INTEGER(retv)[0] == NA_INTEGER ? NA_REAL : INTEGER(retv)[0];
                                else
                                    verror("Evaluation returned neither numeric nor integer");
                            } else
                                verror("Evaluation did not return a scalar");

                            runprotect(1);
                        }
                        runprotect(2);
                        TGStat::itr_idx((eitem - sitem) * (igroup + 1));
                    }
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        SEXP dim, dimnames;

        rprotect(answer = RSaneAllocVector(REALSXP, num_items * num_groups));
        memcpy(REAL(answer), res, res_sizeof);

        rprotect(dim = RSaneAllocVector(INTSXP, 2));
        INTEGER(dim)[0] = num_groups;
        INTEGER(dim)[1] = num_items;
        Rf_setAttrib(answer, R_DimSymbol, dim);

        rprotect(dimnames = RSaneAllocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, rindex_levels);
        SET_VECTOR_ELT(dimnames, 1, !Rf_isNull(_xdimnames) && Rf_xlength(_xdimnames) == 2 ? VECTOR_ELT(_xdimnames, by_rows ? 1 : 0) : R_NilValue);
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
