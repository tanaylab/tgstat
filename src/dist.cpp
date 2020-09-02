#include <algorithm>
#include <cmath>
#include <errno.h>
#include <limits>
#include <sys/mman.h>
#include <unistd.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

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
    void *shm = (double *)MAP_FAILED;
    uint64_t shm_sizeof = 0;

	try {
        TGStat tgstat(_envir);

		if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!isLogical(_tidy) || xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!isReal(_threshold) && !isInteger(_threshold)) || xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = getAttrib(_x, R_DimSymbol);

        if (!isInteger(rdim) || xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));

        uint64_t num_points = nrows(_x);
        uint64_t num_dims = ncols(_x);

        if (num_points < 1 || num_dims < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        uint64_t num_vals = num_points * num_dims;
        vector<double> vals;

        vals.reserve(num_vals);

        for (uint64_t i = 0; i < num_vals; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER))
                vals.push_back(numeric_limits<double>::quiet_NaN());
            else
                vals.push_back(isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i]);
        }

        vdebug("Allocating shared memory for results\n");
        uint64_t res_size = num_points * num_points;
        shm_sizeof = sizeof(double) * res_size;
        shm = mmap(NULL, shm_sizeof, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        if (shm == (double *)MAP_FAILED)
            verror("Failed to allocate shared memory: %s", strerror(errno));

        double *res = (double *)shm;

        int num_processes = (int)min(num_points / 2, (uint64_t)g_tgstat->num_processes());
        double num_row4process = num_points / (double)num_processes;

        ProgressReporter progress;
        progress.init(num_points * num_points / 2 - num_points, 1);

        vdebug("num_processes: %d\n", num_processes);
        TGStat::prepare4multitasking();

        for (int iprocess = 0; iprocess < num_processes; ++iprocess) {
            if (!TGStat::launch_process()) {     // child process
                uint64_t srow[2] = { (uint64_t)(iprocess * num_row4process / 2.), (uint64_t)(num_points - (iprocess + 1) * num_row4process / 2.) };
                uint64_t erow[2] = { (uint64_t)((iprocess + 1) * num_row4process / 2.), (uint64_t)(num_points - iprocess * num_row4process / 2.) };
                uint64_t itr_idx = 0;

                for (int ipart = 0; ipart < 2; ipart++) {
                    for (uint64_t irow1 = 0; irow1 < num_points; ++irow1) {
                        uint64_t start = max(srow[ipart], irow1 + 1);
                        if (start >= erow[ipart])
                            break;

                        uint64_t idx = start * num_points + irow1;

                        for (uint64_t irow2 = start; irow2 < erow[ipart]; ++irow2) {
                            double dist = 0;
                            double dev;
                            uint64_t count = 0;
                            uint64_t idx1 = irow1;
                            uint64_t idx2 = irow2;

                            for (uint64_t icol = 0; icol < num_dims; ++icol) {
                                if (!std::isnan(vals[idx1]) && !std::isnan(vals[idx2])) {
                                    dev = vals[idx1] - vals[idx2];
                                    dist += dev * dev;
                                    ++count;
                                }
                                idx1 += num_points;
                                idx2 += num_points;
                            }

                            if (count) {
                                if (count == num_dims)
                                    res[idx] = sqrt(dist);
                                else
                                    res[idx] = sqrt(dist * num_dims / (double)count);
                            } else
                                res[idx] = NA_REAL;

                            idx += num_points;
                            TGStat::itr_idx(++itr_idx);
                        }
                    }
                }
                rexit();
            }
        }

        while (TGStat::wait_for_kids(3000))
            progress.report(TGStat::itr_idx_sum() - progress.get_elapsed_steps());

        progress.report_last();

        // assemble the answer
        if (tidy) {
            enum { ROW1, ROW2, DIST, num_dims };
            const char *COL_NAMES[num_dims] = { "row1", "row2", "dist" };

            rprotect(answer = RSaneAllocVector(VECSXP, num_dims));

            uint64_t answer_size = 0;

            for (uint64_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                uint64_t idx2 = (ipoint1 + 1) * num_points + ipoint1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (res[idx2] <= threshold)
                        ++answer_size;
                    idx2 += num_points;
                }
            }

            SEXP rcol1, rcol2, rdist, rrownames, rcolnames;

            rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rdist = RSaneAllocVector(REALSXP, answer_size));
            rprotect(rcolnames = RSaneAllocVector(STRSXP, num_dims));
            rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

            for (uint64_t i = 0; i < num_dims; i++)
                SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

            uint64_t idx1 = 0;
            for (uint64_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                uint64_t idx2 = (ipoint1 + 1) * num_points + ipoint1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (res[idx2] <= threshold) {
                        INTEGER(rcol1)[idx1] = ipoint1 + 1;
                        INTEGER(rcol2)[idx1] = ipoint2 + 1;
                        REAL(rdist)[idx1] = res[idx2];
                        INTEGER(rrownames)[idx1] = idx1 + 1;
                        ++idx1;
                    }
                    idx2 += num_points;
                }
            }

            if (_rrownames != R_NilValue) {
                setAttrib(rcol1, R_LevelsSymbol, _rrownames);
                setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
                setAttrib(rcol2, R_LevelsSymbol, _rrownames);
                setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
            }

            SET_VECTOR_ELT(answer, ROW1, rcol1);
            SET_VECTOR_ELT(answer, ROW2, rcol2);
            SET_VECTOR_ELT(answer, DIST, rdist);

            setAttrib(answer, R_NamesSymbol, rcolnames);
            setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
            setAttrib(answer, R_RowNamesSymbol, rrownames);
        } else {
            rprotect(answer = RSaneAllocVector(REALSXP, (uint64_t)num_points * (num_points - 1) / 2));

            uint64_t idx1 = 0;
            double *d = REAL(answer);

            for (uint64_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                uint64_t idx2 = (ipoint1 + 1) * num_points + ipoint1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    d[idx1++] = res[idx2];
                    idx2 += num_points;
                }
            }

            SEXP names = getAttrib(_attrs, R_NamesSymbol);
            for (int i = 0; i < LENGTH(_attrs); i++)
                setAttrib(answer, install(translateChar(STRING_ELT(names, i))), VECTOR_ELT(_attrs, i));
        }
    } catch (TGLException &e) {
        if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
            munmap((char *)shm, shm_sizeof);     // needs to be char * for some versions of Solaris
            shm = (double *)MAP_FAILED;
        }
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    if (!TGStat::is_kid() && shm != (double *)MAP_FAILED) {
        munmap((char *)shm, shm_sizeof);     // needs to be char * for some versions of Solaris
        shm = (double *)MAP_FAILED;
    }
	rreturn(answer);
}

SEXP tgs_dist_blas(SEXP _x, SEXP _attrs, SEXP _tidy, SEXP _threshold, SEXP _rrownames, SEXP _envir)
{
    SEXP answer = R_NilValue;

	try {
        struct Mem {
            double *m{NULL};
            double *mask{NULL};
            double *n{NULL};
            double *res{NULL};
            ~Mem() { free(m); free(mask); free(n); free(res); }
        } mem;

        TGStat tgstat(_envir);

		if ((!isReal(_x) && !isInteger(_x)) || xlength(_x) < 1)
			verror("\"x\" argument must be a matrix of numeric values");

        if (!isLogical(_tidy) || xlength(_tidy) != 1)
            verror("\"tidy\" argument must be a logical value");

        if ((!isReal(_threshold) && !isInteger(_threshold)) || xlength(_threshold) != 1)
            verror("\"threshold\" argument must be a numeric value");

        SEXP rdim = getAttrib(_x, R_DimSymbol);

        if (!isInteger(rdim) || xlength(rdim) != 2)
            verror("\"x\" argument must be a matrix of numeric values");

        bool tidy = asLogical(_tidy);
        double threshold = fabs(asReal(_threshold));

        uint64_t num_points = nrows(_x);
        uint64_t num_dims = ncols(_x);
        int num_points32 = (int)num_points;
        int num_dims32 = (int)num_dims;

        if (num_points < 1 || num_dims < 1)
            verror("\"x\" argument must be a matrix of numeric values");

        uint64_t num_vals = num_points * num_dims;
        uint64_t res_size = num_points * num_points;
        bool nan_in_vals = false;

        // some BLAS implementations ask to align double arrays to 64 for improved efficiency
        if (posix_memalign((void **)&mem.m, 64, sizeof(double) * num_vals))
            verror("%s", strerror(errno));

        if (posix_memalign((void **)&mem.mask, 64, sizeof(double) * num_vals))
            verror("%s", strerror(errno));

        if (posix_memalign((void **)&mem.res, 64, sizeof(double) * res_size))
            verror("%s", strerror(errno));

        for (uint64_t i = 0; i < num_vals; ++i) {
            if ((isReal(_x) && !R_FINITE(REAL(_x)[i])) || (isInteger(_x) && INTEGER(_x)[i] == NA_INTEGER)) {
                mem.m[i] = mem.mask[i] = 0.;
                nan_in_vals = true;
            } else {
                double val = isReal(_x) ? REAL(_x)[i] : INTEGER(_x)[i];
                mem.m[i] = val;
                mem.mask[i] = 1.;
            }
        }

        ProgressReporter progress;
        progress.init(nan_in_vals ? 4 : 2, 1);

        // distance without nans is calculated as:
        //     sqrt(s_x^2 + s_y^2 - 2 * (m %*% t(m)))
        // distance with nans is calculated as:
        //     sqrt((m^2 %*% t(mask) + mask %*% t(m^2) - 2 * (m %*% t(m))) * num_dims / (mask %*% t(mask)))

        // res <- -2 * (m * t(m))
        {
            char uplo = 'L';
            char trans = 'N';
            double alpha = -2;
            double beta = 0;
            F77_NAME(dsyrk)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.m, &num_points32, &beta, mem.res, &num_points32);
            check_interrupt();
            progress.report(1);
        }

        if (nan_in_vals) {
            // m <- m^2
            for (uint64_t i = 0; i < num_vals; ++i)
                mem.m[i] = mem.m[i] * mem.m[i];

            // res <- num_dims * m^2 %*% t(mask) + num_dims * mask %*% t(m^2) + num_dims * res
            {
                char uplo = 'L';
                char trans = 'N';
                double alpha = num_dims;
                double beta = num_dims;
                F77_NAME(dsyr2k)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.m, &num_points32, mem.mask, &num_points32, &beta, mem.res, &num_points32);
                check_interrupt();
                progress.report(1);
            }

            // n <- mask %*% t(mask)
            {
                char uplo = 'L';
                char trans = 'N';
                double alpha = 1;
                double beta = 0;
                if (posix_memalign((void **)&mem.n, 64, sizeof(double) * res_size))
                    verror("%s", strerror(errno));
                F77_NAME(dsyrk)(&uplo, &trans, &num_points32, &num_dims32, &alpha, mem.mask, &num_points32, &beta, mem.n, &num_points32);
                check_interrupt();
                progress.report(1);
            }

            // sqrt(res / n)
            for (uint64_t ipoint1 = 0, idx = 0; ipoint1 < num_points; ++ipoint1) {
                idx += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    mem.res[idx] = sqrt(mem.res[idx] / mem.n[idx]);
                    ++idx;
                }
            }
            check_interrupt();
            progress.report_last();
        } else {
            vector<double> s_v2(num_points, 0.);

            for (uint64_t idim = 0, idx = 0; idim < num_dims; ++idim) {
                for (uint64_t ipoint = 0; ipoint < num_points; ++ipoint) {
                    s_v2[ipoint] += mem.m[idx] * mem.m[idx];
                    ++idx;
                }
            }

            // sqrt(res + s_x^2 + s_y^2)
            for (uint64_t ipoint1 = 0, idx = 0; ipoint1 < num_points; ++ipoint1) {
                idx += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    mem.res[idx] = sqrt(max(mem.res[idx] + s_v2[ipoint1] + s_v2[ipoint2], 0.));
                    ++idx;
                }
            }

            check_interrupt();
            progress.report_last();
        }
        
//{
//SEXP rdims;
//rprotect(answer = RSaneAllocVector(REALSXP, num_points * num_points));
//rprotect(rdims = RSaneAllocVector(INTSXP, 2));
//INTEGER(rdims)[0] = INTEGER(rdims)[1] = num_points;
//setAttrib(answer, R_DimSymbol, rdims);
//memcpy(REAL(answer), mem.res, xlength(answer) * sizeof(REAL(answer)[0]));
//return answer;
//}

        // assemble the answer
        if (tidy) {
            enum { ROW1, ROW2, DIST, num_dims };
            const char *COL_NAMES[num_dims] = { "row1", "row2", "dist" };

            rprotect(answer = RSaneAllocVector(VECSXP, num_dims));

            uint64_t answer_size = 0;

            for (uint64_t ipoint1 = 0, idx = 0; ipoint1 < num_points; ++ipoint1) {
                idx += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (mem.res[idx] <= threshold)
                        ++answer_size;
                    ++idx;
                }
            }

            SEXP rcol1, rcol2, rdist, rrownames, rcolnames;

            rprotect(rcol1 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rcol2 = RSaneAllocVector(INTSXP, answer_size));
            rprotect(rdist = RSaneAllocVector(REALSXP, answer_size));
            rprotect(rcolnames = RSaneAllocVector(STRSXP, num_dims));
            rprotect(rrownames = RSaneAllocVector(INTSXP, answer_size));

            for (uint64_t i = 0; i < num_dims; i++)
                SET_STRING_ELT(rcolnames, i, mkChar(COL_NAMES[i]));

            uint64_t idx1 = 0;
            uint64_t idx2 = 0;

            for (uint64_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx2 += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2) {
                    if (mem.res[idx2] <= threshold) {
                        INTEGER(rcol1)[idx1] = ipoint1 + 1;
                        INTEGER(rcol2)[idx1] = ipoint2 + 1;
                        REAL(rdist)[idx1] = mem.res[idx2];
                        INTEGER(rrownames)[idx1] = idx1 + 1;
                        ++idx1;
                    }
                    ++idx2;
                }
            }

            if (_rrownames != R_NilValue) {
                setAttrib(rcol1, R_LevelsSymbol, _rrownames);
                setAttrib(rcol1, R_ClassSymbol, mkString("factor"));
                setAttrib(rcol2, R_LevelsSymbol, _rrownames);
                setAttrib(rcol2, R_ClassSymbol, mkString("factor"));
            }

            SET_VECTOR_ELT(answer, ROW1, rcol1);
            SET_VECTOR_ELT(answer, ROW2, rcol2);
            SET_VECTOR_ELT(answer, DIST, rdist);

            setAttrib(answer, R_NamesSymbol, rcolnames);
            setAttrib(answer, R_ClassSymbol, mkString("data.frame"));
            setAttrib(answer, R_RowNamesSymbol, rrownames);
        } else {
            rprotect(answer = RSaneAllocVector(REALSXP, (uint64_t)num_points * (num_points - 1) / 2));

            uint64_t idx1 = 0;
            uint64_t idx2 = 0;
            double *d = REAL(answer);

            for (uint64_t ipoint1 = 0; ipoint1 < num_points; ++ipoint1) {
                idx2 += ipoint1 + 1;
                for (uint64_t ipoint2 = ipoint1 + 1; ipoint2 < num_points; ++ipoint2)
                    d[idx1++] = mem.res[idx2++];
            }

            SEXP names = getAttrib(_attrs, R_NamesSymbol);
            for (int i = 0; i < LENGTH(_attrs); i++)
                setAttrib(answer, install(translateChar(STRING_ELT(names, i))), VECTOR_ELT(_attrs, i));
        }
    } catch (TGLException &e) {
		rerror("%s", e.msg());
    } catch (const bad_alloc &e) {
        rerror("Out of memory");
    }

    rreturn(answer);
}

}

