// SPDX-FileCopyrightText: 2026 Jonathan K. Hartlaub <contact@jkhartlaub.com>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file prof_aniso.c
 * @author Jonathan K. Hartlaub
 * @section Description: Profiling driver for epsteinZetaAniso. Two cases:
 * "zero": identity lattice, x = y = 0, alpha = (n, 0, ..., 0).
 * "generic": oblique lattice, generic x and y, alpha spread over all directions.
 * Single mode repeats one evaluation for gprof/perf, sweep mode prints timings
 * and values for d = 2, 3, both cases and a range of |alpha|.
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bits/time.h"
#include "epsteinZeta.h"

enum { DIM_MAX = 6, BATCHES = 3 };

enum prof_case { CASE_ZERO, CASE_GENERIC };

static const char *caseNames[] = {"zero", "generic"};

/**
 * @brief Fills lattice, shift, wave vector and multi-index of a profiling case.
 * @param[in] c: profiling case.
 * @param[in] dim: dimension.
 * @param[in] n: total |alpha|.
 * @param[out] a: lattice matrix, dim x dim.
 * @param[out] x: shift vector.
 * @param[out] y: wave vector.
 * @param[out] alpha: multi-index.
 */
static void setup_case(enum prof_case c, unsigned int dim, unsigned int n, double *a,
                       double *x, double *y, unsigned int *alpha) {
    static const double xGen[DIM_MAX] = {0.31, 0.17, 0.23, 0.11, 0.37, 0.05};
    static const double yGen[DIM_MAX] = {0.13, 0.29, 0.07, 0.41, 0.19, 0.03};
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
            a[(dim * i) + j] = (i == j) ? 1. : 0.;
        }
        if (c == CASE_ZERO) {
            x[i] = 0.;
            y[i] = 0.;
            alpha[i] = (i == 0) ? n : 0;
        } else {
            // upper off-diagonal entries make the lattice oblique
            if (i + 1 < dim) {
                a[(dim * i) + i + 1] = 0.3;
            }
            x[i] = xGen[i];
            y[i] = yGen[i];
            alpha[i] = (n / dim) + ((i < n % dim) ? 1 : 0);
        }
    }
}

/**
 * @brief Wall time in seconds.
 * @return monotonic clock time.
 */
static double now(void) {
    struct timespec ts;
    (void)clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (1e-9 * (double)ts.tv_nsec);
}

/**
 * @brief Time per evaluation of one case, minimum over BATCHES batches of at
 * least minSec / BATCHES seconds each.
 * @param[in] c: profiling case.
 * @param[in] dim: dimension.
 * @param[in] n: total |alpha|.
 * @param[in] nu: exponent.
 * @param[in] minSec: total time budget.
 * @param[out] val: value of epsteinZetaAniso.
 * @return seconds per evaluation.
 */
static double time_case(enum prof_case c, unsigned int dim, unsigned int n,
                        double nu, double minSec, double complex *val) {
    double a[DIM_MAX * DIM_MAX];
    double x[DIM_MAX];
    double y[DIM_MAX];
    unsigned int alpha[DIM_MAX];
    setup_case(c, dim, n, a, x, y, alpha);

    // calibrate repetitions from one evaluation
    double t0 = now();
    *val = epsteinZetaAniso(nu, dim, a, x, y, alpha);
    double tOne = now() - t0;
    double batchSec = minSec / BATCHES;
    long reps = (tOne > 0.) ? (long)(batchSec / tOne) + 1 : 1000000;

    double best = 1e300;
    volatile double complex sink = 0.;
    for (int b = 0; b < BATCHES; b++) {
        t0 = now();
        for (long r = 0; r < reps; r++) {
            sink += epsteinZetaAniso(nu, dim, a, x, y, alpha);
        }
        double t = (now() - t0) / (double)reps;
        best = (t < best) ? t : best;
    }
    (void)sink;
    return best;
}

/**
 * @brief Repeats one case, for gprof or perf.
 * @param[in] argc: argument count.
 * @param[in] argv: dim, n, case, iterations, optionally nu.
 * @return 0 on success.
 */
static int run_single(int argc, char **argv) {
    unsigned int dim = (unsigned int)strtoul(argv[1], NULL, 10);
    unsigned int n = (unsigned int)strtoul(argv[2], NULL, 10);
    enum prof_case c = (strcmp(argv[3], "zero") == 0) ? CASE_ZERO : CASE_GENERIC;
    long iterations = strtol(argv[4], NULL, 10);
    double nu = (argc > 5) ? strtod(argv[5], NULL) : 0.5;
    if (dim < 1 || dim > DIM_MAX) {
        (void)fprintf(stderr, "dim must be in 1..%d\n", DIM_MAX);
        return 1;
    }

    double a[DIM_MAX * DIM_MAX];
    double x[DIM_MAX];
    double y[DIM_MAX];
    unsigned int alpha[DIM_MAX];
    setup_case(c, dim, n, a, x, y, alpha);

    printf("dim = %u, case = %s, nu = %g, iterations = %ld, alpha = (", dim,
           caseNames[c], nu, iterations);
    for (unsigned int i = 0; i < dim; i++) {
        printf((i + 1 < dim) ? "%u, " : "%u)\n", alpha[i]);
    }

    double complex sum = 0.;
    double t0 = now();
    for (long r = 0; r < iterations; r++) {
        sum += epsteinZetaAniso(nu, dim, a, x, y, alpha);
    }
    double t = (now() - t0) / (double)iterations;
    printf("value    = %.16e %+.16e i\n", creal(sum) / (double)iterations,
           cimag(sum) / (double)iterations);
    printf("time     = %.3e s per evaluation\n", t);
    return 0;
}

/**
 * @brief Prints time per evaluation and values for d = 2, 3, both cases and a
 * range of |alpha|.
 * @param[in] minSec: time budget per configuration.
 * @param[in] nu: exponent.
 * @return 0 on success.
 */
static int run_sweep(double minSec, double nu) {
    static const unsigned int ns[] = {0, 1,  2,  3,  4,  5,  6,  7,  8,
                                      9, 10, 11, 12, 14, 16, 20, 24, 30};
    static const unsigned int nMax[] = {0, 0, 30, 20};
    const unsigned int nCount = sizeof(ns) / sizeof(ns[0]);
    double times[4][sizeof(ns) / sizeof(ns[0])];
    double complex vals[4][sizeof(ns) / sizeof(ns[0])];

    for (unsigned int dim = 2; dim <= 3; dim++) {
        for (unsigned int c = 0; c < 2; c++) {
            unsigned int col = (2 * (dim - 2)) + c;
            for (unsigned int j = 0; j < nCount; j++) {
                times[col][j] = -1.;
                if (ns[j] > nMax[dim]) {
                    continue;
                }
                times[col][j] = time_case((enum prof_case)c, dim, ns[j], nu, minSec,
                                          &vals[col][j]);
                (void)fprintf(stderr, "\r%uD %-7s n = %2u", dim, caseNames[c],
                              ns[j]);
            }
        }
    }
    (void)fprintf(stderr, "\n");

    printf("nu = %g, time per evaluation [us] (relative to n = 0)\n", nu);
    printf("%3s  %18s  %18s  %18s  %18s\n", "n", "2D zero", "2D generic", "3D zero",
           "3D generic");
    for (unsigned int j = 0; j < nCount; j++) {
        printf("%3u", ns[j]);
        for (unsigned int col = 0; col < 4; col++) {
            if (times[col][j] < 0.) {
                printf("  %18s", "-");
            } else {
                printf("  %10.1f (%5.1f)", 1e6 * times[col][j],
                       times[col][j] / times[col][0]);
            }
        }
        printf("\n");
    }

    printf("\nvalues (diff against a reference run to catch accuracy changes)\n");
    for (unsigned int col = 0; col < 4; col++) {
        for (unsigned int j = 0; j < nCount; j++) {
            if (times[col][j] >= 0.) {
                printf("%uD %-7s %2u  %+.16e %+.16e\n", 2 + (col / 2),
                       caseNames[col % 2], ns[j], creal(vals[col][j]),
                       cimag(vals[col][j]));
            }
        }
    }
    return 0;
}

/**
 * @brief Dispatches single and sweep mode.
 * @return 0 on success.
 */
int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "sweep") == 0) {
        double minSec = (argc > 2) ? strtod(argv[2], NULL) : 0.3;
        double nu = (argc > 3) ? strtod(argv[3], NULL) : 0.5;
        return run_sweep(minSec, nu);
    }
    if (argc < 5) {
        (void)fprintf(stderr,
                      "usage: %s <dim> <n> <zero|generic> <iterations> [nu]\n"
                      "       %s sweep [seconds per configuration] [nu]\n",
                      argv[0], argv[0]);
        return 1;
    }
    return run_single(argc, argv);
}
