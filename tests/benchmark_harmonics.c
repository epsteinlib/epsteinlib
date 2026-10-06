// SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file benchmark_harmonic.c
 * @author Jonathan Busse
 * @date 06/06/2024
 * @section Description: Benchmark harmonic_h in 1D–4D over a tensor-product
 * y-grid for all k = 0, ..., floor(|alpha|/2), with the adaptive
 * double/double-double evaluation, and in 2D for alpha = (n, 2n) with forced
 * double evaluation. One CSV file per dimension or multi-index. Every point is
 * timed. Every reported time is the median of TIMING_REPEATS independent
 * measurements, each of which is itself an average over an inner repetition
 * loop.
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../src/harmonics.h"
#include "utils.h"

#define BASE_PATH "tests/csv"

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

/** @brief Highest dimension benchmarked in benchmark_harmonic. */
enum { DIM_MAX = 4U };

/** @brief Dimension of the anisotropy stability benchmark. */
enum { DIM_STAB = 2U };

/**
 * @brief Fixed |alpha| used in all dimensions by benchmark_harmonic.
 *
 * Divisible by 1, 2, 3 and 4, so that alpha = (|alpha|/d, ..., |alpha|/d) is a
 * valid multi-index in every dimension.
 */
enum { ALPHA_ABS = 24U };

/** @brief Maximum k = floor(ALPHA_ABS / 2) for benchmark_harmonic. */
#define K_MAX (ALPHA_ABS / 2U)

/**
 * @brief Largest n in the sweep alpha = (n, 2n), n = 0, ..., N_STAB_MAX, of
 * benchmark_harmonic_stab_2D.
 *
 * Here |alpha| = 3n varies with n, so the coefficient table is rebuilt for
 * every multi-index and the maximum k is floor(3n / 2) rather than K_MAX.
 */
enum { N_STAB_MAX = 10U };

/**
 * @brief tol argument of harmonic_h for the adaptive evaluation: double if
 * the computed condition number is at most HARMONIC_DD_COND, else
 * double-double.
 */
#define TOL_ADAPTIVE 0.

/** @brief tol argument of harmonic_h that forces the double evaluation. */
#define TOL_DOUBLE INFINITY

/** @brief Lower bound (inclusive) for each y component. */
#define YMIN (-1.)

/** @brief Upper bound (inclusive) for each y component. */
#define YMAX 1.

/** @brief (Dummy) step size of benchmark_harmonic */
static const double YINC[DIM_MAX] = {1. / 64, 1. / 8, 1. / 2, 1.};
// To replicate the results in the derivatives article, use instead
// static const double YINC[DIM_MAX] = {1. / 512, 1. / 32, 1. / 8, 1. / 4};

/** @brief (Dummy) step size of benchmark_harmonic_stab_2D in both components. */
#define YINCSTAB2D (1. / 4)
// To replicate the results in the derivatives article, use instead
// #define YINCSTAB2D (1. / 64)

/**
 * @brief Inner repetitions per timing sample of benchmark_harmonic, indexed
 * by d-1.
 *
 * Chosen so that one sample takes roughly 10 microseconds at k = 0 (about
 * 2e-8 s, 1.4e-7 s, 1e-6 s and 1e-5 s per call in 1D-4D). Since every point
 * is timed, larger values multiply the total runtime directly.
 */
static const int TIMING_ITERATIONS[DIM_MAX] = {500, 100, 10, 2};

/**
 * @brief Precomputation timing repetitions of benchmark_harmonic per
 * dimension, indexed by d-1.
 */
static const int PRECOMPUTE_ITERATIONS[DIM_MAX] = {1000, 100, 10, 10};

/**
 * @brief Inner repetitions per timing sample of benchmark_harmonic_stab_2D.
 *
 * At most about 2e-7 s per call (alpha = (10, 20), k = 0), so at most about
 * 20 microseconds per sample.
 */
enum { TIMING_ITERATIONS_STAB = 100 };

/**
 * @brief Number of independent timing samples taken per reported time.
 *
 * Each sample is an average over the corresponding inner repetition loop; the
 * value written to the CSV is the median over these samples, taken as the
 * element at index TIMING_REPEATS / 2 of the sorted samples.
 */
enum { TIMING_REPEATS = 10 };

/**
 * @brief Elapsed wall time in seconds between two timespec_get readings.
 *
 * Seconds and nanoseconds are subtracted separately, since a double holding
 * the absolute time since the epoch resolves only about 0.2 microseconds.
 */
static double elapsed_seconds(const struct timespec *t0, const struct timespec *t1) {
    return (double)(t1->tv_sec - t0->tv_sec) +
           ((double)(t1->tv_nsec - t0->tv_nsec) * 1e-9);
}

/**
 * @brief Median over TIMING_REPEATS samples of the time per call of
 * harmonic_h, each sample averaged over iterations calls.
 * @return time per call in seconds.
 */
static double time_harmonic_h(unsigned int k, unsigned int dim, const double *z,
                              unsigned int alphaAbs,
                              const unsigned long long *chunk_offset,
                              const unsigned long long *valid_count,
                              const double *coeffs, const unsigned int *exponents,
                              double tol, int iterations) {
    double samples[TIMING_REPEATS];
    volatile double sink = 0.;
    struct timespec t0;
    struct timespec t1;
    for (int s = 0; s < TIMING_REPEATS; s++) {
        if (timespec_get(&t0, TIME_UTC) != TIME_UTC) {
            return NAN;
        }
        for (int r = 0; r < iterations; r++) {
            sink += harmonic_h(k, dim, z, alphaAbs, chunk_offset, valid_count,
                               coeffs, exponents, tol);
        }
        if (timespec_get(&t1, TIME_UTC) != TIME_UTC) {
            return NAN;
        }
        samples[s] = elapsed_seconds(&t0, &t1) / iterations;
    }
    (void)sink;
    sort(samples, TIMING_REPEATS);
    return samples[TIMING_REPEATS / 2];
}

/**
 * @brief Median over TIMING_REPEATS samples of the time of one full
 * precomputation of chunk sizes and inner sums, each sample averaged over
 * iterations runs.
 * @return time per precomputation in seconds.
 */
static double time_precompute(unsigned int alphaAbs, unsigned int kMax,
                              unsigned int dim, const unsigned int *alpha,
                              unsigned long long *chunk_offset,
                              unsigned long long *valid_count, double *coeffs,
                              unsigned int *exponents, int iterations) {
    double samples[TIMING_REPEATS];
    struct timespec t0;
    struct timespec t1;
    for (int s = 0; s < TIMING_REPEATS; s++) {
        if (timespec_get(&t0, TIME_UTC) != TIME_UTC) {
            return NAN;
        }
        for (int r = 0; r < iterations; r++) {
            precompute_harmonic_h_inner_chunk_size(alphaAbs, kMax, dim, alpha,
                                                   chunk_offset, valid_count);
            precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_offset,
                                            coeffs, exponents);
        }
        if (timespec_get(&t1, TIME_UTC) != TIME_UTC) {
            return NAN;
        }
        samples[s] = elapsed_seconds(&t0, &t1) / iterations;
    }
    sort(samples, TIMING_REPEATS);
    return samples[TIMING_REPEATS / 2];
}

/**
 * @brief Evaluates and times harmonic_h at every point of a tensor-product
 * y-grid in [YMIN, YMAX]^dim with step inc, for k = 0, ..., kMax, and writes
 * one CSV line per (k, y):
 * dim, k, y_1, ..., y_dim, alpha_1, ..., alpha_dim, h_{alpha,k}(y),
 * elapsed_precompute_seconds, elapsed_time_seconds.
 * @return maximum evaluation time over all points.
 */
static double sweep_grid(FILE *file, unsigned int dim, const unsigned int *alpha,
                         unsigned int alphaAbs, unsigned int kMax, double inc,
                         const unsigned long long *chunk_offset,
                         const unsigned long long *valid_count, const double *coeffs,
                         const unsigned int *exponents, double tol, int iterations,
                         double elapsed_precompute) {
    const unsigned int points = (unsigned int)(((YMAX - YMIN) / inc) + 0.5) + 1;
    unsigned long long grid_size = 1;
    for (unsigned int i = 0; i < dim; i++) {
        grid_size *= points;
    }

    double z[DIM_MAX];
    unsigned int index[DIM_MAX];
    double elapsedMax = 0.;

    for (unsigned int k = 0; k <= kMax; k++) {
        for (unsigned int i = 0; i < dim; i++) {
            index[i] = 0;
        }
        for (unsigned long long p = 0; p < grid_size; p++) {
            for (unsigned int i = 0; i < dim; i++) {
                z[i] = YMIN + (index[i] * inc);
            }

            double elapsed =
                time_harmonic_h(k, dim, z, alphaAbs, chunk_offset, valid_count,
                                coeffs, exponents, tol, iterations);
            elapsedMax = (elapsed > elapsedMax) ? elapsed : elapsedMax;
            double result = harmonic_h(k, dim, z, alphaAbs, chunk_offset,
                                       valid_count, coeffs, exponents, tol);

            (void)fprintf(file, "%u,%u", dim, k);
            for (unsigned int i = 0; i < dim; i++) {
                (void)fprintf(file, ",%.17g", z[i]);
            }
            for (unsigned int i = 0; i < dim; i++) {
                (void)fprintf(file, ",%u", alpha[i]);
            }
            (void)fprintf(file, ",%.17g,%.3g,%.3g\n", result, elapsed_precompute,
                          elapsed);

            // advance the tensor-product grid index, last component fastest
            for (int i = (int)dim - 1; i >= 0; i--) {
                index[i]++;
                if (index[i] < points) {
                    break;
                }
                index[i] = 0;
            }
        }
    }
    return elapsedMax;
}

/**
 * @brief Benchmarks the adaptive harmonic_h (tol = TOL_ADAPTIVE) over a
 * tensor-product grid of y values and k = 0..K_MAX in dimension dim.
 *
 * alpha = (ALPHA_ABS/dim, ..., ALPHA_ABS/dim). Precomputed coefficients and
 * exponents are shared across all y evaluations. Results are written to
 * benchmark_harmonic_<dim>D.csv in the layout of sweep_grid.
 *
 * @param[in] dim: dimension, 1 <= dim <= DIM_MAX.
 * @return 0 on success, non-zero on failure.
 */
static int benchmark_harmonic(unsigned int dim) {

    char path[MAX_PATH_LENGTH];
    if (snprintf(path, MAX_PATH_LENGTH, "%s/benchmark_harmonic_%uD.csv", BASE_PATH,
                 dim) >= MAX_PATH_LENGTH) {
        (void)fprintf(stderr, "Error: path too long\n");
        return 1;
    }

    unsigned int alpha[DIM_MAX];
    for (unsigned int i = 0; i < dim; i++) {
        alpha[i] = ALPHA_ABS / dim;
    }

    unsigned long long *chunk_offset =
        malloc((K_MAX + 1) * sizeof(unsigned long long));
    unsigned long long *valid_count =
        malloc((K_MAX + 1) * sizeof(unsigned long long));
    if (chunk_offset == NULL || valid_count == NULL) {
        (void)fprintf(stderr, "Error: allocation failed\n");
        free(chunk_offset);
        free(valid_count);
        return 1;
    }

    unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
        ALPHA_ABS, K_MAX, dim, alpha, chunk_offset, valid_count);

    double *coeffs = malloc(HARMONIC_COEFF_STRIDE * coeffs_size * sizeof(double));
    unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
    if (coeffs == NULL || exponents == NULL) {
        (void)fprintf(stderr, "Error: allocation failed\n");
        free(chunk_offset);
        free(valid_count);
        free(coeffs);
        free(exponents);
        return 1;
    }

    FILE *file = open_file(path, "w");
    if (file == NULL) {
        free(chunk_offset);
        free(valid_count);
        free(coeffs);
        free(exponents);
        return 1;
    }

    double elapsed_precompute =
        time_precompute(ALPHA_ABS, K_MAX, dim, alpha, chunk_offset, valid_count,
                        coeffs, exponents, PRECOMPUTE_ITERATIONS[dim - 1]);

    double elapsedMax =
        sweep_grid(file, dim, alpha, ALPHA_ABS, K_MAX, YINC[dim - 1], chunk_offset,
                   valid_count, coeffs, exponents, TOL_ADAPTIVE,
                   TIMING_ITERATIONS[dim - 1], elapsed_precompute);

    free(chunk_offset);
    free(valid_count);
    free(coeffs);
    free(exponents);

    if (fclose(file) != 0) {
        (void)fprintf(stderr, "Error closing file: %d\n", errno);
        return 1;
    }
    printf("%uD benchmark complete: tpre = %.3g tmax = %.3g\n", dim,
           elapsed_precompute, elapsedMax);
    return 0;
}

/**
 * @brief Benchmarks harmonic_h with forced double evaluation (tol =
 * TOL_DOUBLE) in 2D over a tensor-product y-grid and all
 * k = 0..floor(|alpha|/2), for the anisotropic multi-indices
 * alpha = (n, 2n), n = 0, ..., N_STAB_MAX.
 *
 * Unlike benchmark_harmonic, |alpha| = 3n varies across the sweep, so the
 * coefficient table is reallocated for each n. Results are written to
 * benchmark_harmonic_stab_<alpha_1>_<alpha_2>_2D.csv in the layout of
 * sweep_grid.
 *
 * @return 0 on success, non-zero on failure.
 */
static int benchmark_harmonic_stab_2D(void) {
    const unsigned int dim = DIM_STAB;

    for (unsigned int n = 0; n <= N_STAB_MAX; n++) {
        unsigned int alpha[DIM_STAB] = {n, 2 * n};
        const unsigned int alphaAbs = alpha[0] + alpha[1];
        const unsigned int kMax = alphaAbs / 2U;

        char path[MAX_PATH_LENGTH];
        if (snprintf(path, MAX_PATH_LENGTH,
                     "%s/benchmark_harmonic_stab_%u_%u_2D.csv", BASE_PATH, alpha[0],
                     alpha[1]) >= MAX_PATH_LENGTH) {
            (void)fprintf(stderr, "Error: path too long\n");
            return 1;
        }

        unsigned long long *chunk_offset =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long *valid_count =
            malloc((kMax + 1) * sizeof(unsigned long long));
        if (chunk_offset == NULL || valid_count == NULL) {
            (void)fprintf(stderr, "Error: allocation failed\n");
            free(chunk_offset);
            free(valid_count);
            return 1;
        }

        unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
            alphaAbs, kMax, dim, alpha, chunk_offset, valid_count);

        double *coeffs =
            malloc(HARMONIC_COEFF_STRIDE * coeffs_size * sizeof(double));
        unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
        if (coeffs == NULL || exponents == NULL) {
            (void)fprintf(stderr, "Error: allocation failed\n");
            free(chunk_offset);
            free(valid_count);
            free(coeffs);
            free(exponents);
            return 1;
        }

        FILE *file = open_file(path, "w");
        if (file == NULL) {
            free(chunk_offset);
            free(valid_count);
            free(coeffs);
            free(exponents);
            return 1;
        }

        double elapsed_precompute =
            time_precompute(alphaAbs, kMax, dim, alpha, chunk_offset, valid_count,
                            coeffs, exponents, PRECOMPUTE_ITERATIONS[dim - 1]);

        double elapsedMax =
            sweep_grid(file, dim, alpha, alphaAbs, kMax, YINCSTAB2D, chunk_offset,
                       valid_count, coeffs, exponents, TOL_DOUBLE,
                       TIMING_ITERATIONS_STAB, elapsed_precompute);

        free(chunk_offset);
        free(valid_count);
        free(coeffs);
        free(exponents);

        if (fclose(file) != 0) {
            (void)fprintf(stderr, "Error closing file: %d\n", errno);
            return 1;
        }
        printf("2D stability (double): alpha = (%u,%u), |alpha| = %u tpre = %.3g "
               "tmax = %.3g done.\n",
               alpha[0], alpha[1], alphaAbs, elapsed_precompute, elapsedMax);
    }

    printf("2D stability benchmark complete.\n");
    return 0;
}

/**
 * @brief Main function to run all harmonic polynomial benchmark tests, in the
 * order of the tables in the article: forced double first, adaptive second.
 * @return Number of failed benchmark executions.
 */
int main(void) {
    int failed = benchmark_harmonic_stab_2D();
    for (unsigned int dim = 1; dim <= DIM_MAX; dim++) {
        failed += benchmark_harmonic(dim);
    }
    return failed;
}
