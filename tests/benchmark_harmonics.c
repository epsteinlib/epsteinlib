// SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file benchmark_harmonic.c
 * @author Jonathan Busse
 * @date 06/06/2024
 * @section Description: Benchmark harmonic_h in 1D–4D over a tensor-product
 * y-grid for all k = 0, ..., floor(|alpha|/2). One CSV file per dimension.
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

// Highest dimension benchmarked.
enum { DIM_MAX = 4U };

// Fixed |alpha| used in all dimensions. Divisible by 1, 2, 3 and 4, so that
// alpha = (|alpha|/d, ..., |alpha|/d) is a valid multi-index in every dimension.
enum { ALPHA_ABS = 24U };

// Maximum k = floor(ALPHA_ABS / 2).
#define K_MAX (ALPHA_ABS / 2U)

// Step sizes per dimension, indexed by d-1. Decrease to increase runtime.
// To decrease the runtime and filesizes, choose the following dummy stepsizes:
static const double YINC[DIM_MAX] = {0.05, 0.1, 0.25, 0.5};
// To replicate the results in the derivatives article, choose the following
// stepsizes:
// static const double YINC[DIM_MAX] = {0.001, 0.01, 0.05, 0.25};

// Lower bound (inclusive) for each y component.
#define YMIN (-0.5)

// Upper bound (inclusive) for each y component.
#define YMAX 0.5

// Timing repetitions per dimension, indexed by d-1. Increase for stable measurement.
static const int TIMING_ITERATIONS[DIM_MAX] = {10000, 1000, 100, 10};

// Precomputation timing repetitions per dimension, indexed by d-1.
static const int PRECOMPUTE_ITERATIONS[DIM_MAX] = {1000, 100, 10, 10};

/**
 * @brief Benchmarks harmonic_h over a tensor-product grid of y values and
 * k = 0..K_MAX in dimension dim.
 *
 * alpha = (ALPHA_ABS/dim, ..., ALPHA_ABS/dim). Precomputed coefficients and
 * exponents are shared across all y evaluations. Results are written to
 * benchmark_harmonic_<dim>D.csv, one line per (k, y):
 * dim, k, y_1, ..., y_dim, alpha_1, ..., alpha_dim, h_{alpha,k}(y),
 * elapsed_precompute_seconds, elapsed_time_seconds.
 *
 * @param[in] dim: dimension, 1 <= dim <= DIM_MAX.
 * @return 0 on success, non-zero on failure.
 */
static int benchmark_harmonic(unsigned int dim) { // NOLINT
    char path[MAX_PATH_LENGTH];
    if (snprintf(path, MAX_PATH_LENGTH, "%s/benchmark_harmonic_%uD.csv", BASE_PATH,
                 dim) >= MAX_PATH_LENGTH) {
        (void)fprintf(stderr, "Error: path too long\n");
        return 1;
    }

    const double inc = YINC[dim - 1];
    const unsigned int points = (unsigned int)(((YMAX - YMIN) / inc) + 0.5) + 1;

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

    double *coeffs = malloc(coeffs_size * sizeof(double));
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

    double z[DIM_MAX];
    unsigned int index[DIM_MAX];
    double result = NAN;
    double elapsed = NAN;
    double elapsed_precompute = NAN;
    double sink = 0.;
    clock_t t0;
    clock_t t1;

    t0 = clock();
    for (int r = 0; r < PRECOMPUTE_ITERATIONS[dim - 1]; r++) {
        precompute_harmonic_h_inner_chunk_size(ALPHA_ABS, K_MAX, dim, alpha,
                                               chunk_offset, valid_count);
        precompute_harmonic_h_inner_sum(ALPHA_ABS, dim, alpha, chunk_offset, coeffs,
                                        exponents);
    }
    t1 = clock();
    elapsed_precompute =
        ((double)(t1 - t0)) / CLOCKS_PER_SEC / PRECOMPUTE_ITERATIONS[dim - 1];

    unsigned long long grid_size = 1;
    for (unsigned int i = 0; i < dim; i++) {
        grid_size *= points;
    }

    for (unsigned int k = 0; k <= K_MAX; k++) {
        for (unsigned int i = 0; i < dim; i++) {
            index[i] = 0;
        }
        for (unsigned long long p = 0; p < grid_size; p++) {
            for (unsigned int i = 0; i < dim; i++) {
                z[i] = YMIN + (index[i] * inc);
            }

            sink = 0.;
            t0 = clock();
            for (int r = 0; r < TIMING_ITERATIONS[dim - 1]; r++) {
                sink += harmonic_h(k, dim, z, ALPHA_ABS, chunk_offset, valid_count,
                                   coeffs, exponents);
            }
            t1 = clock();
            (void)sink;

            result = harmonic_h(k, dim, z, ALPHA_ABS, chunk_offset, valid_count,
                                coeffs, exponents);
            elapsed =
                ((double)(t1 - t0)) / CLOCKS_PER_SEC / TIMING_ITERATIONS[dim - 1];

            (void)fprintf(file, "%u,%u", dim, k);
            for (unsigned int i = 0; i < dim; i++) {
                (void)fprintf(file, ",%.16lf", z[i]);
            }
            for (unsigned int i = 0; i < dim; i++) {
                (void)fprintf(file, ",%u", alpha[i]);
            }
            (void)fprintf(file, ",%.16lf,%.16lf,%.16lf\n", result,
                          elapsed_precompute, elapsed);

            // advance the tensor-product grid index, last component fastest
            for (int i = (int)dim - 1; i >= 0; i--) {
                index[i]++;
                if (index[i] < points) {
                    break;
                }
                index[i] = 0;
            }
        }
        printf("%uD: k = %u\n", dim, k);
    }

    free(chunk_offset);
    free(valid_count);
    free(coeffs);
    free(exponents);

    if (fclose(file) != 0) {
        (void)fprintf(stderr, "Error closing file: %d\n", errno);
        return 1;
    }
    printf("%uD benchmark complete.\n", dim);
    return 0;
}

/**
 * @brief Main function to run all harmonic polynomial benchmark tests.
 * @return Number of failed benchmark executions.
 */
int main(void) {
    int failed = 0;
    for (unsigned int dim = 1; dim <= DIM_MAX; dim++) {
        failed += benchmark_harmonic(dim);
    }
    return failed;
}
