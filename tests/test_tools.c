/*
 * SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
 *
 * SPDX-License-Identifier: AGPL-3.0-only
 */

// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "../src/tools.h"
#include "utils.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

#ifndef BASE_PATH
#define BASE_PATH "csv"
#endif

#ifndef MAX_DIM
#define MAX_DIM 12
#endif

/*!
 * @brief Tests matrix invert-transpose against Mathematica reference values.
 * CSV format per line: dim, m[0..dim²-1], ref_transpose_inv[0..dim²-1]
 *
 * @return number of failed tests.
 */
int test_matrix_transpose_inverse(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/svd_ref.csv", BASE_PATH);
    if (result < 0 || result >= (int)sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    char line[4096];
    int testsPassed = 0;
    int totalTests = 0;
    const double tol = 5e-14;
    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        char *ptr = line;
        char *endptr = NULL;

        // parse dimension
        unsigned int dim = (unsigned int)strtoul(ptr, &endptr, 10);
        if (endptr == ptr || *endptr != ',' || dim == 0 || dim > MAX_DIM) {
            printf("Error reading line: %s\n", line);
            printf("Invalid dimension: %u\n", dim);
            continue;
        }
        ptr = endptr + 1;

        unsigned int n = dim * dim;

        // parse input matrix m (row-major, dim x dim)
        double m[dim * dim]; // NOLINT
        int ok = 1;
        for (unsigned int i = 0; i < n; i++) {
            m[i] = strtod(ptr, &endptr);
            if (endptr == ptr || *endptr != ',') {
                ok = 0;
                break;
            }
            ptr = endptr + 1;
        }
        if (!ok) {
            printf("Error reading line: %s\n", line);
            printf("Could not parse input matrix (expected %u values)\n", n);
            continue;
        }

        // parse reference: transpose(inverse(m))
        double ref[dim * dim];
        for (unsigned int i = 0; i < n; i++) {
            ref[i] = strtod(ptr, &endptr);
            if (endptr == ptr) {
                ok = 0;
                break;
            }
            if (i < n - 1) {
                if (*endptr != ',') {
                    ok = 0;
                    break;
                }
                ptr = endptr + 1;
            }
        }
        if (!ok) {
            printf("Error reading line: %s\n", line);
            printf("Could not parse reference matrix (expected %u values)\n", n);
            continue;
        }

        // compute transpose(inverse(m)) via invert()
        double m_copy[dim * dim];
        double m_fourier[dim * dim];
        int p[dim];
        for (unsigned int i = 0; i < n; i++) {
            m_copy[i] = m[i];
        }
        invert(dim, m_copy, p, m_fourier);
        transpose(dim, m_fourier);

        // infinity norm of (m_fourier - ref)
        double diff[dim * dim];
        for (unsigned int i = 0; i < n; i++) {
            diff[i] = m_fourier[i] - ref[i];
        }
        double err = inf_norm(dim, diff);

        errMin = (errMin < err) ? errMin : err;
        errMax = (errMax > err) ? errMax : err;
        errSum += err;

        if (err < tol) {
            testsPassed++;
        } else {
            printf("\n");
            printf("Warning! ");
            printf("matrix transpose inverse (dim=%u):\n", dim);
            printf("inf_norm(m_fourier - ref): %E !< %E  (tolerance)\n", err, tol);
            printf("\n");
            printMatrixUnitTest("m:", m, dim);
            printMatrixUnitTest("reci:", m_fourier, dim);
            printMatrixUnitTest("ref:", ref, dim);
            printf("\n");
        }
        totalTests++;
    }

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / totalTests);
    printf("\n");
    return totalTests - testsPassed;
}

/*!
 * @brief Main function to run all set zeta derivatives function tests.
 *
 * @return number of failed tests.
 */
int main() {
    int failed = 0;
    failed += run_timed_test(test_matrix_transpose_inverse);
    return failed;
}
