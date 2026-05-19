// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file test_harmonics.c
 * @brief Compares the harmonic polynomials to high-precision benchmark values.
 */

#include "../src/crandall.h"
#include "../src/harmonics.h"
#include "../src/tools.h"
#include "utils.h"
#include "wrappers.h"
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

/*!
 * @brief Benchmarks recursive coefficients
 * @return number of failed tests.
 */
int test_coeffs_c_inner(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/coeffs_c_inner_ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double num;
    double ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    double tol = pow(10, -16);

    unsigned int *n = malloc(sizeof(unsigned int));
    unsigned int *i = malloc(sizeof(unsigned int));
    unsigned int *k = malloc(sizeof(unsigned int));
    unsigned int *dim = malloc(sizeof(unsigned int));

    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: n, i, k ,d, result
        scanResult = sscanf(line, "%u,%u,%u,%u,%lf", // NOLINT
                            n, i, k, dim, refRead);

        if (scanResult != 5) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 5\n", scanResult);
            continue;
        }

        ref = refRead[0];
        num = coeffs_c_inner(*n, *i, *k, *dim);

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("coeffs_c    ");
            printf(" %0*.16lf (this implementation) \n\t\t    ≠ "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E ≮ %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("n: %u,  i: %u,  k: %u,  d: %u", *n, *i, *k, *dim);
            printf("\n");
        }
        totalTests++;
    }

    free(n);
    free(i);
    free(k);
    free(dim);
    free(refRead);

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
 * @brief Benchmarks inner harmonic h function call.
 * @return number of failed tests.
 */
int test_harmonic_h_inner(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/harmonic_h_inner_ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double num;
    double ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    double tol = pow(10, -16);

    unsigned int n;
    unsigned int i;

    unsigned int dim = 3;
    unsigned int theta1[dim];
    unsigned int theta2[dim];

    unsigned int *k = malloc(sizeof(unsigned int));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    unsigned int *beta = malloc(dim * sizeof(unsigned int));
    unsigned int *gAmma = malloc(dim * sizeof(unsigned int)); // NOLINT

    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: {alpha, alpha+1, alpha+2},
        //       {beta, beta+1, beta+2}
        //       {gamma, gamma+2, gamma+2}
        //       k, result
        scanResult = sscanf(line, "%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%lf", // NOLINT
                            alpha, alpha + 1, alpha + 2, beta, beta + 1, beta + 2,
                            gAmma, gAmma + 1, gAmma + 2, k, refRead); // NOLINT

        if (scanResult != 11) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 11\n", scanResult);
            continue;
        }

        // set n=|α|, i=|β|, θ₁=α+β−γ, θ₂=γ−β, θ₃=2γ−α−2β.
        n = mult_abs(dim, alpha);
        i = mult_abs(dim, beta);
        for (int j = 0; j < dim; j++) {
            theta1[j] = alpha[j] + beta[j] - gAmma[j]; // NOLINT
            theta2[j] = gAmma[j] - beta[j];            // NOLINT
        }

        ref = refRead[0];
        num = harmonic_h_inner_term(n, i, *k, dim, alpha, beta, theta1, theta2);

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("harmonic_h_inner ");
            printf(" %0*.16lf (this implementation) \n\t\t\t ≠ "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):           %E ≮ %E  (tolerance)\n",
                   errorMaxAbsRel, tol);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printMultiindexUnitTest("beta:\t\t", beta, dim);
            printMultiindexUnitTest("gamma:\t\t", gAmma, dim); // NOLINT
            printf("n: %u,  i: %u,  k: %u,  d: %u", n, i, *k, dim);
        }
        totalTests++;
    }

    free(k);
    free(alpha);
    free(beta);
    free(gAmma);
    free(refRead);

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
 * @brief Benchmarks harmonic polynomials.
 * @return number of failed tests.
 */
int test_harmonic_h_1D(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/harmonic_h_1D_ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double num;
    double ref;
    int scanResult;
    char line[256];
    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;
    int testsPassed = 0;
    int totalTests = 0;
    double tol = pow(10, -15);
    unsigned int alphaAbs;
    unsigned int dim = 1;
    unsigned int *k = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z}, {alpha}, result
        scanResult = sscanf(line, "%u,%lf, %u, %lf", // NOLINT
                            k, z, alpha, refRead);
        if (scanResult != 4) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 4\n", scanResult);
            continue;
        }
        alphaAbs = mult_abs(dim, alpha);
        ref = refRead[0];
        // Precompute coefficients for harmonic polynomials
        unsigned int kMax = alphaAbs / 2;
        unsigned long long *chunk_offset =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long *valid_count =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
            alphaAbs, kMax, dim, alpha, chunk_offset, valid_count);
        double *coeffs = malloc(coeffs_size * sizeof(double));
        unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
        precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_offset, coeffs,
                                        exponents);
        num = harmonic_h(*k, dim, z, alphaAbs, chunk_offset, valid_count, coeffs,
                         exponents);
        free(chunk_offset);
        free(valid_count);
        free(coeffs);
        free(exponents);
        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;
        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;
        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("harmonic_h at k: %u ", *k);
            printf(" %0*.16lf (this implementation) \n\t\t    ≠ "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E ≮ %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printVectorUnitTest("z:    \t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
        }
        totalTests++;
    }
    free(k);
    free(z);
    free(alpha);
    free(refRead);
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
 * @brief Benchmarks harmonic polynomials for arguments z = (1,1,1).
 * This focuses sorely on the precision of the harmonic coefficients.
 * @return number of failed tests.
 */
int test_harmonic_h_3D_unity(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/harmonic_h_3D_unity_ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double num;
    double ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    double tol = pow(10, -15);

    unsigned int alphaAbs;

    unsigned int dim = 3;

    unsigned int *k = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2, z3}, {alpha1, alpha2, alpha3}, result
        scanResult =
            sscanf(line, "%u,%lf, %lf, %lf, %u, %u, %u, %lf", // NOLINT
                   k, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        alphaAbs = mult_abs(dim, alpha);

        ref = refRead[0];

        // Precompute coefficients for harmonic polynomials
        unsigned int kMax = alphaAbs / 2;
        unsigned long long *chunk_offset =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long *valid_count =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
            alphaAbs, kMax, dim, alpha, chunk_offset, valid_count);
        double *coeffs = malloc(coeffs_size * sizeof(double));
        unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
        precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_offset, coeffs,
                                        exponents);
        num = harmonic_h(*k, dim, z, alphaAbs, chunk_offset, valid_count, coeffs,
                         exponents);
        free(chunk_offset);
        free(coeffs);
        free(valid_count);
        free(exponents);

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("harmonic_h at k: %u ", *k);
            printf(" %0*.16lf (this implementation) \n\t\t\t    ≠ "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E ≮ %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printVectorUnitTest("z:    \t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
        }
        totalTests++;
    }

    free(k);
    free(z);
    free(alpha);
    free(refRead);

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    if (totalTests == 0) {
        return fprintf(stderr, "Error: no valid test cases found in %s\n", path);
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
 * @brief Benchmarks harmonic polynomials for random three-dimensional arguments.
 * @return number of failed tests.
 */
int test_harmonic_h_3D_random(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/harmonic_h_3D_random_ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double num;
    double ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    double tol = 5 * pow(10, -13);

    unsigned int alphaAbs;

    unsigned int dim = 3;

    unsigned int *k = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2, z3}, {alpha1, alpha2, alpha3}, result
        scanResult =
            sscanf(line, "%u,%lf, %lf, %lf, %u, %u, %u, %lf", // NOLINT
                   k, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        alphaAbs = mult_abs(dim, alpha);

        ref = refRead[0];

        // Precompute coefficients for harmonic polynomials
        unsigned int kMax = alphaAbs / 2;
        unsigned long long *chunk_offset =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long *valid_count =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
            alphaAbs, kMax, dim, alpha, chunk_offset, valid_count);
        double *coeffs = malloc(coeffs_size * sizeof(double));
        unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
        precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_offset, coeffs,
                                        exponents);
        num = harmonic_h(*k, dim, z, alphaAbs, chunk_offset, valid_count, coeffs,
                         exponents);
        free(chunk_offset);
        free(coeffs);
        free(valid_count);
        free(exponents);

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("harmonic_h at k: %u ", *k);
            printf(" %0*.16lf (this implementation) \n\t\t\t    ≠ "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E ≮ %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printVectorUnitTest("z:    \t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
        }
        totalTests++;
    }

    free(k);
    free(z);
    free(alpha);
    free(refRead);

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    if (totalTests == 0) {
        return fprintf(stderr, "Error: no valid test cases found in %s\n", path);
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
 * @brief Main function to run all harmonic polynomials tests and related function
 * tests.
 *
 * @return number of failed tests.
 */
int main(void) {
    int failed = 0;
    failed += run_timed_test(test_coeffs_c_inner);
    failed += run_timed_test(test_harmonic_h_inner);
    failed += run_timed_test(test_harmonic_h_1D);
    failed += run_timed_test(test_harmonic_h_3D_unity);
    failed += run_timed_test(test_harmonic_h_3D_random);
    return failed;
}
