// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file test_crandall_der.c
 * @brief Compares the derivatives of 3D Crandall functions to high-precision
 * benchmark values.
 */

#include "../crandall.h"
#include "../tools.h"
#include "utils.h"
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

#ifndef BASE_PATH
#define BASE_PATH "csv"
#endif

/** @brief Computes a single summand of h_inner; explicitly
 * (−1/2)^{i} / cn,i,k binom(i+k,k) binom(i,β) binom(α,θ1) θ2! / (θ2 - θ1)!,
 * where n=|α|, i=|β|, θ1=α+β−γ, θ2=γ−β.
 * @param[in] n: index (total of alpha).
 * @param[in] i: index (total of beta).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension of the multi-indices.
 * @param[in] beta: lower multi-index β.
 * @param[in] alpha: upper multi-index α.
 * @param[in] theta1: multi-index α+β−γ.
 * @param[in] theta2: multi-index γ−β.
 * @return value of one summand in of h_inner.
 */
static double harmonic_h_inner_term(unsigned int n, unsigned int i, unsigned int k,
                                    unsigned int dim, const unsigned int *alpha,
                                    const unsigned int *beta,
                                    const unsigned int *theta1,
                                    const unsigned int *theta2) {
    apint_t numerator;

    harmonic_h_inner_term_scalar_apint(n, i, k, dim, &numerator);

    apint_t multi_term;
    harmonic_h_inner_term_multi_apint(dim, alpha, beta, theta1, theta2, &multi_term);

    apint_mul(&numerator, &multi_term, &numerator);

    double denominator = 1.0;

    for (unsigned int l = 1; l <= n / 2 - k; l++) {
        denominator *= coeffs_c_inner(n, l, k, dim);
    }

    double res = apint_to_double(&numerator) / denominator;

    return res;
}

/*!
 * @brief Benchmarks 3D polynomial_p function by comparing to high-precision values
 * over a range of random parameters.
 *
 * @return number of failed tests.
 */
int test_polynomial_p(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/polynomial_p_Ref.csv", // NOLINT
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
    int dim = 3;
    double tol = 5 * pow(10, -15);

    double *y = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    unsigned int *beta = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {z1, z2, z3}, {alpha1, alpha2, alpha3}, {Re[result], Im[result]}
        scanResult = sscanf(line, "%lf,%lf,%lf,%u,%u,%u,%u,%u,%u,%lf", // NOLINT
                            y, y + 1, y + 2, alpha, alpha + 1, alpha + 2, beta,
                            beta + 1, beta + 2, refRead);

        if (scanResult != 10) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 10\n", scanResult);
            continue;
        }

        ref = refRead[0];
        num = polynomial_p(dim, y, alpha, beta);

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
            printf("polynomial_p");
            printf(" %0*.16lf (this implementation) \n\t\t!= "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printVectorUnitTest("y:\t\t", y, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printMultiindexUnitTest("beta:\t\t", beta, dim);
        }
        totalTests++;
    }

    free(y);
    free(alpha);
    free(beta);
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
 * @brief Benchmarks 3D polynomial_l function by comparing to high-precision values
 * over a range of random parameters.
 *
 * @return number of failed tests.
 */
int test_polynomial_l(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/polynomial_l_Ref.csv", // NOLINT
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
    int dim = 3;
    double tol = 5 * pow(10, -15);

    double *y = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    unsigned int *beta = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {z1, z2, z3}, {alpha1, alpha2, alpha3}, {Re[result], Im[result]}
        scanResult = sscanf(line, "%lf,%lf,%lf,%u,%u,%u,%u,%u,%u,%lf", // NOLINT
                            y, y + 1, y + 2, alpha, alpha + 1, alpha + 2, beta,
                            beta + 1, beta + 2, refRead);

        if (scanResult != 10) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 10\n", scanResult);
            continue;
        }

        ref = refRead[0];
        num = polynomial_l(dim, y, alpha, beta);

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
            printf("polynomial_l");
            printf(" %0*.16lf (this implementation) \n\t\t!= "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printVectorUnitTest("y:\t\t", y, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printMultiindexUnitTest("beta:\t\t", beta, dim);
        }
        totalTests++;
    }

    free(y);
    free(alpha);
    free(beta);
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
 * @brief Benchmarks 2D upper Crandall function by computing its taylor series.
 *
 * @return number of failed tests.
 */
int test_crandall_g_der_taylor(void) { // NOLINT
    printf("%s ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double zArgBound;
    double complex valRef;
    double complex valTaylor;
    double tol = pow(10, -14);
    unsigned int dim = 2;
    unsigned int order = 12;
    double zDiff[] = {0.005, 0.01};

    double nu;
    double *z = malloc(dim * sizeof(double));
    double *zPlus = malloc(dim * sizeof(double));

    int testsPassed = 0;
    int totalTests = 0;

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    bool done;
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));

    printf("\n\t ... ");
    printf("generating test values");
    for (int i = 0; i < 100; i++) {

        z[0] = (double)i / 110. + 0.1;
        z[1] = (double)i / 210.;
        nu = -0.500000001 + i / 9.;

        for (int i = 0; i < dim; i++) {
            zPlus[i] = z[i] + zDiff[i];
        }

        zArgBound = assignzArgBound(nu);

        valRef = crandall_g(dim, nu, zPlus, 1, zArgBound);

        // build taylor series around z
        valTaylor = 0;

        // Initialize multi-index
        for (int i = 0; i < dim; i++) {
            alpha[i] = 0;
        }

        // initialize alphaAbs
        unsigned int alphaAbs = 0;

        // Iterate over every multi-index alpha so that every alpha[] < order
        while (true) {

            zArgBound = assignzArgBound(nu);

            valTaylor += mult_pow(dim, alpha, zDiff) / (double)mult_fac(dim, alpha) *
                         crandall_g_der(dim, nu, z, 1., zArgBound, alpha, alphaAbs);

            done = 1;
            for (unsigned int idx = 0; idx < dim; idx++) {
                if (alpha[idx] + 1 <= order) {
                    alpha[idx]++;
                    alphaAbs++;
                    done = 0;
                    break;
                }
                alphaAbs -= alpha[idx];
                alpha[idx] = 0;
            }
            if (done) {
                break;
            }
        }

        errorAbs = errAbs(valRef, valTaylor);
        errorRel = errRel(valRef, valTaylor);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("crandall_g_der: ");
            printf(" %0*.16lf %+.16lf I (as a taylor series) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(valTaylor), cimag(valTaylor), creal(valRef),
                   cimag(valRef));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z0:\t\t", z, dim);
            printVectorUnitTest("zPlus:\t\t", zPlus, dim);
            printVectorUnitTest("zDiff:\t\t", zDiff, dim);
        }
        totalTests++;
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
 * @brief Benchmarks 2D regularized Crandall derivatives by computing its taylor
 * series.
 *
 * @return number of failed tests.
 */
int test_crandall_gReg_der_taylor(void) { // NOLINT
    printf("%s ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex valRef;
    double complex valTaylor;
    double tol = pow(10, -14);
    unsigned int dim = 2;
    unsigned int order = 12;
    double zDiff[] = {0.005, 0.01};

    double nu;
    double *z = malloc(dim * sizeof(double));
    double *zPlus = malloc(dim * sizeof(double));

    int testsPassed = 0;
    int totalTests = 0;

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    bool done;
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));

    printf("\n\t ... ");
    printf("generating test values");
    for (int i = 0; i < 100; i++) {

        z[0] = (double)i / 110. + 0.1;
        z[1] = (double)i / 210.;
        nu = -0.500000001 + i / 9.;

        for (int i = 0; i < dim; i++) {
            zPlus[i] = z[i] + zDiff[i];
        }

        valRef = crandall_gReg(dim, nu, zPlus, 1);

        // build taylor series around z
        valTaylor = 0;

        // Initialize multi-index
        for (int i = 0; i < dim; i++) {
            alpha[i] = 0;
        }

        // initialize alphaAbs
        unsigned int alphaAbs = 0;

        double zArgBound = assignzArgBound(nu);

        // Iterate over every multi-index alpha so that every alpha[] < order
        while (true) {

            valTaylor +=
                mult_pow(dim, alpha, zDiff) / (double)mult_fac(dim, alpha) *
                crandall_gReg_der(dim, nu, z, 1., alpha, alphaAbs, zArgBound);

            done = 1;
            for (unsigned int idx = 0; idx < dim; idx++) {
                if (alpha[idx] + 1 <= order) {
                    alpha[idx]++;
                    alphaAbs++;
                    done = 0;
                    break;
                }
                alphaAbs -= alpha[idx];
                alpha[idx] = 0;
            }
            if (done) {
                break;
            }
        }

        errorAbs = errAbs(valRef, valTaylor);
        errorRel = errRel(valRef, valTaylor);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("crandall_gReg_der: ");
            printf(" %0*.16lf %+.16lf I (as a taylor series) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(valTaylor), cimag(valTaylor), creal(valRef),
                   cimag(valRef));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z0:\t\t", z, dim);
            printVectorUnitTest("zPlus:\t\t", zPlus, dim);
            printVectorUnitTest("zDiff:\t\t", zDiff, dim);
        }
        totalTests++;
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
 * @brief Benchmarks 3D polynomial Y derivatives by comparing to high-precision
 * values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_polynomial_y_der(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/polynomial_y_der_Ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    unsigned int k;
    unsigned int alphaAbs;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double tol = 5 * pow(10, -13);

    unsigned int *kRef = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2, z3}, {alpha1, alpha2i, alpha3}, result
        scanResult =
            sscanf(line, "%u,%lf,%lf,%lf,%u,%u,%u,%lf", // NOLINT
                   kRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        k = kRef[0];

        alphaAbs = mult_abs(dim, alpha);

        num = polynomial_y_der(k, dim, z, alpha, alphaAbs, 1);
        ref = refRead[0] + 0 * I;

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
            printf("polynomial_y_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("k:\t\t %u\n", k);
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(kRef);
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
 * @brief Benchmarks 3D log(pi * y**2) derivatives by comparing to high-precision
 * values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_log_l_der(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/log_l_der_Ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }

    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    unsigned int alphaAbs;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double tol = 5 * pow(10, -12);

    unsigned int *kRef = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2, z3}, {alpha1, alpha2i, alpha3}, result
        scanResult = sscanf(line, "%lf,%lf,%lf,%u,%u,%u,%lf", // NOLINT
                            z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 7) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 7\n", scanResult);
            continue;
        }

        alphaAbs = mult_abs(dim, alpha);

        num = log_l_der(dim, z, alpha, alphaAbs);
        ref = refRead[0] + 0 * I;

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
            printf("log_l_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(kRef);
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
 * @brief Benchmarks recursive coefficients
 * @return number of failed tests.
 */
int test_coeffs_c_inner(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/coeffs_c_inner_Ref.csv", // NOLINT
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
    double tol = 5 * pow(10, -15);

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
        snprintf(path, sizeof(path), "%s/harmonic_h_inner_Ref.csv", // NOLINT
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
    double tol = 5 * pow(10, -15);

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
            printf(" %0*.16lf (this implementation) \n\t\t    ≠ "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E ≮ %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
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
    int result = snprintf(path, sizeof(path), "%s/harmonic_h_1D_Ref.csv", // NOLINT
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
    double tol = 5 * pow(10, -11);

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
        unsigned long long *chunk_size =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long coeffs_size =
            precompute_harmonic_h_inner_chunk_size(alphaAbs, kMax, dim, chunk_size);
        double *coeffs = (double *)malloc(coeffs_size * sizeof(double));
        precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_size, coeffs);

        num = harmonic_h(*k, dim, z, alpha, alphaAbs, chunk_size, coeffs);

        free(coeffs);
        free(chunk_size);

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
 * @brief Benchmarks harmonic polynomials.
 * @return number of failed tests.
 */
int test_harmonic_h_3D(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/harmonic_h_3D_Ref.csv", // NOLINT
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
    double tol = 5 * pow(10, -10);

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
        unsigned long long *chunk_size =
            malloc((kMax + 1) * sizeof(unsigned long long));
        unsigned long long coeffs_size =
            precompute_harmonic_h_inner_chunk_size(alphaAbs, kMax, dim, chunk_size);
        double *coeffs = (double *)malloc(coeffs_size * sizeof(double));
        precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_size, coeffs);

        num = harmonic_h(*k, dim, z, alpha, alphaAbs, chunk_size, coeffs);

        free(coeffs);
        free(chunk_size);

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
 * @brief Benchmarks 3D singularity s derivatives by comparing to high-precision
 * values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_singularity_s_der(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/singularity_s_der_Ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }

    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    unsigned int k;
    unsigned int alphaAbs;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double tol = 5 * pow(10, -12);

    unsigned int *kRef = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2, z3}, {alpha1, alpha2, alpha3}, result
        scanResult =
            sscanf(line, "%u,%lf,%lf,%lf,%u,%u,%u,%lf", // NOLINT
                   kRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        alphaAbs = mult_abs(dim, alpha);

        k = kRef[0];

        num = singularity_s_der(k, dim, z, alpha, alphaAbs);
        ref = refRead[0] + 0 * I;

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
            printf("singularity_s_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t\t  != "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("k:\t\t %u\t\t", k);
            printf("\n");
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(kRef);
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
 * @brief Benchmarks 3D upper Crandall function by comparing to high-precision values
 * over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_crandall_g_der(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/crandall_g_der_Ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double nu;
    double zArgBound;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double prefactor = 1.;
    double tol = 5 * pow(10, -12);

    double *nuRef = malloc(sizeof(double));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(dim * sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {z1, z2}, {alpha1, alpha2}, {Re[result], Im[result]}
        scanResult = sscanf(line, "%lf,%lf,%lf,%lf,%u,%u,%u,%lf,%lf", // NOLINT
                            nuRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2,
                            refRead, refRead + 1);

        if (scanResult != 9) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 9\n", scanResult);
            continue;
        }

        nu = nuRef[0];

        zArgBound = assignzArgBound(nu);

        num = crandall_g_der(dim, nu, z, prefactor, zArgBound, alpha,
                             mult_abs(dim, alpha));
        ref = refRead[0] + refRead[1] * I;

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
            printf("crandall_g_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(nuRef);
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
 * @brief Benchmarks 3D regularized Crandall derivatives by comparing to
 * high-precision values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_crandall_gReg_der(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/crandall_gReg_der_Ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double nu;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double prefactor = 1.;
    double tol = pow(10, -11);

    double *nuRef = malloc(sizeof(double));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(dim * sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {z1, z2}, {alpha1, alpha2}, {Re[result], Im[result]}
        scanResult = sscanf(line, "%lf,%lf,%lf,%lf,%u,%u,%u,%lf,%lf", // NOLINT
                            nuRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2,
                            refRead, refRead + 1);

        if (scanResult != 9) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 9\n", scanResult);
            continue;
        }

        nu = nuRef[0];

        num = crandall_gReg_der(dim, nu, z, prefactor, alpha, mult_abs(dim, alpha),
                                assignzArgBound(nu));
        ref = refRead[0] + refRead[1] * I;

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
            printf("crandall_gReg_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(nuRef);
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
 * @brief Benchmarks 3D regularized Crandall derivatives for nu = dim + 2 * k by
 * comparing to high-precision values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_crandall_gReg_der_d2k_prototype(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/crandall_gReg_der_d2k_prototype_Ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double nu;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 2;
    double prefactor = 1.;
    double tol = 5 * pow(10, -10);

    double *nuRef = malloc(sizeof(double));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(2 * sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2}, {alpha1, alpha2}, Re[result]
        scanResult = sscanf(line, "%lf,%lf,%lf,%u,%u,%lf", // NOLINT
                            nuRef, z, z + 1, alpha, alpha + 1, refRead);

        if (scanResult != 6) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 6\n", scanResult);
            continue;
        }

        nu = nuRef[0];

        num = crandall_gReg_der(dim, nu, z, prefactor, alpha, mult_abs(dim, alpha),
                                assignzArgBound(nu));
        ref = refRead[0] + 0. * I;

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
            printf("%s: ", __func__);
            printf(" %0*.16lf (this implementation) \n\t\t\t\t\t\t\t != "
                   "%.16lf (reference implementation)\n",
                   4, creal(num), creal(ref));
            printf("Min(Emax, Erel): %E !< %E  (tolerance)\n", errorMaxAbsRel, tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(nuRef);
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

int main(void) {

    printf("start ");
    int failed = 0;
    failed += test_polynomial_p();
    failed += test_polynomial_l();
    failed += test_polynomial_y_der();
    failed += test_coeffs_c_inner();
    failed += test_harmonic_h_inner();
    failed += test_harmonic_h_1D();
    failed += test_harmonic_h_3D();
    failed += test_log_l_der();
    failed += test_singularity_s_der();
    failed += test_crandall_g_der();
    failed += test_crandall_gReg_der();
    failed += test_crandall_gReg_der_d2k_prototype();
    failed += test_crandall_g_der_taylor();
    failed += test_crandall_gReg_der_taylor();
    return failed;
}
