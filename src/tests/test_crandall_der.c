// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file test_crandall_der.c
 * @brief Compares the derivatives of 3D Crandall functions to high-precision
 * benchmark values.
 */

#include "../crandall.h"
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
        // Scan: nu, {z1, z2}, {alpha1, alpha2}, {Re[result], Im[result]}
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

        if (errorMaxAbsRel > tol) {
            printf("\n\n");
            printf("Warning! ");
            printf("polynomial_p");
            printf(" %0*.16lf (this implementation) \n\t\t!= "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Emax, Erel):      %E > %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printVectorUnitTest("y:\t\t", y, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printMultiindexUnitTest("beta:\t\t", beta, dim);
        } else {
            testsPassed++;
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
    printf("%d out of %d tests passed with tolerance %E.\n", testsPassed, totalTests,
           tol);

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks 2D upper Crandall function by computing its taylor series.
 *
 * @return number of failed tests.
 */
int test_crandall_g_der_taylor(void) {
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

    bool done;
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));

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

        // Iterate over every multi-index alpha so that every alpha[] < order
        while (true) {

            zArgBound = assignzArgBound(nu);

            valTaylor += mult_pow(dim, alpha, zDiff) / (double)mult_fac(dim, alpha) *
                         crandall_g_der(dim, nu, z, 1., zArgBound, alpha);

            done = 1;
            for (unsigned int idx = 0; idx < dim; idx++) {
                if (alpha[idx] + 1 <= order) {
                    alpha[idx]++;
                    done = 0;
                    break;
                }
                alpha[idx] = 0;
            }
            if (done) {
                break;
            }
        }

        errorAbs = errAbs(valRef, valTaylor);
        errorRel = errRel(valRef, valTaylor);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        if (errorMaxAbsRel > tol) {
            printf("\n\n");
            printf("Warning! ");
            printf("crandall_g: ");
            printf(" %0*.16lf %+.16lf I (as a taylor series) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(valTaylor), cimag(valTaylor), creal(valRef),
                   cimag(valRef));
            printf("Min(Emax, Erel):      %E > %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z0:\t\t", z, dim);
            printVectorUnitTest("zPlus:\t\t", zPlus, dim);
            printVectorUnitTest("zDiff:\t\t", zDiff, dim);
        } else {
            testsPassed++;
        }
        totalTests++;
    }

    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.\n", testsPassed, totalTests,
           tol);

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

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double prefactor = 1.;
    double tol = 5 * pow(10, -13);

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

        num = crandall_g_der(dim, nu, z, prefactor, zArgBound, alpha);
        ref = refRead[0] + refRead[1] * I;

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        if (errorMaxAbsRel > tol) {
            printf("\n\n");
            printf("Warning! ");
            printf("crandall_g_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E > %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        } else {
            testsPassed++;
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
    printf("%d out of %d tests passed with tolerance %E.\n", testsPassed, totalTests,
           tol);

    return totalTests - testsPassed;
}

int main(void) {
    int result1 = test_polynomial_p();
    printf("\n");
    int result2 = test_crandall_g_der();
    printf("\n");
    int result3 = test_crandall_g_der_taylor();
    return result1 + result2 + result3;
}
