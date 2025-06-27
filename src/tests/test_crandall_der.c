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

/*!
 * @brief Benchmarks 3D upper Crandall function by comparing to high-precision values
 * over a range of random parameters.
 *
 * @return 0 if all tests pass, 1 if any test fails.
 */
int test_crandall_g_der(void) {
    printf("%s ... \n", __func__);
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
    double *z = malloc(2 * sizeof(double));
    unsigned int *alpha = malloc(2 * sizeof(unsigned int));
    double *refRead = malloc(2 * sizeof(double));

    printf("\tProcessing file: %s ... ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {z1, z2}, {alpha1, alpha2}, {Re[result], Im[result]}
        scanResult = sscanf(line, "%lf,%lf,%lf,%lf,%u,%u,%u,%lf,%lf", // NOLINT
                            nuRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2,
                            refRead, refRead + 1);

        if (scanResult != 9) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 7\n", scanResult);
            continue;
        }

        nu = nuRef[0];

        zArgBound = assignzArgBound(nu);
        // printf("alpha: (%d, %d)",alpha[0],alpha[1]);
        // printf("alpha/2: (%d, %d)", alpha[0]/2, alpha[1]/2);

        num = crandall_g_der(dim, nu, z, prefactor, zArgBound, alpha);
        ref = refRead[0] + refRead[1] * I;

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        totalTests++;
        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\nWarning! ");
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
        }
    }
    printf("%d out of %d tests passed.\n", testsPassed, totalTests);

    free(nuRef);
    free(z);
    free(alpha);
    free(refRead);

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    return (testsPassed == totalTests) ? 0 : 1;
}

/*!
 * @brief Benchmarks 2D upper Crandall function by computing its taylor series.
 *
 * @return 0 if all tests pass, 1 if any test fails.
 */
int test_crandall_g_der_taylor(void) {
    printf("%s ... ", __func__);
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
        while (1) {

            zArgBound = assignzArgBound(nu);

            valTaylor += mult_pow(dim, alpha, zDiff, 1.) /
                         (double)mult_fac(dim, alpha) *
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

        totalTests++;
        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\nWarning! ");
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
            printf("\n");
        }
    }

    printf("%d out of %d tests passed.\n", testsPassed, totalTests);

    return totalTests - testsPassed;
}

int main(void) {
    bool result1 = test_crandall_g_der();
    bool result2 = test_crandall_g_der_taylor();
    return result1 + result2;
}
