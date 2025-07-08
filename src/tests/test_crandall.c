// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file test_crandall.c
 * @brief Benchmarking of the upper Crandall function.
 */

#include "../crandall.h"
#include "utils.h"
#include <complex.h>
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
 * @brief Test function for crandall_g.
 *
 * @return number of failed tests.
 */
int test_crandall_g(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/crandall_g_Ref.csv", // NOLINT
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
    unsigned int dim = 2;
    double prefactor = 1.;
    double tol = pow(10, -13);

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double *nuRef = malloc(sizeof(double));
    double *z = malloc(2 * sizeof(double));
    double *refRead = malloc(2 * sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf(line, "%lf,%lf,%lf,%lf,%lf", // NOLINT
                            nuRef, z, z + 1, refRead, refRead + 1);

        if (scanResult != 5) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 5\n", scanResult);
            continue;
        }

        nu = nuRef[0];

        zArgBound = assignzArgBound(nu);

        num = crandall_g(dim, nu, z, prefactor, zArgBound);
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
            printf("\n");
            printf("Warning! ");
            printf("crandall_g: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("z:\t\t", z, dim);
        }
        totalTests++;
    }
    free(nuRef);
    free(z);
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
 * @brief Main function to run all Crandall function tests.
 *
 * @return number of failed tests.
 */
int main(void) {
    int failed = test_crandall_g();
    return failed;
}
