// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "epsteinZeta.h"
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
 * @brief Free memory allocated for test resources.
 *
 * This function deallocates memory that was dynamically allocated for various
 * arrays used in the Epstein Zeta function tests.
 *
 * @param[in] a Pointer to the matrix of coefficients.
 * @param[in] nu Pointer to the array of nu values.
 * @param[in] x Pointer to the array of x values.
 * @param[in] y Pointer to the array of y values.
 * @param[in] zetaRef Pointer to the array of reference zeta values.
 */
void freeTestResources(double *a, double *nu, double *x, double *y,
                       double *zetaRef) {
    free(a);
    free(nu);
    free(x);
    free(y);
    free(zetaRef);
}

/*!
 * @brief Test function for Epstein Zeta and Regularized Epstein Zeta functions.
 *
 * This function tests both the epsteinZeta and epsteinZetaReg functions by comparing
 * their outputs with reference values read from data files. It performs multiple
 * test cases for each function and reports the number of passed tests.
 *
 * @return 0 if all tests pass, 1 if any test fails.
 */
int test_epsteinZeta_epsteinZetaReg() {
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/epsteinZeta_Ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *zetaRegRefData = fopen(path, "r");
    if (zetaRegRefData == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    int dim = 2;
    double *a = malloc((int)(dim * dim) * sizeof(double));
    double *nu = malloc(2 * sizeof(double));
    double *x = malloc(dim * sizeof(double));
    double *y = malloc(dim * sizeof(double));
    double *zetaRef = malloc(2 * sizeof(double));
    double complex zetaC;
    double complex zetaM;
    double nuReal;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double tol = pow(10, -13);
    int testsPassed = 0;
    int totalTests = 0;
    int testsPassedOverall = 0;
    int totalTestsOverall = 0;
    printf("Processing file: %s ... ", path);

    int scanResult;
    char line[256];
    while (fgets(line, sizeof(line), zetaRegRefData) != NULL) {
        scanResult =
            sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", // NOLINT
                   nu, nu + 1, a, a + 1, a + 2, a + 3, x, x + 1, y, y + 1, zetaRef,
                   zetaRef + 1);

        if (scanResult != 12) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 12\n", scanResult);
            continue;
        }

        nuReal = nu[0];
        zetaC = epsteinZeta(nuReal, dim, a, x, y);
        zetaM = zetaRef[0] + zetaRef[1] * I;
        errorAbs = errAbs(zetaC, zetaM);
        errorRel = errRel(zetaC, zetaM);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        totalTests++;
        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\nWarning! ");
            { printf("zeta:   "); }
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(zetaC), cimag(zetaC), creal(zetaM), cimag(zetaM));
            printf("Min(Emax, Erel):      %.16lf > %.16lf  (tolerance)\n",
                   errorMaxAbsRel, tol);
            printf("\n");
            printMatrixUnitTest("a:", a, dim);
            printf("nu:\t\t %.16lf + %.16lf I\n", nu[0], nu[1]);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
            printf("\n");
        }
    }
    printf("%d out of %d tests passed.\n", testsPassed, totalTests);

    testsPassedOverall += testsPassed;
    totalTestsOverall += totalTests;

    totalTests = 0;
    testsPassed = 0;

    if (fclose(zetaRegRefData) != 0) {
        free(a);
        free(nu);
        free(x);
        free(y);
        free(zetaRef);
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    result = snprintf(path, sizeof(path), "%s/epsteinZetaReg_Ref.csv", // NOLINT
                      BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        free(a);
        free(nu);
        free(x);
        free(y);
        free(zetaRef);
        return fprintf(stderr, "Error creating file path\n");
    }
    zetaRegRefData = fopen(path, "r");
    if (zetaRegRefData == NULL) {
        free(a);
        free(nu);
        free(x);
        free(y);
        free(zetaRef);
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    printf("Processing file: %s ... ", path);
    while (fgets(line, sizeof(line), zetaRegRefData) != NULL) {
        scanResult =
            sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", // NOLINT
                   nu, nu + 1, a, a + 1, a + 2, a + 3, x, x + 1, y, y + 1, zetaRef,
                   zetaRef + 1);

        if (scanResult != 12) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 12\n", scanResult);
            continue;
        }

        nuReal = nu[0];
        zetaC = epsteinZetaReg(nuReal, dim, a, x, y);
        zetaM = zetaRef[0] + zetaRef[1] * I;
        errorAbs = errAbs(zetaC, zetaM);
        errorRel = errRel(zetaC, zetaM);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        totalTests++;
        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\nWarning! ");
            { printf("zeta reg:"); }
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(zetaC), cimag(zetaC), creal(zetaM), cimag(zetaM));
            printf("Min(Emax, Erel):      %.16lf > %.16lf  (tolerance)\n",
                   errorMaxAbsRel, tol);
            printf("\n");
            printMatrixUnitTest("a:", a, dim);
            printf("nu:\t\t %.16lf + %.16lf I\n", nu[0], nu[1]);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
            printf("\n");
        }
    }

    if (fclose(zetaRegRefData) != 0) {
        freeTestResources(a, nu, x, y, zetaRef);
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }

    freeTestResources(a, nu, x, y, zetaRef);

    printf("%d out of %d tests passed.\n", testsPassed, totalTests);

    testsPassedOverall += testsPassed;
    totalTestsOverall += totalTests;

    return (testsPassedOverall == totalTestsOverall) ? 0 : 1;
}

int main() {
    int result = test_epsteinZeta_epsteinZetaReg();
    return result;
}
