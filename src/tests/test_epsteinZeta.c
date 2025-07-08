// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file test_epsteinZeta.c
 * @brief Benchmarking of the Epstein zeta function.
 */

#include "../tools.h"
#include "epsteinZeta.h"
#include "stdbool.h"
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
 * @brief Calculates the singularity of the Epstein zeta function as y approaches
 * zero.
 *
 * @param[in] nu The parameter nu of the Epstein zeta function.
 * @param[in] dim The dimension of the lattice.
 * @param[in] y Pointer to the array of y values.
 * @return The value of the singularity function.
 */
double sHat(double nu, unsigned int dim, double *y) {
    double ySquared = dot(dim, y, y);
    double k = fmax(0., nearbyint((nu - (double)dim) / 2));
    if (nu == (dim + 2 * k)) {
        return pow(M_PI, (2 * k) + (((double)dim) / 2)) /
               tgamma(k + (((double)dim) / 2)) * pow(-1, k + 1) / tgamma(k + 1) *
               pow(ySquared, k) * log(M_PI * ySquared);
    }
    return pow(M_PI, nu - (((double)dim) / 2)) *
           pow(ySquared, (nu - ((double)dim)) / 2) * tgamma(((double)dim - nu) / 2) /
           tgamma(nu / 2);
}

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
 * @brief Reports an error when the Epstein Zeta function representation test fails.
 *
 * This function prints detailed information about the test case that failed,
 * including the computed values, error margins, and input parameters.
 *
 * @param[in] valZeta The value computed by epsteinZeta.
 * @param[in] valZetaReg The value computed by the epsteinZetaReg representation.
 * @param[in] errorMaxAbsRel The maximum of absolute and relative errors.
 * @param[in] tol The tolerance level for the test.
 * @param[in] m Pointer to the matrix of coefficients.
 * @param[in] dim The dimension of the lattice.
 * @param[in] nu The parameter nu of the Epstein zeta function.
 * @param[in] x Pointer to the array of x values.
 * @param[in] y Pointer to the array of y values.
 */
void reportEpsteinZetaError(double complex valZeta, double complex valZetaReg,
                            double errorMaxAbsRel, double tol, double *m,
                            unsigned int dim, double nu, double *x, double *y) {
    printf("\n");
    printf("Warning! ");
    printf("epsteinZeta:");
    printf(" %0*.16lf %+.16lf I (epsteinZeta) \n\t\t  != "
           "%.16lf "
           "%+.16lf I (epsteinZetaReg representation)\n",
           4, creal(valZeta), cimag(valZeta), creal(valZetaReg), cimag(valZetaReg));
    printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel, tol);
    printf("\n");
    printMatrixUnitTest("m:", m, (int)dim);
    printf("nu:\t\t %.16lf\n", nu);
    printVectorUnitTest("x:\t\t", x, (int)dim);
    printVectorUnitTest("y:\t\t", y, (int)dim);
}

/*!
 * @brief Reports an error when the Epstein Zeta function cutoff test fails.
 *
 * @param[in] testCase A string describing the specific test case that failed.
 * @param[in] zeta1 The first computed value of the Epstein Zeta function.
 * @param[in] zeta2 The second computed value of the Epstein Zeta function for
 * comparison.
 * @param[in] nu The parameter nu of the Epstein zeta function.
 * @param[in] y Pointer to the array of y values.
 * @param[in] dim The dimension of the lattice.
 */
void reportEpsteinZetaCutoffError(const char *testCase, double complex zeta1,
                                  double complex zeta2, double nu, double *y,
                                  unsigned int dim) {
    printf("\n\n");
    printf("Warning! ");
    printf("%s:\n", testCase);
    printf(" %0*.16lf %+.16lf I \n\t\t  != "
           "%.16lf %+.16lf I\n",
           4, creal(zeta1), cimag(zeta1), creal(zeta2), cimag(zeta2));
    printf("nu:\t\t %.16lf\n", nu);
    printf("y:\t\t");
    for (unsigned int i = 0; i < dim; i++) {
        printf("%.32lf", y[i]);
        if (i < dim - 1) {
            printf(", ");
        }
    }
    printf("\n");
}

/*!
 * @brief Test function for Epstein Zeta and Regularized Epstein Zeta functions.
 *
 * This function tests both the epsteinZeta and epsteinZetaReg functions by comparing
 * their outputs with reference values read from data files. It performs multiple
 * test cases for each function and reports the number of passed tests.
 *
 * @return number of failed tests.
 */
int test_epsteinZeta_epsteinZetaReg() { // NOLINT
    printf("%s ", __func__);
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
    double *a = malloc((int)(dim * dim) * sizeof(double)); // NOLINT
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

    double errMin = 0.;
    double errMax = 0.;
    double errSum = 0.;

    int maxIt = 0;

    printf("\n\t ... ");
    printf("processing %s ", path);

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

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n");
            printf("Warning! ");
            printf("zeta:   ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(zetaC), cimag(zetaC), creal(zetaM), cimag(zetaM));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printMatrixUnitTest("a:", a, dim);
            printf("nu:\t\t %.16lf + %.16lf I\n", nu[0], nu[1]);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
        }
        totalTests++;
        maxIt++;
    }
    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / maxIt);

    errMin = 0.;
    errMax = 0.;
    errSum = 0.;

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

    printf("\n\t ... ");
    printf("processing file: %s ", path);

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

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n");
            printf("Warning! ");
            printf("zeta reg:");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(zetaC), cimag(zetaC), creal(zetaM), cimag(zetaM));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printMatrixUnitTest("a:", a, dim);
            printf("nu:\t\t %.16lf + %.16lf I\n", nu[0], nu[1]);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
        }
        totalTests++;
        maxIt++;
    }

    if (fclose(zetaRegRefData) != 0) {
        freeTestResources(a, nu, x, y, zetaRef);
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }

    freeTestResources(a, nu, x, y, zetaRef);

    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / maxIt);
    printf("\n");

    testsPassedOverall += testsPassed;
    totalTestsOverall += totalTests;

    return (testsPassedOverall == totalTestsOverall) ? 0 : 1;
}

/*!
 * @brief Test if the Epstein Zeta function can be represented in terms of the
 * regularized function and the singularity, particularly focusing on cases where nu
 * = dim + 2k, where k is an integer.
 *
 * @return number of failed tests.
 * tests pass.
 */
bool test_epsteinZeta_epsteinZetaReg_represent_as_each_other() {
    printf("%s ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex valZeta;
    double complex valZetaReg;
    double nu;

    int testsPassed = 0;
    int totalTests = 0;

    double errMin = 0.;
    double errMax = 0.;
    double errSum = 0.;

    int maxIt = 200;

    double tol = pow(10, -14);
    unsigned int dim = 2;
    double m[] = {3. / 2, 1. / 5, 1. / 4, 1.};
    double x[] = {0.1, 0.2};
    double y[] = {0, 0.5};
    double vol = 29. / 20;

    for (int i = 0; i < maxIt / 2; i++) {
        nu = -8.5 + (double)i / 5.;

        valZeta = epsteinZeta(nu, dim, m, x, y);
        valZetaReg = cexp(-2 * M_PI * I * dot(dim, x, y)) *
                     (epsteinZetaReg(nu, dim, m, x, y) + sHat(nu, dim, y) / vol);

        errorAbs = errAbs(valZeta, valZetaReg);
        errorRel = errRel(valZeta, valZetaReg);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            reportEpsteinZetaError(valZeta, valZetaReg, errorMaxAbsRel, tol, m, dim,
                                   nu, x, y);
        }
        totalTests++;
    }

    // test identity around zero
    double yZeta[] = {0, pow(10, -16)};
    double yZetaReg[] = {0., 0.};

    for (int i = 0; i < maxIt / 2; i++) {
        nu = -8.5 + (double)i / 5.;

        valZeta = epsteinZeta(nu, dim, m, x, yZeta);
        valZetaReg =
            cexp(-2 * M_PI * I * dot(dim, x, yZeta)) *
            (epsteinZetaReg(nu, dim, m, x, yZetaReg) + sHat(nu, dim, yZeta) / vol);

        errorAbs = errAbs(valZeta, valZetaReg);
        errorRel = errRel(valZeta, valZetaReg);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            reportEpsteinZetaError(valZeta, valZetaReg, errorMaxAbsRel, tol, m, dim,
                                   nu, x, yZetaReg);
        }
        totalTests++;
    }

    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / maxIt);
    printf("\n");

    return totalTests - testsPassed;
}

/*!
 * @brief Test the Epstein Zeta function behavior around the cutoff point for nu <=
 * dim.
 *
 * This function tests the Epstein Zeta function for four cases:
 * 1. At a reference value (y = {0, 0, 0.5})
 * 2. Just before the cutoff (y = {0, 0, 1e-31})
 * 3. Just after the cutoff (y = {0, 0, 1e-33})
 * 4. At zero (y = {0, 0, 0})
 *
 * It then compares these results to check if:
 * - The result after cutoff is the same as at zero
 * - In the case that the result just before the cutoff is different from the
 * reference result: That the result before cutoff is different from after cutoff
 *
 * @return number of failed tests.
 */
bool test_epsteinZeta_cutoff() {
    printf("%s ", __func__);
    double nu;
    unsigned int dim = 3;
    double m[] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Identity matrix 3x3
    double x[] = {0, 0, 0};
    double y_ref[] = {0, 0, 0.5};
    double y_before[] = {0, 0, 1e-31};
    double y_after[] = {0, 0, 1e-33};
    double y_zero[] = {0, 0, 0};

    int testsPassed = 0;
    int totalTests = 0;

    double tol = 1e-15; // Tolerance for comparison

    for (int i = 0; i < 80 + 1; i++) {
        nu = 3 - (double)i / 4;

        double complex zetaRef = epsteinZeta(nu, dim, m, x, y_ref);
        double complex zetaBeforeCutoff = epsteinZeta(nu, dim, m, x, y_before);
        double complex zetaAfterCutoff = epsteinZeta(nu, dim, m, x, y_after);
        double complex zetaZero = epsteinZeta(nu, dim, m, x, y_zero);

        // Check if after cutoff and zero are the same
        if (cabs(zetaAfterCutoff - zetaZero) > tol) {
            reportEpsteinZetaCutoffError(
                "zetaAfterCutoff and zetaZero are not equal", zetaAfterCutoff,
                zetaZero, nu, y_after, dim);
        } else if (cabs(zetaRef - zetaBeforeCutoff) >= tol &&
                   cabs(zetaBeforeCutoff - zetaAfterCutoff) <= tol) {
            // Check if before cutoff and after cutoff are different
            reportEpsteinZetaCutoffError(
                "zetaBeforeCutoff and zetaAfterCutoff are not different",
                zetaBeforeCutoff, zetaAfterCutoff, nu, y_before, dim);
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
 * @brief Main function to run all Epstein Zeta function tests.
 *
 * @return number of failed tests.
 */
int main() {
    int failedQ1 = test_epsteinZeta_epsteinZetaReg();
    int failedQ2 = test_epsteinZeta_epsteinZetaReg_represent_as_each_other();
    int failedQ3 = test_epsteinZeta_cutoff();
    return failedQ1 + failedQ2 + failedQ3;
}
