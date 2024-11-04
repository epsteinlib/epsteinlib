// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "epsteinZeta.h"
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

/**
 * @brief prints vector to terminal, is used in unitTest.
 * @param name: name of the vector as a string.
 * @param vec: vector.
 * @param dim: size of the vector.
 */
void printVectorUnitTest(const char *name, double *vec, int dim) {
    printf("%s[", name);
    for (int i = 0; i < dim; ++i) {
        printf("%.16lf", vec[i]);
        if (i != dim - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

/**
 * @brief prints square matrix to terminal, is used in unitTest.
 * @param name: name of the matrix as a string.
 * @param mat: matrix.
 * @param dim: size of the matrix.
 */
void printMatrixUnitTest(const char *name, const double *mat, int dim) {
    printf("%s", name);
    double matrixEntry;
    for (int i = 0; i < dim; ++i) {
        printf("\t\t [");
        for (int j = 0; j < dim; ++j) {
            matrixEntry = mat[i * dim + j];
            if (pow(matrixEntry, 2) < 100) {
                printf("%.16f", matrixEntry);
            } else {
                printf("%.2e", matrixEntry * pow(10, 9));
            }
            if (j != dim - 1) {
                printf(", ");
            }
        }
        printf("]\n");
    }
}

/*!
 * @brief absolute difference between to complex numbers.
 * @param ref: reference value.
 * @param comp: compparison value.
 * @return sqrt(rev - comp)
 */
double errAbs(double complex ref, double complex comp) {
    double complex diff = ref - comp;
    double norm = creal(diff) * creal(diff) + cimag(diff) * cimag(diff);
    return sqrt(norm);
}

/*!
 * @brief relative difference between to complex numbers.
 * @param ref: reference value.
 * @param comp: compparison value.
 * @return errAbs / norm(ref)
 */
double errRel(double complex ref, double complex comp) {
    if (ref == 0) {
        return errAbs(ref, comp);
    }
    return errAbs(ref, comp) /
           sqrt(creal(ref) * creal(ref) + cimag(ref) * cimag(ref));
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
int epsteinZetaTest() {
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

    while (fscanf(zetaRegRefData, // NOLINT
                  "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", a, a + 1, a + 2,
                  a + 3, nu, nu + 1, x, x + 1, y, y + 1, zetaRef,
                  zetaRef + 1) == 12) {
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
    if (fclose(zetaRegRefData) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    result = snprintf(path, sizeof(path), "%s/epsteinZetaReg_Ref.csv", // NOLINT
                      BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    zetaRegRefData = fopen(path, "r");
    if (zetaRegRefData == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    while (fscanf(zetaRegRefData, // NOLINT
                  "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", a, a + 1, a + 2,
                  a + 3, nu, nu + 1, x, x + 1, y, y + 1, zetaRef,
                  zetaRef + 1) == 12) {
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
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    free(a);
    free(nu);
    free(x);
    free(y);
    free(zetaRef);

    printf("%d out of %d tests passed\n", testsPassed, totalTests);
    return (testsPassed == totalTests) ? 0 : 1;
}

int main() {
    int result = epsteinZetaTest();
    return result;
}
