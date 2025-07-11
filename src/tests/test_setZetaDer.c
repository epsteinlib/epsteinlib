// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "epsteinZeta.h"
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
 * @brief Benchmarks 1D setZetaDer function by comparing to high-precision values
 * from mathematica analytic implementation.
 *
 * @return number of failed tests.
 * */
int test_setZetaDer1D(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/setZetaDer1D_Ref.csv", // NOLINT
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

    int testsPassed = 0;
    int totalTests = 0;
    unsigned int dim = 1;
    double tol = pow(10, -12);

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double *nuRef = malloc(sizeof(double));
    double a[] = {1.};
    double x[] = {0.};
    double *y = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(2 * sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, y, alpha, {Re[result], Im[result]}
        scanResult = sscanf( // NOLINT
            line, "%lf,%lf,%u,%lf,%lf", nuRef, y, alpha, refRead, refRead + 1);

        if (scanResult != 5) {
            printf("\n\t ");
            printf("Error reading line: %s", line);
            printf("\t ");
            printf("Scanned %d values instead of 5", scanResult);
            continue;
        }

        nu = nuRef[0];

        num = setZetaDer(nu, dim, a, x, y, alpha);
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
            printf("setZetaDer: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printMatrixUnitTest("a:", a, dim);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(nuRef);
    free(y);
    free(alpha);
    free(refRead);

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
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
 * @brief Benchmarks 2D setZetaDer function by comparing to high-precision values
 * from mathematica prototype over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_setZetaDer_prototype(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/setZetaDer_prototype_Ref.csv", // NOLINT
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

    int testsPassed = 0;
    int totalTests = 0;
    unsigned int dim = 2;
    double tol = pow(10, -12);

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double *nuRef = malloc(sizeof(double));
    double *a = malloc((unsigned long)dim * (unsigned long)dim * sizeof(double));
    double *x = malloc(dim * sizeof(double));
    double *y = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(2 * sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {a11, a12, a21, a22}, {x1, x2}, {y1, y2}, {alpha1, alpha2},
        // {Re[result], Im[result]}
        scanResult = sscanf( // NOLINT
            line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%u,%u,%lf,%lf", nuRef, a,
            a + 1, a + 2, a + 3, x, x + 1, y, y + 1, alpha, alpha + 1, refRead,
            refRead + 1);

        if (scanResult != 13) {
            printf("\n\t ");
            printf("Error reading line: %s", line);
            printf("\t ");
            printf("Scanned %d values instead of 13", scanResult);
            continue;
        }

        nu = nuRef[0];

        num = setZetaDer(nu, dim, a, x, y, alpha);
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
            printf("setZetaDer: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printMatrixUnitTest("a:", a, dim);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(nuRef);
    free(a);
    free(x);
    free(y);
    free(alpha);
    free(refRead);

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
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
 * @brief Benchmarks 2D set zeta derivatives by computing its taylor series.
 *
 * @return number of failed tests.
 */
int test_setZetaDer_taylor(void) { // NOLINT
    printf("%s ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex valRef;
    double complex valTaylor;

    double tol = 5 * pow(10, -12);
    unsigned int dim = 2;
    unsigned int order = 12;

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double nu = 0.5;
    double m[] = {1., 0.3, 0.3, 1.}; // Non-diagonal matrix with det not 1
    double yDiff[] = {0.01, 0.005};
    unsigned int alpha0[] = {0, 0};
    double *x = malloc(dim * sizeof(double));
    double *y0 = malloc(dim * sizeof(double));
    double *yPlus = malloc(dim * sizeof(double));

    int testsPassed = 0;
    int totalTests = 0;

    bool done;
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));

    for (int i = 0; i < 100; i++) {

        nu = -12.5 + 0.333 * (i + 1);

        x[0] = 0.003 * i;
        x[1] = -0.002 * i;

        y0[0] = -1.1 + 0.1 * i;
        y0[1] = -2.02 + 0.05 * i;

        for (int i = 0; i < dim; i++) {
            yPlus[i] = y0[i] + yDiff[i];
        }

        valRef = setZetaDer(nu, dim, m, x, yPlus, alpha0);

        // build taylor series around z
        valTaylor = 0;

        // Initialize multi-index
        for (int i = 0; i < dim; i++) {
            alpha[i] = 0;
        }

        // Iterate over every multi-index alpha so that every alpha[] < order
        while (true) {

            valTaylor += mult_pow(dim, alpha, yDiff) / (double)mult_fac(dim, alpha) *
                         setZetaDer(nu, dim, m, x, y0, alpha);

            done = true;
            for (unsigned int idx = 0; idx < dim; idx++) {
                if (alpha[idx] + 1 <= order) {
                    alpha[idx]++;
                    done = false;
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

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n");
            printf("Warning! ");
            printf("setZetaDer: ");
            printf(" %0*.16lf %+.16lf I (as a taylor series) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(valTaylor), cimag(valTaylor), creal(valRef),
                   cimag(valRef));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y0:\t\t", y0, dim);
            printVectorUnitTest("yPlus:\t\t", yPlus, dim);
            printVectorUnitTest("yDiff:\t\t", yDiff, dim);
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
 * @brief Benchmarks 2D setZetaDer function by comparing to reference values of the
 * laplacian of set zeta function obtained by finite differences.
 *
 * @return number of failed tests.
 * */
int test_setZetaDer_laplace(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/setZetaDer_laplace_Ref.csv", // NOLINT
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

    int testsPassed = 0;
    int totalTests = 0;
    unsigned int dim = 2;
    double tol = pow(10, -6);

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double *nuRef = malloc(sizeof(double));
    double *a = malloc((unsigned long)dim * (unsigned long)dim * sizeof(double));
    double *x = malloc(dim * sizeof(double));
    double *y = malloc(dim * sizeof(double));
    double *refRead = malloc(2 * sizeof(double));

    unsigned int alpha20[] = {2, 0};
    unsigned int alpha11[] = {1, 1};
    unsigned int alpha02[] = {0, 2};

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {a11, a12, a21, a22}, {x1, x2}, {y1, y2}, {Re[result],
        // Im[result]}
        scanResult = sscanf( // NOLINT
            line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", nuRef, a, a + 1,
            a + 2, a + 3, x, x + 1, y, y + 1, refRead, refRead + 1);

        if (scanResult != 11) {
            printf("\n\t ");
            printf("Error reading line: %s", line);
            printf("\t ");
            printf("Scanned %d values instead of 11", scanResult);
            continue;
        }

        nu = nuRef[0];

        num = setZetaDer(nu, dim, a, x, y, alpha20) +
              2 * setZetaDer(nu, dim, a, x, y, alpha11) +
              setZetaDer(nu, dim, a, x, y, alpha02);
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
            printf("setZetaDer: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printMatrixUnitTest("a:", a, dim);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(nuRef);
    free(a);
    free(x);
    free(y);
    free(refRead);

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
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
 * @brief Benchmarks 2D setZetaDer function by comparing to reference values of the
 * squared laplacian of set zeta function obtained by finite differences.
 *
 * @return number of failed tests.
 * */
int test_setZetaDer_laplace2(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/setZetaDer_laplace2_Ref.csv", // NOLINT
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

    int testsPassed = 0;
    int totalTests = 0;
    unsigned int dim = 2;
    double tol = pow(10, -2);

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double *nuRef = malloc(sizeof(double));
    double *a = malloc((unsigned long)dim * (unsigned long)dim * sizeof(double));
    double *x = malloc(dim * sizeof(double));
    double *y = malloc(dim * sizeof(double));
    double *refRead = malloc(2 * sizeof(double));

    unsigned int alpha40[] = {4, 0};
    unsigned int alpha31[] = {3, 1};
    unsigned int alpha22[] = {2, 2};
    unsigned int alpha13[] = {1, 3};
    unsigned int alpha04[] = {0, 4};

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, {a11, a12, a21, a22}, {x1, x2}, {y1, y2}, {Re[result],
        // Im[result]}
        scanResult = sscanf( // NOLINT
            line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", nuRef, a, a + 1,
            a + 2, a + 3, x, x + 1, y, y + 1, refRead, refRead + 1);

        if (scanResult != 11) {
            printf("\n\t ");
            printf("Error reading line: %s", line);
            printf("\t ");
            printf("Scanned %d values instead of 11", scanResult);
            continue;
        }

        nu = nuRef[0];

        num = setZetaDer(nu, dim, a, x, y, alpha40) +
              4 * setZetaDer(nu, dim, a, x, y, alpha31) +
              6 * setZetaDer(nu, dim, a, x, y, alpha22) +
              4 * setZetaDer(nu, dim, a, x, y, alpha13) +
              setZetaDer(nu, dim, a, x, y, alpha04);

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
            printf("err nb %u\n", totalTests);
            printf("Warning! ");
            printf("setZetaDer: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printMatrixUnitTest("a:", a, dim);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(nuRef);
    free(a);
    free(x);
    free(y);
    free(refRead);

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
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
    int failed1 = test_setZetaDer1D();
    int failed2 = test_setZetaDer_prototype();
    int failed3 = test_setZetaDer_taylor();
    int failed4 = test_setZetaDer_laplace();
    int failed5 = test_setZetaDer_laplace2();
    return failed1 + failed2 + failed3 + failed4 + failed5;
}
