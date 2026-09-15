// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "../src/tools.h"
#include "epsteinZeta.h"
#include "utils.h"
#include "wrappers.h"
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
static int test_setZetaDer_1D(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/setZetaDer_1D_ref.csv", // NOLINT
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
    double *a = malloc((unsigned long)dim * (unsigned long)dim * sizeof(double));
    double *x = malloc(dim * sizeof(double));
    double *y = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(2 * sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: nu, a, x, y, alpha, {Re[result], Im[result]}
        scanResult = sscanf( // NOLINT
            line, "%lf,%lf,%lf,%lf,%u,%lf,%lf", nuRef, a, x, y, alpha, refRead,
            refRead + 1);

        if (scanResult != 7) {
            printf("\n\t ");
            printf("Error reading line: %s", line);
            printf("\t ");
            printf("Scanned %d values instead of 7", scanResult);
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

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks 2D setZetaDer function by comparing to high-precision values
 * from mathematica prototype over a range of random parameters.
 *
 * @return number of failed tests.
 * */
static int test_setZetaDer_2D(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/setZetaDer_2D_ref.csv", // NOLINT
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
    double tol = 5 * pow(10, -12);

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

        // For ref exactly zero, use the epsteinZetaAniso without (potentially big)
        // prefactor! In the future, this file needs to be rewritten for the
        // anisotropic variant to avoid such hacks
        if (cabs(ref) == 0) {
            num = epsteinZetaAniso(nu, dim, a, x, y, alpha);
        }

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

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks 2D epsteinZetaAniso function by comparing to high-precision
 * values from mathematica prototype over a range of random parameters.
 *
 * @return number of failed tests.
 * */
static int test_epsteinZetaAniso_2D_strip(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/epsteinZetaAniso_2D_strip_ref.csv", // NOLINT
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
    int reported = 0;
    unsigned int dim = 2;

    // Known-bad value
    double tol = 100;

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

        num = epsteinZetaAniso(nu, dim, a, x, y, alpha);
        ref = refRead[0] + refRead[1] * I;

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else if (reported < MAX_REPORTS) {
            reported++;
            printf("\n\n");
            printf("Warning! ");
            printf("epsteinZetaAniso: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t         != "
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
            if (reported == MAX_REPORTS) {
                printf("\n\t ... ");
                printf("further failures suppressed");
            }
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

    reportImprovedUnitTest(__func__, errMax, 1E-8);

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks 2D epsteinZetaAniso function by comparing to high-precision
 * values from mathematica prototype over a range of random parameters.
 *
 * @return number of failed tests.
 * */
static int test_epsteinZetaAniso_2D_highorder(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/epsteinZetaAniso_2D_highorder_ref.csv", // NOLINT
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
    int reported = 0;
    unsigned int dim = 2;

    // Known-bad baseline
    double tol = 5E-09;

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

        num = epsteinZetaAniso(nu, dim, a, x, y, alpha);
        ref = refRead[0] + refRead[1] * I;

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);

        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else if (reported < MAX_REPORTS) {
            reported++;
            printf("\n\n");
            printf("Warning! ");
            printf("epsteinZetaAniso: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t         != "
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
            if (reported == MAX_REPORTS) {
                printf("\n\t ... ");
                printf("further failures suppressed");
            }
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

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks 2D set zeta derivatives by computing its taylor series.
 *
 * @return number of failed tests.
 */
static int test_setZetaDer_taylor(void) { // NOLINT
    printf("%s ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex valRef;
    double complex valTaylor;

    double tol = 5 * pow(10, -15);
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

    printf("\n\t ... ");
    printf("generating test values");
    for (int i = 0; i < 10; i++) {

        nu = -12.5 + 3.33 * (i + 1);

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

    free(x);
    free(y0);
    free(yPlus);
    free(alpha);

    printf("\n\t ... ");
    printf("%d  out of %d  tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / totalTests);
    printf("\n");

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return totalTests - testsPassed;
}

/*!
 * @brief Tests the pole structure of epsteinZetaAniso.
 *
 * @return number of failed tests.
 * */
static int test_epsteinZetaAniso_poles(void) { // NOLINT
    printf("%s ", __func__);

    unsigned int dim = 2;
    unsigned int maxPrints = 10;

    double aUnit[] = {1., 0., 0., 1.};
    double aHex[] = {1., 0.5, 0., 0.8660254037844386};
    double aObl[] = {1., 0.3, 0., 1.};  // no z1 -> -z1 mirror symmetry
    double aScal[] = {2., 0., 0., 1.5}; // det(A) = 3
    double *lattice[] = {aUnit, aHex, aObl, aScal};

    // x is irrelevant to the pole, y is not. Both stay in their elementary cell.
    double x[][2] = {{0., 0.}, {0.3, 0.}, {0.3, 0.2}};
    double y[][2] = {{0., 0.}, {0.1, 0.2}};

    unsigned int alpha[][2] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 0},
                               {0, 2}, {2, 1}, {3, 1}, {2, 2}, {4, 0}};
    double nu[] = {1., 2., 2.5, 3., 4., 4.5, 5., 6., 7., 8.};

    unsigned int lattices = sizeof(lattice) / sizeof(lattice[0]);
    unsigned int xs = sizeof(x) / sizeof(x[0]);
    unsigned int ys = sizeof(y) / sizeof(y[0]);
    unsigned int alphas = sizeof(alpha) / sizeof(alpha[0]);
    unsigned int nus = sizeof(nu) / sizeof(nu[0]);

    int testsPassed = 0;
    int totalTests = 0;
    int polesFound = 0;
    unsigned int printed = 0;

    double complex num;
    bool allEven;
    bool yZero;
    bool expectNan;
    bool isNan;
    unsigned int alphaAbs;

    printf("\n\t ... ");
    printf("sweeping %u lattices, %u x, %u y, %u multi-indices and %u exponents",
           lattices, xs, ys, alphas, nus);
    for (unsigned int iA = 0; iA < lattices; iA++) {
        for (unsigned int iX = 0; iX < xs; iX++) {
            for (unsigned int iY = 0; iY < ys; iY++) {
                yZero = (y[iY][0] == 0.) && (y[iY][1] == 0.);
                for (unsigned int iAlpha = 0; iAlpha < alphas; iAlpha++) {

                    allEven = true;
                    alphaAbs = 0;
                    for (unsigned int j = 0; j < dim; j++) {
                        alphaAbs += alpha[iAlpha][j];
                        if (alpha[iAlpha][j] % 2 != 0) {
                            allEven = false;
                        }
                    }

                    for (unsigned int iNu = 0; iNu < nus; iNu++) {
                        expectNan = allEven && yZero &&
                                    (nu[iNu] == (double)(dim + alphaAbs));
                        polesFound += expectNan ? 1 : 0;

                        num = epsteinZetaAniso(nu[iNu], dim, lattice[iA], x[iX],
                                               y[iY], alpha[iAlpha]);
                        isNan = isnan(creal(num)) || isnan(cimag(num));

                        if (isNan == expectNan) {
                            testsPassed++;
                        } else if (printed < maxPrints) {
                            printed++;
                            printf("\n\n");
                            printf("Warning! ");
                            printf("epsteinZetaAniso: ");
                            printf(" %0*.16lf %+.16lf I (this implementation) "
                                   "\n\t\t!= %s (expected)\n",
                                   4, creal(num), cimag(num),
                                   expectNan ? "NaN" : "a finite value");
                            printf("\n");
                            printf("nu:\t\t %.16lf\n", nu[iNu]);
                            printMatrixUnitTest("a:", lattice[iA], dim);
                            printVectorUnitTest("x:\t\t", x[iX], dim);
                            printVectorUnitTest("y:\t\t", y[iY], dim);
                            printMultiindexUnitTest("alpha:\t\t", alpha[iAlpha],
                                                    dim);
                            printf("\n");
                        }
                        totalTests++;
                    }
                }
            }
        }
    }

    if (printed == maxPrints) {
        printf("\n\t ... ");
        printf("further failures suppressed");
    }

    printf("\n\t ... ");
    printf("%d out of %d tests passed.", testsPassed, totalTests);
    printf("\t\t\t\t    ");
    printf("[ %d poles checked\t\t ]", polesFound);
    printf("\n");

    return totalTests - testsPassed;
}

/*!
 * @brief Sweeps one lattice and checks the components that have to vanish.
 *
 * Inversion: conj(Z_alpha(0, y)) = Z_alpha(0, -y) = (-1)^|alpha| Z_alpha(0, y).
 *
 * Half reciprocal wavevector: for 2y in Lambda* the value is real by
 * Lambda*-periodicity, so odd |alpha| vanishes there too. Only checked where
 * checkHalfReci is set, since on lattices whose Fourier window is not symmetric
 * about -y the vanishing holds only up to the truncation error.
 *
 * Mirroring: epsteinZetaAniso vanishes whenever alpha_i is odd, the lattice is
 * mirror symmetric in the i'th component and either 2 x_i e_i is a lattice vector
 * with y_i = 0, or x_i = 0 with 2 y_i e_i a reciprocal lattice vector. The two
 * cases are exchanged by swapping Lambda with Lambda* and x with y, and neither
 * implies the other.
 *
 * @param[in] dim: dimension of the lattice.
 * @param[in] m: lattice basis matrix A, row major, basis vectors in the columns.
 * @param[in] mirror: mirror symmetry of the lattice, one flag per component.
 * @param[in] xs, ys, alphas: flat arrays with stride dim.
 * @param[in] numX, numY, numAlpha: number of entries in xs, ys and alphas.
 * @param[in] nus: exponents nu to sweep over.
 * @param[in] numNu: number of entries in nus.
 * @param[in] checkMirror: also check the mirroring statement.
 * @param[in] checkHalfReci: also check the half reciprocal wavevector statement.
 * @param[in] tol: a component counts as vanishing if its absolute value is
 * below tol.
 * @param[in,out] total: running number of checked components.
 * @param[in,out] reported: running number of printed failures.
 * @param[in,out] errMin, errMax, errSum: running minimum, maximum and sum of
 * the absolute values of the checked components.
 */
static int inversionZeroSweep(unsigned int dim, const double *m, // NOLINT
                              const bool *mirror, const double *xs, int numX,
                              const double *ys, int numY, const unsigned int *alphas,
                              int numAlpha, const double *nus, int numNu,
                              bool checkMirror, bool checkHalfReci, double tol,
                              int *total, unsigned int *reported, double *errMin,
                              double *errMax, double *errSum) {
    const char *partName[] = {"Real", "Imaginary"};
    int failed = 0;

    bool diagonal = true;
    for (unsigned int a = 0; a < dim; a++) {
        for (unsigned int b = 0; b < dim; b++) {
            diagonal = diagonal && ((a == b) || (m[(a * dim) + b] == 0.));
        }
    }

    for (int itx = 0; itx < numX; itx++) {
        const double *x = xs + ((size_t)dim * itx);
        for (int ity = 0; ity < numY; ity++) {
            const double *y = ys + ((size_t)dim * ity);

            // 2y in Lambda* <=> A^T (2y) in Z^d, so Z_alpha(0, -y) = Z_alpha(0, y)
            // by Lambda*-periodicity and the value is real
            bool yHalfReci = true;
            for (unsigned int i = 0; i < dim; i++) {
                double t = 0.;
                for (unsigned int j = 0; j < dim; j++) {
                    t += m[(j * dim) + i] * 2. * y[j];
                }
                yHalfReci = yHalfReci && (fabs(t - nearbyint(t)) < 1e-12);
            }

            for (int ia = 0; ia < numAlpha; ia++) {
                const unsigned int *alpha = alphas + ((size_t)dim * ia);
                unsigned int alphaAbs = mult_abs(dim, alpha);

                bool xZero = true;
                bool yZero = true;
                bool mirrorZero = false;
                for (unsigned int i = 0; i < dim; i++) {
                    xZero = xZero && (x[i] == 0.);
                    yZero = yZero && (y[i] == 0.);
                    if (!checkMirror || (alpha[i] % 2 == 0) || !mirror[i]) {
                        continue;
                    }
                    // 2 x_i e_i in Lambda, so that the reflection about x in the
                    // i'th component is a lattice symmetry. Decided from the
                    // diagonal entry, since the general case needs A^-1.
                    bool xHalfLat = (x[i] == 0.);
                    if (!xHalfLat && diagonal) {
                        double t = 2. * x[i] / m[(i * dim) + i];
                        xHalfLat = (t == nearbyint(t));
                    }
                    // 2 y_i e_i in Lambda*, that is A^T (2 y_i e_i) integer, which
                    // is 2 y_i times the i'th row of m and needs no inverse
                    bool yHalfLat = true;
                    for (unsigned int a = 0; a < dim; a++) {
                        double t = 2. * y[i] * m[(i * dim) + a];
                        yHalfLat = yHalfLat && (fabs(t - nearbyint(t)) < 1e-12);
                    }
                    // the statement and its dual, with Lambda and Lambda*, x and y
                    // exchanged
                    mirrorZero = mirrorZero || (xHalfLat && (y[i] == 0.)) ||
                                 ((x[i] == 0.) && yHalfLat);
                }

                bool odd = (alphaAbs % 2) != 0;

                // index 0 is the real part, index 1 the imaginary part
                bool vanishes[2];
                vanishes[0] = mirrorZero || (xZero && odd);
                vanishes[1] =
                    mirrorZero ||
                    (xZero && (!odd || yZero || (checkHalfReci && yHalfReci)));
                if (!vanishes[0] && !vanishes[1]) {
                    continue;
                }

                for (int in = 0; in < numNu; in++) {
                    double nu = nus[in];

                    // alpha = 0 exercises the isotropic front end
                    double complex num = epsteinZetaAniso(nu, dim, m, x, y, alpha);
                    if (!isfinite(creal(num)) || !isfinite(cimag(num))) {
                        continue;
                    }

                    double part[2] = {creal(num), cimag(num)};
                    for (unsigned int c = 0; c < 2; c++) {
                        if (!vanishes[c]) {
                            continue;
                        }
                        double err = fabs(part[c]);
                        *errMin = (*errMin < err) ? *errMin : err;
                        *errMax = (*errMax > err) ? *errMax : err;
                        *errSum += err;
                        (*total)++;

                        if (err < tol) {
                            continue;
                        }
                        failed++;
                        if (*reported < MAX_REPORTS) {
                            (*reported)++;
                            printf("\n\n");
                            printf("Warning! ");
                            printf("epsteinZetaAniso: ");
                            printf(" %0*.16lf %+.16lf I\n", 4, creal(num),
                                   cimag(num));
                            printf("\t\t\t    %s part should vanish\n", partName[c]);
                            printf("|vanishing part|:           %E !< %E  "
                                   "(tolerance)\n",
                                   err, tol);
                            printf("\n");
                            printf("nu:\t\t %.16lf\n", nu);
                            printMatrixUnitTest("a:", m, dim);
                            printVectorUnitTest("x:\t\t", x, dim);
                            printVectorUnitTest("y:\t\t", y, dim);
                            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
                            printf("\n");
                        }
                    }
                }
            }
        }
    }

    return failed;
}

/*!
 * @brief Checks the components of Epstein zeta that vanish by inversion and by
 * mirroring, over a grid that reaches the strongly anisotropic regime where a
 * cancellation error in the lattice sums grows with the anisotropy order.
 *
 * See inversionZeroSweep for the two identities. alpha = 0 is included and
 * exercises the isotropic front end. Needs no reference values.
 *
 * @return number of failed tests.
 */
static int test_epsteinZetaAniso_inversionZeros(void) { // NOLINT
    printf("%s ", __func__);

    double tol = pow(10, -15);
    int failed = 0;
    int total = 0;
    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;
    unsigned int reported = 0;

    /* ------------------------------- 2D ------------------------------- */
    double a2[5][4] = {{1., 0., 0., 1.},                  // identity
                       {13. / 10, 0., 0., 7. / 10},       // diagonal
                       {1., 0.5, 0., 0.8660254037844386}, // hexagonal
                       {1., 0.3, 0., 1.},                 // oblique
                       {0., 1., 1., 0.}};                 // det(A) < 0
    // the oblique lattice is not mirror symmetric in either component
    bool mirror2[5][2] = {
        {true, true}, {true, true}, {true, true}, {false, false}, {true, true}};

    // Inversion: x = 0, every y, up to |alpha| = 31, where the residue of the
    // vanishing component is unmistakable if the summation order is wrong.
    double xInv2[] = {0., 0.};
    double yInv2[] = {0.,   0.,   // origin, where the whole value vanishes
                      1e-8, 0.,   // near the origin
                      0.5,  0.,   // cell boundary
                      0.,   0.21, // zero in the first component only
                      0.1,  0.2}; // generic
    unsigned int alphaInv2[] = {0, 0, 1, 0, 0, 1, 1,  1, 2,  0, 2,  1,
                                3, 1, 5, 2, 7, 0, 12, 0, 21, 2, 31, 0};
    double nuInv2[] = {2.5, 8., 30.};

    // Mirroring: one component of x and y vanishes and the others do not. Both
    // the anisotropy order and the exponent stay moderate here, since the residue
    // is set by the magnitude of the summands, which grows like the cutoff radius
    // to the power |alpha| and like the distance to the nearest lattice point to
    // the power -nu.
    double xMir2[] = {0.,        37. / 100,  // zero in the first component only
                      29. / 100, 0.,         // zero in the second component only
                      0.5,       0.};        // half lattice point
    double yMir2[] = {0.,        0.,         //
                      0.,        21. / 100,  //
                      13. / 100, 0.,         //
                      0.5,       21. / 100}; // half reciprocal first component
    unsigned int alphaMir2[] = {1, 0, 0, 1, 1, 1, 3, 1, 5, 2, 9, 2};
    double nuMir2[] = {2.5, 8.};

    for (int im = 0; im < 5; im++) {
        failed += inversionZeroSweep(2, a2[im], mirror2[im], xInv2, 1, yInv2, 5,
                                     alphaInv2, 12, nuInv2, 3, false, false, tol,
                                     &total, &reported, &errMin, &errMax, &errSum);
        failed += inversionZeroSweep(2, a2[im], mirror2[im], xMir2, 3, yMir2, 4,
                                     alphaMir2, 6, nuMir2, 2, true, false, tol,
                                     &total, &reported, &errMin, &errMax, &errSum);
    }

    /* ------------------------------- 3D ------------------------------- */
    double a3[3][9] = {{1., 0., 0., 0., 1., 0., 0., 0., 1.},     // identity
                       {2., 0., 0., 0., 3. / 2, 0., 0., 0., 1.}, // diagonal
                       {1., 0.3, 0., 0., 1., 0., 0., 0., 1.}};   // sheared
    // the sheared lattice loses the mirror symmetry in the first two components
    bool mirror3[3][3] = {
        {true, true, true}, {true, true, true}, {false, false, true}};

    double xInv3[] = {0., 0., 0.};
    double yInv3[] = {0., 0., 0., 0., 13. / 100, 29. / 100, 0.1, 0.2, 0.3};
    unsigned int alphaInv3[] = {0, 0, 0, 1, 0, 0, 1, 1, 0,
                                3, 1, 2, 5, 2, 1, 9, 2, 2};
    double nuInv3[] = {2.5, 20.};

    double xMir3[] = {0., 31. / 100, 47. / 100};
    double yMir3[] = {0., 0., 0., 0., 13. / 100, 29. / 100};
    unsigned int alphaMir3[] = {1, 0, 0, 1, 1, 0, 3, 1, 2, 5, 2, 1};
    double nuMir3[] = {2.5, 8.};

    for (int im = 0; im < 3; im++) {
        failed += inversionZeroSweep(3, a3[im], mirror3[im], xInv3, 1, yInv3, 3,
                                     alphaInv3, 6, nuInv3, 2, false, false, tol,
                                     &total, &reported, &errMin, &errMax, &errSum);
        failed += inversionZeroSweep(3, a3[im], mirror3[im], xMir3, 1, yMir3, 2,
                                     alphaMir3, 4, nuMir3, 2, true, false, tol,
                                     &total, &reported, &errMin, &errMax, &errSum);
    }

    /* -------------------- half reciprocal wavevector -------------------- */
    // 2y in Lambda* makes the value real, so odd |alpha| vanishes there too.
    // Only the identity lattice, and only outside the large exponent branch:
    // elsewhere the Fourier window is not symmetric about -y and the vanishing
    // holds only up to the truncation error.
    double yHalf2[] = {0.5, 0.};
    unsigned int alphaHalf2[] = {21, 0, 31, 0};
    double nuHalf2[] = {2.5, 4., 8.};

    failed += inversionZeroSweep(2, a2[0], mirror2[0], xInv2, 1, yHalf2, 1,
                                 alphaHalf2, 2, nuHalf2, 3, false, true, tol, &total,
                                 &reported, &errMin, &errMax, &errSum);

    // 2y = (1, -1/sqrt(3)) is a reciprocal basis vector of the hexagonal lattice,
    // so the value is real and odd |alpha| vanishes, while no single component of
    // 2y spans e_j: the mirror cases cannot catch this one.
    double yHalfHex2[] = {0.5, -0.5 / sqrt(3.)};
    failed += inversionZeroSweep(2, a2[2], mirror2[2], xInv2, 1, yHalfHex2, 1,
                                 alphaHalf2, 2, nuHalf2, 3, false, true, tol, &total,
                                 &reported, &errMin, &errMax, &errSum);

    if (reported >= MAX_REPORTS) {
        printf("\n\t ... ");
        printf("further failures suppressed");
    }

    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", total - failed, total,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / total);
    printf("\n");

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return failed;
}

/*!
 * @brief Tests setZetaDer function for special case nu = dim + Total[alpha] + 2
 * with y = 0, comparing against Mathematica reference values.
 *
 * @return number of failed tests.
 */
static int test_setZetaDer_special_exponents(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/setZetaDer_special_exponents_ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    char line[4096];
    int testsPassed = 0;
    int totalTests = 0;
    double tol = 5 * pow(10, -15);
    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    unsigned int max_test_dim = 4;

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        totalTests++;

        char *ptr = line;
        char *endptr;

        // Parse dimension
        unsigned int dim = (unsigned int)strtoul(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        // Parse nu
        double nu = strtod(ptr, &endptr);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        // Bounds check on dimension
        if (dim < 1 || dim > max_test_dim) {
            (void)fclose(data);
            return fprintf(stderr, "Invalid dimension: %u\n", dim);
        }

        // Use stack allocation with fixed max size
        double a[max_test_dim * max_test_dim];
        double x[max_test_dim];
        double y[max_test_dim];
        unsigned int alpha[max_test_dim];

        // Parse a (dim x dim matrix, flattened)
        int parseOk = 1;
        for (unsigned int i = 0; i < dim * dim && parseOk; i++) {
            a[i] = strtod(ptr, &endptr);
            if (*endptr != ',') {
                parseOk = 0;
            } else {
                ptr = endptr + 1;
            }
        }

        // Parse x (dim values)
        for (unsigned int i = 0; i < dim && parseOk; i++) {
            x[i] = strtod(ptr, &endptr);
            if (*endptr != ',') {
                parseOk = 0;
            } else {
                ptr = endptr + 1;
            }
        }

        // Parse y (dim values)
        for (unsigned int i = 0; i < dim && parseOk; i++) {
            y[i] = strtod(ptr, &endptr);
            if (*endptr != ',') {
                parseOk = 0;
            } else {
                ptr = endptr + 1;
            }
        }

        // Parse alpha (dim unsigned int values)
        for (unsigned int i = 0; i < dim && parseOk; i++) {
            alpha[i] = (unsigned int)strtoul(ptr, &endptr, 10);
            if (*endptr != ',') {
                if (i < dim - 1) {
                    parseOk = 0;
                }
            }
            ptr = endptr + 1;
        }

        if (!parseOk) {
            continue;
        }

        // Parse reference result (Re, Im)
        double refRe = strtod(ptr, &endptr);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;
        double refIm = strtod(ptr, &endptr);

        // Compute result
        double complex num = setZetaDer(nu, dim, a, x, y, alpha);
        double complex ref = refRe + (refIm * I);

        // Compute errors
        double errorAbs = errAbs(ref, num);
        double errorRel = errRel(ref, num);
        double errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

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
                   "%.16lf %+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Eabs, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("dim:\t\t %u\n", dim);
            printf("nu:\t\t %.16lf\n", nu);
            printMatrixUnitTest("a:", a, dim);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
    }

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

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks 6D setZetaDer function by comparing to reference values of the
 * laplacian of set zeta function.
 *
 * @return number of failed tests.
 * */
static int test_setZetaDer_poly_laplace(void) { // NOLINT
    printf("%s ", __func__);

    double nu;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    double tol = pow(10, -13);

    int testsPassed = 0;
    int totalTests = 0;

    int maxn = 5; // For n'th power of the Laplace operator
    unsigned int dim = 4;
    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double *a = malloc((unsigned long)dim * (unsigned long)dim * sizeof(double));
    double *x = malloc(dim * sizeof(double));
    double *y = malloc(dim * sizeof(double));
    unsigned int *alpha0 = malloc(dim * sizeof(unsigned int));

    unsigned int beta[dim];
    for (int i = 0; i < dim; i++) {
        beta[i] = 0;
    }

    unsigned long long betaFact = 1;
    unsigned int allEven;
    unsigned int kMax;

    for (int i = 0; i < dim; ++i) {
        x[i] = 0.1;
        y[i] = 0.2;
        alpha0[i] = 0;
        for (int j = 0; j < dim; ++j) {
            a[(i * dim) + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    int done = 0;
    unsigned int betaAbs = 0;

    printf("\n\t ... ");
    printf("building the %uD PolyLaplace operator Δⁿ, n = 1,2,...,%u", dim, maxn);
    for (int n = 1; n <= maxn; n++) {

        unsigned long long nFact = 1;
        for (int i = 2; i < n + 1; i++) {
            nFact *= i;
        }

        nu = 0.5;

        // Reference value via shifted Epstein zeta
        ref = pow(-1., n) * pow(-2 * M_PI, 2 * n) *
              setZetaDer(nu - (2 * n), dim, a, x, y, alpha0);

        num = 0. + 0. * I;

        while (1) {
            allEven = 1;
            for (int j = 0; j < dim; ++j) {
                if (beta[j] % 2 != 0) {
                    allEven = 0;
                    break;
                }
            }

            if (allEven && betaAbs == 2 * n) {
                // Correct beta factorial: product of (beta[i]/2)!
                betaFact = 1;
                for (unsigned int j = 0; j < dim; j++) {
                    kMax = beta[j] / 2;
                    for (int k = 1; k <= kMax; k++) {
                        betaFact *= k;
                    }
                }

                num += (double)nFact / (double)betaFact *
                       setZetaDer(nu, dim, a, x, y, beta);
            }

            done = 1;
            for (unsigned int idx = 0; idx < dim; idx++) {
                if (beta[idx] + 1 <= 2 * n) {
                    beta[idx]++;
                    betaAbs++;
                    done = 0;
                    break;
                }
                betaAbs -= beta[idx];
                beta[idx] = 0;
            }

            if (done) {
                break;
            }
        }

        errorAbs = errAbs(ref, num);
        errorRel = errRel(ref, num);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;

        } else {
            printf("\n\nWarning! Poly-Laplace mismatch:\n");
            printf("nu: %.16lf\n", nu);
            printMatrixUnitTest("a:", a, dim);
            printVectorUnitTest("x:\t\t ", x, dim);
            printVectorUnitTest("y:\t\t ", y, dim);

            printf("\nComputed (loop-based) poly-Laplace: %0*.16lf %+.16lf I\n", 4,
                   creal(num), cimag(num));
            printf("Reference (shifted Epstein zeta):   %0*.16lf %+.16lf I\n", 4,
                   creal(ref), cimag(ref));
            printf("Min(Emax, Erel):                    %E !< %E  (tolerance)\n",
                   errorMaxAbsRel, tol);
        }
        totalTests++;
    }

    free(a);
    free(x);
    free(y);
    free(alpha0);

    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / totalTests);
    printf("\n");

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks epsteinZetaAniso for second-order anisotropy alpha = 2 e_1
 * and all-equal vector arguments x = x 1, y = y 1 against the poly-Laplacian
 * identity
 *
 *     Z_{Z^d,nu,2 e_j}(x 1, y 1) = Z_{Lambda,nu-2}(x 1, y 1) / d ,
 *
 * where the pole at nu = d + |alpha| = d + 2 for y = 0 is skipped.
 *
 * @return number of failed tests.
 * */
static int test_epsteinZetaAniso_allEqual(void) { // NOLINT
    printf("%s ", __func__);

    unsigned int maxDim = 4;
    unsigned int maxPrints = 10;
    double tol = pow(10, -15);

    // nu = nuMin + i nuStep, x = i xStep, y = i yStep
    double nuMin = -2.;
    double nuStep = 1.; // hits nu = 2 exactly
    unsigned int nus = 9;
    double xStep = 0.2;
    unsigned int xs = 3;
    double yStep = 0.01; // starts at y = 0 exactly
    unsigned int ys = 9;

    double a[16];
    double x[4];
    double y[4];
    unsigned int alpha[4];

    double nu;
    double xVal;
    double yVal;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;

    int testsPassed = 0;
    int totalTests = 0;
    unsigned int printed = 0;

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double errMaxDim;

    printf("\n\t ... ");
    printf("sweeping d ≤ %u,  nu = %g,%g,...,%g, |x| ≤ %g, |y| ≤ %g", maxDim, nuMin,
           nuMin + nuStep, nuMin + ((nus - 1) * nuStep), (xs - 1) * xStep,
           (ys - 1) * yStep);

    for (unsigned int dim = 1; dim <= maxDim; dim++) {

        for (unsigned int i = 0; i < dim * dim; i++) {
            a[i] = (i % (dim + 1) == 0) ? 1. : 0.;
        }
        for (unsigned int i = 0; i < dim; i++) {
            alpha[i] = 0;
        }
        alpha[0] = 2;

        errMaxDim = 0.;

        for (unsigned int iNu = 0; iNu < nus; iNu++) {
            nu = nuMin + (iNu * nuStep);

            for (unsigned int iX = 0; iX < xs; iX++) {
                xVal = iX * xStep;
                for (unsigned int i = 0; i < dim; i++) {
                    x[i] = xVal;
                }

                for (unsigned int iY = 0; iY < ys; iY++) {
                    yVal = iY * yStep;
                    for (unsigned int i = 0; i < dim; i++) {
                        y[i] = yVal;
                    }

                    // pole of both sides
                    if (yVal == 0. && nu == (double)(dim + 2)) {
                        continue;
                    }

                    num = epsteinZetaAniso(nu, dim, a, x, y, alpha);
                    ref = epsteinZeta(nu - 2., dim, a, x, y) / (double)dim;

                    errorAbs = errAbs(ref, num);
                    errorRel = errRel(ref, num);
                    if (cabs(ref) == 0.) {
                        errorRel = errorAbs;
                    }

                    errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

                    errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
                    errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
                    errSum += errorMaxAbsRel;

                    if (!(errorMaxAbsRel < errMaxDim)) {
                        errMaxDim = errorMaxAbsRel;
                    }

                    if (errorMaxAbsRel < tol) {
                        testsPassed++;
                    } else if (printed < maxPrints) {
                        printed++;
                        printf("\n\n");
                        printf("Warning! ");
                        printf("epsteinZetaAniso: ");
                        printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t    "
                               "      != "
                               "%.16lf %+.16lf I (shifted Epstein zeta / d)\n",
                               4, creal(num), cimag(num), creal(ref), cimag(ref));
                        printf(
                            "Min(Eabs, Erel):             %E !< %E  (tolerance)\n",
                            errorMaxAbsRel, tol);
                        printf("\n");
                        printf("dim:\t\t %u\n", dim);
                        printf("nu:\t\t %.16lf\n", nu);
                        printMatrixUnitTest("a:", a, dim);
                        printVectorUnitTest("x:\t\t", x, dim);
                        printVectorUnitTest("y:\t\t", y, dim);
                        printMultiindexUnitTest("alpha:\t\t", alpha, dim);
                        printf("\n");
                    }
                    totalTests++;
                }
            }
        }
    }

    if (printed == maxPrints) {
        printf("\n\t ... ");
        printf("further failures suppressed");
    }

    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / totalTests);
    printf("\n");

    reportImprovedUnitTest(__func__, errMax, tol * REP_IMPR_THRES);

    return totalTests - testsPassed;
}

/*!
 * @brief Main function to run all set zeta derivatives function tests.
 *
 * @return number of failed tests.
 */
int main() {
    int failed = 0;

    failed += run_timed_test(test_setZetaDer_1D);
    failed += run_timed_test(test_setZetaDer_2D);
    failed += run_timed_test(test_setZetaDer_taylor);
    failed += run_timed_test(test_epsteinZetaAniso_poles);
    failed += run_timed_test(test_epsteinZetaAniso_inversionZeros);
    failed += run_timed_test(test_setZetaDer_special_exponents);
    failed += run_timed_test(test_setZetaDer_poly_laplace);
    failed += run_timed_test(test_epsteinZetaAniso_allEqual);

    // tests for functions with known-bad errors
    failed += run_timed_test(test_epsteinZetaAniso_2D_highorder);
    failed += run_timed_test(test_epsteinZetaAniso_2D_strip);

    return failed != 0;
}
