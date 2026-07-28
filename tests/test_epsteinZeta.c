// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file test_epsteinZeta.c
 * @brief Benchmarking of the Epstein zeta function.
 */

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
 * @brief Maximum number of failures reported per test, to keep logs readable.
 */
enum { ANISO_MAX_REPORTS = 20 };

/*!
 * @brief Running statistics for the anisotropic reduction sweeps.
 */
typedef struct {
    int passed;
    int total;
    int reported;
    double errMin;
    double errMax;
    double errSum;
} anisoStats;

/*!
 * @brief Free memory allocated for test resources.
 *
 * This function deallocates memory that was dynamically allocated for various
 * arrays used in the Epstein zeta function tests.
 *
 * @param[in] a Pointer to the matrix of coefficients.
 * @param[in] nu Pointer to the array of nu values.
 * @param[in] x Pointer to the array of x values.
 * @param[in] y Pointer to the array of y values.
 * @param[in] zetaRef Pointer to the array of reference zeta values.
 */
static void freeTestResources(double *a, double *nu, double *x, double *y,
                              double *zetaRef) {
    free(a);
    free(nu);
    free(x);
    free(y);
    free(zetaRef);
}

/*!
 * @brief Reports an error when the Epstein zeta function representation test fails.
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
static void reportEpsteinZetaError(double complex valZeta, double complex valZetaReg,
                                   double errorMaxAbsRel, double tol, double *m,
                                   unsigned int dim, double nu, double *x,
                                   double *y) {
    printf("\n");
    printf("Warning! ");
    printf("epsteinZeta:");
    printf(" %0*.16lf %+.16lf I (epsteinZeta) \n\t\t   != "
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
 * @brief Reports an error when the Epstein zeta function cutoff test fails.
 *
 * @param[in] testCase A string describing the specific test case that failed.
 * @param[in] zeta1 The first computed value of the Epstein zeta function.
 * @param[in] zeta2 The second computed value of the Epstein zeta function for
 * comparison.
 * @param[in] nu The parameter nu of the Epstein zeta function.
 * @param[in] y Pointer to the array of y values.
 * @param[in] dim The dimension of the lattice.
 */
static void reportEpsteinZetaCutoffError(const char *testCase, double complex zeta1,
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
 * @brief Computes out = m . v for a row-major dim x dim matrix.
 */
static void anisoMatVec(unsigned int dim, const double *m, const double *v,
                        double *out) {
    for (unsigned int i = 0; i < dim; i++) {
        out[i] = 0.;
        for (unsigned int j = 0; j < dim; j++) {
            out[i] += m[(dim * i) + j] * v[j];
        }
    }
}

/*!
 * @brief Error measure for the reduction tests. Two NaNs count as agreement (the
 * nu = d pole with y in the reciprocal lattice is NaN on both sides), one NaN does
 * not. Two exact zeros count as agreement (nu a non-positive even integer).
 */
static double anisoErr(double complex a, double complex b) {
    bool aNan = isnan(creal(a)) || isnan(cimag(a));
    bool bNan = isnan(creal(b)) || isnan(cimag(b));
    if (creal(a) == creal(b) && cimag(a) == cimag(b)) {
        return 0.; // includes matching infinities
    }
    if (aNan || bNan) {
        return (aNan && bNan) ? 0. : INFINITY;
    }
    double diff = cabs(a - b);
    double scale = fmax(cabs(a), cabs(b));
    if (scale == 0.) {
        return 0.;
    }
    double rel = diff / scale;
    return (diff < rel) ? diff : rel;
}

/*!
 * @brief Reports a single reduction failure.
 */
static void reportAnisoError(unsigned int dim, const double *m, double nu,
                             const double *x, const double *y, double complex base,
                             double complex aniso, double err, double tol,
                             bool reg) {
    printf("\nWarning! %s: %.16lf %+.16lf I\n",
           reg ? "epsteinZetaReg" : "epsteinZeta", creal(base), cimag(base));
    printf("\t           ");
    if (reg) {
        printf("    ");
    }
    printf("!= %.16lf %+.16lf I (%s, alpha = 0)\n", creal(aniso), cimag(aniso),
           reg ? "epsteinZetaAnisoReg" : "epsteinZetaAniso");
    printf("Min(Eabs, Erel):      %E !< %E  (tolerance)\n\n", err, tol);
    printf("dim:\t\t %u\n", dim);
    printf("nu:\t\t %.16lf\n", nu);
    for (unsigned int i = 0; i < dim; i++) {
        printf("%s\t\t [", (i == 0) ? "A:" : "  ");
        for (unsigned int j = 0; j < dim; j++) {
            printf("%.16lf%s", m[(dim * i) + j], (j + 1 < dim) ? ", " : "");
        }
        printf("]\n");
    }
    printf("x:\t\t[");
    for (unsigned int i = 0; i < dim; i++) {
        printf("%.16lf%s", x[i], (i + 1 < dim) ? ", " : "");
    }
    printf("]\ny:\t\t[");
    for (unsigned int i = 0; i < dim; i++) {
        printf("%.16lf%s", y[i], (i + 1 < dim) ? ", " : "");
    }
    printf("]\n");
}

/*!
 * @brief Guards against a hand-computed reciprocal basis being wrong. The test data
 * supplies both A and A^-T; this checks A^T B = I.
 */
static int anisoCheckReciprocal(unsigned int dim, const double *a,
                                const double *aInvT) {
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 0; j < dim; j++) {
            double s = 0.;
            for (unsigned int k = 0; k < dim; k++) {
                s += a[(dim * k) + i] * aInvT[(dim * k) + j];
            }
            if (fabs(s - ((i == j) ? 1. : 0.)) > 1e-12) {
                printf("\n\t BUG IN TEST DATA: A^T B != I at (%u,%u), got %.3e\n", i,
                       j, s);
                return 1;
            }
        }
    }
    return 0;
}

/*!
 * @brief Sweeps one lattice. x is given in lattice coordinates (x = A xFrac) and y
 * in reciprocal lattice coordinates (y = A^-T yFrac), so that integer entries mean
 * "on the lattice" independently of A.
 */
static void anisoReductionSweep(unsigned int dim, const double *a, // NOLINT
                                const double *aInvT, const double *xFrac, int numX,
                                const double *yFrac, int numY, const double *nus,
                                int numNu, bool reg, double tol, anisoStats *st) {
    unsigned int alpha[3] = {0, 0, 0};
    double x[3];
    double y[3];

    for (int itx = 0; itx < numX; itx++) {
        anisoMatVec(dim, a, xFrac + (dim * itx), x); // NOLINT
        for (int ity = 0; ity < numY; ity++) {
            anisoMatVec(dim, aInvT, yFrac + (dim * ity), y); // NOLINT
            for (int in = 0; in < numNu; in++) {
                double nu = nus[in];
                double complex base;
                double complex aniso;
                if (reg) {
                    base = epsteinZetaReg(nu, dim, a, x, y);
                    aniso = epsteinZetaAnisoReg(nu, dim, a, x, y, alpha);
                } else {
                    base = epsteinZeta(nu, dim, a, x, y);
                    aniso = epsteinZetaAniso(nu, dim, a, x, y, alpha);
                }
                double err = anisoErr(base, aniso);

                st->errMin = (st->errMin < err) ? st->errMin : err;
                st->errMax = (st->errMax > err) ? st->errMax : err;
                if (isfinite(err)) {
                    st->errSum += err;
                }
                st->total++;
                if (err < tol) {
                    st->passed++;
                } else if (st->reported < ANISO_MAX_REPORTS) {
                    reportAnisoError(dim, a, nu, x, y, base, aniso, err, tol, reg);
                    st->reported++;
                }
            }
        }
    }
}

/*!
 * @brief Checks that the anisotropic front ends reduce to the ordinary ones at
 * alpha = 0, over a grid chosen to hit every special case in epsteinZetaInternal:
 * nu a non-positive even integer, nu = 0, nu = d, nu = d + 2k, the nu > 10
 * largeExp branch and its threshold, x = 0, x a nonzero lattice vector, x on the
 * cell boundary, x outside the elementary cell, y = 0, y a nonzero reciprocal
 * lattice vector, y near zero, and lattices that are identity, diagonal, sheared,
 * generic, negative determinant and non-unit volume.
 * @param[in] reg: false compares epsteinZeta / epsteinZetaAniso, true compares
 * epsteinZetaReg / epsteinZetaAnisoReg.
 * @return number of failed tests.
 */
static int anisoReductionAllLattices(bool reg) { // NOLINT
    double tol = 5 * pow(10, -13);
    anisoStats st = {0, 0, 0, NAN, NAN, 0.};
    int dataErrors = 0;

    /* ------------------------------- 2D ------------------------------- */
    double a2[6][4] = {{1., 0., 0., 1.},     // identity
                       {2., 0., 0., 1.},     // diagonal, vol = 2
                       {1., 0.5, 0., 1.},    // unimodular shear, non-diagonal
                       {1.5, 0.2, 0.25, 1.}, // generic, vol = 29/20
                       {0., 1., 1., 0.},     // det < 0, forces pivoting
                       {1., -0.5, 0., sqrt(3.) / 2.}}; // hexagonal, vol = sqrt(3)/2
    double b2[6][4] = {
        {1., 0., 0., 1.},   {0.5, 0., 0., 1.},
        {1., 0., -0.5, 1.}, {20. / 29., -5. / 29., -4. / 29., 30. / 29.},
        {0., 1., 1., 0.},   {1., 0., 1. / sqrt(3.), 2. / sqrt(3.)}};

    double xFrac2[] = {0.,  0.,    // origin
                       1.,  0.,    // nonzero lattice vector
                       2.,  -1.,   // further lattice vector
                       0.3, -0.2,  // generic, inside the cell
                       2.3, -1.2,  // generic, outside the cell
                       0.5, 0.5};  // cell boundary
    double yFrac2[] = {0.,   0.,   // y = 0
                       1.,   0.,   // nonzero reciprocal lattice vector
                       0.3,  0.4,  // generic, inside the cell
                       1.3,  -0.6, // generic, outside the cell
                       0.5,  0.,   // cell boundary
                       1e-8, 0.};  // near zero
    double nu2[] = {-8.5, -8., -7.5, -7.,       -6.25,      -6.,   -5.,        -4.,
                    -3.5, -3., -2.,  -1.,       -0.5,       -0.25, 0.,         0.25,
                    0.5,  1.,  1.5,  2. - 1e-8, 2. - 1e-11, 2.,    2. + 1e-11, 2.5,
                    3.,   3.5, 4.,   4.5,       5.,         6.,    7.,         8.,
                    9.,   9.5, 9.9,  10.,       10.1,       10.5,  11.,        12.};

    for (int im = 0; im < 6; im++) {
        dataErrors += anisoCheckReciprocal(2, a2[im], b2[im]);
        anisoReductionSweep(2, a2[im], b2[im], xFrac2, 6, yFrac2, 6, nu2,
                            (int)(sizeof(nu2) / sizeof(nu2[0])), reg, tol, &st);
    }

    /* ------------------------------- 3D ------------------------------- */
    double a3[4][9] = {{1., 0., 0., 0., 1., 0., 0., 0., 1.},  // identity
                       {2., 0., 0., 0., 1., 0., 0., 0., 0.5}, // diagonal, vol = 1
                       {1., 0.5, 0., 0., 1., 0., 0., 0., 1.}, // shear, non-diagonal
                       {0., 0.5, 0.5, 0.5, 0., 0.5, 0.5, 0.5, 0.}}; // fcc, vol = 1/4
    double b3[4][9] = {{1., 0., 0., 0., 1., 0., 0., 0., 1.},
                       {0.5, 0., 0., 0., 1., 0., 0., 0., 2.},
                       {1., 0., 0., -0.5, 1., 0., 0., 0., 1.},
                       {-1., 1., 1., 1., -1., 1., 1., 1., -1.}};

    double xFrac3[] = {0.,  0.,   0.,    // origin
                       1.,  0.,   0.,    // nonzero lattice vector
                       2.,  -1.,  1.,    // further lattice vector
                       0.3, -0.2, 0.1,   // generic, inside the cell
                       0.5, 0.5,  0.5};  // cell corner
    double yFrac3[] = {0.,   0.,   0.,   // y = 0
                       1.,   0.,   0.,   // nonzero reciprocal lattice vector
                       0.3,  0.4,  -0.2, // generic, inside the cell
                       1.3,  -0.6, 0.2,  // generic, outside the cell
                       1e-8, 0.,   0.};  // near zero
    double nu3[] = {-7.5, -6., -4., -3., -2.,        -1., -0.5,       0.,
                    0.5,  1.,  2.,  2.5, 3. - 1e-11, 3.,  3. + 1e-11, 4.,
                    5.,   6.,  7.,  8.,  9.5,        10., 10.5,       12.};

    for (int im = 0; im < 4; im++) {
        dataErrors += anisoCheckReciprocal(3, a3[im], b3[im]);
        anisoReductionSweep(3, a3[im], b3[im], xFrac3, 5, yFrac3, 5, nu3,
                            (int)(sizeof(nu3) / sizeof(nu3[0])), reg, tol, &st);
    }

    if (st.reported >= ANISO_MAX_REPORTS) {
        printf("\n\t ... further failures suppressed.\n");
    }
    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", st.passed, st.total, tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", st.errMin, st.errMax,
           st.errSum / st.total);
    printf("\n");

    return (st.total - st.passed) + dataErrors;
}

/*!
 * @brief Checks epsteinZetaAniso reduces to epsteinZeta at alpha = 0.
 * @return number of failed tests.
 */
static int test_epsteinZeta_epsteinZetaAniso_reduction() {
    printf("%s ", __func__);
    return anisoReductionAllLattices(false);
}

/*!
 * @brief Checks epsteinZetaAnisoReg reduces to epsteinZetaReg at alpha = 0.
 * @return number of failed tests.
 */
static int test_epsteinZetaReg_epsteinZetaAnisoReg_reduction() {
    printf("%s ", __func__);
    return anisoReductionAllLattices(true);
}

/*!
 * @brief Test function for Epstein zeta function.
 *
 * This function tests the epsteinZeta function by comparing its outputs with
 * reference values read from data files. It performs multiple test cases and
 * reports the number of passed tests.
 *
 * @return number of failed tests.
 */
static int test_epsteinZeta() { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/epsteinZeta_ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *zetaRefData = fopen(path, "r");
    if (zetaRefData == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    int dim = 2;
    double *a = malloc((int)(dim * dim) * sizeof(double)); // NOLINT
    double *nu = malloc(2 * sizeof(double));
    double *x = malloc(dim * sizeof(double));
    double *y = malloc(dim * sizeof(double));
    double *zetaRef = malloc(2 * sizeof(double));
    double complex zetaComputed;
    double complex zetaExpected;
    double nuReal;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double tol = 6 * pow(10, -15);
    int testsPassed = 0;
    int totalTests = 0;

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    printf("\n\t ... ");
    printf("processing %s ", path);

    int scanResult;
    char line[256];
    while (fgets(line, sizeof(line), zetaRefData) != NULL) {
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
        zetaComputed = epsteinZeta(nuReal, dim, a, x, y);
        zetaExpected = zetaRef[0] + zetaRef[1] * I;
        errorAbs = errAbs(zetaComputed, zetaExpected);
        errorRel = errRel(zetaComputed, zetaExpected);
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
                   4, creal(zetaComputed), cimag(zetaComputed), creal(zetaExpected),
                   cimag(zetaExpected));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printMatrixUnitTest("a:", a, dim);
            printf("nu:\t\t %.16lf + %.16lf I\n", nu[0], nu[1]);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
        }
        totalTests++;
    }

    if (fclose(zetaRefData) != 0) {
        freeTestResources(a, nu, x, y, zetaRef);
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }

    freeTestResources(a, nu, x, y, zetaRef);

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
 * @brief Test function for Regularized Epstein zeta function.
 *
 * This function tests the epsteinZetaReg function by comparing its outputs with
 * reference values read from data files. It performs multiple test cases and
 * reports the number of passed tests.
 *
 * @return number of failed tests.
 */
static int test_epsteinZetaReg() { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/epsteinZetaReg_ref.csv", // NOLINT
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
    double complex zetaComputed;
    double complex zetaExpected;
    double nuReal;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double tol = 5 * pow(10, -15);
    int testsPassed = 0;
    int totalTests = 0;

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    printf("\n\t ... ");
    printf("processing file: %s ", path);

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
        zetaComputed = epsteinZetaReg(nuReal, dim, a, x, y);
        zetaExpected = zetaRef[0] + zetaRef[1] * I;
        errorAbs = errAbs(zetaComputed, zetaExpected);
        errorRel = errRel(zetaComputed, zetaExpected);
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
                   4, creal(zetaComputed), cimag(zetaComputed), creal(zetaExpected),
                   cimag(zetaExpected));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printMatrixUnitTest("a:", a, dim);
            printf("nu:\t\t %.16lf + %.16lf I\n", nu[0], nu[1]);
            printVectorUnitTest("x:\t\t", x, dim);
            printVectorUnitTest("y:\t\t", y, dim);
        }
        totalTests++;
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
           errSum / totalTests);
    printf("\n");

    return totalTests - testsPassed;
}

/*!
 * @brief Tests Epstein zeta function for random matrices in multiple dimensions
 * by comparing against Mathematica reference values.
 *
 * @return number of failed tests.
 */
static int test_epsteinZeta_random_matrices(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/epsteinZeta_random_matrices_ref.csv", BASE_PATH);
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
        double complex num = epsteinZeta(nu, dim, a, x, y);
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
            printf("epsteinZeta: ");
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

    return totalTests - testsPassed;
}

/*!
 * @brief Tests Epstein zeta function for random matrices in multiple dimensions
 * by comparing against Mathematica reference values.
 *
 * @return number of failed tests.
 */
static int test_epsteinZetaReg_random_matrices(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/epsteinZetaReg_random_matrices_ref.csv", BASE_PATH);
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
        double complex num = epsteinZetaReg(nu, dim, a, x, y);
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
            printf("epstienZetaReg: ");
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

    return totalTests - testsPassed;
}

/*!
 * @brief Test if the Epstein zeta function can be represented in terms of the
 * regularized function and the singularity, particularly focusing on cases where nu
 * = dim + 2k, where k is an integer.
 *
 * @return number of failed tests.
 * tests pass.
 */
static int test_epsteinZeta_epsteinZetaReg_reduction() { // NOLINT
    printf("%s ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex valZeta;
    double complex valZetaReg;
    double nu;

    int testsPassed = 0;
    int totalTests = 0;

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    double tol = pow(10, -13);
    unsigned int dim = 2;
    double m[] = {3. / 2, 1. / 5, 1. / 4, 1.};
    double x[] = {0.1, 0.2};
    double y[] = {0, 0.5};
    double vol = 29. / 20;

    for (int i = 0; i < 100; i++) {
        nu = -8.5 + (double)i / 5.;

        valZeta = epsteinZeta(nu, dim, m, x, y);
        valZetaReg = cexp(-2 * M_PI * I * dot(dim, x, y)) *
                     (epsteinZetaReg(nu, dim, m, x, y) +
                      singularity_s(nu, dim, dot(dim, y, y)) / vol);

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

    for (int i = 0; i < 100; i++) {
        nu = -8.5 + (double)i / 5.;

        valZeta = epsteinZeta(nu, dim, m, x, yZeta);
        valZetaReg = cexp(-2 * M_PI * I * dot(dim, x, yZeta)) *
                     (epsteinZetaReg(nu, dim, m, x, yZetaReg) +
                      singularity_s(nu, dim, dot(dim, yZeta, yZeta)) / vol);

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

    // test identity for x a nonzero lattice vector (identity lattice, so the
    // projection x_t2 is exactly zero and the nu = 0 special branch is reached)
    double mLat[] = {1., 0., 0., 1.};
    double xLat[] = {1., 0.};
    double yLat[] = {0.3, 0.4};
    double volLat = 1.;
    for (int i = 0; i < 100; i++) {
        nu = -8.45 + (double)i / 5.;
        valZeta = epsteinZeta(nu, dim, mLat, xLat, yLat);
        valZetaReg = cexp(-2 * M_PI * I * dot(dim, xLat, yLat)) *
                     (epsteinZetaReg(nu, dim, mLat, xLat, yLat) +
                      singularity_s(nu, dim, dot(dim, yLat, yLat)) / volLat);
        errorAbs = errAbs(valZeta, valZetaReg);
        errorRel = errRel(valZeta, valZetaReg);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;
        errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
        errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
        errSum += errorMaxAbsRel;
        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            reportEpsteinZetaError(valZeta, valZetaReg, errorMaxAbsRel, tol, mLat,
                                   dim, nu, xLat, yLat);
        }
        totalTests++;
    }

    // nu = 0 exactly, x a nonzero lattice vector: Z_0 = -exp(-2 pi i x.y),
    // singularity_s_0 = 0, hence Z^reg_0 = -1 independently of x and y
    nu = 0.;
    valZeta = epsteinZeta(nu, dim, mLat, xLat, yLat);
    valZetaReg = cexp(-2 * M_PI * I * dot(dim, xLat, yLat)) *
                 (epsteinZetaReg(nu, dim, mLat, xLat, yLat) +
                  singularity_s(nu, dim, dot(dim, yLat, yLat)) / volLat);
    errorAbs = errAbs(valZeta, valZetaReg);
    errorRel = errRel(valZeta, valZetaReg);
    errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;
    errMin = (errMin < errorMaxAbsRel) ? errMin : errorMaxAbsRel;
    errMax = (errMax > errorMaxAbsRel) ? errMax : errorMaxAbsRel;
    errSum += errorMaxAbsRel;
    if (errorMaxAbsRel < tol) {
        testsPassed++;
    } else {
        reportEpsteinZetaError(valZeta, valZetaReg, errorMaxAbsRel, tol, mLat, dim,
                               nu, xLat, yLat);
    }
    totalTests++;

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
 * @brief Test the Epstein zeta function behavior around the cutoff point for nu <=
 * dim.
 *
 * This function tests the Epstein zeta function for four cases:
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
static int test_epsteinZeta_cutoff() {
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

    for (int i = 0; i < 100; i++) {
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
 * @brief Main function to run all Epstein zeta function tests.
 *
 * @return number of failed tests.
 */
int main() {
    int failed = 0;
    failed += run_timed_test(test_epsteinZeta);
    failed += run_timed_test(test_epsteinZetaReg);
    failed += run_timed_test(test_epsteinZeta_random_matrices);
    failed += run_timed_test(test_epsteinZetaReg_random_matrices);
    failed += run_timed_test(test_epsteinZeta_epsteinZetaReg_reduction);
    failed += run_timed_test(test_epsteinZeta_cutoff);
    failed += run_timed_test(test_epsteinZeta_epsteinZetaAniso_reduction);
    failed += run_timed_test(test_epsteinZetaReg_epsteinZetaAnisoReg_reduction);
    return failed;
}
