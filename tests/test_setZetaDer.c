// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

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

    return totalTests - testsPassed;
}

/*!
 * @brief Benchmarks 2D setZetaDer function by comparing to reference values of the
 * laplacian of set zeta function obtained by finite differences.
 *
 * @return number of failed tests.
 * */
static int test_setZetaDer_odd(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/setZetaDer_odd_ref.csv", // NOLINT
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
    unsigned int dim = 3;
    double tol = pow(10, -16);

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
        // Scan: nu, {a11, a12, a13, a21, a22, a23, a21, a32, a33}, {x1, x2, x3},
        // {y1, y2, y3}, {alpha1, alpha2, alpha3}, {Re[result], Im[result]}
        scanResult = sscanf( // NOLINT
            line,
            "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%u,%u,%"
            "u,%lf,%lf",
            nuRef, a, a + 1, a + 2, a + 3, a + 4, a + 5, a + 6, a + 7, a + 8, x,
            x + 1, x + 2, y, y + 1, y + 2, alpha, alpha + 1, alpha + 2, refRead,
            refRead + 1);

        if (scanResult != 21) {
            printf("\n\t ");
            printf("Error reading line: %s", line);
            printf("\t ");
            printf("Scanned %d values instead of 17", scanResult);
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
    int failed = 0;

    failed += run_timed_test(test_setZetaDer_1D);
    failed += run_timed_test(test_setZetaDer_2D);
    failed += run_timed_test(test_setZetaDer_taylor);
    failed += run_timed_test(test_setZetaDer_odd);
    failed += run_timed_test(test_setZetaDer_special_exponents);
    failed += run_timed_test(test_setZetaDer_poly_laplace);

    return failed;
}
