// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file test_crandall_der.c
 * @brief Compares the derivatives of 3D Crandall functions to high-precision
 * benchmark values.
 */

#include "../src/crandall.h"
#include "../src/tools.h"
#include "utils.h"
#include "wrappers.h"
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
 * @brief Benchmarks 3D polynomial Y derivatives by comparing to high-precision
 * values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_polynomial_y_der(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/polynomial_y_der_ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    unsigned int k;
    unsigned int alphaAbs;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double tol = pow(10, -15);

    unsigned int *kRef = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2, z3}, {alpha1, alpha2i, alpha3}, result
        scanResult =
            sscanf(line, "%u,%lf,%lf,%lf,%u,%u,%u,%lf", // NOLINT
                   kRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        k = kRef[0];

        alphaAbs = mult_abs(dim, alpha);

        num = polynomial_y_der(k, dim, z, alpha, alphaAbs, 1);
        ref = refRead[0] + 0 * I;

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
            printf("polynomial_y_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("k:\t\t %u\n", k);
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(kRef);
    free(z);
    free(alpha);
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
 * @brief Benchmarks 3D polynomial Y derivatives by comparing to high-precision
 * values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_polynomial_y_der_harmonic(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/polynomial_y_der_harmonic_ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    unsigned int k;
    unsigned int alphaAbs;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double tol = 5 * pow(10, -15);

    unsigned int *kRef = malloc(sizeof(unsigned int));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        // Scan: k, {z1, z2, z3}, {alpha1, alpha2i, alpha3}, result
        scanResult =
            sscanf(line, "%u,%lf,%lf,%lf,%u,%u,%u,%lf", // NOLINT
                   kRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        k = kRef[0];

        alphaAbs = mult_abs(dim, alpha);

        num = polynomial_y_der_harmonic_wrapper((int)k, dim, z, alpha, alphaAbs);
        ref = refRead[0] + 0 * I;

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
            printf("polynomial_y_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("k:\t\t %u\n", k);
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
        totalTests++;
    }

    free(kRef);
    free(z);
    free(alpha);
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
 * @brief Benchmarks 3D singularity s derivatives by comparing to high-precision
 * values over a range of random parameters.
 *
 * @return number of failed tests.
 * */
int test_singularity_s_der(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/singularity_s_der_ref.csv", // NOLINT
                 BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }

    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double nu;
    unsigned int alphaAbs;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double tol = pow(10, -8);

    double *nuRef = malloc(sizeof(double));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        totalTests++;
        // Scan: nu, {z1, z2, z3}, {alpha1, alpha2, alpha3}, result
        scanResult =
            sscanf(line, "%lf,%lf,%lf,%lf,%u,%u,%u,%lf", // NOLINT
                   nuRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        alphaAbs = mult_abs(dim, alpha);

        nu = nuRef[0];

        num = (double complex)singularity_s_der(nu, dim, z, alpha, alphaAbs);
        ref = refRead[0] + 0 * I;

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
            printf("singularity_s_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t\t  != "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\t\t", nu);
            printf("\n");
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
    }

    free(z);
    free(nuRef);
    free(alpha);
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
 * @brief Benchmarks 3D singularity s derivatives by comparing to high-precision
 * values over a range of random parameters.
 * Test values include the ones for the non-harmonic method as well as higher order
 * derivatives.
 *
 * @return number of failed tests.
 * */
int test_singularity_s_der_harmonic(void) {
    printf("%s ", __func__);

    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path),
                          "%s/singularity_s_der_harmonic_ref.csv", // NOLINT
                          BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }

    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    double nu;
    unsigned int alphaAbs;
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex num;
    double complex ref;
    int scanResult;
    char line[256];

    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    int testsPassed = 0;
    int totalTests = 0;
    int dim = 3;
    double tol = pow(10, -10);

    double *nuRef = malloc(sizeof(double));
    double *z = malloc(dim * sizeof(double));
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));
    double *refRead = malloc(sizeof(double));

    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        totalTests++;
        // Scan: nu, {z1, z2, z3}, {alpha1, alpha2, alpha3}, result
        scanResult =
            sscanf(line, "%lf,%lf,%lf,%lf,%u,%u,%u,%lf", // NOLINT
                   nuRef, z, z + 1, z + 2, alpha, alpha + 1, alpha + 2, refRead);

        if (scanResult != 8) {
            printf("Error reading line: %s\n", line);
            printf("Scanned %d values instead of 8\n", scanResult);
            continue;
        }

        alphaAbs = mult_abs(dim, alpha);

        nu = nuRef[0];

        num = (double complex)singularity_s_der_harmonic_wrapper(nu, dim, z, alpha,
                                                                 alphaAbs);
        ref = refRead[0] + 0 * I;

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
            printf("singularity_s_der: ");
            printf(" %0*.16lf %+.16lf I (this implementation) \n\t\t\t  != "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(num), cimag(num), creal(ref), cimag(ref));
            printf("Min(Emax, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\t\t", nu);
            printf("\n");
            printVectorUnitTest("z:\t\t", z, dim);
            printMultiindexUnitTest("alpha:\t\t", alpha, dim);
            printf("\n");
        }
    }

    free(z);
    free(nuRef);
    free(alpha);
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
 * @brief Tests crandall_gReg_harmonic for random inputs in dimensions 1 to 4
 * by comparing against Mathematica reference values.
 *
 * @return number of failed tests.
 */
int test_crandall_gReg_harmonic(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/crandall_gReg_harmonic_ref.csv", BASE_PATH);
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
    double tol = pow(10, -14);
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

        // Bounds check on dimension
        if (dim < 1 || dim > max_test_dim) {
            (void)fclose(data);
            return fprintf(stderr, "Invalid dimension: %u\n", dim);
        }

        // Parse nu
        double nu = strtod(ptr, &endptr);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        // Use stack allocation with fixed max size
        double z[max_test_dim];

        // Parse z (dim values)
        int parseOk = 1;
        for (unsigned int i = 0; i < dim && parseOk; i++) {
            z[i] = strtod(ptr, &endptr);
            if (*endptr != ',') {
                parseOk = 0;
            } else {
                ptr = endptr + 1;
            }
        }

        // Parse k
        int k = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        // Parse alpha (dim values) and accumulate n = |alpha|
        int n = 0;
        for (unsigned int i = 0; i < dim && parseOk; i++) {
            n += (int)lround(strtod(ptr, &endptr));
            if (*endptr != ',') {
                parseOk = 0;
            } else {
                ptr = endptr + 1;
            }
        }
        if (!parseOk) {
            continue;
        }

        // Parse reference result (real scalar)
        double ref = strtod(ptr, &endptr);

        // Compute result
        double num = crandall_gReg_harmonic(k, n, dim, nu, z, 1.0);

        // Compute errors
        double errorAbs = fabs(ref - num);
        double errorRel = (fabs(ref) > 0.) ? fabs((ref - num) / ref) : errorAbs;
        double errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        errMin =
            (isnan(errMin) || errorMaxAbsRel < errMin) ? errorMaxAbsRel : errMin;
        errMax =
            (isnan(errMax) || errorMaxAbsRel > errMax) ? errorMaxAbsRel : errMax;
        errSum += errorMaxAbsRel;

        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\n\n");
            printf("Warning! ");
            printf("crandall_gReg_harmonic: ");
            printf(" %0*.16lf (this implementation) \n\t\t!= "
                   "%.16lf (reference implementation)\n",
                   4, num, ref);
            printf("Min(Eabs, Erel):      %E !< %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("dim:\t\t %u\n", dim);
            printf("nu:\t\t %.16lf\n", nu);
            printf("k:\t\t %d\n", k);
            printf("n (|alpha|):\t %d\n", n);
            printVectorUnitTest("z:\t\t", z, dim);
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
 * @brief Main function to run all set crandall and related derivatives function
 * tests.
 *
 * @return number of failed tests.
 */
int main(void) {
    int failed = 0;
    failed += run_timed_test(test_crandall_gReg_harmonic);
    failed += run_timed_test(test_polynomial_y_der);

    // tests for functions not part of the EpsteinLib core:
    failed += run_timed_test(test_polynomial_y_der_harmonic);
    failed += run_timed_test(test_singularity_s_der);
    failed += run_timed_test(test_singularity_s_der_harmonic);

    return failed;
}
