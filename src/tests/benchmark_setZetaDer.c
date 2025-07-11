// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file benchmark.c
 * @author Jonathan Busse
 * @date 06/06/2024
 * @section Description: derivative reference values
 */

#include "utils.h"
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "epsteinZeta.h"

#define BASE_PATH "src/tests/csv"

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

/**
 * @brief Benchmarks the Epstein Zeta and regularized Epstein Zeta functions.
 * @param dim Dimension of the lattice.
 * @param a Array representing the lattice matrix.
 * @param x Array of x coordinates.
 * @param y Array of y coordinates.
 * @param zetaDataString Filename for storing Epstein Zeta function data.
 * @param resDataString Filename for storing regularized Epstein Zeta function
 * data.
 * @return  0 on successful execution.
 */
int benchmark(unsigned int dim, double a[], double x[], double y[],
              unsigned int alpha[], char zetaDataString[]) {
    FILE *zetaData = open(zetaDataString, "w");
    if (zetaData == NULL) {
        printf("%s\n", strerror(errno)); // NOLINT
        return 1;
    }
    double nu;
    double complex res;
    double varBenchmark;
    double offset = ldexp(1, -15);
    double elapsedTime;
    int iterations = 25;
    double *elapsedTimes = malloc(iterations * sizeof(double));
    clock_t timeStart;
    clock_t timeEnd;
    for (int i = -250; i < 250 + 1; i++) { // print values to file
        varBenchmark = i * 0.05 + offset;
        nu = varBenchmark;
        // zeta
        for (int n = 0; n < iterations; n++) {
            timeStart = clock();
            res = setZetaDer(nu, dim, a, x, y, alpha);
            timeEnd = clock();
            elapsedTimes[n] = ((double)(timeEnd - timeStart)) / CLOCKS_PER_SEC;
        }
        sort(elapsedTimes, iterations);
        elapsedTime = elapsedTimes[iterations / 2];
        fprintf(zetaData, "%.16lf,%.16lf,%.16lf,%.16lf\n", varBenchmark, // NOLINT
                creal(res), cimag(res), elapsedTime);
        printf("nu:\t %.16lf\t", varBenchmark);
        printf("zeta:\t\t%.16lf %+.16lf I, \t execution time: %.8f seconds\n",
               creal(res), cimag(res), elapsedTime);
    }
    free(elapsedTimes);
    if (fclose(zetaData) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    return 0;
}

/*!
 * @brief Example of evaluating the set zeta derivatives for values where
 * analytic representations exist for benchmarking.
 * @return  0 on successful execution.
 */
int s1() {
    unsigned int dim = 1;
    double a[] = {2. / 3.};
    double x[] = {0};
    double y[] = {0};
    unsigned int alpha[] = {4};
    char zetaDataString[MAX_PATH_LENGTH];

    if (snprintf(zetaDataString, MAX_PATH_LENGTH, "%s/setZetaDer_%s.csv", BASE_PATH,
                 __func__) >= MAX_PATH_LENGTH) {
        return fprintf(stderr, "Error: filename too long\n");
    }

    printf("\n========== Benchmarking %s() ==========\n", __func__);
    return benchmark(dim, a, x, y, alpha, zetaDataString);
}

/*!
 * @brief Example of evaluating the set zeta derivatives for values where
 * analytic representations exist for benchmarking.
 * @return  0 on successful execution.
 */
int s2() {
    unsigned int dim = 2;
    double a[] = {1., 0.5, 0., sqrt(3.) / 2.};
    double x[] = {0., 0.};
    double y[] = {0., 0.};
    unsigned int alpha[] = {2, 0};
    char zetaDataString[MAX_PATH_LENGTH];

    if (snprintf(zetaDataString, MAX_PATH_LENGTH, "%s/setZetaDer_%s.csv", BASE_PATH,
                 __func__) >= MAX_PATH_LENGTH) {
        return fprintf(stderr, "Error: filename too long\n");
    }

    printf("\n========== Benchmarking %s() ==========\n", __func__);
    return benchmark(dim, a, x, y, alpha, zetaDataString);
}

/**
 * @brief Main function to run all benchmark tests.
 * @return number of failed executions.
 */
int main() {
    int failed1 = s1();
    int failed2 = s2();
    return failed1 + failed2;
}
