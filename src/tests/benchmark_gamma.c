// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file benchmark.c
 * @author Jonathan Busse
 * @date 06/06/2024
 * @section Description: Generate epsteinZeta(Reg) reference values
 */

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../gamma.h"

#define BASE_PATH "benchmark/csv"

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

/**
 * @brief Opens a file.
 * @param path Path to the file.
 * @param mode 'r' to read or 'w' to write.
 * @return FILE* Pointer to the opened file.
 * @note Exits the program if the file cannot be opened.
 */
FILE *open(char *path, char *mode) {
    FILE *file = fopen(path, mode);
    if (file) {
        return file;
    }
    printf("File '%s' does not exist.\n", path);
    exit(1); // NOLINT
}

/**
 * @brief Computes and prints reference values of the upper gamma function to a CSV
 * file.
 * @param numin Minimum value for nu.
 * @param numax Maximum value for nu.
 * @param nustep Step size for nu.
 * @param xmin Minimum value for x.
 * @param xmax Maximum value for x.
 * @param xstep Step size for x.
 * @param filename Name of the output CSV file (without extension).
 */
int benchmark_gamma(double xinc, double xbound, const char *filename) {
    char gammaString[MAX_PATH_LENGTH];

    if (snprintf(gammaString, MAX_PATH_LENGTH, "%s/%s.csv", BASE_PATH, filename) >=
        MAX_PATH_LENGTH) {
        return fprintf(stderr, "Error: filename too long\n");
    }

    FILE *gammaData = open(gammaString, "w");
    if (gammaData == NULL) {
        printf("%s\n", strerror(errno)); // NOLINT
        return 1;
    }

    double nu = 0;
    double numin = -12.5;
    double nuinc = ldexp(1, -4);
    double x = 0;
    double xmin = ldexp(1, -12);
    double complex upper_gamma_val;
    for (int i = 0; i < 400 + 1; i++) {
        for (int j = 0; j < 400 + 1; j++) {
            nu = numin + i * nuinc;
            x = xmin + j * xinc;
            if (x > xbound) {
                break;
            }
            upper_gamma_val = egf_ugamma(nu, x);
            fprintf(gammaData, "%.16lf,%.16lf,%.16lf,%.16lf\n", nu, x, // NOLINT
                    creal(upper_gamma_val), cimag(upper_gamma_val));
            printf("nu: %.16lf, x: %.16lf, upper gamma: %.16lf + %.16lfi\n", nu, x,
                   creal(upper_gamma_val), cimag(upper_gamma_val));
        }
    }
    if (fclose(gammaData) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    return 0;
}

/**
 * @brief Calls benchmark_gamma with parameters suitable for
 * testing the upper gamma function with larger x values (up to 2.01).
 *
 * @return  0 on successful execution.
 */
int gamma_big() {
    // Parameters from "Computation and Properties of the Epstein Zeta Function"
    // Using larger stepsize to reduce evaluation time and file size
    //    double xinc = ldexp(1, -4);
    double xinc = 20. * ldexp(1, -4);
    double xbound = 20.1;
    return benchmark_gamma(xinc, xbound, __func__);
}

/**
 * @brief Calls benchmark_gamma with parameters suitable for
 * testing the upper gamma function with smaller x values (up to 2.01).
 *
 * @return  0 on successful execution.
 */
int gamma_small() {
    // Parameters from "Computation and Properties of the Epstein Zeta Function"
    // Using larger stepsize to reduce evaluation time and file size
    //    double xinc = ldexp(1, -7);
    double xinc = 20. * ldexp(1, -7);
    double xbound = 2.01;
    return benchmark_gamma(xinc, xbound, __func__);
}

/**
 * @brief Main function to run all benchmark tests.
 * @return number of failed executions.
 */
int main() {
    int failed1 = gamma_big();
    int failed2 = gamma_small();
    return failed1 + failed2;
}
