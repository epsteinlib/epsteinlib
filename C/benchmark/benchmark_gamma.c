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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "gamma.h"

#define BASE_PATH "benchmark/csv"

#define MAX_PATH_LENGTH 100

#define STEPS 100

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
    exit(1);
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
void benchmark_gamma(double xinc, double xbound, const char *filename) {
    char gammaString[MAX_PATH_LENGTH];
    snprintf(gammaString, MAX_PATH_LENGTH, "%s/%s.csv", BASE_PATH, filename);
    FILE *gammaData = open(gammaString, "w");
    if (gammaData == NULL) {
        printf("Error opening file: %s\n", strerror(errno));
        return;
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
            fprintf(gammaData, "%.16lf,%.16lf,%.16lf,%.16lf\n", nu, x,
                    creal(upper_gamma_val), cimag(upper_gamma_val));
            printf("nu: %.16lf, x: %.16lf, upper gamma: %.16lf + %.16lfi\n", nu, x,
                   creal(upper_gamma_val), cimag(upper_gamma_val));
        }
    }
    fclose(gammaData);
}

/**
 * @brief Benchmarks the upper gamma function for larger x values.
 *
 * This function calls benchmark_gamma with parameters suitable for
 * testing the upper gamma function with larger x values (up to 20.1).
 */
void gamma_big() {
    double xinc = ldexp(1, -4);
    double xbound = 20.1;
    benchmark_gamma(xinc, xbound, __func__);
}

/**
 * @brief Benchmarks the upper gamma function for smaller x values.
 *
 * This function calls benchmark_gamma with parameters suitable for
 * testing the upper gamma function with smaller x values (up to 2.01).
 */
void gamma_small() {
    double xinc = ldexp(1, -7);
    double xbound = 2.01;
    benchmark_gamma(xinc, xbound, __func__);
}

/**
 * @brief Main function to run all benchmark tests.
 * @return int Returns 0 on successful execution.
 */
int main() {
    gamma_big();
    gamma_small();
    return 0;
}
