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

#include "epsteinZeta.h"

#define BASE_PATH "benchmark/csv"

#define MAX_PATH_LENGHT 100

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
 * @brief Sorts an array of doubles in ascending order using bubble sort.
 * @param arr: The array to be sorted.
 * @param size: The size of the array.
 */
void sort(double *arr, int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (*(arr + j) > *(arr + j + 1)) {
                double temp = *(arr + j);
                *(arr + j) = *(arr + j + 1);
                *(arr + j + 1) = temp;
            }
        }
    }
}

/**
 * @brief Benchmarks the Epstein Zeta and regularized Epstein Zeta functions.
 * @param dim Dimension of the lattice.
 * @param a Array representing the lattice matrix.
 * @param x Array of x coordinates.
 * @param y Array of y coordinates.
 * @param zetaDataString Filename for storing Epstein Zeta function data.
 * @param zetaRegDataString Filename for storing regularized Epstein Zeta function
 * data.
 */
void benchmark(int dim, double a[], double x[], double y[], char zetaDataString[],
               char zetaRegDataString[]) {
    FILE *zetaData = open(zetaDataString, "w");
    FILE *zetaRegData = open(zetaRegDataString, "w");
    if (zetaData == NULL)
        printf("%s\n", strerror(errno));
    double nu;
    double complex zetaReg;
    double varBenchmark;
    double offset = ldexp(1, -15);
    double elapsedTime;
    int iterations = 25;
    int iterationsReg = 25;
    double *elapsedTimes = malloc(iterations * sizeof(double));
    clock_t timeStart, timeEnd;
    for (int i = -250; i < 250 + 1; i++) { // print values to file
        varBenchmark = i * 0.05 + offset;
        nu = varBenchmark;
        // zeta
        elapsedTime = 0;
        for (int n = 0; n < iterations; n++) {
            timeStart = clock();
            zetaReg = epsteinZeta(nu, dim, a, x, y);
            timeEnd = clock();
            elapsedTimes[n] = ((double)(timeEnd - timeStart)) / CLOCKS_PER_SEC;
        }
        sort(elapsedTimes, iterations);
        elapsedTime = elapsedTimes[(int)iterations / 2];
        fprintf(zetaData, "%.16lf,%.16lf,%.16lf,%.16lf\n", varBenchmark,
                creal(zetaReg), cimag(zetaReg), elapsedTime);
        printf("nu:\t %.16lf\t", varBenchmark);
        printf("zeta:\t\t%.16lf %+.16lf I, \t execution time: %.8f seconds\n",
               creal(zetaReg), cimag(zetaReg), elapsedTime);
        //       // zeta Reg
        elapsedTime = 0;
        for (int n = 0; n < iterationsReg; n++) {
            timeStart = clock();
            zetaReg = epsteinZetaReg(nu, dim, a, x, y);
            timeEnd = clock();
            elapsedTimes[n] = ((double)(timeEnd - timeStart)) / CLOCKS_PER_SEC;
        }
        sort(elapsedTimes, iterations);
        elapsedTime = elapsedTimes[(int)(iterations / 2)];
        fprintf(zetaRegData, "%.16lf,%.16lf,%.16lf,%.16lf\n", varBenchmark,
                creal(zetaReg), cimag(zetaReg), elapsedTime);
        printf("nu:\t %.16lf\t", varBenchmark);
        printf("zetaReg:\t%.16lf %+.16lf I, \t execution time: %.8f "
               "seconds\n",
               creal(zetaReg), cimag(zetaReg), elapsedTime);
    }
    free(elapsedTimes);
    fclose(zetaData);
    fclose(zetaRegData);
}

/*!
 * @brief Example of evaluating the Epstein Zeta function for values where
 * analytic representations exist for benchmarking.
 */
void s1() {
    int dim = 1;
    double a[] = {1};
    double x[] = {-0.5};
    double y[] = {0};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking %s() ==========\n", __func__);
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/*!
 * @brief Example of evaluating the Epstein Zeta function for values where
 * analytic representations exist for benchmarking.
 */
void s21() {
    int dim = 2;
    double a[] = {1, 0, 0, 2};
    double x[] = {-1, -2};
    double y[] = {0, 0};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking diag12_m1m2_00() ==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Benchmarks the Epstein Zeta function for a 2D lattice with specific
 * parameters.
 */
void s22() {
    int dim = 2;
    double a[] = {1, 1. / 2, 0, sqrt(3.) / 2};
    double x[] = {0.0, 0.0};
    double y[] = {-0, 0};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking onehalf0sqrt3haf_00_00() ==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Benchmarks the Epstein Zeta function for a 3D lattice with specific
 * parameters.
 */
void s31() {
    int dim = 3;
    double a[] = {1, 0, 0, 0, 1, 0, 0, 0, 2};
    double x[] = {0, 0, -0.5};
    double y[] = {0.5, 0, 0};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking diag112_00mhalf_half00() ==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Benchmarks the Epstein Zeta function for another 3D lattice configuration.
 */
void s32() {
    int dim = 3;
    double a[] = {6, 0, 0, 0, 6, 0, 0, 0, 6};
    double x[] = {-1, -1, -1};
    double y[] = {1. / 12, 1. / 12, 1. / 12};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking diag666_m1m1m1_twelthtwelthtwelth() "
           "==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Benchmarks the Epstein Zeta function for a third 3D lattice configuration.
 */
void s33() {
    int dim = 3;
    double a[] = {2 * sqrt(2.), 0, 0, 0, 4, 0, 0, 0, 2};
    double x[] = {0, -1, -1};
    double y[] = {1. / (4 * sqrt(2.)), 0, 0};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking diag2sqrt242_0m1m1_4sqrt2th00() ==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Benchmarks the Epstein Zeta function for a 4D lattice with specific
 * parameters.
 */
void s4() {
    int dim = 4;
    double a[] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    double x[] = {0.5, 0.0, 0.0, 0.0};
    double y[] = {0.0, 0.0, 0.0, 0.0};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking Id_half000_0000() ==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Benchmarks the Epstein Zeta function for a 6D lattice with specific
 * parameters.
 */
void s6() {
    int dim = 6;
    double a[dim * dim];
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            a[dim * i + j] = (i == j) ? 1 : 0;
        }
    }
    double x[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double y[] = {0.5, 0.5, 0.0, 0.0, 0.0, 0.0};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking Id_000000_halfhalf0000() ==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Benchmarks the Epstein Zeta function for an 8D lattice with specific
 * parameters.
 */
void s8() {
    int dim = 8;
    double a[dim * dim];
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            a[dim * i + j] = (i == j) ? 1 : 0;
        }
    }
    double x[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double y[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    char zetaDataString[MAX_PATH_LENGHT];
    char zetaRegDataString[MAX_PATH_LENGHT];
    sprintf(zetaDataString, "%s/epsteinZeta_%s.csv", BASE_PATH, __func__);
    sprintf(zetaRegDataString, "%s/epsteinZetaReg_%s.csv", BASE_PATH, __func__);
    printf("\n========== Benchmarking "
           "Id_00000000_halfhalfhalfhalfhalfhalfhalfhalf() ==========\n");
    benchmark(dim, a, x, y, zetaDataString, zetaRegDataString);
}

/**
 * @brief Main function to run all benchmark tests.
 * @return int Returns 0 on successful execution.
 */
int main() {
    s1();
    s21();
    s22();
    s31();
    s32();
    s33();
    s4();
    s6();
    s8();
    return 0;
}
