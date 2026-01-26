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
 * @brief Opens a file.
 * @param path Path to the file.
 * @param mode 'r' to read or 'w' to write.
 * @return FILE* Pointer to the opened file.
 * @note Exits the program if the file cannot be opened.
 */
static FILE *open(char *path, char *mode) {
    FILE *file = fopen(path, mode);
    if (file) {
        return file;
    }
    printf("File '%s' does not exist.\n", path);
    exit(1); // NOLINT
}

/**
 * @brief Benchmarks the evaluation time of the setZetaDer function
 *
 * @return  0 on successful execution.
 */
int benchmark_setZetaDer_timing() { // NOLINT

    printf("\n========== Benchmarking %s ==========\n", __func__);

    char zetaDataString[MAX_PATH_LENGTH];
    unsigned int dimMax = 6;
    double nu = 0.5;
    double complex res;
    double elapsedTime;
    double elapsedTimeMax;
    int iterations = 10;
    clock_t timeStart;
    clock_t timeEnd;

    FILE *zetaData = NULL;
    double *elapsedTimes = NULL;
    double *a = NULL;
    double *x = NULL;
    double *y = NULL;
    unsigned int *alpha = NULL;

    for (unsigned int dim = 1; dim < dimMax + 1; dim++) {

        if (snprintf(zetaDataString, MAX_PATH_LENGTH, "%s/setZetaDer_timing_%uD.csv",
                     BASE_PATH, dim) >= MAX_PATH_LENGTH) {
            return fprintf(stderr, "Error: filename too long\n");
        }
        zetaData = open(zetaDataString, "w");
        if (zetaData == NULL) {
            printf("%s\n", strerror(errno)); // NOLINT
            return 1;
        }
        elapsedTimes = malloc(iterations * sizeof(double));
        a = malloc((unsigned long)dim * (unsigned long)dim * sizeof(double));
        x = malloc(dim * sizeof(double));
        y = malloc(dim * sizeof(double));
        alpha = malloc(dim * sizeof(unsigned int));

        for (int i = 0; i < dim; i++) {
            x[i] = 0;
            alpha[i] = 0;
            for (int j = 0; j < dim; j++) {
                a[(dim * i) + j] = (i == j) ? 1 : 0;
            }
        }

        printf("\n");
        printf("d: %u, t: ", dim);

        for (int j = 0; j < 10 + 1; j++) {

            if (dim > 4 && j > 5) {
                continue;
            }

            elapsedTimeMax = 0.;

            alpha[0] = j;

            for (int i = -5; i < 5 + 1; i++) {
                y[0] = (double)i * 0.1;
                for (int n = 0; n < iterations; n++) {
                    timeStart = clock();
                    res = setZetaDer(nu, dim, a, x, y, alpha);
                    timeEnd = clock();
                    elapsedTimes[n] =
                        ((double)(timeEnd - timeStart)) / CLOCKS_PER_SEC;
                }
                sort(elapsedTimes, iterations);
                elapsedTime = elapsedTimes[iterations / 2];
                elapsedTimeMax =
                    (elapsedTime > elapsedTimeMax) ? elapsedTime : elapsedTimeMax;
            }

            fprintf(zetaData, "%u,%.16lf\n", alpha[0], elapsedTimeMax); // NOLINT
            printf("%.2E ", elapsedTimeMax);
        }

        if (fclose(zetaData) != 0) {
            return fprintf(stderr, "Error closing file: %d\n", errno);
        }
        free(a);
        free(x);
        free(y);
        free(alpha);
        free(elapsedTimes);
    }

    printf("\n");

    (void)res;

    return 0;
}

/**
 * @brief Main function to run all benchmark tests.
 * @return number of failed executions.
 */
int main() {
    int failed = 0;
    failed += benchmark_setZetaDer_timing();
    return failed;
}
