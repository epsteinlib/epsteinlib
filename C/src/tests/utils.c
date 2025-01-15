// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include <complex.h>
#include <math.h>
#include <stdio.h>

#define EPSILON_REF (double)pow(10, -62)

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

#ifndef BASE_PATH
#define BASE_PATH "csv"
#endif

/*!
 * @brief absolute difference between to complex numbers.
 * @param ref: reference value.
 * @param comp: compparison value.
 * @return sqrt(rev - comp)
 */
double errAbs(double complex ref, double complex comp) {
    double complex diff = ref - comp;
    double norm = creal(diff) * creal(diff) + cimag(diff) * cimag(diff);
    return sqrt(norm);
}

/*!
 * @brief relative difference between to complex numbers.
 * @param ref: reference value.
 * @param comp: compparison value.
 * @return errAbs / norm(ref)
 */
double errRel(double complex ref, double complex comp) {
    double absRef = sqrt(creal(ref) * creal(ref) + cimag(ref) * cimag(ref));
    if (absRef < EPSILON_REF) {
        return errAbs(ref, comp);
    }
    return errAbs(ref, comp) / absRef;
}

/**
 * @brief prints vector to terminal, is used in unitTest.
 * @param name: name of the vector as a string.
 * @param vec: vector.
 * @param dim: size of the vector.
 */
void printVectorUnitTest(const char *name, double *vec, int dim) {
    printf("%s[", name);
    for (int i = 0; i < dim; ++i) {
        printf("%.16lf", vec[i]);
        if (i != dim - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

/**
 * @brief prints square matrix to terminal, is used in unitTest.
 * @param name: name of the matrix as a string.
 * @param mat: matrix.
 * @param dim: size of the matrix.
 */
void printMatrixUnitTest(const char *name, const double *mat, int dim) {
    printf("%s", name);
    double matrixEntry;
    for (int i = 0; i < dim; ++i) {
        printf("\t\t [");
        for (int j = 0; j < dim; ++j) {
            matrixEntry = mat[i * dim + j];
            if (pow(matrixEntry, 2) < 100) {
                printf("%.16f", matrixEntry);
            } else {
                printf("%.2e", matrixEntry * pow(10, 9));
            }
            if (j != dim - 1) {
                printf(", ");
            }
        }
        printf("]\n");
    }
}
