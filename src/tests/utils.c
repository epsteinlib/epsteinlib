// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
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
 * @param[in] ref: reference value.
 * @param[in] comp: compparison value.
 * @return sqrt(rev - comp)
 */
double errAbs(double complex ref, double complex comp) {
    double complex diff = ref - comp;
    double norm = (creal(diff) * creal(diff)) + (cimag(diff) * cimag(diff));
    return sqrt(norm);
}

/*!
 * @brief relative difference between to complex numbers.
 * @param[in] ref: reference value.
 * @param[in] comp: compparison value.
 * @return errAbs / norm(ref)
 */
double errRel(double complex ref, double complex comp) {
    double absRef = sqrt((creal(ref) * creal(ref)) + (cimag(ref) * cimag(ref)));
    if (absRef < EPSILON_REF) {
        return errAbs(ref, comp);
    }
    return errAbs(ref, comp) / absRef;
}

/**
 * @brief prints vector to terminal, is used in unitTest.
 * @param[in] name: name of the vector as a string.
 * @param[in] vec: vector.
 * @param[in] dim: size of the vector.
 */
void printVectorUnitTest(const char *name, double *vec, unsigned int dim) {
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
 * @brief prints vector to terminal, is used in unitTest.
 * @param[in] name: name of the vector as a string.
 * @param[in] vec: vector.
 * @param[in] dim: size of the vector.
 */
void printMultiindexUnitTest(const char *name, unsigned int *vec, unsigned int dim) {
    printf("%s[", name);
    for (int i = 0; i < dim; ++i) {
        printf("%u", vec[i]);
        if (i != dim - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

/**
 * @brief prints square matrix to terminal, is used in unitTest.
 * @param[in] name: name of the matrix as a string.
 * @param[in] mat: matrix.
 * @param[in] dim: size of the matrix.
 */
void printMatrixUnitTest(const char *name, const double *mat, unsigned int dim) {
    printf("%s", name);
    double matrixEntry;
    for (int i = 0; i < dim; ++i) {
        printf("\t\t [");
        for (int j = 0; j < dim; ++j) {
            matrixEntry = mat[(i * dim) + j];
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
