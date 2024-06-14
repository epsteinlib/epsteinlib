// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file tools.c
 * @brief  Minimal linear algebra for matrix vector operations.
 */

#include <math.h>
#include <stdbool.h>

/*!
 * @brief minimal distance of two vector elements considered unequal.
 */
#define EPS ldexp(1, -32)

/**
 * @brief euclidean dot product.
 * @param dim: dimension of the input vectors
 * @param v1: first vector.
 * @param v2: second vector.
 * @return dot product of v1 and v2.
 */
double dot(int dim, double *v1, double *v2) {
    double r = 0;
    for (int i = 0; i < dim; i++)
        r += v1[i] * v2[i];
    return r;
}

/**
 * @brief matrix - (integer) vector multiplication.
 * @param dim: dimension of the square matrix and the integer vector.
 * @param m: square matrix.
 * @param v: integer vector.
 * @param res: solution vector of the vector matrix multiplication.
 */
void matrix_intVector(int dim, double *m, int *v, double *res) {
    for (int i = 0; i < dim; i++) {
        res[i] = 0;
        for (int j = 0; j < dim; j++)
            res[i] += m[i * dim + j] * v[j];
    }
    return;
}

/**
 * @brief square matrix transpose.
 * @param dim: dimension of the square matrix.
 * @param m: square matrix.
 */
void transpose(int dim, double *m) {
    double swap;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < i; j++) {
            swap = m[dim * i + j];
            m[dim * i + j] = m[dim * j + i];
            m[dim * j + i] = swap;
        }
    }
}

/**
 * @brief check if two vectors are equal.
 * @param dim: dimension of the vectors.
 * @param v1: first vector.
 * @param v2: second vector.
 * @return true if the vectors are equal, false if the vectors are not equal.
 */
bool equals(int dim, double *v1, double *v2) {
    bool eq = true;
    for (int i = 0; i < dim && eq; i++)
        eq = eq && fabs(v1[i] - v2[i]) < EPS;
    return eq;
}

/**
 * @brief check if vector is zero.
 * @param dim: dimension of the vectors.
 * @param v: vector.
 * @return true if the vector is zero.
 */
bool equalsZero(int dim, double *v) {
    bool eq = true;
    for (int i = 0; i < dim && eq; i++)
        eq = eq && fabs(v[i]) < EPS;
    return eq;
}
#undef EPS
