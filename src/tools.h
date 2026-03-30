// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024-2026 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file tools.h
 * @brief  Minimal linear algebra for matrix vector operations.
 */

#ifndef EPSTEIN_TOOLS
#define EPSTEIN_TOOLS
#include <complex.h>
#include <stdbool.h>

/**
 * @brief Kahan compensated summation step.
 * @param[in,out] sum: running compensated sum.
 * @param[in,out] epsilon: running compensation term.
 * @param[in] x: summand to add.
 */
static inline void kahan_add(double complex *restrict sum,
                             double complex *restrict epsilon, double complex x) {
    double complex y = x - *epsilon;
    double complex t = *sum + y;
    *epsilon = (t - *sum) - y;
    *sum = t;
}

/**
 * @brief euclidean dot product.
 * @param[in] dim: dimension of the input vectors
 * @param[in] v1: first vector.
 * @param[in] v2: second vector.
 * @return dot product of v1 and v2.
 */
double dot(unsigned int dim, const double *v1, const double *v2);

/**
 * @brief matrix - (integer) vector multiplication.
 * @param[in] dim: dimension of the square matrix and the integer vector.
 * @param[in] m: square matrix.
 * @param[in] v: integer vector.
 * @param[in,out] res: solution vector of the vector matrix multiplication.
 */
void matrix_intVector(unsigned int dim, const double *m, const int *v, double *res);

/**
 * @brief square matrix transpose.
 * @param[in] dim: dimension of the square matrix.
 * @param[in,out] m: square matrix.
 */
void transpose(unsigned int dim, double *m);

/**
 * @brief check if two vectors are equal.
 * @param[in] dim: dimension of the vectors.
 * @param[in] v1: first vector.
 * @param[in] v2: second vector.
 * @return true if the vectors are equal, false if the vectors are not equal.
 */
bool equals(unsigned int dim, const double *v1, const double *v2);

/**
 * @brief Invert matrix.
 * @param[in] dim: dimension of the vectors.
 * @param[in, out] m: matrix to invert. overwritten bei LU-decomposition.
 * @param[out] p: permutation vector.
 * @param[out] r: where inverse matrix is stored.
 */
void invert(unsigned int dim, double *m, int *p, double *r);

/**
 * @brief Compute infinity norm (maximum sum row norm).
 * @param[in] dim: dimension of the vectors.
 * @param[in] m: matrix to compute infinity norm of.
 */
double inf_norm(unsigned int dim, const double *m);
#endif
