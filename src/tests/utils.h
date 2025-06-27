// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include <complex.h>

#ifndef UTILS_H
#define UTILS_H

/*!
 * @brief absolute difference between to complex numbers.
 * @param[in] ref: reference value.
 * @param[in] comp: compparison value.
 * @return sqrt(rev - comp)
 */
double errAbs(double complex ref, double complex comp);

/*!
 * @brief relative difference between to complex numbers.
 * @param[in] ref: reference value.
 * @param[in] comp: compparison value.
 * @return errAbs / norm(ref)
 */
double errRel(double complex ref, double complex comp);

/**
 * @brief prints double vector to terminal, is used in unitTest.
 * @param[in] name: name of the vector as a string.
 * @param[in] vec: vector.
 * @param[in] dim: size of the vector.
 */
void printVectorUnitTest(const char *name, const double *vec, unsigned int dim);

/**
 * @brief prints unsigned integer vector to terminal, is used in unitTest.
 * @param[in] name: name of the vector as a string.
 * @param[in] vec: vector.
 * @param[in] dim: size of the vector.
 */
void printMultiindexUnitTest(const char *name, const unsigned int *vec,
                             unsigned int dim);

/**
 * @brief prints square matrix to terminal, is used in unitTest.
 * @param[in] name: name of the matrix as a string.
 * @param[in] mat: matrix.
 * @param[in] dim: size of the matrix.
 */
void printMatrixUnitTest(const char *name, const double *mat, unsigned int dim);

/**
 * @brief Compute the factorial of a multi-index.
 * @param[in] alpha: multi-index.
 * @return factorial of alpha.
 */
unsigned int mult_fac(unsigned int dim, const unsigned int *alpha);

/**
 * @brief Compute a vector to the power of a multi-index.
 * @param[in] dim: dimension of alpha end vec.
 * @param[in] alpha: multi-index.
 * @param[in] vec: base vector.
 * @param[in] prefactor: prefactor of the vector.
 * @return (prefactor * vec) ** alpha.
 */
double mult_pow(unsigned int dim, const unsigned int *alpha, const double *vec);

#endif
