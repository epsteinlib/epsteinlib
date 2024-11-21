// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include <complex.h>

#ifndef UTILS_H
#define UTILS_H

/*!
 * @brief absolute difference between to complex numbers.
 * @param ref: reference value.
 * @param comp: compparison value.
 * @return sqrt(rev - comp)
 */
double errAbs(double complex ref, double complex comp);

/*!
 * @brief relative difference between to complex numbers.
 * @param ref: reference value.
 * @param comp: compparison value.
 * @return errAbs / norm(ref)
 */
double errRel(double complex ref, double complex comp);

/**
 * @brief prints vector to terminal, is used in unitTest.
 * @param name: name of the vector as a string.
 * @param vec: vector.
 * @param dim: size of the vector.
 */
void printVectorUnitTest(const char *name, double *vec, int dim);

/**
 * @brief prints square matrix to terminal, is used in unitTest.
 * @param name: name of the matrix as a string.
 * @param mat: matrix.
 * @param dim: size of the matrix.
 */
void printMatrixUnitTest(const char *name, const double *mat, int dim);
#endif
