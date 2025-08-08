// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file zeta.h
 * @brief Calculates the (regularized) Epstein zeta function and the deriavtives of
 * the set zeta function.
 */

#ifndef ZETA_H
#define ZETA_H
#include <complex.h>

/**
 * @brief calculates the (regularized) Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 * @param[in] lambda: relative weight of the sums in Crandall's formula.
 * @param[in] variant: 0 for no regularization
 *                    1 for the regularization
 *                    2 for set Zeta (including derivatives)
 *                    3 for the regularization (including Derivatives).
 * @return function value of the regularized Epstein zeta.
 */
double complex epsteinZetaInternal(double nu, unsigned int dim, const double *m,
                                   const double *x, const double *y, double lambda,
                                   unsigned int variant, const unsigned int *alpha);
#endif
