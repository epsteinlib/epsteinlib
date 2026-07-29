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
#include <stdbool.h>

/**
 * @brief calculates the (regularized) Epstein zeta function as well es the
 * derivatives of the (regularized) set zeta function with a prefactor.
 * @param[in] nu: exponent.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice .
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] lambda: relative weight of the sums in Crandall's formula.
 * @param[in] reg: false for no regularization, true for the regularization.
 * @param[in] aniso: false for no anisotropy, true for the anisotropic variant.
 * @param[in] alpha: multiindex for the derivatives.
 * @return function value.
 */
double complex epsteinZetaInternal(double nu, unsigned int dim, const double *m,
                                   const double *x, const double *y, double lambda,
                                   bool reg, bool aniso, const unsigned int *alpha);
#endif
