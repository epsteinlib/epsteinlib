// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file zeta.h
 * @brief Calculates the (regularized) Epstein zeta function.
 */

#ifndef _ZETA_H_
#define _ZETA_H_
#include <complex.h>

/**
 * @brief calculates the (regularized) Epstein Zeta function.
 * @param nu: exponent for the Epstein zeta function.
 * @param dim: dimension of the input vectors.
 * @param m: matrix that transforms the lattice in the Epstein Zeta function.
 * @param x: x vector of the Epstein Zeta function.
 * @param y: y vector of the Epstein Zeta function.
 * @param lambda: relative weight of the sums in Crandall's formula.
 * @param regBool: 0 for no regularization, > 0 for the regularization.
 * @return function value of the regularized Epstein zeta.
 */
double complex __epsteinZeta(double nu, int dim, double *m, double *x, double *y,
                             double lambda, int regBool);
#endif
