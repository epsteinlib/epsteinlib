// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file epsteinZeta.c
 * @brief Calculates the (regularized) Epstein zeta function.
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @see Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta
 * variants. Algorithmic Reflections: Selected Works. PSIpress (2012).
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @date 06/13/2024
 */

#include <complex.h>
#include <stdbool.h>

#include "epsteinZeta.h"
#include "zeta.h"

/**
 * @brief calculates the Epstein Zeta function.
 * @param nu: exponent for the Epstein zeta function.
 * @param dim: dimension of the input vectors.
 * @param a: matrix that transforms the lattice in the Epstein Zeta function.
 * @param x: x vector of the Epstein Zeta function.
 * @param y: y vector of the Epstein Zeta function.
 * @return function value of the Epstein zeta.
 */
double complex epsteinZeta(double nu, int dim, double *a, double *x, double *y) {
    return __epsteinZeta(nu, dim, a, x, y, 1, false);
}

/**
 * @brief calculates the regularized Epstein Zeta function.
 * @param nu: exponent for the Epstein zeta function.
 * @param dim: dimension of the input vectors.
 * @param a: matrix that transforms the lattice in the Epstein Zeta function.
 * @param x: x vector of the Epstein Zeta function.
 * @param y: y vector of the Epstein Zeta function.
 * @return function value of the regularized Epstein zeta.
 */
double complex epsteinZetaReg(double nu, int dim, double *a, double *x, double *y) {
    return __epsteinZeta(nu, dim, a, x, y, 1, true);
}
