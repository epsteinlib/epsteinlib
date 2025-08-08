// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file epsteinZeta.h
 * @brief Calculates the (regularized) Epstein zeta function.
 *
 * Main header file to include when using epsteinZeta library.
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @see Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta
 * variants. Algorithmic Reflections: Selected Works. PSIpress (2012).
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @date 06/13/2024
 */

#include <complex.h>

#ifndef EPSTEIN_H
#define EPSTEIN_H

/**
 * @brief calculates the Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta function.
 * @param[in] x: x vector of the Epstein Zeta function.
 * @param[in] y: y vector of the Epstein Zeta function.
 * @return function value of the regularized Epstein zeta.
 */
double complex epsteinZeta(double nu, unsigned int dim, const double *a,
                           const double *x, const double *y);

/**
 * @brief calculates the regularized Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta function.
 * @param[in] x: x vector of the Epstein Zeta function.
 * @param[in] y: y vector of the Epstein Zeta function.
 * @return function value of the regularized Epstein zeta.
 */
double complex epsteinZetaReg(double nu, unsigned int dim, const double *a,
                              const double *x, const double *y);

/**
 * @brief calculates the derivatives of the set zeta function for lattices.
 * @param[in] nu: exponent for the set zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the set zeta
 * function.
 * @param[in] x: x vector of the set zeta function.
 * @param[in] y: y vector of the set zeta function.
 * @param[in] alpha: multiindex for the derivative of the set zeta function.
 * @return function value of the set zeta.
 */
double complex setZetaDer(double nu, unsigned int dim, const double *a,
                          const double *x, const double *y,
                          const unsigned int *alpha);
/**
 * @brief calculates the derivatives of the regularized Epstein zeta function for
 * lattices.
 * @param[in] nu: exponent for the set zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the set zeta
 * function.
 * @param[in] x: x vector of the set zeta function.
 * @param[in] y: y vector of the set zeta function.
 * @param[in] alpha: multiindex for the derivative of the set zeta function.
 * @return function value of the set zeta.
 */
double complex epsteinZetaRegDer(double nu, unsigned int dim, const double *a,
                                 const double *x, const double *y,
                                 const unsigned int *alpha);

/**
 * @brief Calculates the incomplete bessel function.
 * @param[in] nu: exponent of the function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] k: input vector of the function.
 * @param[in] r: input vector of the function.
 * @return 2 int_0**1 t**(-nu - 1) exp(-pi k**2 / t**2) exp(-pi r**2 t**2) dt
 */
double incomplete_bessel_g(double nu, unsigned int dim, const double *k,
                           const double *r);
#endif
