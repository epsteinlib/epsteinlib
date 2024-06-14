// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
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
 * @param nu: exponent for the Epstein zeta function.
 * @param dim: dimension of the input vectors.
 * @param a: matrix that transforms the lattice in the Epstein Zeta function.
 * @param x: x vector of the Epstein Zeta function.
 * @param y: y vector of the Epstein Zeta function.
 * @return function value of the regularized Epstein zeta.
 */
double complex epsteinZeta(double nu, int dim, double *a, double *x, double *y);

#ifndef EPSTEIN_CRANDALL

/**
 * @brief Assumes x and y to be in the respective elementary lattice cell.
 * Multiply with exp(2 * PI * i * x * y) to get the second sum in Crandall's
 * @param dim: dimension of the input vectors.
 * @param nu: exponent of the regularized Epstein zeta function.
 * @param z: input vector of the function
 * @param prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @return upperGamma(nu/2,pi prefactor * z**2)
 *      / (pi * prefactor z**2)^(nu / 2) in
 */
double crandall_g(int dim, double nu, double *z, double prefactor,
                  double zArgBound);

/**
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula. That is, the summand for k = 0.
 * @param dim: dimension of the input vectors
 * @param nu: exponent of the regularized Epstein zeta function.
 * @param z: input vector of the function
 * @param prefactor: prefactor of the vector, e. g. lambda
 * @return gamma(nu/2) * gammaStar(nu/2, pi * prefactor * z**2), where
 * gammaStar is the twice regularized lower incomplete gamma function
 * gamma(s,x) / (gamma(s) * x ** s)
 */
double crandall_gReg(int dim, double nu, double *z, double prefactor);
#endif
#endif
