// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file crandall.h
 * @brief Calculates the summand function G and related functions
 * in Crandall's formula.
 */

#include <complex.h>

#ifndef EPSTEIN_CRANDALL
#define EPSTEIN_CRANDALL
/**
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula. That is, the summand for k = 0.
 * @param[in] dim: dimension of the input vectors
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function
 * @param[in] prefactor: prefactor of the vector, e. g. lambda
 * @return gamma(nu/2) * gammaStar(nu/2, pi * prefactor * z**2), where
 * gammaStar is the twice regularized lower incomplete gamma function
 * gamma(s,x) / (gamma(s) * x ** s)
 */
double complex crandall_gReg(unsigned int dim, double nu, const double *z,
                             double prefactor);

/**
 * @brief calculates bounds on when to use asymptotic expansion of the
 * upper incomplete gamma function, depending on the value of nu.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @return minimum value of z, when to use the fast asymptotic expansion in the
 * calculation of the incomplete upper gamma function upperGamma(nu, z).
 */
double assignzArgBound(double nu);

/**
 * @brief Assumes x and y to be in the respective elementary lattice cell.
 * Multiply with exp(2 * PI * i * x * y) to get the second sum in Crandall's
 * @param[in] dim: dimension of the input vectors.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function
 * @param[in] prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @return upperGamma(nu/2,pi prefactor * z**2)
 *      / (pi * prefactor z**2)^(nu / 2) in
 */
double complex crandall_g(unsigned int dim, double nu, const double *z,
                          double prefactor, double zArgBound);
#endif
