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
 * sum in Crandall's formula in the special case of
 * nu = dim + 2k for some natural number k.
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function.
 * @param[in] arg: input of the function
 * @param[in] k: k = - s / 2 = (nu - d) / 2 as an integer
 * @param[in] lambda: scaling parameter of crandalls formula
 * @return arg ** (- s / 2) * (gamma(s / 2, arg) + ((-1)^k / k! ) * (log(arg) -
 * log(lambda ** 2))
 */
double complex crandall_gReg_nuequalsdim(double s, double arg, double k,
                                         double lambda);

/**
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula.
 * @param[in] dim: dimension of the input vectors
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function,
 * that is d - nu
 * @param[in] z: input vector of the function
 * @param[in] prefactor: prefactor of the vector, e. g. lambda
 * @return - gamma(s/2) * gammaStar(s/2, pi * prefactor * z**2),
 * where gammaStar is the twice regularized lower incomplete gamma function if s is
 * not equal to - 2k and (pi * prefactor * y ** 2) ** (- s / 2)
 * (gamma(s / 2, pi * prefactor * z ** 2) + ((-1)^k / k! ) * (log(pi * y ** 2) -
 * log(prefactor ** 2))) if s is  equal to - 2k for non negative natural number k
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
