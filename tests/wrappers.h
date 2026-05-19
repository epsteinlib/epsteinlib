// SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file wrappers.h
 * @brief Wrappers of functions that are included in unit tests but do not appear as
 * is in the core library.
 */

#include <complex.h>

#ifndef WRAPPERS_H
#define WRAPPERS_H

/** @brief Computes a single summand of h_inner; explicitly
 * (−1/2)^{i} / cn,i,k binom(i+k,k) binom(i,β) binom(α,θ1) θ2! / (θ2 - θ1)!,
 * where n=|α|, i=|β|, θ1=α+β−γ, θ2=γ−β.
 * @param[in] n: index (total of alpha).
 * @param[in] i: index (total of beta).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension of the multi-indices.
 * @param[in] beta: lower multi-index β.
 * @param[in] alpha: upper multi-index α.
 * @param[in] theta1: multi-index α+β−γ.
 * @param[in] theta2: multi-index γ−β.
 * @return value of one summand in of h_inner.
 */
double harmonic_h_inner_term(unsigned int n, unsigned int i, unsigned int k,
                             unsigned int dim, const unsigned int *alpha,
                             const unsigned int *beta, const unsigned int *theta1,
                             const unsigned int *theta2);

/** @brief Computes p_(alpha,beta)(y) = (-pi)^(alpha - beta) * (alpha choose beta)
 * ((alpha - beta)! / (alpha - 2*beta)!) * (2*y)^(alpha - 2*beta)
 * where 2 beta =< alpha
 * @param[in] dim: dimension of alpha, beta and y.
 * @param[in] y: vector of the polynomial.
 * @parma[in] alpha: upper multi-index.
 * @parma[in] beta: lower multi-index.
 * @return p(y).
 */
double polynomial_p(unsigned int dim, const double *z, const unsigned int *alpha,
                    const unsigned int *beta);

/**
 * @brief Calculates the upper Crandall function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast
 * @param[in] alpha: multi-index of the partial derivatives
 * @param[in] alphaAbs: sum of the elements of alpha
 * * @return upperGamma(nu / 2,pi prefactor * z**2) / (pi * prefactor z**2)^(nu /
 * 2).
 */
double complex crandall_g_der(unsigned int dim, double nu, const double *z,
                              double prefactor, double zArgBound,
                              const unsigned int *alpha, unsigned int alphaAbs);

/** @brief Calculates the polynomial l_(alpha,beta)(y) = - (-1)**|alpha - beta| *
 * binom(alpha,beta) * (alpha-beta)! / (alpha - 2 beta)! |alpha - beta|! / |alpha -
 * beta| * (2 * y)**(alpha - 2 beta) where 2 beta =< alpha
 * @param[in] dim: dimension of alpha, beta and y.
 * @param[in] y: vector of the polynomial.
 * @parma[in] alpha: upper multi-index.
 * @parma[in] beta: lower multi-index.
 * @return p(y).
 */
double polynomial_l(unsigned int dim, const double *z, const unsigned int *alpha,
                    const unsigned int *beta);

/** @brief Calculates the derivatives of L(z) = log(pi * z**2)
 * @param[in] dim: dimension of z.
 * @param[in] z: vector of the polynomial.
 * @parma[in] alpha: multi-index for the derivative.
 * @parma[in] alphaAbs: absolute value of the multi-index alpha.
 * @return partial derivative of L(z).
 */
double complex log_l_der(unsigned int dim, const double *z,
                         const unsigned int *alpha, unsigned int alphaAbs);

/** @brief Calculates the singularity s_{d+2k}(z) = pi**(k + d / 2) / gamma(k + d /
 * 2) * (-1)**(k+1) / k! * (pi * z**2)**k * log(pi * z**2)
 * @param[in] k: non-negative integer so that d + 2k is the argument of s.
 * @param[in] dim: dimension of z.
 * @param[in] z: vector of the singularity.
 * @parma[in] alpha: multi-index for the derivative.
 * @parma[in] alphaAbs: absolute value of the multi-index alpha.
 * @return partial derivative of s_{d+2k}(z).
 */
double singularity_s_der(unsigned int k, unsigned int dim, const double *z,
                         const unsigned int *alpha, unsigned int alphaAbs);

/**
 * @brief Calculates the derivatives of regularization of the zero summand in the
 * second sum in Crandall's formula in the special case of nu = dim + 2k for some
 * natural number k.
 * @param[in] s: s = (d - nu).
 * @param[in] k: k = (nu - d) / 2 as an integer.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] z: input vector of the function.
 * @param[in] lambda: scaling parameter of crandalls formula.
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast asymptotic
 * @return partial derivative of the regularized Crandall function. */
double crandall_gReg_nuequalsdimplus2k_der(double s, unsigned int k,
                                           unsigned int dim, const double *z,
                                           double lambda, const unsigned int *alpha,
                                           unsigned int alphaAbs, double zArgBound);

/**
 * @brief Calculates the derivatives of the regularization of the zero summand in the
 * second sum in Crandall's formula.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function,
 * that is d - nu.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda.
 * @param[in] alpha: multi-index of the partial derivatives
 * @param[in] alphaAbs: sum of the elements of alpha
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast asymptotic
 * @return partial derivatives of - gamma(s/2) * gammaStar(s/2, pi * prefactor *
 * z**2), where gammaStar is the twice regularized lower incomplete gamma function if
 * s is not equal to - 2k and partial derivatives of (pi * prefactor * y ** 2) ** (-
 * s / 2) (gamma(s / 2, pi * prefactor * z ** 2) + ((-1)^k / k! ) * (log(pi * y ** 2)
 * - log(prefactor ** 2))) if s is  equal to - 2k for non negative natural number k.
 */
double complex crandall_gReg_der(unsigned int dim, double s, const double *z,
                                 double prefactor, const unsigned int *alpha,
                                 unsigned int alphaAbs, double zArgBound);

#endif
