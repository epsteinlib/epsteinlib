// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file crandall.h
 * @brief Calculates the summand function G and related functions
 * in Crandall's formula.
 */

#include "tools.h"
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
 * @brief Calculates the upper Crandall function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast asymptotic
 * expansion in the calculation of the Crandall function.
 * @return upperGamma(nu / 2,pi prefactor * z**2) / (pi * prefactor z**2)^(nu / 2) if
 |z| > 0 and - 2 / nu otherwise.

 */
double complex crandall_g(unsigned int dim, double nu, const double *z,
                          double prefactor, double zArgBound);

/**
 * @brief Calculates the lower Crandall function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @return lowerGamma(nu / 2,pi prefactor * z**2) / (pi * prefactor z**2)^(nu / 2) if
 * |z| > 0 and 2 / nu otherwise.
 */
double complex crandall_g_lower(unsigned int dim, double nu, const double *z,
                                double prefactor);

/** @brief Calculates the polynomial p_(alpha,beta)(y) = (-pi)^(alpha - beta) *
 * (alpha choose beta) *
 * ((alpha - beta)! / (alpha - 2*beta)!) * (2*y)^(alpha - 2*beta)
 * where 2 beta =< alpha
 * @param[in] dim: dimension of alpha, beta and y.
 * @param[in] z: vector of the polynomial.
 * @parma[in] alpha: upper multi-index.
 * @parma[in] beta: lower multi-index.
 * @return p(y).
 */
double polynomial_p(unsigned int dim, const double *z, const unsigned int *alpha,
                    const unsigned int *beta);

/**
 * @brief Calculates the partiall deritavies of the upper Crandall function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast asymptotic
 * expansion in the calculation of the upper Crandall function.
 * @param[in] alpha: multi-index for the partiall derivative.
 * @return upperGamma(nu / 2,pi prefactor * z**2) / (pi * prefactor z**2)^(nu / 2).
 */
double complex crandall_g_der(unsigned int dim, double nu, const double *z,
                              double prefactor, double zArgBound,
                              const unsigned int *alpha, unsigned int alphaAbs);

/** @brief Computes cₙ,ᵢ,ₖ＝ ∏_{j=i+1}^{⌊n/2⌋-k} (2n＋d−2−4k−2j).
 * @param[in] n: index (total of alpha).
 * @param[in] start: lower bound of product (exclusive).
 * @param[in] end: upper bound of product (inclusive).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension of the zeta function inputs.
 * @return cₙ,ᵢ,ₖ.
 */
double coeffs_c_inner(long long n, long long i, long long k, long long dim);

/** @brief Computes cₙ,ₖ＝2⁻ᵏ/ cₙ,₀,ₖ ∏ⱼ₌₁ᵏ(2n＋d−2j)∕((2n＋d＋2−4j)(2n＋d−4j)).
 * @param[in] n: index (total of alpha).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension of the zeta function inputs.
 * @return cₙ,ₖ.
 */
double coeffs_c_outer(long long n, long long k, long long dim);

/** @brief Computes cₙ,ᵢ,ₖ＝ ∏_{j=i+1}^{⌊n/2⌋-k} (2n＋d−2−4k−2j) as apint.
 * @param[in] n: index (total of alpha).
 * @param[in] i: lower bound of product (exclusive).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result apint.
 */
void coeffs_c_inner_apint(unsigned int n, unsigned int i, unsigned int k,
                          unsigned int dim, apint_t *out);

/** @brief Computes (-2)**(-i) · c_{n,i,k} · binom(i+k,k) as apint.
 * @param[in] n: index (total of alpha).
 * @param[in] i: index (total of beta).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result apint.
 */
void harmonic_h_inner_term_scalar_apint(unsigned int n, unsigned int i,
                                        unsigned int k, unsigned int dim,
                                        apint_t *out);

/** @brief Computes binom(i,β) × binom(α,θ₁) × θ₂!/(θ₂−θ₁)! as apint.
 * @param[in] dim: dimension of the multi-indices.
 * @param[in] alpha: upper multi-index α.
 * @param[in] beta: lower multi-index β.
 * @param[in] theta1: multi-index α+β−γ.
 * @param[in] theta2: multi-index γ−β.
 * @param[out] out: result apint (always positive).
 */
void harmonic_h_inner_term_multi_apint(unsigned int dim, const unsigned int *alpha,
                                       const unsigned int *beta,
                                       const unsigned int *theta1,
                                       const unsigned int *theta2, apint_t *out);

/** @brief Computes the inner sum h_inner(α,γ,k) using exact apint arithmetic.
 * @param[in] k: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of alpha, beta and gamma.
 * @param[in] alpha: upper multi-index.
 * @param[in] gamma: fixed multi-index gamma.
 * @param[in] alphaAbs: total of alpha.
 * @return h_inner(α,γ,k).
 */
double harmonic_h_inner_sum(unsigned int k, // NOLINT
                            unsigned int dim, const unsigned int *alpha,
                            const unsigned int *gamma, unsigned int alphaAbs);

/** @brief Computes chunk offsets and valid entry counts for precomputed
 * inner harmonic sums corresponding to k = 0, ..., floor(|alpha|/2).
 * For each k, chunk_offset[k] stores the starting offset in the coeffs array
 * where values corresponding to |gamma| = |alpha| - k are stored.
 * valid_count[k] stores the number of valid gamma for that k.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] kMax: floor(|alpha|/2).
 * @param[in] dim: dimension of alpha and gamma.
 * @param[in] alpha: upper multi-index.
 * @param[out] chunk_offset: array of length kMax+1 storing offsets.
 * @param[out] valid_count: array of length kMax+1 storing valid entry counts.
 * @return total length required for coeffs array.
 */
unsigned long long
precompute_harmonic_h_inner_chunk_size(unsigned int alphaAbs, unsigned int kMax,
                                       unsigned int dim, const unsigned int *alpha,
                                       unsigned long long *chunk_offset,
                                       unsigned long long *valid_count);

/** @brief Precomputes and stores inner harmonic sums h_inner(α,γ,k)
 * and exponents (2γ-α) for all k = 0, ..., floor(|alpha|/2) and all gamma
 * with |gamma| = |alpha| - k satisfying 2γ[i] >= α[i].
 * Coefficients are stored in coeffs starting at offsets given by chunk_offset[k].
 * Exponents are stored in exponents with stride dim, i.e., exponent i for
 * entry n of chunk k is at exponents[(chunk_offset[k] + n) * dim + i].
 * Gamma is ordered identically to harmonic_h_update_outer_index.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] dim: dimension of alpha and gamma.
 * @param[in] alpha: upper multi-index.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[out] coeffs: array storing precomputed inner harmonic sums.
 * @param[out] exponents: array storing precomputed exponents (2γ-α), size =
 * totalSize * dim.
 */
void precompute_harmonic_h_inner_sum(unsigned int alphaAbs, // NOLINT
                                     unsigned int dim, const unsigned int *alpha,
                                     const unsigned long long *chunk_offset,
                                     double *coeffs, unsigned int *exponents);

/** @brief Calculates the homogeneous harmonic polynomial h₍α,k₎
 * of degree |α|−2k such that y^α = ∑ₖ (y·y)^k h₍α,k₎(y) for d = 1 and k =
 * floor(|alpha|/2), explicitly, h₍α,k₎(y)=c_{|α|,k} ∑{|γ|=|α|−k} y^{2γ−α}
 * h_inner(α,γ,k). In 1D, h₍α,k₎(y) is y ** (alphaAbs mod 2) and zero otherwise. Uses
 * precomputed coefficients and exponents to avoid multi-index iteration.
 * @param[in] z: vector of the polynomial.
 * @param[in] alphaAbs: total of alpha.
 * @return h₍α,k₎(z).
 */
double harmonic_h_1D_kMax(const double *z, unsigned int alphaAbs);

/** @brief Calculates the homogeneous harmonic polynomial h₍α,k₎
 * of degree |α|−2k such that y^α = ∑ₖ (y·y)^k h₍α,k₎(y);
 * explicitly, h₍α,k₎(y)=c_{|α|,k} ∑{|γ|=|α|−k} y^{2γ−α} h_inner(α,γ,k).
 * Uses precomputed coefficients and exponents to avoid multi-index iteration.
 * @param[in] k: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of alpha, gamma and y.
 * @param[in] z: vector of the polynomial.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid entries for each k.
 * @param[in] coeffs: array storing precomputed inner harmonic sums.
 * @param[in] exponents: array storing precomputed exponents (2γ-α).
 * @return h₍α,k₎(z).
 */
double harmonic_h(unsigned int k, unsigned int dim, const double *z,
                  unsigned int alphaAbs, const unsigned long long *chunk_offset,
                  const unsigned long long *valid_count, const double *coeffs,
                  const unsigned int *exponents);

/** @brief Calculates the polynomial l_(alpha,beta)(z) = - (-1)**|alpha - beta| *
 * binom(alpha,beta) * (alpha-beta)! / (alpha - 2 beta)! |alpha - beta|! / |alpha -
 * beta| * (2 * z)**(alpha - 2 beta) where 2 beta =< alpha
 * @param[in] dim: dimension of alpha, beta and z.
 * @param[in] z: vector of the polynomial.
 * @parma[in] alpha: upper multi-index.
 * @parma[in] beta: lower multi-index.
 * @return l(z).
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

/** @brief Calculates the derivatives of Y_k(z) / n! = (pi * z**2)**k / n! where n <=
 * k.
 * @param[in] k: integer power.
 * @param[in] dim: dimension of z.
 * @param[in] y: vector of the polynomial.
 * @parma[in] alpha: multi-index for the derivative.
 * @param[in] n: factorial divisor smaller than k.
 * @return partial derivative of Y_k(z) / n!.
 */
double polynomial_y_der(unsigned int k, unsigned int dim, const double *z, // NOLINT
                        const unsigned int *alpha, unsigned int alphaAbs,
                        unsigned int n);

/** @brief Calculates the singularity s_{d+2k}(z) = pi**(k + d / 2) / gamma(k + d /
 * 2) * (-1)**(k+1) / k! * (pi * z**2)**k * log(pi * z**2)
 * @param[in] k: non-negative integer so that d + 2k is the argument of s.
 * @param[in] dim: dimension of z.
 * @param[in] z: vector of the singularity.
 * @parma[in] alpha: multi-index for the derivative.
 * @parma[in] alphaAbs: absolute value of the multi-index alpha.
 * @return partial derivative of s_{d+2k}(z).
 */
double complex singularity_s_der(unsigned int k, unsigned int dim, const double *z,
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
 * @param[in] alpha: multi-index of the partial derivatives
 * @param[in] alphaAbs: sum of the elements of alpha
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast asymptotic
 * @return partial derivative of the regularized Crandall function.
 */
double complex crandall_gReg_nuequalsdimplus2k_der(double s, unsigned int k,
                                                   unsigned int dim, const double *z,
                                                   double lambda,
                                                   const unsigned int *alpha,
                                                   unsigned int alphaAbs);

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
