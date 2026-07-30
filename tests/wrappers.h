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

/**
 * @brief Computes the singularity s_ν(z).
 * For the generic case (ν ≠ d + 2ℓ for any non-negative integer ℓ):
 * s_ν(z) = π^(ν/2) / Γ(ν/2) · Γ((d−ν)/2) · (π z·z)^((ν−d)/2)
 * and for the logarithmic case (ν = d + 2ℓ, ℓ ≥ 0):
 * s_{d+2ℓ}(z) = π^(ℓ+d/2) / Γ(ℓ+d/2) · (−1)^(ℓ+1) / ℓ! · (π z·z)^ℓ · log(π z·z)
 * @param[in] nu: epxonent of singularity.
 * @param[in] dim: dimension of z.
 * @param[in] zSquared: |z|^2 for vector z of the singularity.
 * @return s_ν(z).
 */
double singularity_s(double nu, unsigned int dim, double zSquared);

/** @brief Calculates the derivatives of Y_ell(z) = (pi * z**2)**ell.
 * @param[in] ell: integer power.
 * @param[in] dim: dimension of z.
 * @param[in] z: vector z of the polynomial.
 * @parma[in] alpha: multi-index for the derivative.
 * @param[in] alphaAbs: |alpha|.
 * @return partial derivative of Y_ell(z).
 */
double polynomial_y_der_harmonic_wrapper(int ell, unsigned int dim,
                                         const double *z, // NOLINT
                                         const unsigned int *alpha,
                                         unsigned int alphaAbs);

/** @brief Wraps the derivatives of the singularity s_ν(z) by the harmonic method.
 *
 * For the generic case (ν ≠ d + 2ℓ for any
 * non-negative integer ℓ): s_ν(z) = π^(ν/2) / Γ(ν/2) · Γ((d−ν)/2) · (π
 * z·z)^((ν−d)/2) and for the logarithmic case (ν = d + 2ℓ, ℓ ≥ 0): s_{d+2ℓ}(z) =
 * π^(ℓ+d/2) / Γ(ℓ+d/2) · (−1)^(ℓ+1) / ℓ! · (π z·z)^ℓ · log(π z·z)
 * @param[in] nu: epxonent of singularity.
 * @param[in] dim: dimension of z.
 * @param[in] z: vector of the singularity.
 * @parma[in] alphaAbs: absolute value of the multi-index alpha.
 * @return partial derivative of s_ν(z).
 */
double singularity_s_der_harmonic_wrapper(double nu, unsigned int dim,
                                          const double *z, const unsigned int *alpha,
                                          unsigned int alphaAbs);

/** @brief Computes the derivatives singularity s_ν(z) by the harmonic method.
 * For the generic case (ν ≠ d + 2ℓ for any
 * non-negative integer ℓ): s_ν(z) = π^(ν/2) / Γ(ν/2) · Γ((d−ν)/2) · (π
 * z·z)^((ν−d)/2) and for the logarithmic case (ν = d + 2ℓ, ℓ ≥ 0): s_{d+2ℓ}(z) =
 * π^(ℓ+d/2) / Γ(ℓ+d/2) · (−1)^(ℓ+1) / ℓ! · (π z·z)^ℓ · log(π z·z)
 * @param[in] nu: epxonent of singularity.
 * @param[in] dim: dimension of z.
 * @param[in] z: vector of the singularity.
 * @parma[in] alphaAbs: absolute value of the multi-index alpha.
 * @return partial derivative of s_ν(z).
 */
double singularity_s_der(double nu, unsigned int dim, const double *z, // NOLINT
                         const unsigned int *alpha, unsigned int alphaAbs);

#endif
