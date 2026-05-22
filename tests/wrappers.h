// SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file wrappers.h
 * @brief Wrappers of functions that are included in unit tests but do not appear as
 * is in the core library.
 */

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

/** @brief Computes the derivatives singularity s_ν(z).
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
double singularity_s_der(double nu, unsigned int dim, const double *z,
                         const unsigned int *alpha, unsigned int alphaAbs);

#endif
