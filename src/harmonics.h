// SPDX-FileCopyrightText: 2025-2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file harmonics.c
 * @brief Calculates the harmonic polynomials.
 */

#include "hpdyad.h"

#ifndef EPSTEIN_HARMONICS
#define EPSTEIN_HARMONICS

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

/** @brief Computes cₙ,ᵢ,ₖ＝ ∏_{j=i+1}^{⌊n/2⌋-k} (2n＋d−2−4k−2j) as hpdyad.
 * @param[in] n: index (total of alpha).
 * @param[in] i: lower bound of product (exclusive).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result hpdyad.
 */
void coeffs_c_inner_hpdyad(unsigned int n, unsigned int i, unsigned int k,
                           unsigned int dim, hpdyad_t *out);

/** @brief Computes (-2)**(-i) · c_{n,i,k} · binom(i+k,k) as hpdyad.
 * @param[in] n: index (total of alpha).
 * @param[in] i: index (total of beta).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result hpdyad.
 */
void harmonic_h_inner_term_scalar_hpdyad(unsigned int n, unsigned int i,
                                         unsigned int k, unsigned int dim,
                                         hpdyad_t *out);

/** @brief Computes binom(i,β) × binom(α,θ₁) × θ₂!/(θ₂−θ₁)! as hpdyad.
 * @param[in] dim: dimension of the multi-indices.
 * @param[in] alpha: upper multi-index α.
 * @param[in] beta: lower multi-index β.
 * @param[in] theta1: multi-index α+β−γ.
 * @param[in] theta2: multi-index γ−β.
 * @param[out] out: result hpdyad (always positive).
 */
void harmonic_h_inner_term_multi_hpdyad(unsigned int dim, const unsigned int *alpha,
                                        const unsigned int *beta,
                                        const unsigned int *theta1,
                                        const unsigned int *theta2, hpdyad_t *out);

/** @brief Computes the inner sum h_inner(α,γ,k) using exact hpdyad arithmetic.
 * @pre k ≤ ⌊alphaAbs/2⌋
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

#endif
