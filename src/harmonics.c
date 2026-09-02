// SPDX-FileCopyrightText: 2025-2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file harmonics.c
 * @brief Calculates the harmonic polynomials.
 */

#include "harmonics.h"
#include "hpdyad.h"
#include "stdbool.h"
#include "tools.h"
#include <math.h>
#include <stddef.h>

/** @brief Computes cₙ,ᵢ,ₖ＝ ∏_{j=i+1}^{⌊n/2⌋-k} (2n＋d−2−4k−2j).
 * @param[in] n: index (total of alpha).
 * @param[in] i: lower bound of product (exclusive).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension of the zeta function inputs.
 * @return cₙ,ᵢ,ₖ.
 */
double coeffs_c_inner(long long n, long long i, long long k, long long dim) {
    long long end = (n / 2) - k;

    double res = 1.0;
    long long term = (2 * n) + dim - 4 - (4 * k) - (2 * i); // value at j = i+1

    for (long long j = i + 1; j <= end; j++) {
        res *= (double)term;
        term -= 2;
    }
    return res;
}

/** @brief Computes cₙ,ₖ＝2⁻ᵏ/ cₙ,₀,ₖ ∏ⱼ₌₁ᵏ(2n＋d−2j)∕((2n＋d＋2−4j)(2n＋d−4j)).
 * @param[in] n: index (total of alpha).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension of the zeta function inputs.
 * @return cₙ,ₖ.
 */
double coeffs_c_outer(long long n, long long k, long long dim) {
    long long a = (2 * n) + dim;

    double num_prod = 1.0;
    double den_prod = 1.0;

    long long num_term = a - 2;  // a - 2j  at j=1
    long long den1_term = a - 2; // a + 2 - 4j at j=1
    long long den2_term = a - 4; // a - 4j at j=1

    for (long long j = 0; j < k; j++) {
        num_prod *= (double)num_term;
        den_prod *= (double)den1_term * (double)den2_term;
        num_term -= 2;
        den1_term -= 4;
        den2_term -= 4;
    }

    double inner = coeffs_c_inner(n, 0, k, dim);

    return ldexp(num_prod / (den_prod * inner), -(int)k);
}

/** @brief Computes cₙ,ᵢ,ₖ＝ ∏_{j=i+1}^{⌊n/2⌋-k} (2n＋d−2−4k−2j) as hpdyad.
 * @param[in] n: index (total of alpha).
 * @param[in] i: lower bound of product (exclusive).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result hpdyad.
 */
void coeffs_c_inner_hpdyad(unsigned int n, unsigned int i, unsigned int k,
                           unsigned int dim, hpdyad_t *out) {

    long long end = (long long)(n / 2) - (long long)k;

    hpdyad_set_ull(out, 1, 1);
    if (i >= end) {
        hpdyad_normalize(out);
        return;
    }
    int base = (2 * (int)n) + (int)dim - 2 - (4 * (int)k);
    for (unsigned int j = i + 1; j <= end; j++) {
        int factor = base - (2 * (int)j);
        hpdyad_t tmp;
        hpdyad_set_ull(&tmp, (unsigned long long)factor, 1);
        hpdyad_mul(out, out, &tmp);
    }
}

/** @brief Computes (-2)**(-i) · c_{n,i,k} · binom(i+k,k) as hpdyad.
 * @param[in] n: index (total of alpha).
 * @param[in] i: index (total of beta).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result hpdyad.
 */
void harmonic_h_inner_term_scalar_hpdyad(unsigned int n, unsigned int i,
                                         unsigned int k, unsigned int dim,
                                         hpdyad_t *out) {
    unsigned long long b = binom((long long)(i + k), (long long)k);
    hpdyad_set_ull(out, b, 1);

    // c_{n,i,k}
    hpdyad_t coeff;
    coeffs_c_inner_hpdyad(n, i, k, dim, &coeff);
    hpdyad_mul(out, out, &coeff);

    hpdyad_normalize(out);

    // (-1/2)^i
    if (i & 1) {
        out->sign = -1;
    }
    out->exp2 -= (int)i;
}

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
                                        const unsigned int *theta2, hpdyad_t *out) {
    hpdyad_set_ull(out, 1, 1);

    // prod_j binom(betaAbsDim, beta[j])
    unsigned long long betaAbsDim = 0;
    for (unsigned int j = 0; j < dim; j++) {
        betaAbsDim += beta[j];
        unsigned long long b = binom((long long)betaAbsDim, (long long)beta[j]);
        hpdyad_t tmp;
        hpdyad_set_ull(&tmp, b, 1);
        hpdyad_mul(out, out, &tmp);
    }

    // prod_j binom(alpha[j], theta1[j])
    for (unsigned int j = 0; j < dim; j++) {
        unsigned long long b = binom((long long)alpha[j], (long long)theta1[j]);
        hpdyad_t tmp;
        hpdyad_set_ull(&tmp, b, 1);
        hpdyad_mul(out, out, &tmp);
    }

    // prod_j falling factorial: theta2[j]! / (theta2[j] - theta1[j])!
    for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int l = theta2[j] - theta1[j] + 1; l <= theta2[j]; l++) {
            hpdyad_t tmp;
            hpdyad_set_ull(&tmp, (unsigned long long)l, 1);
            hpdyad_mul(out, out, &tmp);
        }
    }

    hpdyad_normalize(out);
}

/** @brief Computes the inner sum h_inner(α,γ,k) using exact hpdyad arithmetic.
 * @pre k ≤ ⌊alphaAbs/2⌋
 * @param[in] k: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of alpha, beta and gamma.
 * @param[in] alpha: upper multi-index.
 * @param[in] gamma: fixed multi-index gamma.
 * @param[in] alphaAbs: total of alpha.
 * @return h_inner(α,γ,k).
 */
// fast paths for incremental multi-index enumeration increases nesting
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
/** @brief Minimal double-double arithmetic used to evaluate the outer sum of
 * harmonic_h. The outer sum is ill-conditioned: at |alpha| = 41 the ratio
 * sum|term| / |sum| reaches 1e6, so each coefficient's stored rounding error of
 * eps/2 is amplified to ~1e-10 relative. Carrying a residual limb for both the
 * coefficients and the accumulator reduces that to eps^2 * cond, i.e. below the
 * final rounding.
 */
typedef struct {
    double hi;
    double lo;
} dd_t;

/** @brief Exact sum of two doubles as a double-double (Knuth TwoSum). */
static inline dd_t dd_two_sum(double a, double b) {
    dd_t r;
    double bb;
    r.hi = a + b;
    bb = r.hi - a;
    r.lo = (a - (r.hi - bb)) + (b - bb);
    return r;
}

/** @brief Exact product of two doubles as a double-double (TwoProduct via fma).
 */
static inline dd_t dd_two_prod(double a, double b) {
    dd_t r;
    r.hi = a * b;
    r.lo = fma(a, b, -r.hi);
    return r;
}

/** @brief Sum of two double-doubles, renormalized. */
static inline dd_t dd_add(dd_t a, dd_t b) {
    dd_t s = dd_two_sum(a.hi, b.hi);
    s.lo += a.lo + b.lo;
    return dd_two_sum(s.hi, s.lo);
}

/** @brief Product of two double-doubles, renormalized. */
static inline dd_t dd_mul(dd_t a, dd_t b) {
    dd_t p = dd_two_prod(a.hi, b.hi);
    p.lo += (a.hi * b.lo) + (a.lo * b.hi);
    return dd_two_sum(p.hi, p.lo);
}

/** @brief Represent a double exactly as an hpdyad_t.
 * @param[out] a: receives the exact value of x.
 * @param[in] x: finite double.
 * @return void
 */
static void hpdyad_set_double_exact(hpdyad_t *a, double x) {
    signed char sg;
    int e;
    double f;
    unsigned long long m;

    if (x == 0.) {
        hpdyad_set_ull(a, 0, 1);
        return;
    }
    sg = (x < 0.) ? -1 : 1;
    f = frexp(fabs(x), &e);
    m = (unsigned long long)ldexp(f, 53);
    hpdyad_set_ull(a, m, sg);
    a->exp2 += e - 53;
}

/** @brief Split an exact hpdyad value into the nearest double and the exact
 * residual. The residual is exact because hpdyad_add takes the subtraction path
 * on a sign mismatch.
 * @param[in] v: exact value.
 * @param[out] hi: correctly rounded double nearest to v.
 * @param[out] lo: double nearest to v - hi.
 * @return void
 */
static void hpdyad_split_double(const hpdyad_t *v, double *hi, double *lo) {
    hpdyad_t h;
    hpdyad_t neg;
    hpdyad_t rest;

    *hi = hpdyad_to_double(v);
    hpdyad_set_double_exact(&h, *hi);
    neg = h;
    neg.sign = (signed char)(-neg.sign);
    hpdyad_add(&rest, v, &neg);
    *lo = hpdyad_to_double(&rest);
}

/** @brief Computes the inner sum h_inner(α,γ,k) using exact hpdyad arithmetic.
 * @pre k ≤ ⌊alphaAbs/2⌋
 * @param[in] k: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of alpha, beta and gamma.
 * @param[in] alpha: upper multi-index.
 * @param[in] gamma: fixed multi-index gamma.
 * @param[in] alphaAbs: total of alpha.
 * @param[out] residual: if non-NULL, receives the exact remainder
 * h_inner(α,γ,k) minus the returned double. The outer sum in harmonic_h is
 * ill-conditioned, so this second limb is what keeps the result correctly
 * rounded at high |alpha|.
 * @return h_inner(α,γ,k), correctly rounded.
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
double harmonic_h_inner_sum(unsigned int k, unsigned int dim,
                            const unsigned int *alpha, const unsigned int *gamma,
                            unsigned int alphaAbs, double *residual) {

    unsigned int beta[dim];
    unsigned int theta1[dim];
    unsigned int theta2[dim];
    for (unsigned int j = 0; j < dim; j++) {
        beta[j] = 0;
        theta1[j] = alpha[j] + beta[j] - gamma[j];
        theta2[j] = gamma[j] - beta[j];
    }

    bool redotheta1 = false;
    unsigned int betaAbs = 0;

    unsigned long long totalCount = 1;
    for (unsigned int j = 0; j < dim; j++) {
        totalCount *= gamma[j] - (alpha[j] + 1) / 2 + 1;
    }

    unsigned int lastScalarIndex = (alphaAbs / 2) - k;

    // scalar coefficients as hpdyad
    hpdyad_t scalar_coeffs[lastScalarIndex + 1];
    for (unsigned int i = 0; i <= lastScalarIndex; i++) {
        harmonic_h_inner_term_scalar_hpdyad(alphaAbs, i, k, dim, &scalar_coeffs[i]);
    }

    // multi-index coefficients as hpdyad (accumulated via hpdyad_add)
    hpdyad_t multi_coeffs[lastScalarIndex + 1];
    for (unsigned int i = 0; i <= lastScalarIndex; i++) {
        hpdyad_set_ull(&multi_coeffs[i], 0, 1);
    }

    // loop over beta multi-indices
    for (unsigned long long i = 0; i < totalCount; i++) {

        bool skip = false;
        for (unsigned int j = 0; j < dim; j++) {
            if (gamma[j] > alpha[j] && gamma[j] - alpha[j] > beta[j]) {
                skip = true;
                break;
            }
        }

        if (!skip) {

            if (redotheta1) {
                for (unsigned int j = 0; j < dim; j++) {
                    theta1[j] = alpha[j] + beta[j] - gamma[j];
                }
            }
            redotheta1 = false;

            // unconditionally update theta2
            for (unsigned int j = 0; j < dim; j++) {
                theta2[j] = gamma[j] - beta[j];
            }

            hpdyad_t term;
            harmonic_h_inner_term_multi_hpdyad(dim, alpha, beta, theta1, theta2,
                                               &term);
            hpdyad_add(&multi_coeffs[betaAbs], &multi_coeffs[betaAbs], &term);
        }

        for (unsigned int idx = 0; idx < dim; idx++) {
            if (2 * beta[idx] + 2 <= 2 * gamma[idx] - alpha[idx]) {
                beta[idx]++;
                theta1[idx]++;
                betaAbs++;
                break;
            }
            betaAbs -= beta[idx];
            beta[idx] = 0;
            redotheta1 = true;
        }
    }

    // combine: sum_i scalar_coeffs[i] * multi_coeffs[i]
    hpdyad_t sumInner;
    hpdyad_set_ull(&sumInner, 0, 1);

    for (unsigned int i = 0; i <= lastScalarIndex; i++) {
        hpdyad_t term;
        hpdyad_mul(&term, &scalar_coeffs[i], &multi_coeffs[i]);
        hpdyad_add(&sumInner, &sumInner, &term);
    }

    // convert to double, keeping the residual so that the ill-conditioned outer
    // sum in harmonic_h can be evaluated in double-double arithmetic
    double result;
    double lo;
    hpdyad_split_double(&sumInner, &result, &lo);
    if (residual != NULL) {
        *residual = lo;
    }

    return result;
}

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
                                       unsigned long long *valid_count) {
    unsigned int k;
    unsigned long long totalSize = 0;
    unsigned int ceilHalfAlphaAbs = 0;
    unsigned int i;
    long long remainder;

    // Compute sum of ceil(alpha[i]/2)
    for (i = 0; i < dim; i++) {
        ceilHalfAlphaAbs += (alpha[i] + 1) / 2;
    }

    for (k = 0; k <= kMax; k++) {
        chunk_offset[k] = totalSize;
        remainder = (long long)(alphaAbs - k) - (long long)ceilHalfAlphaAbs;
        if (remainder < 0) {
            valid_count[k] = 0;
        } else {
            valid_count[k] = binom(remainder + dim - 1, dim - 1);
        }
        totalSize += valid_count[k];
    }

    return totalSize;
}

/** @brief Precomputes and stores inner harmonic sums h_inner(α,γ,k)
 * and exponents (2γ-α) for all k = 0, ..., floor(|alpha|/2) and all gamma
 * with |gamma| = |alpha| - k satisfying 2γ[i] >= α[i].
 * Coefficients are stored in coeffs starting at offsets given by chunk_offset[k].
 * Exponents are stored in exponents with stride dim, i.e., exponent i for
 * entry n of chunk k is at exponents[(chunk_offset[k] + n) * dim + i].
 * Gamma is iterated in graded lexicographic order with total degree |alpha| - k.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] dim: dimension of alpha and gamma.
 * @param[in] alpha: upper multi-index.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[out] coeffs: array storing precomputed inner harmonic sums.
 * @param[out] exponents: array storing precomputed exponents (2γ-α), size =
 * totalSize * dim.
 */
// fast paths for incremental multi-index enumeration increases nesting
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void precompute_harmonic_h_inner_sum(unsigned int alphaAbs, unsigned int dim,
                                     const unsigned int *alpha,
                                     const unsigned long long *chunk_offset,
                                     double *coeffs, unsigned int *exponents) {
    unsigned int kMax = alphaAbs / 2;

    for (unsigned int k = 0; k <= kMax; k++) {

        unsigned int gamma[dim];
        for (unsigned int j = 0; j < dim; j++) {
            gamma[j] = 0;
        }
        unsigned int gammaAbs = 0;
        unsigned int maxAbs = alphaAbs - k;

        // totalCount = (maxAbs + 1)^dim
        unsigned long long totalCount = 1;
        for (unsigned int j = 0; j < dim; j++) {
            totalCount *= maxAbs + 1;
        }

        unsigned long long n = 0;
        for (unsigned long long i = 0; i < totalCount; i++) {

            bool skip = (gammaAbs != maxAbs);
            if (!skip) {
                for (unsigned int j = 0; j < dim; j++) {
                    if (2 * gamma[j] < alpha[j]) {
                        skip = true;
                        break;
                    }
                }
            }

            if (!skip) {
                // coeffs holds interleaved (hi, lo) pairs: the exact inner sum
                // does not fit a double, and the discarded limb dominates the
                // error of the ill-conditioned outer sum in harmonic_h
                coeffs[2 * (chunk_offset[k] + n)] =
                    harmonic_h_inner_sum(k, dim, alpha, gamma, alphaAbs,
                                         &coeffs[(2 * (chunk_offset[k] + n)) + 1]);
                unsigned long long expIdx = (chunk_offset[k] + n) * dim;
                for (unsigned int j = 0; j < dim; j++) {
                    exponents[expIdx + j] = (2 * gamma[j]) - alpha[j];
                }
                n++;
            }

            for (unsigned int idx = 0; idx < dim; idx++) {
                if (gamma[idx] + 1 <= maxAbs) {
                    gamma[idx]++;
                    gammaAbs++;
                    break;
                }
                gammaAbs -= gamma[idx];
                gamma[idx] = 0;
            }
        }
    }
}

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
                  const unsigned int *exponents) {

    double sumOuter;
    double maxTerm = 0.0;
    unsigned long long n;
    unsigned long long count = valid_count[k];
    unsigned long long baseIdx = chunk_offset[k];
    unsigned long long expIdx;
    unsigned int i;
    unsigned int e;
    // exponents sum to |alpha| - 2k, so no single one exceeds that bound;
    // sizing the table to alphaAbs would roughly double the setup cost
    unsigned int maxExp = alphaAbs - (2 * k);
    dd_t acc;
    // z[i]^e for every exponent the chunk can contain. Tabulating costs one
    // multiplication per entry and removes repeated squaring from the inner
    // loop, so double-double evaluation stays close to the cost of the plain
    // power chain it replaces.
    dd_t zPowTab[dim][maxExp + 1];

    for (i = 0; i < dim; i++) {
        dd_t zi;
        zi.hi = z[i];
        zi.lo = 0.0;
        zPowTab[i][0].hi = 1.0;
        zPowTab[i][0].lo = 0.0;
        for (e = 1; e <= maxExp; e++) {
            zPowTab[i][e] = dd_mul(zPowTab[i][e - 1], zi);
        }
    }

    acc.hi = 0.0;
    acc.lo = 0.0;

    for (n = 0; n < count; n++) {
        dd_t coeff;
        dd_t zPow;
        dd_t summand;
        coeff.hi = coeffs[2 * (baseIdx + n)];
        coeff.lo = coeffs[(2 * (baseIdx + n)) + 1];
        expIdx = (baseIdx + n) * dim;
        zPow = zPowTab[0][exponents[expIdx]];
        for (i = 1; i < dim; i++) {
            zPow = dd_mul(zPow, zPowTab[i][exponents[expIdx + i]]);
        }
        summand = dd_mul(zPow, coeff);
        maxTerm = fmax(maxTerm, fabs(summand.hi));
        acc = dd_add(acc, summand);
    }

    sumOuter = acc.hi + acc.lo;

    // avoid accumulation error that is pure rounding noise to be amplified by
    // divergent Crandall factor upon return
    if (fabs(sumOuter) <= EPS_CANCELLATION_DD * (double)count * maxTerm) {
        return 0.0;
    }

    sumOuter *= coeffs_c_outer(alphaAbs, k, dim);

    return sumOuter;
}
