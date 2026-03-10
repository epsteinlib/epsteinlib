// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file crandall.c
 * @brief Calculates the summand function G and related functions
 * in Crandall's formula.
 */

#include "gamma.h"
#include "tools.h"
#include <complex.h>
#include <float.h>
#include <math.h>

/*!
 * @brief epsilon for the cutoff around nu = dimension.
 */
#define EPS ldexp(1, -30)

/*!
 * @brief epsilon for the cutoff around arguments of the crandall function around
 * zero.
 */
#define EPS_ZERO_PIY (M_PI * 1e-64)

/**
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula in the special case of
 * nu = dim + 2k for some natural number k.
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function.
 * @param[in] arg: input of the function.
 * @param[in] k: k = - s / 2 = (nu - d) / 2 as an integer.
 * @param[in] lambda: scaling parameter of crandalls formula.
 * @return arg ** (- s / 2) * (gamma(s / 2, arg) + ((-1)^k / k! ) * (log(arg) -
 * log(lambda ** 2)).
 */
double complex crandall_gReg_nuequalsdimplus2k(double s, double arg, double k,
                                               double lambda) {
    double complex gReg = 0;
    // Taylor expansion if nu = dim and y close to zero.
    double taylorCutoff = 0.1 * 0.1 * M_PI;
    if (s == 0 && arg < taylorCutoff) {
        double eulerGamma = 0.57721566490153286555;
        double taylorCoeffs[10] = {-eulerGamma,
                                   1,
                                   -0.25,
                                   0.05555555555555555,
                                   -0.010416666666666666,
                                   0.0016666666666666668,
                                   -0.0002314814814814815,
                                   0.00002834467120181406,
                                   -3.1001984126984127e-6,
                                   3.0619243582206544e-7};
        for (int i = 0; i < 10; i++) {
            gReg += taylorCoeffs[i] * pow(arg, i);
        }
    } else if (arg == 0) {
        gReg = 1 / k;
    } else {
        gReg = pow(arg, k) *
               (egf_ugamma(-k, arg) + (pow(-1, k) / tgamma(k + 1)) * log(arg));
    }
    // subtract polynomial of order k due to free parameter lambda
    gReg -= pow(arg, k) * log(lambda * lambda);
    return gReg;
}

/**
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function,
 * that is d - nu.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda.
 * @return - gamma(s/2) * gammaStar(s/2, pi * prefactor * z**2),
 * where gammaStar is the twice regularized lower incomplete gamma function if s is
 * not equal to - 2k and (pi * prefactor * y ** 2) ** (- s / 2)
 * (gamma(s / 2, pi * prefactor * z ** 2) + ((-1)^k / k! ) * (log(pi * y ** 2) -
 * log(prefactor ** 2))) if s is  equal to - 2k for non negative natural number k.
 */
double complex crandall_gReg(unsigned int dim, double s, const double *z,
                             double prefactor) {
    double zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;
    double k = -nearbyint(s / 2.);
    if (s < 1 && (s == -2 * k)) {
        return crandall_gReg_nuequalsdimplus2k(s, zArgument, k, prefactor);
    }
    return -tgamma(s / 2) * egf_gammaStar(s / 2, zArgument);
}

/**
 * @brief Calculates bounds on when to use asymptotic expansion of the
 * upper incomplete gamma function, depending on the value of nu.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @return minimum value of z, when to use the fast asymptotic expansion in the
 * calculation of the incomplete upper gamma function upperGamma(nu, z).
 */
double assignzArgBound(double nu) {
    if (nu > -1 && nu < 5) {
        return M_PI * 3.15 * 3.15;
    }
    if (nu > -20 && nu < 20) {
        return M_PI * 3.35 * 3.35;
    }
    if (nu > -150 && nu < 60) {
        return M_PI * 3.5 * 3.5;
    }
    if (nu > -400 && nu < 100) {
        return M_PI * 3.65 * 3.65;
    }
    return DBL_MAX; // do not use expansion if nu is to big
}

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
 * |z| > 0 and - 2 / nu otherwise.
 */
double complex crandall_g(unsigned int dim, double nu, const double *z,
                          double prefactor, double zArgBound) {
    double zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;

    if (zArgument < EPS_ZERO_PIY) {
        return -2. / nu;
    }
    if (zArgument > zArgBound) {
        return exp(-zArgument) * (-2 + 2 * zArgument + nu) /
               (2 * zArgument * zArgument);
    }
    return egf_ugamma(nu / 2, zArgument) / pow(zArgument, nu / 2);
}

/**
 * @brief Calculates the lower Crandall function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast asymptotic
 * expansion in the calculation of the Crandall function.
 * @return lowerGamma(nu / 2,pi prefactor * z**2) / (pi * prefactor z**2)^(nu / 2) if
 * |z| > 0 and 2 / nu otherwise.
 */
double complex crandall_g_lower(unsigned int dim, double nu, const double *z,
                                double prefactor) {
    double zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;

    if (zArgument < EPS_ZERO_PIY) {
        return 2. / nu;
    }

    return tgamma(nu / 2) * egf_gammaStar(nu / 2, zArgument);
}

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
                    const unsigned int *beta) {
    double res = 1;
    unsigned int ai = 0;
    unsigned int bi = 0;
    unsigned long long factFrac;
    unsigned int aMinusb = 0;
    for (int i = 0; i < dim; i++) {
        ai = alpha[i];
        bi = beta[i];
        aMinusb += ai - bi;
        factFrac = 1; // Calculate (alpha - beta)!/(alpha - 2beta)!
        for (unsigned int j = ai - (2 * bi) + 1; j <= ai - bi; j++) {
            factFrac *= j;
        }
        res *= (double)binom((long long)ai, (long long)bi) * (double)factFrac *
               real_int_pow(2 * z[i], ai - (2 * bi));
    }

    res *= real_int_pow(-M_PI, aMinusb);

    return res;
}

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

/** @brief Computes cₙ,ᵢ,ₖ＝ ∏_{j=i+1}^{⌊n/2⌋-k} (2n＋d−2−4k−2j) as apint.
 * @param[in] n: index (total of alpha).
 * @param[in] i: lower bound of product (exclusive).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result apint.
 */
void coeffs_c_inner_apint(unsigned int n, unsigned int i, unsigned int k,
                          unsigned int dim, apint_t *out) {

    long long end = (n / 2) - k;

    apint_set_ull(out, 1, 1);
    if (i >= end) {
        apint_normalize(out);
        return;
    }
    int base = (2 * (int)n) + (int)dim - 2 - (4 * (int)k);
    for (unsigned int j = i + 1; j <= end; j++) {
        int factor = base - (2 * (int)j);
        apint_t tmp;
        apint_set_ull(&tmp, (unsigned long long)factor, 1);
        apint_mul(out, out, &tmp);
    }
}

/** @brief Computes (-2)**(-i) · c_{n,i,k} · binom(i+k,k) as apint.
 * @param[in] n: index (total of alpha).
 * @param[in] i: index (total of beta).
 * @param[in] k: index (specifies degree |alpha| - 2k).
 * @param[in] dim: dimension.
 * @param[out] out: result apint.
 */
void harmonic_h_inner_term_scalar_apint(unsigned int n, unsigned int i,
                                        unsigned int k, unsigned int dim,
                                        apint_t *out) {
    unsigned long long b = binom((long long)(i + k), (long long)k);
    apint_set_ull(out, b, 1);

    // c_{n,i,k}
    apint_t coeff;
    coeffs_c_inner_apint(n, i, k, dim, &coeff);
    apint_mul(out, out, &coeff);

    apint_normalize(out);

    // (-1/2)^i
    if (i & 1) {
        out->sign = -1;
    }
    out->exp2 -= (int)i;
}

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
                                       const unsigned int *theta2, apint_t *out) {
    apint_set_ull(out, 1, 1);

    // prod_j binom(betaAbsDim, beta[j])
    unsigned long long betaAbsDim = 0;
    for (unsigned int j = 0; j < dim; j++) {
        betaAbsDim += beta[j];
        unsigned long long b = binom((long long)betaAbsDim, (long long)beta[j]);
        apint_t tmp;
        apint_set_ull(&tmp, b, 1);
        apint_mul(out, out, &tmp);
    }

    // prod_j binom(alpha[j], theta1[j])
    for (unsigned int j = 0; j < dim; j++) {
        unsigned long long b = binom((long long)alpha[j], (long long)theta1[j]);
        apint_t tmp;
        apint_set_ull(&tmp, b, 1);
        apint_mul(out, out, &tmp);
    }

    // prod_j falling factorial: theta2[j]! / (theta2[j] - theta1[j])!
    for (unsigned int j = 0; j < dim; j++) {
        for (unsigned int l = theta2[j] - theta1[j] + 1; l <= theta2[j]; l++) {
            apint_t tmp;
            apint_set_ull(&tmp, (unsigned long long)l, 1);
            apint_mul(out, out, &tmp);
        }
    }

    apint_normalize(out);
}

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
                            const unsigned int *gamma, unsigned int alphaAbs) {

    unsigned int beta[dim];
    unsigned int theta1[dim];
    unsigned int theta2[dim];
    for (unsigned int j = 0; j < dim; j++) {
        beta[j] = 0;
        theta1[j] = alpha[j] + beta[j] - gamma[j];
        theta2[j] = gamma[j] - beta[j];
    }

    unsigned int redotheta1 = 0;
    unsigned int redotheta2 = 0;
    unsigned int betaAbs = 0;
    int done;
    int skip;

    unsigned int lastScalarIndex = (alphaAbs / 2) - k;

    // scalar coefficients as apint
    apint_t scalar_coeffs[lastScalarIndex + 1];
    for (unsigned int i = 0; i <= lastScalarIndex; i++) {
        harmonic_h_inner_term_scalar_apint(alphaAbs, i, k, dim, &scalar_coeffs[i]);
    }

    // multi-index coefficients as apint (accumulated via apint_add)
    apint_t multi_coeffs[lastScalarIndex + 1];
    for (unsigned int i = 0; i <= lastScalarIndex; i++) {
        apint_set_ull(&multi_coeffs[i], 0, 1);
    }

    // loop over beta multi-indices
    while (1) {

        skip = 0;
        for (unsigned int j = 0; j < dim; j++) {
            if (gamma[j] > alpha[j] && gamma[j] - alpha[j] > beta[j]) {
                skip = 1;
            }
        }

        if (!skip) {

            if (redotheta1) {
                for (unsigned int j = 0; j < dim; j++) {
                    theta1[j] = alpha[j] + beta[j] - gamma[j];
                }
            }
            redotheta1 = 0;

            if (redotheta2) {
                for (unsigned int j = 0; j < dim; j++) {
                    theta2[j] = gamma[j] - beta[j];
                }
            }
            redotheta2 = 0;

            apint_t term;
            harmonic_h_inner_term_multi_apint(dim, alpha, beta, theta1, theta2,
                                              &term);
            apint_add(&multi_coeffs[betaAbs], &multi_coeffs[betaAbs], &term);
        }

        done = 1;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (2 * beta[idx] + 2 <= 2 * gamma[idx] - alpha[idx]) {
                beta[idx]++;
                theta1[idx]++;
                betaAbs++;
                redotheta2 = 1;
                done = 0;
                break;
            }
            betaAbs -= beta[idx];
            theta2[idx] = beta[idx];
            beta[idx] = 0;
            redotheta1 = 1;
        }
        if (done) {
            break;
        }
    }

    // combine: sum_i scalar_coeffs[i] * multi_coeffs[i]
    apint_t sumInner;
    apint_set_ull(&sumInner, 0, 1);

    for (unsigned int i = 0; i <= lastScalarIndex; i++) {
        apint_t term;
        apint_mul(&term, &scalar_coeffs[i], &multi_coeffs[i]);
        apint_add(&sumInner, &sumInner, &term);
    }

    // convert to double
    double result = apint_to_double(&sumInner);

    return result;
}

/** @brief Updates the multi-index gamma in graded lexicographic order
 * while maintaining its total degree.
 * @param[in] maxAbs: maximum allowed total degree.
 * @param[in] dim: dimension of gamma.
 * @param[in,out] gamma: current multi-index to be updated.
 * @param[in,out] gammaAbs: total of gamma.
 * @return 1 if iteration is finished, 0 otherwise.
 */
static int harmonic_h_update_outer_index(unsigned int maxAbs, unsigned int dim,
                                         unsigned int *gamma,
                                         unsigned int *gammaAbs) {
    for (unsigned int idx = 0; idx < dim; idx++) {
        if (gamma[idx] + 1 <= maxAbs) {
            gamma[idx]++;
            (*gammaAbs)++;
            return 0;
        }
        (*gammaAbs) -= gamma[idx];
        gamma[idx] = 0;
    }
    return 1;
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
                                     double *coeffs, unsigned int *exponents) {
    unsigned int gamma[dim];
    unsigned int gammaAbs;
    int skip;
    int done;
    unsigned int k;
    unsigned int kMax = alphaAbs / 2;
    long long n;
    unsigned long long expIdx;
    unsigned int i;

    for (k = 0; k <= kMax; k++) {
        /* Initialize gamma and count n */
        for (i = 0; i < dim; i++) {
            gamma[i] = 0;
        }
        gammaAbs = 0;
        n = 0;
        while (1) {
            skip = 0;
            if (gammaAbs != alphaAbs - k) {
                skip = 1;
            }
            if (!skip) {
                for (i = 0; i < dim; i++) {
                    if (2 * gamma[i] < alpha[i]) {
                        skip = 1;
                        break;
                    }
                }
            }
            if (!skip) {
                coeffs[chunk_offset[k] + n] =
                    harmonic_h_inner_sum(k, dim, alpha, gamma, alphaAbs);
                expIdx = (chunk_offset[k] + n) * dim;
                for (i = 0; i < dim; i++) {
                    exponents[expIdx + i] = (2 * gamma[i]) - alpha[i];
                }
                n++;
            }
            done =
                harmonic_h_update_outer_index(alphaAbs - k, dim, gamma, &gammaAbs);
            if (done) {
                break;
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
    double zPow;
    double sumOuter = 0.0;
    double epsilonOuter = 0.0;
    double auxtOuter;
    double auxyOuter;
    unsigned long long n;
    unsigned long long count = valid_count[k];
    unsigned long long baseIdx = chunk_offset[k];
    unsigned long long expIdx;
    unsigned int i;

    for (n = 0; n < count; n++) {
        double sumInner = coeffs[baseIdx + n];
        expIdx = (baseIdx + n) * dim;
        zPow = 1.0;
        for (i = 0; i < dim; i++) {
            zPow *= real_int_pow(z[i], exponents[expIdx + i]);
        }
        auxyOuter = zPow * sumInner - epsilonOuter;
        auxtOuter = sumOuter + auxyOuter;
        epsilonOuter = (auxtOuter - sumOuter) - auxyOuter;
        sumOuter = auxtOuter;
    }
    sumOuter *= coeffs_c_outer(alphaAbs, k, dim);
    return sumOuter;
}

/** @brief Calculates the homogeneous harmonic polynomial h₍α,k₎
 * of degree |α|−2k such that y^α = ∑ₖ (y·y)^k h₍α,k₎(y) for d = 1 and k =
 * floor(|alpha|/2), explicitly, h₍α,k₎(y)=c_{|α|,k} ∑{|γ|=|α|−k} y^{2γ−α}
 * h_inner(α,γ,k). In 1D, h₍α,k₎(y) is y ** (alphaAbs mod 2) and zero otherwise. Uses
 * precomputed coefficients and exponents to avoid multi-index iteration.
 * @param[in] z: vector of the polynomial.
 * @param[in] alphaAbs: total of alpha.
 * @return h₍α,k₎(z).
 */
double harmonic_h_1D_kMax(const double *z, unsigned int alphaAbs) {
    if (alphaAbs % 2) {
        return z[0];
    }
    return 1.;
}
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
                              const unsigned int *alpha, unsigned int alphaAbs) {
    if (mult_abs(dim, alpha) == 0) {
        return crandall_g(dim, nu, z, prefactor, zArgBound);
    }

    unsigned int beta[dim];
    for (int i = 0; i < dim; i++) {
        beta[i] = 0;
    }

    double nuIt;
    double zArgBoundIt;

    int done = 0;
    unsigned int betaAbs = 0;

    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;

    // zArgBounds[i] is the zArgBound for nu = nu + 2 * alphaAbs - 2 * i;
    double zArgBounds[(alphaAbs / 2) + 1];
    for (int i = 0; i < (alphaAbs / 2) + 1; i++) {
        nuIt = nu + 2 * (double)alphaAbs - 2 * (double)i;
        zArgBounds[i] = assignzArgBound(nuIt);
    }

    // Iterate over every multi-index beta so that 2 beta <= alpha
    while (1) {

        nuIt = nu + 2 * alphaAbs - 2 * betaAbs;
        zArgBoundIt = zArgBounds[betaAbs]; // NOLINT

        // summing using Kahan's method
        // catch vanishing polynomials
        double p = polynomial_p(dim, z, alpha, beta);
        if (p) {
            auxy = p * crandall_g(dim, nuIt, z, prefactor, zArgBoundIt) - epsilon;
            auxt = sum + auxy;
            epsilon = (auxt - sum) - auxy;
            sum = auxt;
        }

        done = 1;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] + 1 <= alpha[idx] / 2) {
                beta[idx]++;
                betaAbs++;
                done = 0;
                break;
            }
            betaAbs -= beta[idx];
            beta[idx] = 0;
        }
        if (done) {
            break;
        }
    }

    return sum;
}

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
                    const unsigned int *beta) {

    double res = 1;
    unsigned int ai = 0;
    unsigned int bi = 0;
    unsigned long long factFrac;
    unsigned int aMinusb = 0;

    for (int i = 0; i < dim; i++) {
        ai = alpha[i];
        bi = beta[i];
        aMinusb += ai - bi;
        factFrac = 1; // Calculate (alpha - beta)!/(alpha - 2beta)!
        for (unsigned int j = ai - (2 * bi) + 1; j <= ai - bi; j++) {
            factFrac *= j;
        }
        res *= (double)binom((long long)ai, (long long)bi) * (double)factFrac *
               real_int_pow(2 * z[i], ai - (2 * bi));
    }

    unsigned long long factorial = 1;
    for (int j = 1; j < aMinusb; j++) {
        factorial *= j;
    }

    res *= (double)factorial;

    if (!(aMinusb % 2)) {
        res *= -1.;
    }

    return res;
}

/** @brief Calculates the derivatives of L(z) = log(pi * z**2)
 * @param[in] dim: dimension of z.
 * @param[in] z: vector of the polynomial.
 * @parma[in] alpha: multi-index for the derivative.
 * @parma[in] alphaAbs: absolute value of the multi-index alpha.
 * @return partial derivative of L(z).
 */
double complex log_l_der(unsigned int dim, const double *z,
                         const unsigned int *alpha, unsigned int alphaAbs) {

    double zArg = dot(dim, z, z);

    // Return function if there is no derivative
    if (!alphaAbs) {
        return log(M_PI * zArg);
    }

    unsigned int beta[dim];
    for (int i = 0; i < dim; i++) {
        beta[i] = 0;
    }

    unsigned int aMinusb = alphaAbs;

    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;

    int done = 0;

    // Iterate over every multi-index beta so that 2 beta <= alpha
    while (1) {

        // summing using Kahan's method
        auxy = polynomial_l(dim, z, alpha, beta) / real_int_pow(zArg, aMinusb) -
               epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;

        done = 1;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] + 1 <= alpha[idx] / 2) {
                beta[idx]++;
                aMinusb--;
                done = 0;
                break;
            }
            aMinusb += beta[idx];
            beta[idx] = 0;
        }
        if (done) {
            break;
        }
    }

    return sum;
}

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
                        unsigned int n) {

    // Return function if there is no derivative
    if (!alphaAbs) {
        return real_int_pow(M_PI * dot(dim, z, z), k);
    }

    unsigned int betaMin[dim];
    for (int i = 0; i < dim; i++) {
        betaMin[i] = (alpha[i] + 1) / 2;
    }

    unsigned int absMin = 0;
    for (int i = 0; i < dim; i++) {
        absMin += betaMin[i];
    }

    // Higher derivatives vanish
    if (absMin > k) {
        return 0.;
    }

    unsigned long long factMin = 1;
    for (int i = 0; i < dim; i++) {
        for (int j = 1; j < betaMin[i] + 1; j++) {
            factMin *= j;
        }
    }

    unsigned int beta[dim];
    for (int i = 0; i < dim; i++) {
        beta[i] = betaMin[i];
    }

    unsigned int betaAbs = absMin;
    unsigned long long betaFact = factMin;

    double sum = 0.;
    double epsilon = 0.;
    double auxt;
    double auxy;

    double summand;
    double res;

    unsigned int redoFact = 0;
    unsigned int done = 0;
    while (1) {

        if (!(betaAbs - k)) {

            // Recalculate factorial (expensive but stable)
            if (redoFact) {
                betaFact = factMin;
                for (unsigned int i = 0; i < dim; i++) {
                    for (unsigned int j = betaMin[i] + 1; j < beta[i] + 1; j++) {
                        betaFact *= j;
                    }
                }
                redoFact = 0;
            }

            summand = 1.;
            for (int i = 0; i < dim; i++) {
                summand *=
                    (double)binom((long long)(2 * beta[i]), (long long)alpha[i]) *
                    real_int_pow(z[i], (2 * beta[i]) - alpha[i]);
            }

            summand = summand / (double)betaFact;

            // summing using Kahan's method
            auxy = summand - epsilon;
            auxt = sum + auxy;
            epsilon = (auxt - sum) - auxy;
            sum = auxt;
        }

        done = 1;
        // Loop over multi-indexes beta so that 2 beta >= alpha and |beta| <= k
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] < k - absMin + betaMin[idx]) {
                // Fast path: increment current dimension
                beta[idx]++;
                betaAbs++;
                betaFact *= beta[idx];
                done = 0;
                break;
            }

            // Slow path: reset this dimension and continue
            betaAbs -= beta[idx] - betaMin[idx];
            beta[idx] = betaMin[idx];
            redoFact = 1;
        }
        if (done) {
            break;
        }
    }

    // factorial = alpha! k! / n!
    unsigned long long factorial = 1;
    for (unsigned int i = 0; i < dim; i++) {
        for (int j = 1; j < alpha[i] + 1; j++) {
            factorial *= j;
        }
    }
    for (unsigned int j = n + 1; j < k + 1; j++) {
        factorial *= j;
    }

    res = (double)factorial * real_int_pow(M_PI, k) * sum;

    return res;
}

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
                                 const unsigned int *alpha, unsigned int alphaAbs) {

    double res;

    unsigned long long kFact = 1;
    for (int j = 1; j < k + 1; j++) {
        kFact *= j;
    }

    double prefactor = pow(M_PI, (double)k + ((double)dim / 2.)) /
                       tgamma((double)k + ((double)dim / 2.)) *
                       ((k % 2) ? 1. : -1.) / (double)kFact;

    unsigned int beta[dim];
    unsigned int gamma[dim]; // gamma = alpha - beta
    for (int i = 0; i < dim; i++) {
        beta[i] = 0;
        gamma[i] = alpha[i];
    }

    // Return function if there is no derivative
    if (!alphaAbs) {
        double zArg = M_PI * dot(dim, z, z);
        res = prefactor * real_int_pow(zArg, k) * log(zArg);
        return res;
    }

    unsigned int betaAbs = 0;
    unsigned int gammaAbs = alphaAbs;

    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;

    unsigned int multBinom;

    int done = 0;
    // Iterate over every multi-index beta so that beta + gamma = alpha
    while (1) {

        if (betaAbs / 2 < k + dim) {
            multBinom = 1;
            for (int i = 0; i < dim; i++) {
                multBinom *= binom((long long)alpha[i], (long long)beta[i]);
            }

            // summing using Kahan's method
            auxy = (double)multBinom *
                       polynomial_y_der(k, dim, z, beta, betaAbs, 1) *
                       log_l_der(dim, z, gamma, gammaAbs) -
                   epsilon;
            auxt = sum + auxy;
            epsilon = (auxt - sum) - auxy;
            sum = auxt;
        }

        done = 1;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] < alpha[idx]) {
                betaAbs++;
                gammaAbs--;
                beta[idx]++;
                gamma[idx]--;
                done = 0;
                break;
            }
            betaAbs -= beta[idx];
            gammaAbs += beta[idx];
            gamma[idx] += beta[idx];
            beta[idx] = 0;
            done = 1;
        }
        if (done) {
            break;
        }
    }

    res = prefactor * sum;

    return res;
}

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
 * @return partial derivative of the regularized Crandall function.
 */
double complex crandall_gReg_nuequalsdimplus2k_der(
    double s, unsigned int k, unsigned int dim, const double *z, double lambda,
    const unsigned int *alpha, unsigned int alphaAbs, double zArgBound) {
    double res;
    double zArg = M_PI * dot(dim, z, z);
    double taylorBranchBound = M_PI * 0.65 * 0.65;
    unsigned int taylorSeriesLimit = 25;
    if (!alphaAbs) {
        // No derivative
        res = crandall_gReg_nuequalsdimplus2k(s, zArg, (double)k, lambda);
    } else if (zArg < taylorBranchBound) {
        // Taylor expansion if arg close to zero.

        double complex sum = 0.0;
        double complex epsilon = 0.0;
        double complex auxt;
        double complex auxy;

        double eulerGamma = 0.57721566490153286555;

        double harmonic = 0;
        for (int i = 1; i < k + 1; i++) {
            harmonic += 1. / (double)i;
        }

        res = ((k % 2) ? -1. : 1.) *
              (harmonic - eulerGamma -
               real_int_pow(lambda, 2 * k) * log(lambda * lambda)) *
              polynomial_y_der(k, dim, z, alpha, alphaAbs, k);

        // summand n = 0
        if (k) {
            res -= polynomial_y_der(0, dim, z, alpha, alphaAbs, 0) / (double)(-k);
        }

        // summands 0 < n < k
        for (int n = 1; n < taylorSeriesLimit + 1; n++) {
            if (n - k) {
                // Summing using Kahan's method
                auxy = ((n % 2) ? -1. : 1.) *
                           polynomial_y_der(n, dim, z, alpha, alphaAbs, n) /
                           (double)(n - (int)k) -
                       epsilon;
                auxt = sum + auxy;
                epsilon = (auxt - sum) - auxy;
                sum = auxt;
            }
        }

        res -= sum;

    } else {
        // Evaluate difference of
        res = crandall_g_der(dim, s, z, lambda, zArgBound, alpha, alphaAbs) -
              pow(M_PI, -((double)dim - s) / 2.) * tgamma(((double)dim - s) / 2.) *
                  singularity_s_der(k, dim, z, alpha, alphaAbs);
    }

    return res;
}

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
                                 unsigned int alphaAbs, double zArgBound) {
    double zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;
    unsigned int k = (unsigned int)(-nearbyint(s / 2.));
    if (s < 1 && (-s == 2 * k)) {
        return crandall_gReg_nuequalsdimplus2k_der(s, k, dim, z, prefactor, alpha,
                                                   alphaAbs, zArgBound);
    }

    unsigned int beta[dim];
    for (int i = 0; i < dim; i++) {
        beta[i] = 0;
    }

    double sIt;

    int done = 0;
    unsigned int betaAbs = 0;

    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;

    // Iterate over every multi-index beta so that 2 beta <= alpha
    while (1) {

        sIt = s + 2 * alphaAbs - 2 * betaAbs;

        // Summing using Kahan's method
        auxy = -polynomial_p(dim, z, alpha, beta) * tgamma(sIt / 2) *
                   egf_gammaStar(sIt / 2, zArgument) -
               epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;

        done = 1;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] + 1 <= alpha[idx] / 2) {
                beta[idx]++;
                betaAbs++;
                done = 0;
                break;
            }
            betaAbs -= beta[idx];
            beta[idx] = 0;
        }
        if (done) {
            break;
        }
    }

    return sum;
}

#undef EPS
#undef G_CUTOFF
