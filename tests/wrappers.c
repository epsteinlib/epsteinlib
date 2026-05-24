// SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file wrappers.c
 * @brief Wrappers of functions that are included in unit tests but do not appear as
 * is in the core library.
 */

#include "../src/crandall.h"
#include "../src/harmonics.h"
#include "../src/hpdyad.h"
#include "../src/tools.h"
#include <math.h>
#include <stdlib.h>

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
                             const unsigned int *theta2) {
    hpdyad_t numerator;
    harmonic_h_inner_term_scalar_hpdyad(n, i, k, dim, &numerator);

    hpdyad_t multi_term;
    harmonic_h_inner_term_multi_hpdyad(dim, alpha, beta, theta1, theta2,
                                       &multi_term);

    hpdyad_mul(&numerator, &multi_term, &numerator);

    double res = hpdyad_to_double(&numerator);

    return res;
}

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
static double singularity_s(double nu, unsigned int dim, double zSquared) {

    double zArg = M_PI * zSquared;
    int ell = (int)nearbyint((nu - (double)dim) / 2.);

    if (ell >= 0 && nu == (double)(dim + (2 * ell))) {
        double fact = 1;
        for (int i = 2; i <= ell; i++) {
            fact *= (double)i;
        }

        return pow(M_PI, nu / 2) / tgamma(nu / 2) * negative_one_pow(ell + 1) /
               fact * real_int_pow(zArg, ell) * log(zArg);
    }

    return pow(M_PI, nu / 2) / tgamma(nu / 2) * tgamma(((double)dim - nu) / 2.) *
           pow(zArg, (nu - (double)dim) / 2);
}

/** @brief Computes the derivatives singularity s_ν(z) where the coefficients for the
 * harmonic polynomials are provided as arguments. For the generic case (ν ≠ d + 2ℓ
 * for any non-negative integer ℓ): s_ν(z) = π^(ν/2) / Γ(ν/2) · Γ((d−ν)/2) · (π
 * z·z)^((ν−d)/2) and for the logarithmic case (ν = d + 2ℓ, ℓ ≥ 0): s_{d+2ℓ}(z) =
 * π^(ℓ+d/2) / Γ(ℓ+d/2) · (−1)^(ℓ+1) / ℓ! · (π z·z)^ℓ · log(π z·z)
 * @param[in] nu: epxonent of singularity.
 * @param[in] dim: dimension of z.
 * @param[in] z: vector of the singularity.
 * @parma[in] alphaAbs: absolute value of the multi-index alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return partial derivative of s_nu(z).
 */
static double singularity_s_der_harmonic(double nu, unsigned int dim,
                                         const double *z, unsigned int alphaAbs,
                                         const unsigned long long *chunk_offset,
                                         const unsigned long long *valid_count,
                                         const double *coeffs,
                                         const unsigned int *exponents) {

    double res = 0;
    double zSquared = dot(dim, z, z);

    if (zSquared == 0) {
        return 0;
    }

    // Compute set zeta derivatives by harmonic method
    for (int k = 0; k <= alphaAbs / 2; k++) {

        double h = harmonic_h(k, dim, z, alphaAbs, chunk_offset, valid_count, coeffs,
                              exponents);

        int ell = (int)nearbyint((nu - (double)dim) / 2.);
        if (h) {

            int m = (int)alphaAbs - (2 * k);

            if (ell >= m + k && nu == (double)(dim + (2 * ell))) {
                double poch = 1;
                for (int i = 1; i <= m; i++) {
                    poch *= nu / 2 - k - i;
                }
                // difference of digamma functions ψ(ℓ+d/2)−ψ(ℓ+d/2−k)
                double digamma = 0;
                for (int i = 1; i <= k; i++) {
                    digamma += 1. / ((double)(ell - i) + ((double)dim) / 2.);
                }

                // difference of harmonic numbers Hₗ−Hₗ₊ₖ₋ₙ for n = alphaAbs
                double harmonic = 0;
                for (int i = 1; i <= (int)alphaAbs - k; i++) {
                    harmonic += 1. / ((double)(ell + k - (int)alphaAbs) + (double)i);
                }

                res += h * ldexp(1., k) * real_int_pow(-2 * M_PI * M_PI, m + k) /
                       poch *
                       singularity_s(nu - (2 * (double)(m + k)), dim, zSquared) *
                       (1 + (harmonic + digamma) / log(M_PI * zSquared));

                // skip zeros 1/gamma(s) of the reciprocal gamma function, where s =
                // nu/2 -k is an even, non-positive integer
            } else if (!((nu / 2.) - (double)k <= 0. &&
                         nearbyint((nu / 2.) - (double)k) ==
                             (nu / 2.) - (double)k)) {

                res += h * tgamma((((double)dim - nu) / 2.) + k + m) /
                       tgamma((nu / 2.) - (double)k) * negative_one_pow(k + m) *
                       ldexp(1., (2 * k) + m) * pow(M_PI, nu - ((double)dim / 2)) *
                       pow(dot(dim, z, z), ((nu - (double)dim) / 2.) - k - m);
            }
        }
    }

    return res;
}

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
double singularity_s_der_harmonic_wrapper(double nu, unsigned int dim,
                                          const double *z, const unsigned int *alpha,
                                          unsigned int alphaAbs) {

    // Precompute coefficients for harmonic polynomials
    unsigned int kMax = alphaAbs / 2;
    unsigned long long *chunk_offset =
        malloc((kMax + 1) * sizeof(unsigned long long));
    unsigned long long *valid_count =
        malloc((kMax + 1) * sizeof(unsigned long long));
    unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
        alphaAbs, kMax, dim, alpha, chunk_offset, valid_count);
    double *coeffs = malloc(coeffs_size * sizeof(double));
    unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
    precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_offset, coeffs,
                                    exponents);

    return singularity_s_der_harmonic(nu, dim, z, alphaAbs, chunk_offset,
                                      valid_count, coeffs, exponents);
}
