// SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file wrappers.c
 * @brief Wrappers of functions that are included in unit tests but do not appear as
 * is in the core library.
 */

#include "wrappers.h"
#include "../src/crandall.h"
#include "../src/harmonics.h"
#include "../src/hpdyad.h"
#include "../src/tools.h"
#include "../src/zeta.h"
#include "complex.h"
#include "epsteinZeta.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

/**
 * @brief calculates the derivatives of the set zeta function for lattices.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 * @param[in] alpha: multiindex for the derivative of the set zeta function.
 * @return function value of the Epstein zeta.
 */
double complex setZetaDer(double nu, unsigned int dim, const double *a,
                          const double *x, const double *y,
                          const unsigned int *alpha) {
    unsigned int alphaAbs = mult_abs(dim, alpha);
    double complex prefactor = imaginary_int_pow(alphaAbs) *
                               real_int_pow(-2 * M_PI, alphaAbs) *
                               cexp(2 * M_PI * I * dot(dim, x, y));
    return prefactor * epsteinZetaInternal(nu, dim, a, x, y, 1, 2, alpha);
}

/**
 * @brief calculates the derivatives of the regularized Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 * @param[in] alpha: multiindex for the derivative of the set zeta function.
 * @return function value of the Epstein zeta.
 */
double complex epsteinZetaRegDer(double nu, unsigned int dim, const double *a,
                                 const double *x, const double *y,
                                 const unsigned int *alpha) {
    unsigned int alphaAbs = mult_abs(dim, alpha);
    double complex prefactor =
        imaginary_int_pow(alphaAbs) * real_int_pow(-2 * M_PI, alphaAbs);
    return prefactor * epsteinZetaAnisoReg(nu, dim, a, x, y, alpha);
}

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
double singularity_s(double nu, unsigned int dim, double zSquared) {

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

/** @brief Calculates the derivatives of Y_ell(z) = (pi * z**2)**ell.
 * @param[in] ell: integer power.
 * @param[in] dim: dimension of z.
 * @param[in] z: vector z of the polynomial.
 * @param[in] n: |alpha| .
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return partial derivative of Y_ell(z).
 */
static double polynomial_y_der_harmonic(int ell, unsigned int dim,
                                        const double *z, // NOLINT
                                        int n,
                                        const unsigned long long *chunk_offset,
                                        const unsigned long long *valid_count,
                                        const double *coeffs,
                                        const unsigned int *exponents) {

    double zSquared = dot(dim, z, z);
    double res = 0;

    for (int k = 0; k <= n / 2; k++) {
        double h =
            harmonic_h(k, dim, z, n, chunk_offset, valid_count, coeffs, exponents);

        // falling pochhammer symbol (ell)_(n-k)(ell + dim/2 -1)_k
        double poch = 1;
        for (int i = 0; i < n - k; i++) {
            poch *= ell - i;
        }
        for (int i = 1; i <= k; i++) {
            poch *= ell + ((double)dim / 2.) - i;
        }

        if (h && poch) {
            res += h * poch * real_int_pow(zSquared, ell + k - n);
        }
    }

    res *= real_int_pow(M_PI, ell) * ldexp(1., n);

    return res;
}

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

    double res = polynomial_y_der_harmonic(ell, dim, z, (int)alphaAbs, chunk_offset,
                                           valid_count, coeffs, exponents);

    free(chunk_offset);
    free(valid_count);
    free(exponents);

    return res;
}

/** @brief Computes the derivatives of the singularity s_ν(z) where the coefficients
 * for the harmonic polynomials are provided as arguments. For the generic case (ν ≠
 * d + 2ℓ for any non-negative integer ℓ): s_ν(z) = π^(ν/2) / Γ(ν/2) · Γ((d−ν)/2) ·
 * (π z·z)^((ν−d)/2) and for the logarithmic case (ν = d + 2ℓ, ℓ ≥ 0): s_{d+2ℓ}(z) =
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
                // pochhamer symbol (nu/2-k-1)_m
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

    double res = singularity_s_der_harmonic(nu, dim, z, alphaAbs, chunk_offset,
                                            valid_count, coeffs, exponents);

    free(chunk_offset);
    free(valid_count);
    free(exponents);

    return res;
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
static double polynomial_l(unsigned int dim, const double *z,
                           const unsigned int *alpha, const unsigned int *beta) {

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
static double log_l_der(unsigned int dim, const double *z, const unsigned int *alpha,
                        unsigned int alphaAbs) {
    double zArg = dot(dim, z, z);
    // Return function if there is no derivative
    if (!alphaAbs) {
        return log(M_PI * zArg);
    }

    unsigned int beta[dim];
    for (unsigned int j = 0; j < dim; j++) {
        beta[j] = 0;
    }
    unsigned int aMinusb = alphaAbs;

    unsigned long long totalCount = 1;
    for (unsigned int j = 0; j < dim; j++) {
        totalCount *= alpha[j] / 2 + 1;
    }

    double sum = 0.0;
    double epsilon = 0.0;

    // Iterate over every multi-index beta so that 2 * beta <= alpha
    for (unsigned long long i = 0; i < totalCount; i++) {

        double complex summand =
            polynomial_l(dim, z, alpha, beta) / real_int_pow(zArg, aMinusb);
        kahan_add_r(&sum, &epsilon, summand);

        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] + 1 <= alpha[idx] / 2) {
                beta[idx]++;
                aMinusb--;
                break;
            }
            aMinusb += beta[idx];
            beta[idx] = 0;
        }
    }

    return sum;
}

/** @brief Computes p_(alpha,beta)(y) = (alpha choose beta)
 * ((alpha - beta)! / (alpha - 2*beta)!) * (2*y)^(alpha - 2*beta)
 * where 2 beta =< alpha
 * @param[in] dim: dimension of alpha, beta and y.
 * @param[in] y: vector of the polynomial.
 * @parma[in] alpha: upper multi-index.
 * @parma[in] beta: lower multi-index.
 * @return p(y).
 */
static double polynomial_p(unsigned int dim, const double *z,
                           const unsigned int *alpha, const unsigned int *beta) {

    double res = 1;
    unsigned int ai = 0;
    unsigned int bi = 0;
    unsigned long long factFrac;

    for (int i = 0; i < dim; i++) {
        ai = alpha[i];
        bi = beta[i];
        factFrac = 1; // Calculate (alpha - beta)!/(alpha - 2beta)!
        for (unsigned int j = ai - (2 * bi) + 1; j <= ai - bi; j++) {
            factFrac *= j;
        }
        res *= (double)binom((long long)ai, (long long)bi) * (double)factFrac *
               real_int_pow(2 * z[i], ai - (2 * bi));
    }

    return res;
}

/** @brief Calculates the derivatives of Y_nu(z) = (pi * z**2)**nu for real nu
 * where nu != dim + 2*k for any non-negative integer k.
 * @param[in] k: integer power.
 * @param[in] dim: dimension of z.
 * @param[in] y: vector of the polynomial.
 * @parma[in] alpha: multi-index for the derivative.
 * @param[in] n: factorial divisor smaller than k.
 * @return partial derivative of Y_k(z) / n!.
 */
static double polynomial_y_der_general(double nu, unsigned int dim,
                                       const double *z, // NOLINT
                                       const unsigned int *alpha,
                                       unsigned int alphaAbs) {

    double zArg = dot(dim, z, z);

    unsigned int beta[dim];
    for (unsigned int j = 0; j < dim; j++) {
        beta[j] = 0;
    }
    unsigned int aMinusb = alphaAbs;

    unsigned long long totalCount = 1;
    for (unsigned int j = 0; j < dim; j++) {
        totalCount *= alpha[j] / 2 + 1;
    }

    double sum = 0.0;
    double epsilon = 0.0;

    // Iterate over every multi-index beta so that 2 * beta <= alpha
    for (unsigned long long i = 0; i < totalCount; i++) {

        double poch = 1;
        for (int j = 0; j < aMinusb; j++) {
            poch *= nu - j;
        }

        double complex summand = poch * polynomial_p(dim, z, alpha, beta) /
                                 pow(zArg, (double)aMinusb - nu);
        kahan_add_r(&sum, &epsilon, summand);

        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] + 1 <= alpha[idx] / 2) {
                beta[idx]++;
                aMinusb--;
                break;
            }
            aMinusb += beta[idx];
            beta[idx] = 0;
        }
    }

    return pow(M_PI, nu) * sum;
}

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
                         const unsigned int *alpha, unsigned int alphaAbs) {

    // only works in this case
    unsigned int k = (unsigned int)fmax(0., nearbyint((nu - (double)dim) / 2));

    if (!(nu == (dim + 2 * (double)k))) {
        return pow(M_PI, nu / 2.) * tgamma(((double)dim - nu) / 2.) /
               tgamma(nu / 2.) *
               polynomial_y_der_general((nu - (double)dim) / 2., dim, z, alpha,
                                        alphaAbs);
    }

    // derivatives for special case nu = dim + 2* k
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
        return prefactor * real_int_pow(zArg, k) * log(zArg);
    }

    unsigned int betaAbs = 0;
    unsigned int gammaAbs = alphaAbs;

    double sum = 0.0;
    double epsilon = 0.0;

    unsigned int multBinom;

    bool done = false;
    // Iterate over every multi-index beta so that beta + gamma = alpha
    while (1) {

        if (betaAbs / 2 < k + dim) {
            multBinom = 1;
            for (int i = 0; i < dim; i++) {
                multBinom *= binom((long long)alpha[i], (long long)beta[i]);
            }
            double summand = (double)multBinom *
                             polynomial_y_der(k, dim, z, beta, betaAbs, 1) *
                             log_l_der(dim, z, gamma, gammaAbs);
            kahan_add_r(&sum, &epsilon, summand);
        }

        done = true;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] < alpha[idx]) {
                betaAbs++;
                gammaAbs--;
                beta[idx]++;
                gamma[idx]--;
                done = false;
                break;
            }
            betaAbs -= beta[idx];
            gammaAbs += beta[idx];
            gamma[idx] += beta[idx];
            beta[idx] = 0;
            done = true;
        }

        if (done) {
            break;
        }
    }

    return prefactor * sum;
}
