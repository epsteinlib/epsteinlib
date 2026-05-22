// SPDX-FileCopyrightText: 2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file wrappers.c
 * @brief Wrappers of functions that are included in unit tests but do not appear as
 * is in the core library.
 */

#include "../src/crandall.h"
#include "../src/gamma.h"
#include "../src/harmonics.h"
#include "../src/hpdyad.h"
#include "../src/tools.h"
#include <complex.h>
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
    for (unsigned int j = 0; j < dim; j++) {
        beta[j] = 0;
    }
    unsigned int betaAbs = 0;

    // zArgBounds[i] is the zArgBound for nu = nu + 2 * alphaAbs - 2 * i
    double zArgBounds[(alphaAbs / 2) + 1];
    for (unsigned int i = 0; i <= alphaAbs / 2; i++) {
        zArgBounds[i] =
            assignzArgBound(nu + (2 * (double)alphaAbs) - (2 * (double)i));
    }

    unsigned long long totalCount = 1;
    for (unsigned int j = 0; j < dim; j++) {
        totalCount *= alpha[j] / 2 + 1;
    }

    double complex sum = 0.0;
    double complex epsilon = 0.0;

    // Iterate over every multi-index beta so that 2 * beta <= alpha
    for (unsigned long long i = 0; i < totalCount; i++) {

        double p = polynomial_p(dim, z, alpha, beta);
        if (p) {
            double nuIt = nu + (2 * alphaAbs) - (2 * betaAbs);
            double zArgBoundIt = zArgBounds[betaAbs]; // NOLINT
            double complex summand =
                p * crandall_g(dim, nuIt, z, prefactor, zArgBoundIt);
            kahan_add_c(&sum, &epsilon, summand);
        }

        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] + 1 <= alpha[idx] / 2) {
                beta[idx]++;
                betaAbs++;
                break;
            }
            betaAbs -= beta[idx];
            beta[idx] = 0;
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
    for (unsigned int j = 0; j < dim; j++) {
        beta[j] = 0;
    }
    unsigned int aMinusb = alphaAbs;

    unsigned long long totalCount = 1;
    for (unsigned int j = 0; j < dim; j++) {
        totalCount *= alpha[j] / 2 + 1;
    }

    double complex sum = 0.0;
    double complex epsilon = 0.0;

    // Iterate over every multi-index beta so that 2 * beta <= alpha
    for (unsigned long long i = 0; i < totalCount; i++) {

        double complex summand =
            polynomial_l(dim, z, alpha, beta) / real_int_pow(zArg, aMinusb);
        kahan_add_c(&sum, &epsilon, summand);

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
double singularity_s_der(double nu, unsigned int dim, const double *z,
                         const unsigned int *alpha, unsigned int alphaAbs) {

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

/**
 * @brief Calculates the derivatives of regularization of the zero summand in the
 * second sum in Crandall's formula in the special case of nu = dim + 2k for some
 * natural number k.
 * @param[in] s: s = (d - nu).
 * @param[in] k: k = (nu - d) / 2 as an integer.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] z: input vector of the function.
 * @param[in] lambda: scaling parameter of crandalls formula.
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast
 * asymptotic
 * @return partial derivative of the regularized Crandall function. */
double crandall_gReg_nuequalsdimplus2k_der(double s, unsigned int k,
                                           unsigned int dim, const double *z,
                                           double lambda, const unsigned int *alpha,
                                           unsigned int alphaAbs, double zArgBound) {
    double res;
    double zArg = M_PI * dot(dim, z, z);
    double taylorBranchBound = M_PI * 0.65 * 0.65;
    unsigned int taylorSeriesLimit = 25;
    if (!alphaAbs) {
        // No derivative
        res = crandall_gReg_nuequalsdimplus2k(s, zArg, (double)k, lambda);
    } else if (zArg < taylorBranchBound) {
        // Taylor expansion if arg close to zero.

        double sum = 0.0;
        double epsilon = 0.0;

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
                double summand = ((n % 2) ? -1. : 1.) *
                                 polynomial_y_der(n, dim, z, alpha, alphaAbs, n) /
                                 (double)(n - (int)k);
                kahan_add_r(&sum, &epsilon, summand);
            }
        }

        res -= sum;

    } else {
        // Evaluate difference of
        res = crandall_g_der(dim, s, z, lambda, zArgBound, alpha, alphaAbs) -
              pow(M_PI, -((double)dim - s) / 2.) * tgamma(((double)dim - s) / 2.) *
                  singularity_s_der(dim + (2 * k), dim, z, alpha, alphaAbs);
    }

    return res;
}

/**
 * @brief Calculates the derivatives of the regularization of the zero summand in
 * the second sum in Crandall's formula.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta
 * function, that is d - nu.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda.
 * @param[in] alpha: multi-index of the partial derivatives
 * @param[in] alphaAbs: sum of the elements of alpha
 * @param[in] zArgBound: minimum value of pi * z**2, when to use the fast
 * asymptotic
 * @return partial derivatives of - gamma(s/2) * gammaStar(s/2, pi * prefactor *
 * z**2), where gammaStar is the twice regularized lower incomplete gamma
 * function if s is not equal to - 2k and partial derivatives of (pi * prefactor
 * * y ** 2) ** (- s / 2) (gamma(s / 2, pi * prefactor * z ** 2) + ((-1)^k / k! )
 * * (log(pi * y ** 2)
 * - log(prefactor ** 2))) if s is  equal to - 2k for non negative natural number
 * k.
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
    for (unsigned int j = 0; j < dim; j++) {
        beta[j] = 0;
    }
    unsigned int betaAbs = 0;

    unsigned long long totalCount = 1;
    for (unsigned int j = 0; j < dim; j++) {
        totalCount *= alpha[j] / 2 + 1;
    }

    double complex sum = 0.0;
    double complex epsilon = 0.0;

    // Iterate over every multi-index beta so that 2 * beta <= alpha
    for (unsigned long long i = 0; i < totalCount; i++) {

        double sIt = s + (2 * alphaAbs) - (2 * betaAbs);

        double complex summand = -polynomial_p(dim, z, alpha, beta) *
                                 tgamma(sIt / 2) * egf_gammaStar(sIt / 2, zArgument);

        kahan_add_c(&sum, &epsilon, summand);

        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] + 1 <= alpha[idx] / 2) {
                beta[idx]++;
                betaAbs++;
                break;
            }
            betaAbs -= beta[idx];
            beta[idx] = 0;
        }
    }

    return sum;
}
