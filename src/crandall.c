// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024-2026 Jonathan Busse <jonathan@jbusse.de>
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
#include "stdbool.h"
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
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula in the special case of
 * nu = dim + 2l for some natural number l >= k.
 * @param[in] l: integer so that s = d-nu = -2 l for non-negative integer l >= k.
 * @param[in] k: specifies degree |alpha| - 2k.
 * @param[in] n: |alpha|.
 * @param[in] arg: of the function, same as pi * z.z.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] z: input vector of the function.
 * @param[in] lambda: scaling parameter of crandalls formula.
 * @return Regularized upper Crandall function for logarithm special case.
 */
double complex crandall_gReg_nuequalsdimplus2l_harmonic(int l, int k, int n,
                                                        double arg, unsigned int dim,
                                                        const double *z,
                                                        double lambda) {
    double complex res = NAN;

    // difference of digamma functions ψ(ℓ+d/2)−ψ(ℓ+d/2−k)
    double digamma = 0;
    for (int i = 1; i <= k; i++) {
        digamma += 1. / ((double)(l - i) + .5 * ((double)dim));
    }

    // Employ series expansion only in (unstable) base case
    double taylorCutoff = 0.1 * 0.1 * M_PI;
    if (l == n - k && arg < taylorCutoff) {

        // harmonic number H
        double harmonic = 0;
        for (int i = 1; i <= n - k; i++) {
            harmonic += 1. / (double)i;
        }

        double eulerGamma = 0.57721566490153286555;

        // taylorCoeffs[i-1] = (-1)^i / (i * i!)
        double taylorCoeffs[9] = {1,
                                  -0.25,
                                  0.05555555555555555,
                                  -0.010416666666666666,
                                  0.0016666666666666668,
                                  -0.0002314814814814815,
                                  0.00002834467120181406,
                                  -3.1001984126984127e-6,
                                  3.0619243582206544e-7};

        res = harmonic + digamma - eulerGamma -
              2 * pow(lambda, 2 * (l - k)) * log(lambda);

        for (int i = 1; i < 10; i++) {
            res += taylorCoeffs[i - 1] * pow(arg, i);
        }

        return res;
    }

    double exp = l + k - n;

    double zArgBound = assignzArgBound(-2 * exp);
    res = crandall_g(dim, -2 * exp, z, lambda, zArgBound);

    // skip correction term in the origin
    if (arg) {

        // difference of harmonic numbers Hₗ−Hₗ₊ₖ₋ₙ
        double harmonic = 0;
        for (int i = 1; i <= n - k; i++) {
            harmonic += 1. / (exp + (double)i);
        }

        // factorial exp!
        double fact = 1;
        for (int i = 2; i <= exp; i++) {
            fact *= i;
        }

        res += negative_one_pow((unsigned int)exp) * pow(arg, exp) / fact *
               (log(arg) + harmonic + digamma);
    }

    return res;
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

/**
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula for the harmonic method.
 * @param[in] k: specifies degree |alpha| - 2k.
 * @param[in] n: |alpha|.
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function,
 * that is d - nu.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] z: input vector of the function.
 * @param[in] prefactor: prefactor of the vector, e. g. lambda.
 * @return - gamma(s/2) * gammaStar(s/2, pi * prefactor * z**2),
 * where gammaStar is the twice regularized lower incomplete gamma function if s is
 * and the special definition if s is  equal to - 2k for non negative natural number
 * k <= l.
 */
double complex crandall_gReg_harmonic(int k, int n, unsigned int dim, double s,
                                      const double *z, double prefactor) {

    double zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;

    // catch log term special treatment
    double l = -nearbyint(s / 2.);
    if (-s == 2 * l && l >= n - k) {
        return crandall_gReg_nuequalsdimplus2l_harmonic((int)l, k, n, zArgument, dim,
                                                        z, prefactor);
    }

    double exp = s + (2 * (n - k));

    return -crandall_g_lower(dim, exp, z, prefactor);
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
    for (unsigned int i = 0; i < dim; i++) {
        betaMin[i] = (alpha[i] + 1) / 2;
    }

    unsigned int absMin = 0;
    for (unsigned int i = 0; i < dim; i++) {
        absMin += betaMin[i];
    }

    // Higher derivatives vanish
    if (absMin > k) {
        return 0.;
    }

    unsigned long long factMin = 1;
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 1; j <= betaMin[i]; j++) {
            factMin *= j;
        }
    }

    unsigned int beta[dim];
    for (unsigned int i = 0; i < dim; i++) {
        beta[i] = betaMin[i];
    }
    unsigned int betaAbs = absMin;
    unsigned long long betaFact = factMin;

    // cutoff is k - absMin per dimension, totalCount = (k - absMin + 1)^dim
    unsigned int range = k - absMin;
    unsigned long long totalCount = 1;
    for (unsigned int j = 0; j < dim; j++) {
        totalCount *= range + 1;
    }

    double sum = 0.;
    double epsilon = 0.;
    bool redoFact = false;

    // Iterate over multi-indices beta so that 2 * beta >= alpha and |beta| <= k
    for (unsigned long long i = 0; i < totalCount; i++) {

        if (betaAbs == k) {
            // Recalculate factorial (expensive but stable)
            if (redoFact) {
                betaFact = factMin;
                for (unsigned int j = 0; j < dim; j++) {
                    for (unsigned int l = betaMin[j] + 1; l <= beta[j]; l++) {
                        betaFact *= l;
                    }
                }
                redoFact = false;
            }

            double summand = 1.;
            for (unsigned int j = 0; j < dim; j++) {
                summand *=
                    (double)binom((long long)(2 * beta[j]), (long long)alpha[j]) *
                    real_int_pow(z[j], (2 * beta[j]) - alpha[j]);
            }
            summand /= (double)betaFact;
            kahan_add_r(&sum, &epsilon, summand);
        }

        for (unsigned int idx = 0; idx < dim; idx++) {
            if (beta[idx] < betaMin[idx] + range) {
                // Fast path: increment current dimension
                beta[idx]++;
                betaAbs++;
                betaFact *= beta[idx];
                break;
            }
            // Slow path: reset this dimension and carry
            betaAbs -= beta[idx] - betaMin[idx];
            beta[idx] = betaMin[idx];
            redoFact = true;
        }
    }

    // factorial = alpha! * k! / n!
    unsigned long long factorial = 1;
    for (unsigned int i = 0; i < dim; i++) {
        for (unsigned int j = 1; j <= alpha[i]; j++) {
            factorial *= j;
        }
    }
    for (unsigned int j = n + 1; j <= k; j++) {
        factorial *= j;
    }

    return (double)factorial * real_int_pow(M_PI, k) * sum;
}

#undef EPS
#undef G_CUTOFF
