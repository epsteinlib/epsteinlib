// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
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
#include <math.h>

/*!
 * @brief epsilon for the cutoff around nu = dimension.
 */
#define EPS ldexp(1, -30)

/**
 * @brief Calculates the regularization of the zero summand in the second
 * sum in Crandall's formula. That is, the summand for k = 0.
 * @param[in] dim: dimension of the input vectors
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function
 * @param[in] prefactor: prefactor of the vector, e. g. lambda
 * @return gamma(nu/2) * gammaStar(nu/2, pi * prefactor * z**2), where
 * gammaStar is the twice regularized lower incomplete gamma function
 * gamma(s,x) / (gamma(s) * x ** s)
 */
double complex crandall_gReg(int dim, double nu, const double *z, double prefactor) {
    double complex zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;
    return -tgamma(nu / 2) * egf_gammaStar(nu / 2, zArgument);
}

/**
 * @brief calculates bounds on when to use asymptotic expansion of the
 * upper incomplete gamma function, depending on the value of nu.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @return minimum value of z, when to use the fast asymptotic expansion in the
 * calculation of the incomplete upper gamma function upperGamma(nu, z).
 */
double assignzArgBound(double nu) {
    if ((nu > 2 - EPS && nu < 2 + EPS) || (nu > 4 - EPS && nu < 4 + EPS)) {
        return M_PI * 2.6 * 2.6;
    }
    if (nu > 1.6 && nu < 4.4) {
        return M_PI * 2.99 * 2.99;
    }
    if (nu > -3 && nu < 8) {
        return M_PI * 3.15 * 3.15;
    }
    if (nu > -70 && nu < 40) {
        return M_PI * 3.35 * 3.35;
    }
    if (nu > -600 && nu < 80) {
        return M_PI * 3.5 * 3.5;
    }
    return pow(10, 16); // do not use expansion if nu is to big
}

/**
 * @brief Assumes x and y to be in the respective elementary lattice cell.
 * Multiply with exp(2 * PI * i * x * y) to get the second sum in Crandall's
 * @param[in] dim: dimension of the input vectors.
 * @param[in] nu: exponent of the regularized Epstein zeta function.
 * @param[in] z: input vector of the function
 * @param[in] prefactor: prefactor of the vector, e. g. lambda or 1/lambda in
 *      Crandall's formula
 * @return upperGamma(nu/2,pi prefactor * z**2)
 *      / (pi * prefactor z**2)^(nu / 2) in
 */
double complex crandall_g(int dim, double nu, const double *z, double prefactor,
                          double zArgBound) {
    double zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;

    if (zArgument < ldexp(1, -62)) {
        return -2. / nu;
    }
    if (zArgument > zArgBound) {
        return exp(-zArgument) * (-2 + 2 * zArgument + nu) /
               (2 * zArgument * zArgument);
    }
    return egf_ugamma(nu / 2, zArgument) / pow(zArgument, nu / 2);
}
#undef EPS
#undef G_CUTOFF
