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
 * sum in Crandall's formula in the special case of
 * nu = dim + 2k for some natural number k.
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function.
 * @param[in] arg: input of the function
 * @param[in] k: k = - s / 2 = (nu - d) / 2 as an integer
 * @param[in] lambda: scaling parameter of crandalls formula
 * @return arg ** (- s / 2) * (gamma(s / 2, arg) + ((-1)^k / k! ) * (log(arg) -
 * log(lambda ** 2))
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
 * @param[in] dim: dimension of the input vectors
 * @param[in] s: dimension minus exponent of the regularized Epstein zeta function,
 * that is d - nu
 * @param[in] z: input vector of the function
 * @param[in] prefactor: prefactor of the vector, e. g. lambda
 * @return - gamma(s/2) * gammaStar(s/2, pi * prefactor * z**2),
 * where gammaStar is the twice regularized lower incomplete gamma function if s is
 * not equal to - 2k and (pi * prefactor * y ** 2) ** (- s / 2)
 * (gamma(s / 2, pi * prefactor * z ** 2) + ((-1)^k / k! ) * (log(pi * y ** 2) -
 * log(prefactor ** 2))) if s is  equal to - 2k for non negative natural number k
 */
double complex crandall_gReg(unsigned int dim, double s, const double *z,
                             double prefactor) {
    double zArgument = dot(dim, z, z);
    zArgument *= M_PI * prefactor * prefactor;
    double k = -(double)nearbyint(s / 2.);
    if (s < 1 && (s == -2 * k)) {
        return crandall_gReg_nuequalsdimplus2k(s, zArgument, k, prefactor);
    }
    return -tgamma(s / 2) * egf_gammaStar(s / 2, zArgument);
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
double complex crandall_g(unsigned int dim, double nu, const double *z,
                          double prefactor, double zArgBound) {
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
