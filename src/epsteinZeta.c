// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file epsteinZeta.c
 * @brief Calculates the (regularized) Epstein zeta function.
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @see Crandall, R., Unified algorithms for polylogarithm, L-series, and zeta
 * variants. Algorithmic Reflections: Selected Works. PSIpress (2012).
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @date 06/13/2024
 */

#include <complex.h>
#include <math.h>

#include "crandall.h"
#include "epsteinZeta.h"
#include "gamma.h"
#include "tools.h"
#include "zeta.h"

/**
 * @brief calculates the Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 * @return function value of the Epstein zeta.
 */
double complex epsteinZeta(double nu, unsigned int dim, const double *a,
                           const double *x, const double *y) {
    return epsteinZetaInternal(nu, dim, a, x, y, 1, 0, (unsigned int[]){0});
}

/**
 * @brief calculates the regularized Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 * @return function value of the regularized Epstein zeta.
 */
double complex epsteinZetaReg(double nu, unsigned int dim, const double *a,
                              const double *x, const double *y) {
    return epsteinZetaInternal(nu, dim, a, x, y, 1, 1, (unsigned int[]){0});
}

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

    return epsteinZetaInternal(nu, dim, a, x, y, 1, 2, alpha);
}

/**
 * @brief calculates the derivatives of the regularized Epstein zeta function for
 * lattices.
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

    return epsteinZetaInternal(nu, dim, a, x, y, 1, 3, alpha);
}

/**
 * @brief Calculates the incomplete bessel function.
 * @param[in] nu: exponent of the function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] k: input vector of the function.
 * @param[in] r: input vector of the function.
 * @return 2 int_0**1 t**(-nu - 1) exp(-pi k**2 / t**2) exp(-pi r**2 t**2) dt
 */
double incomplete_bessel_g(double nu, unsigned int dim, const double *k,
                           const double *r) {

    double eps = 1e-32;

    double s = -nu / 2.;
    double x = M_PI * dot(dim, k, k);
    double y = M_PI * dot(dim, r, r);

    // Vanishing arguments
    if (x + y < eps) {
        return s;
    }

    // Vanishing first argument
    if (x < eps) {
        // Lower Crandall function
        return tgamma(s) * egf_gammaStar(s, y);
    }

    // Vanishing second argument
    if (y < eps) {
        return crandall_g(dim, nu, k, 1., assignzArgBound(nu));
    }

    // Swap not implemented yet, since K is missing
    //    // Reflect parameters for upper half-plane
    //    bool swap = (x + 0.1 < y);
    //    if (swap) {
    //        s = -s;
    //        double z = x;
    //        x = y;
    //        y = z;
    //    }

    double result = 0.0;

    // Choose series expansion close to origin
    if (x + y < 1.5) {
        result = pow(x, s) * egf_ugamma(-s, x);
        unsigned long long fact = 1;
        for (int j = 1; j <= 20; j++) {
            fact *= j;
            result += pow(x, s + j) * egf_ugamma(-s - j, x) * int_pow(-y, j) /
                      (double)fact;
        }
    } else {
        // Recursive algorithm away from origin

        // Initialize numerators
        double n1 = 0.0;
        double n2 = 0.0;
        double n3 = 1.0;

        // Initialize denominators
        double d1 = 0.0;
        double d2 = exp(x + y);
        double d3 = (x - y + s + 1.0) * d2;

        // Initialize final fraction
        double N = 0.;
        double D = 0.;

        for (int j = 2; j <= 100; j++) {
            N = ((x - y + s + 1 + 2 * (j - 1)) * n3 + (2 * y - s - (j - 1)) * n2 -
                 y * n1) /
                (double)j;
            D = ((x - y + s + 1 + 2 * (j - 1)) * d3 +
                 ((2 * y - s - (j - 1)) * d2 - y * d1)) /
                (double)j;

            n1 = n2;
            n2 = n3;
            n3 = N;

            d1 = d2;
            d2 = d3;
            d3 = D;
        }
        result = N / D;
    }

    //    // Reflect result for upper half-plane
    //    if (swap) {
    //        result = 2.0 * pow(x / y, s / 2.0) * bessel_K(-s, 2.0 * sqrt(x * y)) -
    //        result;
    //    }

    return result;
}
