// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/*
 * @file zeta_mathematica.c
 * @brief Epstein Zeta function call for Mathematica interface
 */

#include "zeta_mathematica.h"
#include "epsteinZeta.h"
#include <complex.h>

/*
 * @brief Saves real and imaginary parts of the Epstein Zeta function and
 * in an array.
 * @param[in, out] out: 2D array, where real and imaginary parts of the Epstein Zeta
 * function will be stored.
 *  @param[in] nu: Exponent for the Epstein zeta function.
 *  @param[in] dim: Dimension of the input vectors.
 *  @param[in] a: Matrix the transforms the lattice in the Epstein Zeta function.
 *  @param[in] x: x vector of the Epstein Zeta function.
 *  @param[in] y: y vector of the Epstein Zeta function.
 *  @return 0
 */
int epstein_zeta_mathematica_call(double *out, double nu, int dim, double *a,
                                  double *x, double *y) {
    double complex asg = epsteinZeta(nu, dim, a, x, y);
    out[0] = creal(asg);
    out[1] = cimag(asg);
    return 0;
}
/*
 * @brief Saves real and imaginary parts of the regularized Epstein Zeta function
 * in an array.
 * @param[in, out] out: 2D array, where real and imaginary parts of the regularized
 * Epstein Zeta function will be stored.
 * @param[in] nu: Exponent for the regularized Epstein zeta function.
 * @param[in] dim: Dimension of the input vectors.
 * @param[in] a: Matrix that transforms the lattice in the Epstein Zeta function.
 * @param[in] x: x vector of the Epstein Zeta function.
 * @param[in] y: y vector of the Epstein Zeta function.
 * @return 0
 */
int epstein_zeta_reg_mathematica_call(double *out, double nu, int dim, double *a,
                                      double *x, double *y) {
    double complex asg = epsteinZetaReg(nu, dim, a, x, y);
    out[0] = creal(asg);
    out[1] = cimag(asg);
    return 0;
}
