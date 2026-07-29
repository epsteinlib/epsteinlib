// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file epsteinZeta.c
 * @brief Calculates the (regularized) Epstein zeta function and the (regularized)
 * anisotropic Epstein zeta function.
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @see Buchheit, A. A., Busse, J. K., & Gutendorf, R. (2026). "Computation and
 * properties of the Epstein zeta function with applications to quantum systems." IMA
 * Journal of Numerical Analysis, drag057. DOI: 10.1093/imanum/drag057 *
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @date 06/13/2024
 */

#include <complex.h>
#include <stdbool.h>

#include "epsteinZeta.h"
#include "zeta.h"

/**
 * @brief calculates the Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice.
 * @param[in] x: shift vector.
 * @param[in] y: wavevector.
 * @return function value of the Epstein zeta function.
 */
double complex epsteinZeta(double nu, unsigned int dim, const double *a,
                           const double *x, const double *y) {
    return epsteinZetaInternal(nu, dim, a, x, y, 1, false, false,
                               (unsigned int[]){0});
}

/**
 * @brief calculates the regularized Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice.
 * @param[in] x: shift vector.
 * @param[in] y: wavevector.
 * @return function value of the regularized Epstein zeta function.
 */
double complex epsteinZetaReg(double nu, unsigned int dim, const double *a,
                              const double *x, const double *y) {
    return epsteinZetaInternal(nu, dim, a, x, y, 1, true, false,
                               (unsigned int[]){0});
}

/**
 * @brief calculates the anisotropic Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice.
 * @param[in] x: shift vector.
 * @param[in] y: wavevector.
 * @param[in] alpha: multiindex for the anisotropy.
 * @return function value of the anisotropic Epstein zeta function.
 */
double complex epsteinZetaAniso(double nu, unsigned int dim, const double *a,
                                const double *x, const double *y,
                                const unsigned int *alpha) {
    return epsteinZetaInternal(nu, dim, a, x, y, 1, false, true, alpha);
}

/**
 * @brief calculates the regularized anisotropic Epstein zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice.
 * @param[in] x: shift vector.
 * @param[in] y: wavevector.
 * @param[in] alpha: multiindex for the anisotropy.
 * @return function value of the regularized anisotropic Epstein zeta function.
 */
double complex epsteinZetaAnisoReg(double nu, unsigned int dim, const double *a,
                                   const double *x, const double *y,
                                   const unsigned int *alpha) {
    return epsteinZetaInternal(nu, dim, a, x, y, 1, true, true, alpha);
}
