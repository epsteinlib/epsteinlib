// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#ifndef EPSTEIN_MATHEMATICA
#define EPSTEIN_MATHEMATICA

/*
 * @brief Saves real and imaginary parts of the Epstein zeta function and
 * in an array.
 * @param[in, out] out: 2D array, where real and imaginary parts of the Epstein zeta
 * function will be stored.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 *  @return 0
 */
int epstein_zeta_mathematica_call(double *out, double nu, int dim, const double *a,
                                  const double *x, const double *y);

/*
 * @brief Saves real and imaginary parts of the regularized Epstein zeta function
 * in an array.
 * @param[in, out] out: 2D array, where real and imaginary parts of the regularized
 * Epstein zeta function will be stored.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] a: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 *  @return 0
 */
int epstein_zeta_reg_mathematica_call(double *out, double nu, int dim,
                                      const double *a, const double *x,
                                      const double *y);
#endif
