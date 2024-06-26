// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file gamma.h
 * @brief Gamma functions.
 *
 * Calculates the gamma function, the incomplete upper gamma function and the
regularized lower incomplete gamma function for evaluations of Crandall's
formula.
 * @see Walter Gautschi. “A Computational Procedure for
Incomplete Gamma Func-257 tions”. In: ACM Trans. Math. Softw. 5 (1979), pp.
466–481
 */

#ifndef GAMMA_H
#define GAMMA_H
/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param a: exponent of the upper incomplete gamma function.
 * @param x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_ugamma(double a, double x);
/**
 * @brief calculate the upper incomplete gamma function as in Gautschi.
 * @param a: exponent of the upper incomplete gamma function.
 * @param x: lower integral boundary of the upper incomplete gamma function.
 * @return function value of the upper incomplete gamma function.
 */
double egf_gammaStar(double a, double x);
#endif
