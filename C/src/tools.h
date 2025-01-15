// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file tools.h
 * @brief  Minimal linear algebra for matrix vector operations.
 */

#ifndef EPSTEIN_TOOLS
#define EPSTEIN_TOOLS
#include <stdbool.h>
double dot(unsigned int dim, const double *v1, const double *v2);
void matrix_intVector(unsigned int dim, const double *m, const int *v, double *res);
void transpose(unsigned int dim, double *m);
bool equals(unsigned int dim, const double *v1, const double *v2);
bool equalsZero(unsigned int dim, const double *v);
void invert(unsigned int dim, double *m, int *p, double *r);
double inf_norm(unsigned int dim, const double *m);
#endif
