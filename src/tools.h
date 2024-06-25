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
double dot(int dim, double *v1, double *v2);
void matrix_intVector(int dim, double *m, int *v, double *res);
void transpose(int dim, double *m);
bool equals(int dim, double *v1, double *v2);
bool equalsZero(int dim, double *v);
#endif
