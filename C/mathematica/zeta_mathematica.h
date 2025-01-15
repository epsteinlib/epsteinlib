// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#ifndef EPSTEIN_MATHEMATICA
#define EPSTEIN_MATHEMATICA
int epstein_zeta_mathematica_call(double *out, double nu, int dim, double *a,
                                  double *x, double *y);
int epstein_zeta_reg_mathematica_call(double *out, double nu, int dim, double *a,
                                      double *x, double *y);
#endif
