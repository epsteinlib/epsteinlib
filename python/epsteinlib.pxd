# SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
# SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
# SPDX-FileCopyrightText: 2024 Ruben Gutendorf <ruben.gutendorf@uni-saarland.de>
#
# SPDX-License-Identifier: AGPL-3.0-only

cdef extern from "../include/epsteinZeta.h":
    double complex epsteinZeta(double nu, int dim, const double *a, const double *x, const double *y)
    double complex epsteinZetaReg(double nu, int dim, const double *a, const double *x, const double *y)
