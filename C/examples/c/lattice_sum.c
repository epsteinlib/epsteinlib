// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file lattice_sum.c
 * @brief Calculates the Madelung constant.
 *
 * Minimal working example for the Epstein Zeta Library.
 * If the library is installed, compile with `gcc -o lattice_sum lattice_sum.c
 * -lm -lepsteinZeta` If the library is not installed, compile with `gcc -o
 * lattice_sum lattice_sum.c -lm -L/path/to/library -lepsteinZeta
 * -I/path/to/include`
 * @author Andreas Buchheit, Jonathan Busse and Ruben Gutendorf.
 * @date 06/06/2024
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
// If library is not installed, use #include "/path/to/epsteinZeta.h"
#include "epsteinZeta.h"

/**
 * @brief calculate madelung constant and compare to precomputed value.
 *
 * Madelung constant:
 * sum_{i, j, k in Z} (-1)**(i + j + k) / sqrt(i**2 + j**2 + k**2)
 * @return true, if the difference to precomputed value is smaller than
 * 10**(-14)
 */
int main() {

    // Madelung constant found in literature
    double madelungRef = -1.7475645946331821906362120355443974;
    int dim = 3;
    double m[] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // identity matrix for whole numbers
    double x[] = {0, 0, 0};                   // no shift
    double y[] = {0.5, 0.5, 0.5};             // alternating sum
    double nu = 1.0;
    double madelung = creal(epsteinZeta(nu, dim, m, x, y));
    printf("Madelung sum in 3 dimensions:\t %.16lf\n", creal(madelung));
    printf("Reference value:\t\t %.16lf\n", madelungRef);
    printf("Relative error:\t\t\t +%.2e\n",
           fabs(madelungRef - madelung) / fabs(madelungRef));

    return fabs(madelung - madelungRef) > pow(10, -14);
}
