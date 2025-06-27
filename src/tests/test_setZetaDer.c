// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "../tools.h"
#include "epsteinZeta.h"
#include "utils.h"
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

#ifndef BASE_PATH
#define BASE_PATH "csv"
#endif

/*!
 * @brief Benchmarks 2D set zeta derivatives by computing its taylor series.
 *
 * @return 0 if all tests pass, 1 if any test fails.
 */
int test_setZetaDer_taylor(void) {
    printf("%s ... ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex valRef;
    double complex valTaylor;

    double tol = pow(10, -14);
    unsigned int dim = 2;
    unsigned int order = 12;
    double zDiff[] = {0.005, 0.01};

    double nu = 0.5;
    double m[] = {1., 0., 0., 1.};
    double x[] = {0., 0.};
    double y0[] = {0., 0.1};
    double yDiff[] = {0., 0.01};
    unsigned int alpha0[] = {0, 0};
    double *yPlus = malloc(dim * sizeof(double));

    int testsPassed = 0;
    int totalTests = 0;

    bool done;
    unsigned int *alpha = malloc(dim * sizeof(unsigned int));

    for (int i = 0; i < 1; i++) {

        for (int i = 0; i < dim; i++) {
            yPlus[i] = y0[i] + yDiff[i];
        }

        valRef = setZetaDer(nu, dim, m, x, yPlus, alpha0);

        // build taylor series around z
        valTaylor = 0;

        // Initialize multi-index
        for (int i = 0; i < dim; i++) {
            alpha[i] = 0;
        }

        // Iterate over every multi-index alpha so that every alpha[] < order
        while (1) {

            valTaylor += mult_pow(dim, alpha, zDiff, 1.) /
                         (double)mult_fac(dim, alpha) *
                         setZetaDer(nu, dim, m, x, y0, alpha);

            done = true;
            for (unsigned int idx = 0; idx < dim; idx++) {
                if (alpha[idx] + 1 <= order) {
                    alpha[idx]++;
                    done = false;
                    break;
                }
                alpha[idx] = 0;
            }
            if (done) {
                break;
            }
        }

        errorAbs = errAbs(valRef, valTaylor);
        errorRel = errRel(valRef, valTaylor);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        totalTests++;
        if (errorMaxAbsRel < tol) {
            testsPassed++;
        } else {
            printf("\nWarning! ");
            printf("setZetaDer: ");
            printf(" %0*.16lf %+.16lf I (as a taylor series) \n\t\t!= "
                   "%.16lf "
                   "%+.16lf I (reference implementation)\n",
                   4, creal(valTaylor), cimag(valTaylor), creal(valRef),
                   cimag(valRef));
            printf("Min(Emax, Erel):      %E > %E  (tolerance)\n", errorMaxAbsRel,
                   tol);
            printf("\n");
            printf("nu:\t\t %.16lf\n", nu);
            printVectorUnitTest("y0:\t\t", y0, dim);
            printVectorUnitTest("yPlus:\t\t", yPlus, dim);
            printVectorUnitTest("yDiff:\t\t", yDiff, dim);
            printf("\n");
        }
    }

    printf("%d out of %d tests passed.\n", testsPassed, totalTests);

    return totalTests - testsPassed;
}

/*!
 * @brief Main function to run all set zeta derivatives function tests.
 *
 * @return 0 if all tests pass, non-zero if any test fails.
 */
int main() {
    bool failed = test_setZetaDer_taylor();
    return failed;
}
