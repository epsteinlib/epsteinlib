// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

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
 * @brief Free memory allocated for test resources.
 *
 * This function deallocates memory that was dynamically allocated for various
 * arrays used in the set zeta function tests.
 *
 * @param[in] a Pointer to the matrix of coefficients.
 * @param[in] nu Pointer to the array of nu values.
 * @param[in] x Pointer to the array of x values.
 * @param[in] y Pointer to the array of y values.
 * @param[in] zetaRef Pointer to the array of reference zeta values.
 */
void freeTestResources(double *a, double *nu, double *x, double *y,
                       double *zetaRef) {
    free(a);
    free(nu);
    free(x);
    free(y);
    free(zetaRef);
}

/*!
 * @brief Reports an error when the set zeta function representation test fails.
 *
 * This function prints detailed information about the test case that failed,
 * including the computed values, error margins, and input parameters.
 *
 * @param[in] valZeta: value computed by setZetaDer.
 * @param[in] valZetaTaylor: value computed by the setZetaDerReg representation.
 * @param[in] errorMaxAbsRel: maximum of absolute and relative errors.
 * @param[in] tol: tolerance level for the test.
 * @param[in] m: pointer to the matrix of coefficients.
 * @param[in] dim: dimension of the lattice.
 * @param[in] nu: parameter nu of the Epstein zeta function.
 * @param[in] x: pointer to the array of x values.
 * @param[in] y: pointer to the array of y values.
 * @param[in] alpha: multiindex for the deriavtives.
 */
void reportSetZetaDerError(double complex valZeta, double complex valZetaTaylor,
                           double errorMaxAbsRel, double tol, const double *m,
                           unsigned int dim, const double nu, const double *x,
                           const double *y, const unsigned int *alpha) {
    printf("\nWarning! ");
    printf("setZetaDer:");
    printf(" %0*.16lf %+.16lf I (setZetaDer) \n\t\t  != "
           "%.16lf "
           "%+.16lf I (setZetaDerReg representation)\n",
           4, creal(valZeta), cimag(valZeta), creal(valZetaTaylor),
           cimag(valZetaTaylor));
    printf("Min(Emax, Erel):      %.16lf > %.16lf  (tolerance)\n", errorMaxAbsRel,
           tol);
    printf("\n");
    printMatrixUnitTest("m:", m, (int)dim);
    printf("nu:\t\t %.16lf\n", nu);
    printVectorUnitTest("x:\t\t", x, (int)dim);
    printVectorUnitTest("y:\t\t", y, (int)dim);
    printMultiindexUnitTest("alpha:\t\t", alpha, (int)dim);
}

/*!
 * @brief Test if the set zeta function can be represented in terms of the
 * taylor series.
 * @return true if any test case fails (difference exceeds tolerance), false if all
 * tests pass.
 */
bool test_setZetaDer() {
    printf("%s ... ", __func__);
    double errorAbs;
    double errorRel;
    double errorMaxAbsRel;
    double complex valZeta;
    double complex valZetaTaylor;
    double nu;

    double tol = pow(10, -14);
    unsigned int dim = 2;
    double m[] = {1, 0, 0, 1};
    double x[] = {0.1, 0.2};
    double y[] = {0, 0.5};
    double yPlus[] = {0, 0.5 + 0.0};
    unsigned int alpha0[] = {0, 0};
    unsigned int alpha[] = {0, 0};

    for (int i = 0; i < 100 + 1; i++) {
        nu = -12.5 + (double)i / 4;
        valZeta = setZetaDer(nu, dim, m, x, y, alpha0);
        valZetaTaylor = setZetaDer(nu, dim, m, x, yPlus, alpha);

        errorAbs = errAbs(valZeta, valZetaTaylor);
        errorRel = errRel(valZeta, valZetaTaylor);
        errorMaxAbsRel = (errorAbs < errorRel) ? errorAbs : errorRel;

        if (errorMaxAbsRel > tol) {
            reportSetZetaDerError(valZeta, valZetaTaylor, errorMaxAbsRel, tol, m,
                                  dim, nu, x, y, alpha);
            return true;
        }
    }

    printf("passed.\n");
    return false;
}

/*!
 * @brief Main function to run all set zeta derivatives function tests.
 *
 * @return 0 if all tests pass, non-zero if any test fails.
 */
int main() {
    bool failed = test_setZetaDer();
    return failed;
}
