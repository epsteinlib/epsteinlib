// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file zeta.c
 * @brief Calculates the (regularized) Epstein zeta function.
 */

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "crandall.h"
#include "tools.h"

#include "zeta.h"

/*!
   @brief Smallest value z such that G(nu, z) is negligible for
   nu < 10.
*/
#define G_BOUND 3.2

/*!
 * @brief epsilon for the cutoff around nu = dimension.
 */
#define EPS ldexp(1, -30)

/*!
 * @brief epsilon for the cutoff around x = 0 and y = 0
 */
#define EPS_ZERO_Y 1e-64

/**
 * @brief calculates the first sum in Crandall's formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameters that decides the weight of each sum.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @return helper function for the first sum in crandalls formula. Calculates
 * sum_{z in m whole_numbers ** dim} G_{nu}((z - x) / lambda))
 * X exp(-2 * PI * I * z * y)
 */
double complex sum_real(double nu, unsigned int dim, double lambda, const double *m,
                        const double *x, const double *y, const int cutoffs[],
                        double zArgBound) {
    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    long totalCutoffs[dim + 1];
    for (int k = 0; k < dim; k++) {
        totalCutoffs[k] = totalSummands;
        totalSummands *= 2 * cutoffs[k] + 1;
    }
    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;
    double complex rot;
    // First Sum (in real space)
    for (long n = 0; n < totalSummands; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] =
                (((int)(n / totalCutoffs[k])) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m, zv, lv);
        rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] - x[i];
        }
        // summing using Kahan's method
        auxy = rot * crandall_g(dim, nu, lv, 1. / lambda, zArgBound) - epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    return sum;
}

/**
 * @brief calculates the derivative of the first sum in Crandall's formula including
 * a phase-factor.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameters that decides the weight of each sum.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alpha: multi-index for the derivative.
 * @return helper function for the first sum in crandalls formula. Calculates
 * sum_{z in m whole_numbers ** dim} (-2 * PI * I * (z-x) ) ** alpha
 * * G_{nu}^{(alpha)}((z - x) / lambda)) X exp(-2 * PI * I * (z * y)
 */
double complex sum_real_der(double nu, unsigned int dim, double lambda,
                            const double *m, const double *x, const double *y,
                            const int cutoffs[], double zArgBound,
                            const unsigned int *alpha) {
    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    long totalCutoffs[dim + 1];
    for (int k = 0; k < dim; k++) {
        totalCutoffs[k] = totalSummands;
        totalSummands *= 2 * cutoffs[k] + 1;
    }
    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;
    double complex rot;
    double complex mon;
    // First Sum (in real space)
    for (long n = 0; n < totalSummands; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] =
                (((int)(n / totalCutoffs[k])) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m, zv, lv);
        rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] - x[i];
        }
        // Calculate (- 2 PI I (z - x)) ** alpha
        mon = 1.;
        for (int i = 0; i < dim; i++) {
            if (alpha[i]) {
                mon *= int_pow(-2 * M_PI * I * lv[i], alpha[i]);
            }
        }
        // summing using Kahan's method
        auxy = rot * mon * crandall_g(dim, nu, lv, 1. / lambda, zArgBound) - epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    return sum;
}

/**
 * @brief calculates the second sum in Crandall's formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameters that decides the weight of each sum.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @return helper function for the second sum in crandalls formula. Calculates
 * sum_{k in m_invt whole_numbers ** dim without zero} G_{dim - nu}(lambda * (k + y))
 * X exp(-2 * PI * I * x * (k + y))
 */
double complex sum_fourier(double nu, unsigned int dim, double lambda,
                           const double *m_invt, const double *x, const double *y,
                           const int cutoffs[], double zArgBound) {
    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    long totalCutoffs[dim + 1];
    for (int k = 0; k < dim; k++) {
        totalCutoffs[k] = totalSummands;
        totalSummands *= 2 * cutoffs[k] + 1;
    };
    long zeroIndex = (totalSummands - 1) / 2;
    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;
    // second sum (in fourier space)
    for (long n = 0; n < zeroIndex; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] =
                (((int)(n / totalCutoffs[k])) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m_invt, zv, lv);
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] + y[i];
        }
        double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot * crandall_g(dim, dim - nu, lv, lambda, zArgBound) - epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    // skips zero
    for (long n = zeroIndex + 1; n < totalSummands; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] =
                (((int)(n / totalCutoffs[k])) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m_invt, zv, lv);
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] + y[i];
        }
        double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot * crandall_g(dim, dim - nu, lv, lambda, zArgBound) - epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    return sum;
}

/**
 * @brief calculates the second sum  for the derivatives in Crandall's formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameters that decides the weight of each sum.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @return helper function for the second sum in crandalls formula. Calculates
 * sum_{k in m_invt whole_numbers ** dim without zero}
 * * G_{dim - nu}^{(alpha)}(lambda * (k + y)) * exp(-2 * PI * I * x * (k + y))
 */
double complex sum_fourier_der(double nu, unsigned int dim, double lambda,
                               const double *m_invt, const double *x,
                               const double *y, const int cutoffs[],
                               double zArgBound, const unsigned int *alpha) {
    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    long totalCutoffs[dim + 1];
    for (int k = 0; k < dim; k++) {
        totalCutoffs[k] = totalSummands;
        totalSummands *= 2 * cutoffs[k] + 1;
    };
    long zeroIndex = (totalSummands - 1) / 2;
    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;
    // second sum (in fourier space)
    for (long n = 0; n < zeroIndex; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] =
                (((int)(n / totalCutoffs[k])) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m_invt, zv, lv);
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] + y[i];
        }
        double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot * crandall_g_der(dim, dim - nu, lv, lambda, zArgBound, alpha) -
               epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    // skips zero
    for (long n = zeroIndex + 1; n < totalSummands; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] =
                (((int)(n / totalCutoffs[k])) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m_invt, zv, lv);
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] + y[i];
        }
        double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot * crandall_g_der(dim, dim - nu, lv, lambda, zArgBound, alpha) -
               epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    return sum;
}

/**
 * @brief calculate projection of vector to elementary lattice cell.
 * @param[in] dim: dimension of the input vectors
 * @param[in] m: matrix that transforms the lattice in the function.
 * @param[in] m_invt: inverse of m.
 * @param[in] v: vector for which the projection to the elementary lattice cell
 * is needet.
 * @return projection of v to the elementary lattice cell.
 */
double *vectorProj(unsigned int dim, const double *m, const double *m_invt,
                   const double *v) {
    bool todo = false;
    double *vt = malloc(dim * sizeof(double));
    for (int i = 0; i < dim; i++) {
        vt[i] = 0;
        for (int j = 0; j < dim; j++) {
            vt[i] += m_invt[(dim * j) + i] * v[j];
        }
    }
    // check if projection is needed, else copy
    for (int i = 0; i < dim && !todo; i++) {
        todo = todo || (vt[i] <= -0.5 || vt[i] >= 0.5);
    }
    if (todo) {
        for (int i = 0; i < dim; i++) {
            vt[i] = remainder(vt[i], 1);
        }
        double *vres = malloc(dim * sizeof(double));
        for (int i = 0; i < dim; i++) {
            vres[i] = 0;
            for (int j = 0; j < dim; j++) {
                vres[i] += m[(dim * i) + j] * vt[j];
            }
        }
        free(vt);
        return vres;
    }
    for (int i = 0; i < dim; i++) {
        vt[i] = v[i];
    }
    return vt;
}

/**
 * @brief calculates the (regularized) Epstein zeta function as well es the
 * derivatives of the set zeta function.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: x vector of the Epstein zeta function.
 * @param[in] y: y vector of the Epstein zeta function.
 * @param[in] lambda: relative weight of the sums in Crandall's formula.
 * @param[in] variant: 0 for no regularization
 *                    1 for the regularization
 *                    2 for set Zeta (including derivatives)
 *                    3 for the regularization (including Derivatives).
 * @param[in] alpha: multiindex for the derivatives of the set zeta function. *
 * @return function value of the regularized Epstein zeta.
 */
double complex epsteinZetaInternal(double nu, unsigned int dim, // NOLINT
                                   const double *m, const double *x, const double *y,
                                   double lambda, unsigned int variant,
                                   const unsigned int *alpha) {
    unsigned int alphaAbs = mult_abs(dim, alpha);
    // Early return for 0th derivative special cases
    if (variant == 2 && !alphaAbs) {
        return cexp(2 * M_PI * I * dot(dim, x, y)) *
               epsteinZetaInternal(nu, dim, m, x, y, 1, 0, alpha);
    }

    if (variant == 3 && !alphaAbs) {
        return epsteinZetaInternal(nu, dim, m, x, y, 1, 1, alpha);
    }
    // 1. Transform: Compute determinant and fourier transformed matrix, scale
    // both of them
    double m_fourier[dim * dim];
    double m_copy[dim * dim];
    double m_real[dim * dim];
    double x_t1[dim];
    double y_t1[dim];
    int p[dim];
    bool isDiagonal = 1;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            m_copy[(dim * i) + j] = m[(dim * i) + j];
            m_real[(dim * i) + j] = m[(dim * i) + j];
            isDiagonal = isDiagonal && ((i == j) || (m[(dim * i) + j] == 0));
        }
    }
    invert(dim, m_copy, p, m_fourier);
    double vol = 1;
    for (int k = 0; k < dim; k++) {
        vol *= m_copy[(dim * k) + k];
    }
    transpose(dim, m_fourier);
    vol = fabs(vol);
    double ms = pow(vol, -1. / dim);
    for (int i = 0; i < dim * dim; i++) {
        m_real[i] *= ms;
        m_fourier[i] /= ms;
    }
    for (int i = 0; i < dim; i++) {
        x_t1[i] = x[i] * ms;
        y_t1[i] = y[i] / ms;
    }
    // 2. transform: get x and y in their respective elementary cells
    double *x_t2 = vectorProj(dim, m_real, m_fourier, x_t1);
    double *y_t2 = vectorProj(dim, m_fourier, m_real, y_t1);
    // set cutoffs
    int cutoffsReal[dim];
    int cutoffsFourier[dim];
    double cutoff_id = G_BOUND + 0.5;
    if (isDiagonal) {
        // Chose absolute diag. entries for cutoff
        for (int k = 0; k < dim; k++) {
            cutoffsReal[k] = floor(cutoff_id / fabs(m_real[(dim * k) + k]));
            cutoffsFourier[k] = floor(cutoff_id * fabs(m_real[(dim * k) + k]));
        }
    } else {
        // choose cutoff depending on smallest and biggest abs eigenvalue
        double ev_abs_max = inf_norm(dim, m_real);
        double ev_abs_min_r = inf_norm(dim, m_fourier);
        for (int k = 0; k < dim; k++) {
            cutoffsReal[k] = floor(cutoff_id * ev_abs_min_r);
            cutoffsFourier[k] = floor(cutoff_id * ev_abs_max);
        }
    }
    // handle special case of non-positive integer values nu.
    double complex res;
    if (nu < 1 && fabs((nu / 2.) - nearbyint(nu / 2.)) < EPS &&
        (variant == 0 || variant == 1)) {
        if (dot(dim, x_t2, x_t2) == 0 && nu == 0) {
            res = -1 * cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2));
        } else {
            res = 0;
        }
    } else if (fabs(nu - dim) < EPS && dot(dim, y_t2, y_t2) < EPS_ZERO_Y &&
               variant == 0) {
        res = NAN;
    } else {
        double zArgBound = assignzArgBound(nu);
        double zArgBoundReci = assignzArgBound(dim - nu);
        double complex s1 = 0;
        double complex s2 = 0;
        double complex nc = 0;
        double complex rot = 1;
        double complex xfactor = 1;
        double vx[dim];
        for (int i = 0; i < dim; i++) {
            vx[i] = x_t1[i] - x_t2[i];
        }
        xfactor = cexp(-2 * M_PI * I * dot(dim, vx, y_t1));
        if (variant == 0) {
            // calculate non regularized Epstein zeta function values.
            nc = crandall_g(dim, dim - nu, y_t2, lambda, zArgBound) *
                 cexp(-2 * M_PI * I * dot(dim, x_t2, y_t2));
            s1 = sum_real(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                          zArgBound);
            s2 = sum_fourier(nu, dim, lambda, m_fourier, x_t2, y_t2, cutoffsFourier,
                             zArgBoundReci) +
                 nc;
        } else if (variant == 1) {
            // calculate regularized Epstein zeta function values.
            nc = crandall_gReg(dim, dim - nu, y_t1, lambda);
            rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));
            s2 = sum_fourier(nu, dim, lambda, m_fourier, x_t1, y_t2, cutoffsFourier,
                             zArgBoundReci);
            // correct wrong zero summand in regularized fourier sum.
            if (!equals(dim, y_t1, y_t2)) {
                s2 += crandall_g(dim, dim - nu, y_t2, lambda, zArgBoundReci) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2)) -
                      crandall_g(dim, dim - nu, y_t1, lambda, zArgBoundReci) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t1));
            }
            s2 = s2 * rot + nc;
            s1 = sum_real(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                          zArgBound) *
                 rot * xfactor;
            xfactor = 1;
        } else if (variant == 2) {
            // calculate non set zeta  derivative function values.
            rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));
            if (equals(dim, y_t1, y_t2)) {
                nc = crandall_g_der(dim, dim - nu, y_t1, lambda, zArgBoundReci,
                                    alpha);
            } else {
                nc = crandall_g_der(dim, dim - nu, y_t2, lambda, zArgBoundReci,
                                    alpha) *
                     cexp(-2 * M_PI * I * dot(dim, y_t2, x_t1)) * rot;
            }
            s2 = sum_fourier_der(nu, dim, lambda, m_fourier, x_t1, y_t2,
                                 cutoffsFourier, zArgBoundReci, alpha);
            s2 = int_pow(lambda, alphaAbs) * (s2 * rot + nc);
            s1 = sum_real_der(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                              zArgBound, alpha) *
                 rot * xfactor;
            xfactor = 1. / int_pow(ms, alphaAbs);
        }
        res = xfactor * pow(lambda * lambda / M_PI, -nu / 2.) / tgamma(nu / 2.) *
              (s1 + pow(lambda, dim) * s2);
    }
    free(x_t2);
    free(y_t2);
    res *= pow(ms, nu);
    // apply correction to matrix scaling if nu = d + 2k
    double k = fmax(0., nearbyint((nu - (double)dim) / 2));
    if ((variant == 1) && (nu == (dim + 2 * k))) {
        if (k == 0) {
            res += pow(M_PI, (double)dim / 2) / tgamma((double)dim / 2) *
                   log(ms * ms) / vol;
        } else {
            double ySquared = 0;
            for (int i = 0; i < dim; i++) {
                ySquared += y[i] * y[i];
            }
            res -= pow(M_PI, (2 * k) + ((double)dim / 2)) /
                   tgamma(k + ((double)dim / 2)) * pow(-1, k + 1) / tgamma(k + 1) *
                   pow(ySquared, k) * log(ms * ms) / vol;
        }
    }
    return res;
}
#undef G_BOUND
