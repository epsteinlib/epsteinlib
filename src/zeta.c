// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file zeta.c
 * @brief Calculates the (regularized) Epstein zeta function.
 */

#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "crandall.h"
#include "gamma.h"
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

/**
 * @brief calculates the first sum in Crandall's formula.
 * @param nu: exponent for the Epstein zeta function.
 * @param dim: dimension of the input vectors.
 * @param lambda: parameters that decides the weight of each sum.
 * @param m: matrix that transforms the lattice in the Epstein Zeta function.
 * @param x: projection of x vector to elementary lattice cell.
 * @param y: projection of y vector to elementary lattice cell.
 * @param cutoffs: how many summands in each direction are considered.
 * @param zArgBound: global bound on when to use the asymptotic expansion in the
 * incomplete gamma evaluation.
 * @return the first sum in crandalls formula, if
 * x, y are in the elementary lattice cells. Multiply with exp(2 * Pi * I * x *
 * y) otherwise.
 */
double complex sum_real(double nu, short dim, double lambda, double *m, double *x,
                        double *y, int cutoffs[], double zArgBound) {
    // 1. Transform: Compute determinant and fourier transformed matrix,
    // scale both of them
    int *zv = malloc(dim * sizeof(int));       // counting vector in Z^dim
    double *lv = malloc(dim * sizeof(double)); // lattice vector
    double complex s1 = 0;
    // cuboid cutoffs
    long totalSummands = 1;
    long totalCutoffs[dim + 1];
    for (int k = 0; k < dim; k++) {
        totalCutoffs[k] = totalSummands;
        totalSummands *= 2 * cutoffs[k] + 1;
    };
    // First Sum (in real space)
    for (long n = 0; n < totalSummands; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] = ((n / totalCutoffs[k]) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m, zv, lv);
        double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
        for (int i = 0; i < dim; i++)
            lv[i] = (lv[i] - x[i]) / lambda;
        s1 += rot * crandall_g(dim, nu, lv, 1, zArgBound);
        // s1 += rot * crandall_g(dim, nu, lv, 1);
    }
    free(zv);
    free(lv);
    return s1;
}

/**
 * @brief calculates the second sum in Crandall's formula.
 * @param nu: exponent for the Epstein zeta function.
 * @param dim: dimension of the input vectors.
 * @param lambda: parameters that decides the weight of each sum.
 * @param m: matrix that transforms the lattice in the Epstein Zeta function.
 * @param x: projection of x vector to elementary lattice cell.
 * @param y: projection of y vector to elementary lattice cell.
 * @param cutoffs: how many summands in each direction are considered.
 * @param zArgBound: global bound on when to use the asymptotic expansion in the
 * incomplete gamma evaluation.
 * @return the second sum in crandalls formula, if
 * x, y are in the elementary lattice cells. Multiply with exp(2 * Pi * I * x *
 * y) otherwise. Add the zero summand in the regularization case.
 */
double complex sum_fourier(double nu, short dim, double lambda, double *m, double *x,
                           double *y, int cutoffs[], double zArgBound) {
    int *zv = malloc(dim * sizeof(int));       // counting vector in Z^dim
    double *lv = malloc(dim * sizeof(double)); // lattice vector
    double complex s2 = 0;
    // cuboid cutoffs
    long totalSummands = 1;
    long totalCutoffs[dim + 1];
    for (int k = 0; k < dim; k++) {
        totalCutoffs[k] = totalSummands;
        totalSummands *= 2 * cutoffs[k] + 1;
    };
    long zeroIndex = (totalSummands - 1) / 2;
    // First Sum (in fourier space)
    for (long n = 0; n < zeroIndex; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] = ((n / totalCutoffs[k]) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m, zv, lv);
        for (int i = 0; i < dim; i++)
            lv[i] = lv[i] + y[i];
        double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        s2 += rot * crandall_g(dim, dim - nu, lv, lambda, zArgBound);
    }
    // skips zero
    for (long n = zeroIndex + 1; n < totalSummands; n++) {
        for (int k = 0; k < dim; k++) {
            zv[k] = ((n / totalCutoffs[k]) % (2 * cutoffs[k] + 1)) - cutoffs[k];
        }
        matrix_intVector(dim, m, zv, lv);
        for (int i = 0; i < dim; i++)
            lv[i] = lv[i] + y[i];
        double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        s2 += rot * crandall_g(dim, dim - nu, lv, lambda, zArgBound);
    }
    free(zv);
    free(lv);
    return s2;
}

/**
 * @brief calculate projection of vector to elementary lattice cell.
 * @param dim: dimension of the input vectors
 * @param m: matrix that transforms the lattice in the function.
 * @param m_invt: inverse of m.
 * @param v: vector for which the projection to the elementary lattice cell is
 * needet.
 * @return projection of v to the elementary lattice cell.
 */
double *vectorProj(short dim, double *m, double *m_invt, double *v) {
    bool todo = false;
    double *vt = malloc(dim * sizeof(double));
    for (int i = 0; i < dim; i++) {
        vt[i] = 0;
        for (int j = 0; j < dim; j++) {
            vt[i] += m_invt[dim * j + i] * v[j];
        }
    }
    for (int i = 0; i < dim && !todo; i++)
        todo = todo || !(vt[i] >= -0.5 && vt[i] <= 0.5);
    if (todo) {
        for (int i = 0; i < dim; i++)
            vt[i] = remainder(vt[i], 1);
        double *vres = malloc(dim * sizeof(double));
        for (int i = 0; i < dim; i++) {
            vres[i] = 0;
            for (int j = 0; j < dim; j++)
                vres[i] += m[dim * i + j] * vt[j];
        }
        free(vt);
        return vres;
    } else {
        for (int i = 0; i < dim; i++)
            vt[i] = v[i];
        return vt;
    }
}

/**
 * @brief calculates the (regularized) Epstein Zeta function.
 * @param nu: exponent for the Epstein zeta function.
 * @param dim: dimension of the input vectors.
 * @param m: matrix that transforms the lattice in the Epstein Zeta function.
 * @param x: x vector of the Epstein Zeta function.
 * @param y: y vector of the Epstein Zeta function.
 * @param lambda: relative weight of the sums in Crandall's formula.
 * @param regBool: 0 for no regularization, > 0 for the regularization.
 * @return function value of the regularized Epstein zeta.
 */
double complex __epsteinZeta(double nu, int dim, double *m, double *x, double *y,
                             double lambda, int reg) {
    // 1. Transform: Compute determinant and fourier transformed matrix, scale
    // both of them
    double *m_fourier = malloc(dim * dim * sizeof(double));
    double *m_copy = malloc(dim * dim * sizeof(double));
    double *m_real = malloc(dim * dim * sizeof(double));
    double *x_t1 = malloc(dim * sizeof(double));
    double *y_t1 = malloc(dim * sizeof(double));
    int *p = malloc(dim * sizeof(int));
    bool isDiagonal = 1;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            m_copy[dim * i + j] = m[dim * i + j];
            m_real[dim * i + j] = m[dim * i + j];
            m_fourier[dim * i + j] = (i == j) ? 1 : 0;
            isDiagonal = isDiagonal && ((i == j) || (m[dim * i + j] == 0));
        }
    }
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, dim, dim, m_copy, dim, p, m_fourier, dim);
    double vol = 1;
    for (int k = 0; k < dim; k++)
        vol *= m_copy[dim * k + k];
    transpose(dim, m_fourier);
    free(p);
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
            cutoffsReal[k] = floor(cutoff_id / fabs(m_real[dim * k + k]));
            cutoffsFourier[k] = floor(cutoff_id * fabs(m_real[dim * k + k]));
        }
    } else {
        // choose cutoff depending on smallest and biggest abs eigenvalue
        double *m_ev = malloc(dim * sizeof(double));
        double *m_aux = malloc(dim * sizeof(double));
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                m_copy[dim * i + j] = m_real[dim * i + j];
        }
        LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', dim, m_copy, dim, m_ev, m_aux,
                      NULL, dim, NULL, dim);
        double ev_abs_min = fabs(m_ev[0]);
        double ev_abs_max = fabs(m_ev[0]);
        for (int k = 1; k < dim; k++) {
            ev_abs_min = (fabs(m_ev[k]) < ev_abs_min) ? fabs(m_ev[k]) : ev_abs_min;
            ev_abs_max = (fabs(m_ev[k]) > ev_abs_max) ? fabs(m_ev[k]) : ev_abs_max;
        }
        for (int k = 0; k < dim; k++) {
            cutoffsReal[k] = floor(cutoff_id / ev_abs_min);
            cutoffsFourier[k] = floor(cutoff_id * ev_abs_max);
        }
        free(m_ev);
        free(m_aux);
    }
    free(m_copy);
    // handle special case of non-positive integer values nu.
    double complex res;
    if (nu < 1 && fabs(nu / 2. - nearbyint(nu / 2.)) < EPS) {
        if (dot(dim, x_t2, x_t2) == 0 && nu == 0) {
            res = -1 * cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2));
        } else {
            res = 0;
        }
    } else if (fabs(nu - dim) < EPS && equalsZero(dim, y_t2)) {
        res = NAN;

    } else {
        double zArgBound = assignzArgBound(nu);
        double complex s1, s2, nc;
        double complex rot = 1;
        double complex xfactor = 1;
        double *vx = malloc(dim * sizeof(double));
        for (int i = 0; i < dim; i++)
            vx[i] = x_t1[i] - x_t2[i];
        xfactor = cexp(-2 * M_PI * I * dot(dim, vx, y_t1));
        free(vx);
        if (reg) {
            // calculate regularized Epstein Zeta function values.
            nc = crandall_gReg(dim, dim - nu, y_t1, lambda);
            rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));
            s2 = sum_fourier(nu, dim, lambda, m_fourier, x_t1, y_t2, cutoffsFourier,
                             zArgBound);
            // correct wrong zero summand in regularized fourier sum.
            if (!equals(dim, y_t1, y_t2)) {
                s2 += crandall_g(dim, dim - nu, y_t2, lambda, zArgBound) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2)) -
                      crandall_g(dim, dim - nu, y_t1, lambda, zArgBound) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t1));
            }
            s2 = s2 * rot + nc;
            s1 = sum_real(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                          zArgBound) *
                 rot * xfactor;
            xfactor = 1;
        } else {
            // calculate non regularized Epstein Zeta function values.
            nc = crandall_g(dim, dim - nu, y_t2, lambda, zArgBound) *
                 cexp(-2 * M_PI * I * dot(dim, x_t2, y_t2));
            s1 = sum_real(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                          zArgBound);
            s2 = sum_fourier(nu, dim, lambda, m_fourier, x_t2, y_t2, cutoffsFourier,
                             zArgBound) +
                 nc;
        }
        res = xfactor * pow(lambda * lambda / M_PI, -nu / 2.) / tgamma(nu / 2.) *
              (s1 + pow(lambda, dim) * s2);
    }
    free(x_t1);
    free(y_t1);
    free(x_t2);
    free(y_t2);
    free(m_real);
    free(m_fourier);
    return pow(ms, nu) * res;
}
#undef G_BOUND
