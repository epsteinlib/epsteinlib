// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file zeta.c
 * @brief Calculates the (regularized) Epstein zeta function and derivavies of the
 * set zeta function.
 */

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "crandall.h"
#include "tools.h"

#include "zeta.h"
#include <assert.h>

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
 * @brief epsilon for the cutoff around x = 0 and y = 0.
 * It holds that M_PI * EPS_ZERO_Y = EPS_ZERO_PIY
 */
#define EPS_ZERO_Y 1e-64

/*!
 * @brief Epsilon to catch exact cancellation to zero in inner sum of the singular
 * sum in real space for the set Zeta derivatives.
 */
#define EPS_CANCELLATION 1e-16

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
static double complex sum_real(double nu, unsigned int dim, double lambda,
                               const double *m, const double *x, const double *y,
                               const int cutoffs[], double zArgBound) {
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
static double complex sum_real_der(double nu, unsigned int dim, double lambda,
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
            unsigned int alphai = alpha[i];
            if (alphai) {
                mon *= imaginary_int_pow(alphai) *
                       real_int_pow(-2 * M_PI * lv[i], alphai);
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
 * @brief calculates the an harmonic polynomial applied to the summands in real
 * space.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return helper function for the first sum in crandalls formula. Calculates
 * sum_{z in m whole_numbers ** dim} h₍α,kIndex₎(y) * G_{nu}(z - x) *
 * exp(-2 * PI * I * (z * y)
 */
static double complex sum_real_harmonic(
    double nu, unsigned int kIndex, unsigned int dim, const double *m,
    const double *x, const double *y, const int cutoffs[], double zArgBound,
    unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {

    double lambda = 1.; // parameter that decides the weight of each sum

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

    double h;

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

        h = harmonic_h(kIndex, dim, lv, alphaAbs, chunk_offset, valid_count, coeffs,
                       exponents);

        // summing using Kahan's method
        auxy = rot * h * crandall_g(dim, nu, lv, 1. / lambda, zArgBound) - epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    return sum;
}

/**
 * @brief calculates the an harmonic polynomial applied to the summands in real
 * space. Computes (- lower Crandall) instead of upper Crandall for |z| = 0 and |z|
 * = 1.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return helper function for the first sum in crandalls formula. Calculates
 * sum_{z in m whole_numbers ** dim} I ** (|α| - 2k) h₍α,kIndex₎(y) * G_{nu}(z - x) *
 * exp(-2 * PI * I * (z * y)
 */
static double complex sum_real_harmonic_large_exp(
    double nu, unsigned int kIndex, unsigned int dim, const double *m,
    const double *x, const double *y, const int cutoffs[], double zArgBound,
    unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {

    double lambda = 1.; // parameter that decides the weight of each sum

    int zv[dim]; // counting vector in Z^dim
    unsigned int zv_1_norm = 0;
    for (int i = 0; i < dim; i++) {
        zv[i] = -cutoffs[i];
        zv_1_norm += cutoffs[i];
    }

    double lv[dim]; // lattice vector

    double complex crandall;

    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;
    double complex rot;

    double h;

    unsigned int specialCase;
    int done = 0;

    // Iterate over every over zv in the whole numbers, so that
    // - cutoff[i] <= zv[i] <= cutoffs[i]
    while (1) {

        matrix_intVector(dim, m, zv, lv);
        rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] - x[i];
        }

        h = harmonic_h(kIndex, dim, lv, alphaAbs, chunk_offset, valid_count, coeffs,
                       exponents);

        if (h) {

            // use lower Crandall for origin and its the nearest neighbors
            specialCase = (zv_1_norm <= 1);
            if (specialCase) {
                crandall = -crandall_g_lower(dim, nu, lv, 1. / lambda);
            } else {
                crandall = crandall_g(dim, nu, lv, 1. / lambda, zArgBound);
            }
            // summing using Kahan's method
            auxy = rot * h * crandall - epsilon;
            auxt = sum + auxy;
            epsilon = (auxt - sum) - auxy;
            sum = auxt;
        }

        done = 1;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (zv[idx] + 1 <= cutoffs[idx]) {
                zv[idx]++;
                zv_1_norm += (zv[idx] > 0) ? 1 : -1;
                done = 0;
                break;
            }
            zv[idx] = -cutoffs[idx];
        }
        if (done) {
            break;
        }
    }

    return sum;
}

/**
 * @brief calculates the sum over the harmonic polynomials applied to the singularity
 * summands around the origin in real space.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kMax: k sum limit ⌊|α|/2⌋.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return (-2πi)^{|α|} ∑_{z ∈ Aℤᵈ − x, |z| = 0 or |z| = 1} e^(−2πi y·z)
 * ∑_{k=0}^{⌊|α|/2⌋} h_{α,k}(−2πi z) (z·z)^{−(ν−2k)/2}
 */
static double complex sum_real_harmonic_large_exp_singularity_sum( // NOLINT
    double nu, unsigned int kMax, unsigned int dim, const double *m, const double *x,
    const double *y, unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {

    int zv[dim];
    for (int i = 0; i < dim; i++) {
        zv[i] = -1;
    }
    double lv[dim];
    double complex sum = 0.0;
    double complex epsilon = 0.0;
    double complex auxt;
    double complex auxy;
    double complex rot;
    double h;
    double nuIt;
    int done = 0;
    unsigned int zv_1_norm;
    unsigned int specialCase;

    while (1) {

        // compute 1-norm
        zv_1_norm = 0;
        for (unsigned int i = 0; i < dim; i++) {
            zv_1_norm += abs(zv[i]);
        }

        // only iterate over origin and its the nearest neighbors
        specialCase = (zv_1_norm <= 1);
        if (specialCase) {

            matrix_intVector(dim, m, zv, lv);
            rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
            for (int i = 0; i < dim; i++) {
                lv[i] = lv[i] - x[i];
            }

            double lvSquared = dot(dim, lv, lv);

            double sumInner = 0.0;
            double epsilonInner = 0.0;
            double auxtInner;
            double auxyInner;

            for (unsigned int k = 0; k <= kMax; k++) {

                nuIt = nu - (2 * k);

                // These special summands are already handled in the real sum
                if (nuIt < -1 && fabs((nuIt / 2.) - nearbyint(nuIt / 2.)) < EPS) {
                    continue;
                }

                if (fabs(nuIt / 2) < EPS) {
                    continue;
                }

                h = harmonic_h(k, dim, lv, alphaAbs, chunk_offset, valid_count,
                               coeffs, exponents);

                // summing using Kahan's method
                if (h && lvSquared > EPS_ZERO_Y) {
                    auxyInner = h * real_int_pow(lvSquared, k) - epsilonInner;
                    auxtInner = sumInner + auxyInner;
                    epsilonInner = (auxtInner - sumInner) - auxyInner;
                    sumInner = auxtInner;
                }
            }

            if (fabs(sumInner) > EPS_CANCELLATION) {

                auxy = rot * sumInner * pow(lvSquared, -nu / 2) - epsilon;
                auxt = sum + auxy;
                epsilon = (auxt - sum) - auxy;
                sum = auxt;
            }
        }

        done = 1;
        for (unsigned int idx = 0; idx < dim; idx++) {
            if (zv[idx] + 1 <= 1) {
                zv[idx]++;
                done = 0;
                break;
            }
            zv[idx] = -1;
        }
        if (done) {
            break;
        }
    }

    return sum;
}

/**
 * @brief calculates the an harmonic polynomial applied to the summands in real
 * space.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @return helper function for the first sum in crandalls formula. Calculates
 * sum_{z in m whole_numbers ** dim} I ** (|α| - 2k) h₍α,kIndex₎(y) * G_{nu}(z - x) *
 * exp(-2 * PI * I * (z * y)
 */
static double complex sum_real_harmonic_1D(double nu, unsigned int dim,
                                           const double *m, const double *x,
                                           const double *y, const int cutoffs[],
                                           double zArgBound, unsigned int alphaAbs) {

    double lambda = 1.; // parameter that decides the weight of each sum

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

    double h;

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

        h = harmonic_h_1D_kMax(lv, alphaAbs);

        auxy = rot * h * crandall_g(dim, nu, lv, 1. / lambda, zArgBound) - epsilon;
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
static double complex sum_fourier(double nu, unsigned int dim, double lambda,
                                  const double *m_invt, const double *x,
                                  const double *y, const int cutoffs[],
                                  double zArgBound) {
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
    double complex rot;
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
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
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
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
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
static double complex sum_fourier_der(double nu, unsigned int dim, double lambda,
                                      const double *m_invt, const double *x,
                                      const double *y, const int cutoffs[],
                                      double zArgBound, const unsigned int *alpha,
                                      unsigned int alphaAbs) {
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
    double complex rot;
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
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot * crandall_g_der(dim, dim - nu, lv, lambda, zArgBound, alpha,
                                    alphaAbs) -
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
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot * crandall_g_der(dim, dim - nu, lv, lambda, zArgBound, alpha,
                                    alphaAbs) -
               epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    return sum;
}

/**
 * @brief calculates the harmonic polynomials applied to the sum in fourier space in
 * 1D.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @return helper function for the second sum in crandalls formula. Calculates
 * sum_{k in m_invt whole_numbers ** dim without zero} h₍α,kIndex₎(y + k) G_{dim - nu
 * + 2 * |α| - 4 * kIndex}(k + y) * exp(-2 * PI * I * x * (k + y))
 */
static double complex sum_fourier_harmonic_1D(double nu, unsigned int dim,
                                              const double *m_invt, const double *x,
                                              const double *y, const int cutoffs[],
                                              double zArgBound,
                                              unsigned int alphaAbs) {
    double lambda = 1.; // parameter that decides the weight of each sum

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
    double complex rot;
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
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        // symmetric addition of terms to catch same - same cancellation errors
        auxy = rot * harmonic_h_1D_kMax(lv, alphaAbs) *
               crandall_g(dim, dim - nu, lv, lambda, zArgBound);
        for (int k = 0; k < dim; k++) {
            zv[k] = (((int)((totalSummands - 1 - n) / totalCutoffs[k])) %
                     (2 * cutoffs[k] + 1)) -
                    cutoffs[k];
        }
        matrix_intVector(dim, m_invt, zv, lv);
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] + y[i];
        }
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy += rot * harmonic_h_1D_kMax(lv, alphaAbs) *
                crandall_g(dim, dim - nu, lv, lambda, zArgBound);
        auxy -= epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
    }
    return sum;
}

/**
 * @brief calculates the harmonic polynomials applied to the sum in fourier space.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: projection of x vector to elementary lattice cell.
 * @param[in] y: projection of y vector to elementary lattice cell.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return helper function for the second sum in crandalls formula. Calculates
 * sum_{k in m_invt whole_numbers ** dim without zero} h₍α,kIndex₎(y + k) G_{dim - nu
 * + 2 * |α| - 4 * kIndex}(k + y) * exp(-2 * PI * I * x * (k + y))
 */
static double complex sum_fourier_harmonic(
    double nu, unsigned int kIndex, unsigned int dim, const double *m_invt,
    const double *x, const double *y, const int cutoffs[], double zArgBound,
    unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {
    double lambda = 1.; // parameter that decides the weight of each sum

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
    double complex rot;
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
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot *
                   harmonic_h(kIndex, dim, lv, alphaAbs, chunk_offset, valid_count,
                              coeffs, exponents) *
                   crandall_g(dim, dim - nu, lv, lambda, zArgBound) -
               epsilon;
        auxt = sum + auxy;
        epsilon = (auxt - sum) - auxy;
        sum = auxt;
        // symmerict addition to catch +- identical terms in y = 0
        for (int k = 0; k < dim; k++) {
            zv[k] = (((int)((totalSummands - 1 - n) / totalCutoffs[k])) %
                     (2 * cutoffs[k] + 1)) -
                    cutoffs[k];
        }
        matrix_intVector(dim, m_invt, zv, lv);
        for (int i = 0; i < dim; i++) {
            lv[i] = lv[i] + y[i];
        }
        rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
        auxy = rot *
                   harmonic_h(kIndex, dim, lv, alphaAbs, chunk_offset, valid_count,
                              coeffs, exponents) *
                   crandall_g(dim, dim - nu, lv, lambda, zArgBound) -
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
static double *vectorProj(unsigned int dim, const double *m, const double *m_invt,
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
 * @brief calculates set zeta derivatives by the 1-D harmonic method.
 * @param[in] res: accumulated result from prior computation.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] alphaAbs: total |α| of the multi-index.
 * @param[in] lambda: splitting parameter.
 * @param[in] ms: lattice scaling factor, pow(vol, -1./dim).
 * @param[in] m_real: matrix that transforms the lattice (real sum).
 * @param[in] m_fourier: matrix that transforms the lattice (Fourier sum).
 * @param[in] x_t1: first projection of x vector to elementary lattice cell.
 * @param[in] x_t2: second projection of x vector to elementary lattice cell.
 * @param[in] y_t1: first projection of y vector to elementary lattice cell.
 * @param[in] y_t2: second projection of y vector to elementary lattice cell.
 * @param[in] x_t2_squared: precomputed dot(dim, x_t2, x_t2).
 * @param[in] cutoffsReal: how many summands in each direction (real sum).
 * @param[in] cutoffsFourier: how many summands in each direction (Fourier sum).
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] xfactor: precomputed prefactor for the real sum.
 * @return updated res after adding the 1-D harmonic contribution.
 */
static double complex summation_harmonic_1D(
    double complex res, double nu, unsigned int dim, unsigned int alphaAbs,
    double lambda, double ms, const double *m_real, const double *m_fourier,
    const double *x_t1, const double *x_t2, const double *y_t1, const double *y_t2,
    double x_t2_squared, const int cutoffsReal[], const int cutoffsFourier[],
    double zArgBound, double complex xfactor) {

    double complex s1;
    double complex s2;
    double complex nc;

    // Compute set zeta derivatives by harmonic method

    unsigned int k = alphaAbs / 2;
    double complex resIt;
    double nuIt = nu - (2 * k);

    // skip iterartions where nuIt is a negative even integer, as
    // 1/gamma(nIt) = 0
    if (nuIt < -1 && fabs((nuIt / 2.) - nearbyint(nuIt / 2.)) < EPS) {
        resIt = 0;
        // if nu = 0, everything except the zero summand in the real sum
        // vanishes
    } else if (fabs(nuIt / 2) < EPS) {
        if (x_t2_squared > EPS_ZERO_Y) {
            resIt = 0.;
        } else {
            resIt =
                -imaginary_int_pow(alphaAbs) * harmonic_h_1D_kMax(x_t2, alphaAbs);
        }
    } else {

        double nuReci = nuIt - (2 * alphaAbs) + (4 * k);
        double zArgBoundReci = assignzArgBound(dim - nuReci);

        // calculate set zeta derivative function values.
        double complex rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));

        if (equals(dim, y_t1, y_t2)) {
            nc = harmonic_h_1D_kMax(y_t1, alphaAbs) *
                 crandall_g(dim, dim - nuReci, y_t1, lambda, zArgBoundReci);
        } else {
            nc = harmonic_h_1D_kMax(y_t2, alphaAbs) *
                 crandall_g(dim, dim - nuReci, y_t2, lambda, zArgBoundReci) *
                 cexp(-2 * M_PI * I * dot(dim, y_t2, x_t1)) * rot;
        }

        s2 = sum_fourier_harmonic_1D(nuReci, dim, m_fourier, x_t1, y_t2,
                                     cutoffsFourier, zArgBoundReci, alphaAbs);

        s2 = negative_one_pow(k) * (s2 * rot + nc);

        s1 = sum_real_harmonic_1D(nuIt, dim, m_real, x_t2, y_t2, cutoffsReal,
                                  zArgBound, alphaAbs) *
             imaginary_int_pow(alphaAbs) * rot * xfactor;
        resIt = (s1 + pow(lambda, dim) * s2) / tgamma(nuIt / 2.);
    }
    resIt *= pow(lambda * lambda / M_PI, -nuIt / 2.);

    res += resIt;
    res *= real_int_pow(-2 * M_PI, alphaAbs) / real_int_pow(ms, alphaAbs);

    return res;
}

/**
 * @brief calculates set zeta derivatives by the harmonic method for large nu,
 * isolating singularities in the Fourier sum.
 * @param[in] res: accumulated result from prior computation.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] alphaAbs: total |α| of the multi-index.
 * @param[in] alpha: multi-index α (length dim).
 * @param[in] lambda: splitting parameter.
 * @param[in] ms: lattice scaling factor, pow(vol, -1./dim).
 * @param[in] m_real: matrix that transforms the lattice (real sum).
 * @param[in] m_fourier: matrix that transforms the lattice (Fourier sum).
 * @param[in] x_t1: first projection of x vector to elementary lattice cell.
 * @param[in] x_t2: second projection of x vector to elementary lattice cell.
 * @param[in] y_t1: first projection of y vector to elementary lattice cell.
 * @param[in] y_t2: second projection of y vector to elementary lattice cell.
 * @param[in] x_t2_squared: precomputed dot(dim, x_t2, x_t2).
 * @param[in] cutoffsReal: how many summands in each direction (real sum).
 * @param[in] cutoffsFourier: how many summands in each direction (Fourier sum).
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] xfactor: precomputed prefactor for the real sum.
 * @return updated res after adding the high-exponent harmonic contribution.
 */
static double complex summation_harmonic_high_exp(
    double complex res, double nu, unsigned int dim, unsigned int alphaAbs,
    const unsigned int *alpha, double lambda, double ms, const double *m_real,
    const double *m_fourier, const double *x_t1, const double *x_t2,
    const double *y_t1, const double *y_t2, double x_t2_squared,
    const int cutoffsReal[], const int cutoffsFourier[], double zArgBound,
    double complex xfactor) {
    // Compute set zeta derivatives by harmonic method for large nu,
    // Isolate singularities in fourier sum

    double complex s1;
    double complex s2;
    double complex nc;

    // Precompute coefficients for harmonic polynomials once
    unsigned int kMax = alphaAbs / 2;
    unsigned long long *chunk_offset =
        malloc((kMax + 1) * sizeof(unsigned long long));
    unsigned long long *valid_count =
        malloc((kMax + 1) * sizeof(unsigned long long));
    unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
        alphaAbs, kMax, dim, alpha, chunk_offset, valid_count);
    double *coeffs = malloc(coeffs_size * sizeof(double));
    unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
    precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_offset, coeffs,
                                    exponents);

    double nuIt;
    double nuReci;
    double complex resIt;

    double complex rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));

    // Main loop without singularities around the origin in real sum
    for (unsigned int k = 0; k <= kMax; k++) {

        nuIt = nu - (2 * k);

        // skip iterartions where nuIt is a negative even integer, as
        // 1/gamma(nIt) = 0
        if (nuIt < -1 && fabs((nuIt / 2.) - nearbyint(nuIt / 2.)) < EPS) {
            continue;
        }

        // if nu = 0, everything except the zero summand in the real sum
        // vanishes
        if (fabs(nuIt / 2) < EPS) {
            if (x_t2_squared > EPS_ZERO_Y) {
                resIt = 0.;
            } else {
                resIt = -imaginary_int_pow(alphaAbs) *
                        harmonic_h(k, dim, x_t2, alphaAbs, chunk_offset, valid_count,
                                   coeffs, exponents);
            }
        } else {
            nuReci = nuIt - (2 * alphaAbs) + (4 * k);
            double zArgBoundReci = assignzArgBound(dim - nuReci);

            // calculate set zeta derivative function values.

            // skip zero summand if harmonic polynomial vanishes
            double h = harmonic_h(k, dim, y_t2, alphaAbs, chunk_offset, valid_count,
                                  coeffs, exponents);
            if (h) {
                if (equals(dim, y_t1, y_t2)) {
                    nc = h *
                         crandall_g(dim, dim - nuReci, y_t1, lambda, zArgBoundReci);
                } else {
                    nc = h *
                         crandall_g(dim, dim - nuReci, y_t2, lambda, zArgBoundReci) *
                         cexp(-2 * M_PI * I * dot(dim, y_t2, x_t1)) * rot;
                }
            } else {
                nc = 0;
            }

            s2 = sum_fourier_harmonic(nuReci, k, dim, m_fourier, x_t1, y_t2,
                                      cutoffsFourier, zArgBoundReci, alphaAbs,
                                      chunk_offset, valid_count, coeffs, exponents);

            s2 = negative_one_pow(k) * (s2 * rot + nc);

            s1 = sum_real_harmonic_large_exp(
                     nuIt, k, dim, m_real, x_t2, y_t2, cutoffsReal, zArgBound,
                     alphaAbs, chunk_offset, valid_count, coeffs, exponents) *
                 imaginary_int_pow(alphaAbs) * rot * xfactor;

            resIt = (s1 + pow(lambda, dim) * s2) / tgamma(nuIt / 2.);
        }
        resIt *= pow(lambda * lambda / M_PI, -nuIt / 2.);
        res += resIt;
    }

    double complex s1_singularity =
        xfactor * rot * imaginary_int_pow(alphaAbs) *
        sum_real_harmonic_large_exp_singularity_sum(nu, kMax, dim, m_real, x_t2,
                                                    y_t2, alphaAbs, chunk_offset,
                                                    valid_count, coeffs, exponents);

    res += s1_singularity;

    res *= real_int_pow(-2 * M_PI, alphaAbs) / real_int_pow(ms, alphaAbs);

    free(chunk_offset);
    free(valid_count);
    free(coeffs);
    free(exponents);

    return res;
}

/**
 * @brief calculates set zeta derivatives by the general harmonic method.
 * @param[in] res: accumulated result from prior computation.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] alphaAbs: total |α| of the multi-index.
 * @param[in] alpha: multi-index α (length dim).
 * @param[in] lambda: splitting parameter.
 * @param[in] ms: lattice scaling factor, pow(vol, -1./dim).
 * @param[in] m_real: matrix that transforms the lattice (real sum).
 * @param[in] m_fourier: matrix that transforms the lattice (Fourier sum).
 * @param[in] x_t1: first projection of x vector to elementary lattice cell.
 * @param[in] x_t2: second projection of x vector to elementary lattice cell.
 * @param[in] y_t1: first projection of y vector to elementary lattice cell.
 * @param[in] y_t2: second projection of y vector to elementary lattice cell.
 * @param[in] cutoffsReal: how many summands in each direction (real sum).
 * @param[in] cutoffsFourier: how many summands in each direction (Fourier sum).
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] xfactor: precomputed prefactor for the real sum.
 * @return updated res after adding the general harmonic contribution.
 */
static double complex summation_harmonic(
    double complex res, double nu, unsigned int dim, unsigned int alphaAbs,
    const unsigned int *alpha, double lambda, double ms, const double *m_real,
    const double *m_fourier, const double *x_t1, const double *x_t2,
    const double *y_t1, const double *y_t2, const int cutoffsReal[],
    const int cutoffsFourier[], double zArgBound, double complex xfactor) {

    // Compute set zeta derivatives by harmonic method

    double complex s1;
    double complex s2;
    double complex nc;

    // Precompute coefficients for harmonic polynomials once
    unsigned int kMax = alphaAbs / 2;
    unsigned long long *chunk_offset =
        malloc((kMax + 1) * sizeof(unsigned long long));
    unsigned long long *valid_count =
        malloc((kMax + 1) * sizeof(unsigned long long));
    unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
        alphaAbs, kMax, dim, alpha, chunk_offset, valid_count);
    double *coeffs = malloc(coeffs_size * sizeof(double));
    unsigned int *exponents = malloc(coeffs_size * dim * sizeof(unsigned int));
    precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, chunk_offset, coeffs,
                                    exponents);

    double nuIt;
    double nuReci;
    double complex resIt;

    for (unsigned int k = 0; k <= kMax; k++) {

        nuIt = nu - (2 * k);

        // skip iterartions where nuIt is a negative even integer, as
        // 1/gamma(nIt) = 0
        if (nuIt < -1 && fabs((nuIt / 2.) - nearbyint(nuIt / 2.)) < EPS) {
            continue;
        }

        // if nu = 0, everything except the zero summand in the real sum
        // vanishes
        if (fabs(nuIt / 2) < EPS) {
            if (dot(dim, x_t2, x_t2) > EPS_ZERO_Y) {
                resIt = 0.;
            } else {
                resIt = -imaginary_int_pow(alphaAbs) *
                        harmonic_h(k, dim, x_t2, alphaAbs, chunk_offset, valid_count,
                                   coeffs, exponents);
            }
        } else {
            nuReci = nuIt - (2 * alphaAbs) + (4 * k);
            double zArgBoundReci = assignzArgBound(dim - nuReci);

            // calculate set zeta derivative function values.
            double complex rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));

            // skip zero summand if harmonic polynomial vanishes
            double h = harmonic_h(k, dim, y_t2, alphaAbs, chunk_offset, valid_count,
                                  coeffs, exponents);
            if (h) {
                if (equals(dim, y_t1, y_t2)) {
                    nc = h *
                         crandall_g(dim, dim - nuReci, y_t1, lambda, zArgBoundReci);
                } else {
                    nc = h *
                         crandall_g(dim, dim - nuReci, y_t2, lambda, zArgBoundReci) *
                         cexp(-2 * M_PI * I * dot(dim, y_t2, x_t1)) * rot;
                }
            } else {
                nc = 0;
            }

            s2 = sum_fourier_harmonic(nuReci, k, dim, m_fourier, x_t1, y_t2,
                                      cutoffsFourier, zArgBoundReci, alphaAbs,
                                      chunk_offset, valid_count, coeffs, exponents);

            s2 = negative_one_pow(k) * (s2 * rot + nc);

            s1 = sum_real_harmonic(nuIt, k, dim, m_real, x_t2, y_t2, cutoffsReal,
                                   zArgBound, alphaAbs, chunk_offset, valid_count,
                                   coeffs, exponents) *
                 imaginary_int_pow(alphaAbs) * rot * xfactor;

            resIt = (s1 + pow(lambda, dim) * s2) / tgamma(nuIt / 2.);
        }
        resIt *= pow(lambda * lambda / M_PI, -nuIt / 2.);
        res += resIt;
    }

    res *= real_int_pow(-2 * M_PI, alphaAbs) / real_int_pow(ms, alphaAbs);

    free(chunk_offset);
    free(valid_count);
    free(coeffs);
    free(exponents);

    return res;
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
    // Early return for 0th derivative special cases
    unsigned int alphaAbs = (variant > 1) ? mult_abs(dim, alpha) : 0;
    if (variant == 2 && !alphaAbs) {
        return cexp(2 * M_PI * I * dot(dim, x, y)) *
               epsteinZetaInternal(nu, dim, m, x, y, 1, 0, alpha);
    }
    if (variant == 3 && !alphaAbs) {
        return epsteinZetaInternal(nu, dim, m, x, y, 1, 1, alpha);
    }
    // Odd derivatives in zero
    if (variant == 2) {
        for (int i = 0; i < dim; i++) {
            if (alpha[i] % 2 && fabs(x[i]) < EPS_ZERO_Y && fabs(y[i]) < EPS_ZERO_Y) {
                return 0.;
            }
        }
    }
    // 1. Transform: Compute determinant and fourier transformed matrix, scale
    // both of them
    assert(dim > 0);
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
    // choose harmonic method for setZetaDer in most cases
    bool harmonicMethod = (variant == 2);
    bool harmonicMethod1D = harmonicMethod && (dim == 1);
    bool harmonicMethodHighExp = harmonicMethod && (nu > 10);
    // if variant == 2 and all of the above false, polynomial p method is used
    // works for low derivatives and arbitrary nu, slightly
    // faster than harmonic Methods in d > 1

    // handle special case of non-positive integer values nu.
    double complex res = 0.;
    double x_t2_squared = dot(dim, x_t2, x_t2);
    double y_t2_squared = dot(dim, y_t2, y_t2);
    if (variant < 3 && nu < 1 && fabs((nu / 2.) - nearbyint(nu / 2.)) < EPS) {
        if (variant < 2 && x_t2_squared == 0 && nu == 0) {
            res = -1 * cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2));
        } else {
            res = 0;
        }
    } else if ((variant == 0 || variant == 2) && fabs(nu - dim - alphaAbs) < EPS &&
               y_t2_squared < EPS_ZERO_Y) {
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
            if (harmonicMethod1D) {
                res = summation_harmonic_1D(res, nu, dim, alphaAbs, lambda, ms,
                                            m_real, m_fourier, x_t1, x_t2, y_t1,
                                            y_t2, x_t2_squared, cutoffsReal,
                                            cutoffsFourier, zArgBound, xfactor);
            } else if (harmonicMethodHighExp) {
                res = summation_harmonic_high_exp(
                    res, nu, dim, alphaAbs, alpha, lambda, ms, m_real, m_fourier,
                    x_t1, x_t2, y_t1, y_t2, x_t2_squared, cutoffsReal,
                    cutoffsFourier, zArgBound, xfactor);
            } else if (harmonicMethod) {
                res = summation_harmonic(res, nu, dim, alphaAbs, alpha, lambda, ms,
                                         m_real, m_fourier, x_t1, x_t2, y_t1, y_t2,
                                         cutoffsReal, cutoffsFourier, zArgBound,
                                         xfactor);
            } else {
                // Compute set zeta derivatives by polynomial p method
                rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));
                if (equals(dim, y_t1, y_t2)) {
                    nc = crandall_g_der(dim, dim - nu, y_t1, lambda, zArgBoundReci,
                                        alpha, alphaAbs);
                } else {
                    nc = crandall_g_der(dim, dim - nu, y_t2, lambda, zArgBoundReci,
                                        alpha, alphaAbs) *
                         cexp(-2 * M_PI * I * dot(dim, y_t2, x_t1)) * rot;
                }
                s2 = sum_fourier_der(nu, dim, lambda, m_fourier, x_t1, y_t2,
                                     cutoffsFourier, zArgBoundReci, alpha, alphaAbs);
                s2 = real_int_pow(lambda, alphaAbs) * (s2 * rot + nc);
                s1 = sum_real_der(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                                  zArgBound, alpha) *
                     rot * xfactor;
                xfactor = 1. / real_int_pow(ms, alphaAbs);
            }
        } else if (variant == 3) {
            // calculate Epstein zeta reg derivative function values.
            nc = crandall_gReg_der(dim, dim - nu, y_t1, lambda, alpha, alphaAbs,
                                   zArgBoundReci);
            rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));
            s2 = sum_fourier_der(nu, dim, lambda, m_fourier, x_t1, y_t2,
                                 cutoffsFourier, zArgBoundReci, alpha, alphaAbs);
            // correct wrong zero summand in regularized fourier sum.
            if (!equals(dim, y_t1, y_t2)) {
                s2 += crandall_g_der(dim, dim - nu, y_t2, lambda, zArgBoundReci,
                                     alpha, alphaAbs) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2)) -
                      crandall_g_der(dim, dim - nu, y_t1, lambda, zArgBoundReci,
                                     alpha, alphaAbs) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t1));
            }
            s2 = real_int_pow(lambda, alphaAbs) * (s2 * rot + nc);
            s1 = sum_real_der(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                              zArgBound, alpha) *
                 rot * xfactor;
            xfactor = 1. / real_int_pow(ms, alphaAbs);
        }
        // In the harmonic method, the res is already set as there is no global
        // nu-dependent coefficient there
        if (!(harmonicMethod)) {
            res = xfactor * pow(lambda * lambda / M_PI, -nu / 2.) / tgamma(nu / 2.) *
                  (s1 + pow(lambda, dim) * s2);
        }
    }
    free(x_t2);
    free(y_t2);
    res *= pow(ms, nu);
    // apply correction to matrix scaling if nu = d + 2k
    unsigned int k = (unsigned int)fmax(0., nearbyint((nu - (double)dim) / 2));
    if ((variant == 1 || variant == 3) && (nu == (dim + 2 * (double)k))) {
        if (alphaAbs) {
            res -= pow(M_PI, (double)k + ((double)dim / 2)) /
                   tgamma((double)k + ((double)dim / 2)) * negative_one_pow(k + 1) *
                   polynomial_y_der(k, dim, y, alpha, alphaAbs, k) * log(ms * ms) /
                   vol;
        } else {
            if (k) {
                double ySquared = 0;
                for (int i = 0; i < dim; i++) {
                    ySquared += y[i] * y[i];
                }

                res -= pow(M_PI, (double)(2 * k) + ((double)dim / 2)) *
                       tgamma((double)k + ((double)dim / 2)) *
                       negative_one_pow(k + 1) / tgamma((double)k + 1) *
                       real_int_pow(ySquared, k) * log(ms * ms) / vol;
            } else {
                res += pow(M_PI, (double)dim / 2) / tgamma((double)dim / 2) *
                       log(ms * ms) / vol;
            }
        }
    }
    return res;
}
#undef G_BOUND
