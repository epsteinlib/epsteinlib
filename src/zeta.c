// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024-2026 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file zeta.c
 * @brief Calculates the (regularized) Epstein zeta function and the (regularized)
 * anisotropic Epstein zeta function.
 */

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "crandall.h"
#include "harmonics.h"
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
 * @brief epsilon for the cutoff around x = 0 and y = 0.
 * It holds that M_PI * EPS_ZERO_Y = EPS_ZERO_PIY
 */
#define EPS_ZERO_Y 1e-64

/*!
 * @brief Epsilon to catch exact cancellation to zero in inner sum of the singular
 * sum in real space for the set Zeta derivatives.
 */
#define EPS_CANCELLATION 4e-16

/**
 * @brief Increments the integer lattice vector to the next lattice point.
 * @param[in] dim: dimension of the lattice.
 * @param[in] cutoffs: number of summands in each direction.
 * @param[in,out] zv: integer lattice vector, incremented in-place.
 */
static inline void lattice_vector_increment(unsigned int dim, const int cutoffs[],
                                            int zv[]) {
    for (int k = 0; k < dim; k++) {
        if (++zv[k] <= cutoffs[k]) {
            break;
        }
        zv[k] = -cutoffs[k];
    }
}

/**
 * @brief Increments the integer lattice vector to the next lattice point and
 * tracks the L1 norm.
 * @param[in] dim: dimension of the lattice.
 * @param[in] cutoffs: number of summands in each direction.
 * @param[in,out] zv: integer lattice vector, incremented in-place.
 * @param[in,out] zv_1_norm: L1 norm of zv, updated in-place.
 */
static inline void lattice_vector_increment_norm(unsigned int dim,
                                                 const int cutoffs[], int zv[],
                                                 unsigned int *restrict zv_1_norm) {
    for (int k = 0; k < dim; k++) {
        if (++zv[k] <= cutoffs[k]) {
            *zv_1_norm += (zv[k] > 0) ? 1 : -1;
            break;
        }
        zv[k] = -cutoffs[k];
    }
}

/**
 * @brief Computes one summand of the first sum in Crandall's formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameter that decides the weight of each sum.
 * @param[in,out] lv: lattice vector, shifted by x in-place.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @return rot * G_{nu}(lv - x), the real space summand.
 */
static inline double complex summand_real(double nu, unsigned int dim, double lambda,
                                          double lv[], const double *x,
                                          const double *y, double zArgBound) {
    double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
    for (int i = 0; i < dim; i++) {
        lv[i] -= x[i];
    }
    return rot * crandall_g(dim, nu, lv, 1. / lambda, zArgBound);
}

/**
 * @brief calculates the first sum in Crandall's formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameters that decides the weight of each sum.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * @param[in] diag: 1 iff the lattice matrix is diagonal.
 * the incomplete gamma evaluation.
 * @return helper function for the first sum in crandalls formula. Calculates
 * sum_{z in m whole_numbers ** dim} G_{nu}((z - x) / lambda))
 * X exp(-2 * PI * I * z * y)
 */
static double complex sum_real(double nu, unsigned int dim, double lambda,
                               const double *m, const double *x, const double *y,
                               const int cutoffs[], double zArgBound, bool diag) {
    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    for (int k = 0; k < dim; k++) {
        zv[k] = -cutoffs[k]; // lattice vector initialized
        totalSummands *= 2 * cutoffs[k] + 1;
    }
    double complex sum = 0.0;
    double complex epsilon = 0.0;

    // Sum in real space
    for (long n = 0; n < totalSummands; n++) {
        matrix_intVector(dim, m, zv, lv, diag);
        double complex summand = summand_real(nu, dim, lambda, lv, x, y, zArgBound);
        kahan_add_c(&sum, &epsilon, summand);
        lattice_vector_increment(dim, cutoffs, zv);
    }

    return sum;
}

/**
 * @brief Computes one summand of the harmonic method of the first sum in Crandall's
 * formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameter that decides the weight of each sum.
 * @param[in,out] lv: lattice vector, shifted by x in-place.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return h₍α,kIndex₎(y) * G_{nu}(z - x) * exp(-2 * PI * I * (z * y), the real space
 * summand of the harmonic method.
 */
static inline double complex summand_real_harmonic(
    double nu, unsigned int kIndex, unsigned int dim, double lambda, double lv[],
    const double *x, const double *y, double zArgBound, unsigned int alphaAbs,
    const unsigned long long *chunk_offset, const unsigned long long *valid_count,
    const double *coeffs, const unsigned int *exponents) {

    double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
    for (int i = 0; i < dim; i++) {
        lv[i] -= x[i];
    }
    double h = harmonic_h(kIndex, dim, lv, alphaAbs, chunk_offset, valid_count,
                          coeffs, exponents);

    return rot * h * crandall_g(dim, nu, lv, 1. / lambda, zArgBound);
}

/**
 * @brief calculates the an harmonic polynomial applied to the summands in real
 * space.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] diag: true, if the lattice matrix is diagonal.
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
    bool diag, unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {

    double lambda = 1.; // parameter that decides the weight of each sum

    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    for (int k = 0; k < dim; k++) {
        zv[k] = -cutoffs[k]; // lattice vector initialized
        totalSummands *= 2 * cutoffs[k] + 1;
    }
    double complex sum = 0.0;
    double complex epsilon = 0.0;

    // First Sum (in real space)
    for (long n = 0; n < totalSummands; n++) {
        matrix_intVector(dim, m, zv, lv, diag);
        double complex summand = summand_real_harmonic(
            nu, kIndex, dim, lambda, lv, x, y, zArgBound, alphaAbs, chunk_offset,
            valid_count, coeffs, exponents);
        kahan_add_c(&sum, &epsilon, summand);
        lattice_vector_increment(dim, cutoffs, zv);
    }

    return sum;
}

/**
 * @brief Computes one summand of the derivative of the first sum in Crandall's
 * formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameter that decides the weight of each sum.
 * @param[in,out] lv: lattice vector, shifted by x in-place.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @param[in] nearOrigin: true if zv is zero or a nearest neighbor.
 * @return I ** (|α| - 2k) h₍α,kIndex₎(y) * G_{nu}(z - x) * exp(-2 * PI * I * (z *
 * y).
 */
static inline double complex summand_real_harmonic_large_exp(
    double nu, unsigned int kIndex, unsigned int dim, double lambda, double lv[],
    const double *x, const double *y, double zArgBound, unsigned int alphaAbs,
    const unsigned long long *chunk_offset, const unsigned long long *valid_count,
    const double *coeffs, const unsigned int *exponents, bool nearOrigin) {

    double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
    for (int i = 0; i < dim; i++) {
        lv[i] -= x[i];
    }

    double h = harmonic_h(kIndex, dim, lv, alphaAbs, chunk_offset, valid_count,
                          coeffs, exponents);

    // skip for h = 0 for optimization only, as lower crandall smooth near the origin
    if (h) {
        // use lower Crandall for the origin and its the nearest neighbors
        double complex crandall;
        if (nearOrigin) {
            crandall = -crandall_g_lower(dim, nu, lv, 1. / lambda);
        } else {
            crandall = crandall_g(dim, nu, lv, 1. / lambda, zArgBound);
        }
        return rot * h * crandall;
    }

    return 0;
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
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] diag: true, if the lattice matrix is diagonal.
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
    bool diag, unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {

    double lambda = 1.; // parameter that decides the weight of each sum

    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    unsigned int zv_1_norm = 0;
    for (int k = 0; k < dim; k++) {
        zv[k] = -cutoffs[k]; // lattice vector initialized
        zv_1_norm += cutoffs[k];
        totalSummands *= 2 * cutoffs[k] + 1;
    }
    double complex sum = 0.0;
    double complex epsilon = 0.0;

    for (long n = 0; n < totalSummands; n++) {
        bool nearOrigin = (zv_1_norm <= 1);
        matrix_intVector(dim, m, zv, lv, diag);
        double complex summand = summand_real_harmonic_large_exp(
            nu, kIndex, dim, lambda, lv, x, y, zArgBound, alphaAbs, chunk_offset,
            valid_count, coeffs, exponents, nearOrigin);
        kahan_add_c(&sum, &epsilon, summand);
        lattice_vector_increment_norm(dim, cutoffs, zv, &zv_1_norm);
    }

    return sum;
}

/**
 * @brief Computes one summand of the derivative of the first sum in Crandall's
 * formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kMax: k sum limit ⌊|α|/2⌋.
 * @param[in] dim: dimension of the input vectors.
 * @param[in,out] lv: lattice vector, shifted by x in-place.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return I ** (|α| - 2k) h₍α,kIndex₎(y) * G_{nu}(z - x) * exp(-2 * PI * I * (z *
 * y).
 */
static inline double complex summand_real_harmonic_large_exp_singularity_sum(
    double nu, unsigned int kMax, unsigned int dim, double lv[], const double *x,
    const double *y, unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {

    double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, y));
    for (int i = 0; i < dim; i++) {
        lv[i] -= x[i];
    }

    double lvSquared = dot(dim, lv, lv);

    double sumInner = 0.0;
    double epsilonInner = 0.0;
    double maxInner = 0.0;

    for (unsigned int k = 0; k <= kMax; k++) {

        double nuIt = nu - (2 * k);

        // These special summands are already handled in the real sum
        if (nuIt < .5 && fabs((nuIt / 2.) - nearbyint(nuIt / 2.)) < EPS) {
            continue;
        }

        double h = harmonic_h(k, dim, lv, alphaAbs, chunk_offset, valid_count,
                              coeffs, exponents);

        if (h && lvSquared > EPS_ZERO_Y) {
            double summand = h * real_int_pow(lvSquared, k);
            kahan_add_r(&sumInner, &epsilonInner, summand);
            maxInner = fmax(maxInner, fabs(summand));
        }
    }

    if (fabs(sumInner) > EPS_CANCELLATION * (kMax + 1) * maxInner) {
        return rot * sumInner * pow(lvSquared, -nu / 2);
    }
    return 0;
}

/**
 * @brief calculates the sum over the harmonic polynomials applied to the singularity
 * summands around the origin in real space.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kMax: k sum limit ⌊|α|/2⌋.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] diag: true, if the lattice matrix is diagonal.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return (-2πi)^{|α|} ∑_{z ∈ Aℤᵈ − x, |z| = 0 or |z| = 1} e^(−2πi y·z)
 * ∑_{k=0}^{⌊|α|/2⌋} h_{α,k}(−2πi z) (z·z)^{−(ν−2k)/2}
 */
static double complex sum_real_harmonic_large_exp_singularity_sum(
    double nu, unsigned int kMax, unsigned int dim, const double *m, const double *x,
    const double *y, bool diag, unsigned int alphaAbs,
    const unsigned long long *chunk_offset, const unsigned long long *valid_count,
    const double *coeffs, const unsigned int *exponents) {

    int zv[dim];
    double lv[dim];
    int cutoffs[dim];
    long totalSummands = 1;
    // only iterate over origin and its the nearest neighbors
    for (int k = 0; k < dim; k++) {
        cutoffs[k] = 1;
        zv[k] = -cutoffs[k];
        totalSummands *= 2 * cutoffs[k] + 1;
    }
    unsigned int zv_1_norm = dim;
    double complex sum = 0.0;
    double complex epsilon = 0.0;

    for (long n = 0; n < totalSummands; n++) {
        bool nearOrigin = (zv_1_norm <= 1);
        if (nearOrigin) {
            matrix_intVector(dim, m, zv, lv, diag);
            double complex summand = summand_real_harmonic_large_exp_singularity_sum(
                nu, kMax, dim, lv, x, y, alphaAbs, chunk_offset, valid_count, coeffs,
                exponents);
            kahan_add_c(&sum, &epsilon, summand);
        }
        lattice_vector_increment_norm(dim, cutoffs, zv, &zv_1_norm);
    }

    return sum;
}

/**
 * @brief Computes one summand of the second sum in Crandall's formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameter that decides the weight of each sum.
 * @param[in,out] lv: lattice vector, shifted by y in-place.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @return rot * G_{dim - nu}(lv + y), the Fourier space summand.
 */
static inline double complex summand_fourier(double nu, unsigned int dim,
                                             double lambda, double lv[],
                                             const double *x, const double *y,
                                             double zArgBound) {
    for (int i = 0; i < dim; i++) {
        lv[i] += y[i];
    }
    double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
    return rot * crandall_g(dim, dim - nu, lv, lambda, zArgBound);
}

/**
 * @brief calculates the second sum in Crandall's formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameters that decides the weight of each sum.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * @param[in] diag: 1 iff the lattice matrix is diagonal.
 * the incomplete gamma evaluation.
 * @return helper function for the second sum in crandalls formula. Calculates
 * sum_{k in m_invt whole_numbers ** dim without zero} G_{dim - nu}(lambda * (k + y))
 * X exp(-2 * PI * I * x * (k + y))
 */
static double complex sum_fourier(double nu, unsigned int dim, double lambda,
                                  const double *m_invt, const double *x,
                                  const double *y, const int cutoffs[],
                                  double zArgBound, bool diag) {
    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    for (int k = 0; k < dim; k++) {
        zv[k] = -cutoffs[k]; // lattice vector initialized
        totalSummands *= 2 * cutoffs[k] + 1;
    };
    long zeroIndex = (totalSummands - 1) / 2;
    double complex sum = 0.0;
    double complex epsilon = 0.0;
    // second sum (in fourier space)
    for (long n = 0; n < zeroIndex; n++) {
        matrix_intVector(dim, m_invt, zv, lv, diag);
        double complex summand =
            summand_fourier(nu, dim, lambda, lv, x, y, zArgBound);
        kahan_add_c(&sum, &epsilon, summand);
        lattice_vector_increment(dim, cutoffs, zv);
    }
    lattice_vector_increment(dim, cutoffs, zv); // skips zero
    for (long n = zeroIndex + 1; n < totalSummands; n++) {
        matrix_intVector(dim, m_invt, zv, lv, diag);
        double complex summand =
            summand_fourier(nu, dim, lambda, lv, x, y, zArgBound);
        kahan_add_c(&sum, &epsilon, summand);
        lattice_vector_increment(dim, cutoffs, zv);
    }
    return sum;
}

/**
 * @brief Computes one summand of the harmonic method of the second sum in Crandall's
 * formula.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] lambda: parameter that decides the weight of each sum.
 * @param[in,out] lv: lattice vector, shifted by y in-place.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return h₍α,kIndex₎(y + k) G_{dim - nu + 2 * |α| - 4 * kIndex}(k + y) *
 * exp(-2 * PI * I * x * (k + y))
 */
static inline double complex summand_fourier_harmonic(
    double nu, unsigned int kIndex, unsigned int dim, double lambda, double lv[],
    const double *x, const double *y, double zArgBound, unsigned int alphaAbs,
    const unsigned long long *chunk_offset, const unsigned long long *valid_count,
    const double *coeffs, const unsigned int *exponents) {
    for (int i = 0; i < dim; i++) {
        lv[i] += y[i];
    }
    double complex rot = cexp(-2 * M_PI * I * dot(dim, lv, x));
    return rot *
           harmonic_h(kIndex, dim, lv, alphaAbs, chunk_offset, valid_count, coeffs,
                      exponents) *
           crandall_g(dim, dim - nu, lv, lambda, zArgBound);
}

/**
 * @brief calculates the harmonic polynomials applied to the sum in fourier space.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] kIndex: specifies degree |alpha| - 2k.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice in the Epstein zeta
 * function.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] cutoffs: how many summands in each direction are considered.
 * @param[in] zArgBound: global bound on when to use the asymptotic expansion in
 * the incomplete gamma evaluation.
 * @param[in] diag: true, if the lattice matrix is diagonal.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] chunk_offset: starting offsets for each k.
 * @param[in] valid_count: number of valid gamma entries for each k.
 * @param[in] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[in] exponents: precomputed exponents (2γ-α), stride dim per entry.
 * @return helper function for the second sum in crandalls formula. Calculates
 * sum_{k in m_invt whole_numbers ** dim without zero} h₍α,kIndex₎(y + k) *
 * G_{dim - nu + 2 * |α| - 4 * kIndex}(k + y) * exp(-2 * PI * I * x * (k + y))
 */
static double complex sum_fourier_harmonic(
    double nu, unsigned int kIndex, unsigned int dim, const double *m_invt,
    const double *x, const double *y, const int cutoffs[], double zArgBound,
    bool diag, unsigned int alphaAbs, const unsigned long long *chunk_offset,
    const unsigned long long *valid_count, const double *coeffs,
    const unsigned int *exponents) {
    double lambda = 1.; // parameter that decides the weight of each sum

    int zv[dim];    // counting vector in Z^dim
    double lv[dim]; // lattice vector
    // cuboid cutoffs
    long totalSummands = 1;
    for (int k = 0; k < dim; k++) {
        zv[k] = -cutoffs[k]; // lattice vector initialized
        totalSummands *= 2 * cutoffs[k] + 1;
    }
    double complex sum = 0.0;
    double complex epsilon = 0.0;

    long zeroIndex = (totalSummands - 1) / 2;

    // second sum (in fourier space)
    for (long n = 0; n < zeroIndex; n++) {
        matrix_intVector(dim, m_invt, zv, lv, diag);
        double complex summand = summand_fourier_harmonic(
            nu, kIndex, dim, lambda, lv, x, y, zArgBound, alphaAbs, chunk_offset,
            valid_count, coeffs, exponents);
        kahan_add_c(&sum, &epsilon, summand);
        lattice_vector_increment(dim, cutoffs, zv);
    }
    lattice_vector_increment(dim, cutoffs, zv); // skips zero
    for (long n = zeroIndex + 1; n < totalSummands; n++) {
        matrix_intVector(dim, m_invt, zv, lv, diag);
        double complex summand = summand_fourier_harmonic(
            nu, kIndex, dim, lambda, lv, x, y, zArgBound, alphaAbs, chunk_offset,
            valid_count, coeffs, exponents);
        kahan_add_c(&sum, &epsilon, summand);
        lattice_vector_increment(dim, cutoffs, zv);
    }

    return sum;
}

/**
 * @brief allocates and fills the scratch buffers for the harmonic polynomial
 * coefficients.
 * @param[in] alphaAbs: total of alpha.
 * @param[in] kMax: largest degree index, alphaAbs / 2.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] alpha: multi-index alpha, length dim.
 * @param[out] chunk_offset: starting offsets for each k.
 * @param[out] valid_count: number of valid gamma entries for each k.
 * @param[out] coeffs: precomputed inner harmonic sums h_inner(α,γ,k).
 * @param[out] exponents: precomputed exponents (2γ-α).
 * @return 0 on success, 1 on allocation failure.
 */
static int harmonic_coeffs_alloc(unsigned int alphaAbs, unsigned int kMax,
                                 unsigned int dim, const unsigned int *alpha,
                                 unsigned long long **chunk_offset,
                                 unsigned long long **valid_count, double **coeffs,
                                 unsigned int **exponents) {
    *chunk_offset = malloc((kMax + 1) * sizeof **chunk_offset);
    *valid_count = malloc((kMax + 1) * sizeof **valid_count);
    *coeffs = NULL;
    *exponents = NULL;

    if (!*chunk_offset || !*valid_count) {
        free(*chunk_offset);
        free(*valid_count);
        *chunk_offset = NULL;
        *valid_count = NULL;
        return 1;
    }

    unsigned long long coeffs_size = precompute_harmonic_h_inner_chunk_size(
        alphaAbs, kMax, dim, alpha, *chunk_offset, *valid_count);

    // overflow unreachable for |alpha| < 200 in 2D and |alpha| < 80 in 3D
    // rather, the hpdyad arithmetic is the bottleneck
    *coeffs = malloc(coeffs_size * sizeof **coeffs);
    *exponents = malloc(coeffs_size * dim * sizeof **exponents);

    if (!*coeffs || !*exponents) {
        free(*chunk_offset);
        free(*valid_count);
        free(*coeffs);
        free(*exponents);
        *chunk_offset = NULL;
        *valid_count = NULL;
        *coeffs = NULL;
        *exponents = NULL;
        return 1;
    }

    precompute_harmonic_h_inner_sum(alphaAbs, dim, alpha, *chunk_offset, *coeffs,
                                    *exponents);
    return 0;
}

/**
 * @brief calculates the regularized anisotropic Epstein zeta function by
 * the harmonic method. Sibling of summation_harmoonic, the regularized version
 * carries the oscillating factor rot = exp(2πi x·y), evaluates the harmonic
 * polynomial in the zero summand at y_t1 rather than y_t2, calls
 * crandall_gReg_harmonic in place of crandall_g, and corrects the zero summand of
 * the Fourier sum when y_t1 != y_t2.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] alphaAbs: total |α| of the multi-index.
 * @param[in] alpha: multi-index α (length dim).
 * @param[in] lambda: splitting parameter.
 * @param[in] ms: lattice scaling factor, pow(vol, -1./dim).
 * @param[in] m_real: matrix that transforms the lattice (real sum).
 * @param[in] m_fourier: matrix that transforms the lattice (Fourier sum).
 * @param[in] x_t1: rescaled shift vector.
 * @param[in] x_t2: projection of shift vector to the elementary lattice cell.
 * @param[in] y_t1: rescaled wave vector.
 * @param[in] y_t2: projection of wave vector to reciprocal elementary lattice cell.
 * @param[in] cutoffsReal: how many summands in each direction (real sum).
 * @param[in] cutoffsFourier: how many summands in each direction (Fourier sum).
 * @param[in] diag: true, if the lattice matrix is diagonal.
 * the incomplete gamma evaluation.
 * @param[in] xfactor: precomputed prefactor for the real sum.
 * @return regularized anisotropic Epstein zeta function value without correction
 * terms.
 */
static double complex summation_harmonic_reg(
    double nu, unsigned int dim, unsigned int alphaAbs, const unsigned int *alpha,
    double lambda, double ms, const double *m_real, const double *m_fourier,
    const double *x_t1, const double *x_t2, const double *y_t1, const double *y_t2,
    const int cutoffsReal[], const int cutoffsFourier[], bool diag,
    double complex xfactor) {

    // Precompute coefficients for harmonic polynomials once
    unsigned int kMax = alphaAbs / 2;
    unsigned long long *chunk_offset;
    unsigned long long *valid_count;
    double *coeffs;
    unsigned int *exponents;
    if (harmonic_coeffs_alloc(alphaAbs, kMax, dim, alpha, &chunk_offset,
                              &valid_count, &coeffs, &exponents)) {
        return NAN;
    }

    // Decide when to isolate singularies in the real sum for large exponents
    bool largeExp = nu > 10 && alphaAbs > 1;

    double complex res = 0.;
    double complex rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));

    // Outer loop of the harmonic method
    for (unsigned int k = 0; k <= kMax; k++) {

        double nuIt = nu - (2 * k);
        double zArgBoundIt = assignzArgBound(nuIt);

        // skip iterartions where nuIt is a negative even integer, as
        // 1/gamma(nIt) = 0
        if (nuIt < -1 && fabs((nuIt / 2.) - nearbyint(nuIt / 2.)) < EPS) {
            continue;
        }

        double complex resIt;

        // if nu = 0, everything except the zero summand in the real sum
        // vanishes
        if (fabs(nuIt / 2) < EPS) {
            if (dot(dim, x_t2, x_t2) > EPS_ZERO_Y) {
                resIt = 0.;
            } else {
                resIt = -harmonic_h(k, dim, x_t2, alphaAbs, chunk_offset,
                                    valid_count, coeffs, exponents);
            }
        } else {

            double complex s1;
            double complex s2;
            double complex nc = 0.;

            double nuReci = nuIt - (2 * alphaAbs) + (4 * k);
            double zArgBoundReci = assignzArgBound(dim - nuReci);

            // skip zero summand if harmonic polynomial vanishes
            double h = harmonic_h(k, dim, y_t1, alphaAbs, chunk_offset, valid_count,
                                  coeffs, exponents);

            // guards 0 * inf where correct value is h(y) g(y) -> 0 for y -> 0
            // note that h(0) = 0 exactly without cancellation
            if (h) {
                nc += h * crandall_gReg_harmonic((int)k, (int)alphaAbs, dim,
                                                 dim - nu, y_t1, lambda);
            }

            s2 = sum_fourier_harmonic(nuReci, k, dim, m_fourier, x_t1, y_t2,
                                      cutoffsFourier, zArgBoundReci, diag, alphaAbs,
                                      chunk_offset, valid_count, coeffs, exponents);

            // correct wrong zero summand in regularized fourier sum.
            if (!equals(dim, y_t1, y_t2)) {
                s2 += harmonic_h(k, dim, y_t2, alphaAbs, chunk_offset, valid_count,
                                 coeffs, exponents) *
                      crandall_g(dim, dim - nuReci, y_t2, lambda, zArgBoundReci) *
                      cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2));
                s2 -= harmonic_h(k, dim, y_t1, alphaAbs, chunk_offset, valid_count,
                                 coeffs, exponents) *
                      crandall_g(dim, dim - nuReci, y_t1, lambda, zArgBoundReci) *
                      cexp(-2 * M_PI * I * dot(dim, x_t1, y_t1));
            }

            s2 = negative_one_pow(k) * inverse_imaginary_int_pow(alphaAbs) *
                 (rot * s2 + nc);

            if (largeExp) {
                s1 = sum_real_harmonic_large_exp(nuIt, k, dim, m_real, x_t2, y_t2,
                                                 cutoffsReal, zArgBoundIt, diag,
                                                 alphaAbs, chunk_offset, valid_count,
                                                 coeffs, exponents) *
                     rot * xfactor;
            } else {
                s1 = sum_real_harmonic(nuIt, k, dim, m_real, x_t2, y_t2, cutoffsReal,
                                       zArgBoundIt, diag, alphaAbs, chunk_offset,
                                       valid_count, coeffs, exponents) *
                     rot * xfactor;
            }

            resIt = (s1 + pow(lambda, dim) * s2) / tgamma(nuIt / 2.);
        }
        resIt *= pow(lambda * lambda / M_PI, -nuIt / 2.);
        res += resIt;
    }

    if (largeExp) {
        double complex s1_singularity =
            xfactor * rot *
            sum_real_harmonic_large_exp_singularity_sum(
                nu, kMax, dim, m_real, x_t2, y_t2, diag, alphaAbs, chunk_offset,
                valid_count, coeffs, exponents);
        res += s1_singularity;
    }

    res /= real_int_pow(ms, alphaAbs);

    free(chunk_offset);
    free(valid_count);
    free(coeffs);
    free(exponents);

    return res;
}

/**
 * @brief calculates the anisotropic Epstein zeta function by
 * the harmonic method.
 * @param[in] nu: exponent for the Epstein zeta function.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] alphaAbs: total |α| of the multi-index.
 * @param[in] alpha: multi-index α (length dim).
 * @param[in] lambda: splitting parameter.
 * @param[in] ms: lattice scaling factor, pow(vol, -1./dim).
 * @param[in] m_real: matrix that transforms the lattice (real sum).
 * @param[in] m_fourier: matrix that transforms the lattice (Fourier sum).
 * @param[in] x_t1: rescaled shift vector.
 * @param[in] x_t2: projection of shift vector to the elementary lattice cell.
 * @param[in] y_t1: rescaled wave vector.
 * @param[in] y_t2: projection of wave vector to reciprocal elementary lattice cell.
 * @param[in] cutoffsReal: how many summands in each direction (real sum).
 * @param[in] cutoffsFourier: how many summands in each direction (Fourier sum).
 * @param[in] diag: true, if the lattice matrix is diagonal.
 * the incomplete gamma evaluation.
 * @param[in] xfactor: precomputed prefactor for the real sum.
 * @return updated res after adding the general harmonic contribution.
 */
static double complex
summation_harmonic(double nu, unsigned int dim, unsigned int alphaAbs,
                   const unsigned int *alpha, double lambda, double ms,
                   const double *m_real, const double *m_fourier, const double *x_t1,
                   const double *x_t2, const double *y_t2, const int cutoffsReal[],
                   const int cutoffsFourier[], bool diag, double complex xfactor) {

    unsigned int kMax = alphaAbs / 2;
    unsigned long long *chunk_offset;
    unsigned long long *valid_count;
    double *coeffs;
    unsigned int *exponents;
    if (harmonic_coeffs_alloc(alphaAbs, kMax, dim, alpha, &chunk_offset,
                              &valid_count, &coeffs, &exponents)) {
        return NAN;
    }

    // Decide when to isolate singularies in the real sum for large exponents
    bool largeExp = nu > 10 && alphaAbs > 1;

    double complex res = 0.;

    // Outer loop of the harmonic method
    for (unsigned int k = 0; k <= kMax; k++) {

        double nuIt = nu - (2 * k);
        double zArgBoundIt = assignzArgBound(nuIt);

        // skip iterartions where nuIt is a negative even integer, as
        // 1/gamma(nIt) = 0
        if (nuIt < -1 && fabs((nuIt / 2.) - nearbyint(nuIt / 2.)) < EPS) {
            continue;
        }

        double complex resIt;

        // if nu = 0, everything except the zero summand in the real sum
        // vanishes
        if (fabs(nuIt / 2) < EPS) {
            if (dot(dim, x_t2, x_t2) > EPS_ZERO_Y) {
                resIt = 0.;
            } else {
                resIt = -xfactor * harmonic_h(k, dim, x_t2, alphaAbs, chunk_offset,
                                              valid_count, coeffs, exponents);
            }
        } else {

            double complex s1;
            double complex s2;
            double complex nc = 0.;

            double nuReci = nuIt - (2 * alphaAbs) + (4 * k);
            double zArgBoundReci = assignzArgBound(dim - nuReci);

            // skip zero summand if harmonic polynomial vanishes
            double h = harmonic_h(k, dim, y_t2, alphaAbs, chunk_offset, valid_count,
                                  coeffs, exponents);

            // guards 0 * inf where correct value is h(y) g(y) -> 0 for y -> 0
            // note that h(0) = 0 exactly without cancellation
            if (h) {
                nc += h *
                      crandall_g(dim, dim - nuReci, y_t2, lambda, zArgBoundReci) *
                      cexp(-2 * M_PI * I * dot(dim, y_t2, x_t1));
            }

            s2 = sum_fourier_harmonic(nuReci, k, dim, m_fourier, x_t1, y_t2,
                                      cutoffsFourier, zArgBoundReci, diag, alphaAbs,
                                      chunk_offset, valid_count, coeffs, exponents);
            s2 = negative_one_pow(k) * inverse_imaginary_int_pow(alphaAbs) *
                 (s2 + nc);

            if (largeExp) {
                s1 = sum_real_harmonic_large_exp(nuIt, k, dim, m_real, x_t2, y_t2,
                                                 cutoffsReal, zArgBoundIt, diag,
                                                 alphaAbs, chunk_offset, valid_count,
                                                 coeffs, exponents) *
                     xfactor;
            } else {
                s1 = sum_real_harmonic(nuIt, k, dim, m_real, x_t2, y_t2, cutoffsReal,
                                       zArgBoundIt, diag, alphaAbs, chunk_offset,
                                       valid_count, coeffs, exponents) *
                     xfactor;
            }

            resIt = (s1 + pow(lambda, dim) * s2) / tgamma(nuIt / 2.);
        }
        resIt *= pow(lambda * lambda / M_PI, -nuIt / 2.);
        res += resIt;
    }

    if (largeExp) {
        double complex s1_singularity =
            xfactor * sum_real_harmonic_large_exp_singularity_sum(
                          nu, kMax, dim, m_real, x_t2, y_t2, diag, alphaAbs,
                          chunk_offset, valid_count, coeffs, exponents);
        res += s1_singularity;
    }

    res /= real_int_pow(ms, alphaAbs);

    free(chunk_offset);
    free(valid_count);
    free(coeffs);
    free(exponents);

    return res;
}

/**
 * @brief calculates the (regularized) Epstein zeta function as well as the
 * (regularized) anisotropic Epstein zeta function with a prefactor.
 * @param[in] nu: exponent.
 * @param[in] dim: dimension of the input vectors.
 * @param[in] m: matrix that transforms the lattice.
 * @param[in] x: shift vector.
 * @param[in] y: wave vector.
 * @param[in] lambda: relative weight of the sums in Crandall's formula.
 * @param[in] reg: false for no regularization, true for the regularization.
 * @param[in] aniso: false for no anisotropy, true for the anisotropic variant.
 * @param[in] alpha: multiindex for the derivatives.
 * @return function value.
 */
// main logic for all four variants, splitting would duplicate shared setup
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
double complex epsteinZetaInternal(double nu, unsigned int dim, const double *m,
                                   const double *x, const double *y, double lambda,
                                   bool reg, bool aniso, const unsigned int *alpha) {
    // 1. Transform: Compute determinant and fourier transformed matrix, scale
    // both of them
    double m_fourier[dim * dim];
    double m_copy[dim * dim];
    double m_real[dim * dim];
    double x_t1[dim];
    double y_t1[dim];
    int p[dim];
    bool diag = 1;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            m_copy[(dim * i) + j] = m[(dim * i) + j];
            m_real[(dim * i) + j] = m[(dim * i) + j];
            diag = diag && ((i == j) || (m[(dim * i) + j] == 0));
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
    if (diag) {
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
    double complex res = NAN;
    double x_t2_squared = dot(dim, x_t2, x_t2);
    double y_t2_squared = dot(dim, y_t2, y_t2);
    unsigned int alphaAbs = aniso ? mult_abs(dim, alpha) : 0;
    if (!aniso && nu < 1 && fabs((nu / 2.) - nearbyint(nu / 2.)) < EPS) {
        if (x_t2_squared < EPS_ZERO_Y && nu == 0) {
            if (reg) {
                res = -1; // reg already carries phase by definition
            } else {
                res = -cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2));
            }
        } else {
            res = 0;
        }
    } else if (!reg && fabs(nu - dim - alphaAbs) < EPS &&
               y_t2_squared < EPS_ZERO_Y) {
        res = NAN;
    } else {
        double zArgBound = assignzArgBound(nu);
        double zArgBoundReci = assignzArgBound(dim - nu);
        double complex s1 = NAN;
        double complex s2 = NAN;
        double complex nc;
        double complex rot = 1;
        double vx[dim];
        for (int i = 0; i < dim; i++) {
            vx[i] = x_t1[i] - x_t2[i];
        }
        double complex xfactor = cexp(-2 * M_PI * I * dot(dim, vx, y_t1));
        if (!reg && !aniso) {
            // calculate non regularized Epstein zeta function values.
            nc = crandall_g(dim, dim - nu, y_t2, lambda, zArgBoundReci) *
                 cexp(-2 * M_PI * I * dot(dim, x_t2, y_t2));
            s1 = sum_real(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                          zArgBound, diag);
            s2 = sum_fourier(nu, dim, lambda, m_fourier, x_t2, y_t2, cutoffsFourier,
                             zArgBoundReci, diag) +
                 nc;
        } else if (reg && !aniso) {
            // calculate regularized Epstein zeta function values.
            nc = crandall_gReg(dim, dim - nu, y_t1, lambda);
            rot = cexp(2 * M_PI * I * dot(dim, x_t1, y_t1));
            s2 = sum_fourier(nu, dim, lambda, m_fourier, x_t1, y_t2, cutoffsFourier,
                             zArgBoundReci, diag);
            // correct wrong zero summand in regularized fourier sum.
            if (!equals(dim, y_t1, y_t2)) {
                s2 += crandall_g(dim, dim - nu, y_t2, lambda, zArgBoundReci) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t2)) -
                      crandall_g(dim, dim - nu, y_t1, lambda, zArgBoundReci) *
                          cexp(-2 * M_PI * I * dot(dim, x_t1, y_t1));
            }
            s2 = s2 * rot + nc;
            s1 = sum_real(nu, dim, lambda, m_real, x_t2, y_t2, cutoffsReal,
                          zArgBound, diag) *
                 rot * xfactor;
            xfactor = 1;
        } else if (!reg && aniso) {
            // Catch vanishing function values where alpha_i odd an x_i = y_i = 0
            bool oddDersInZero = false;
            for (int i = 0; i < dim; i++) {
                if (alpha[i] % 2 && fabs(x[i]) < EPS_ZERO_Y &&
                    fabs(y[i]) < EPS_ZERO_Y) {
                    oddDersInZero = true;
                }
            }
            if (oddDersInZero) {
                res = 0.;
            } else {
                res = summation_harmonic(nu, dim, alphaAbs, alpha, lambda, ms,
                                         m_real, m_fourier, x_t1, x_t2, y_t2,
                                         cutoffsReal, cutoffsFourier, diag, xfactor);
            }
        } else if (reg && aniso) {
            res = summation_harmonic_reg(nu, dim, alphaAbs, alpha, lambda, ms,
                                         m_real, m_fourier, x_t1, x_t2, y_t1, y_t2,
                                         cutoffsReal, cutoffsFourier, diag, xfactor);
        }

        // In the harmonic method, the res is already set as there is no global
        // nu-dependent coefficient there
        if (!aniso) {
            res = xfactor * pow(lambda * lambda / M_PI, -nu / 2.) / tgamma(nu / 2.) *
                  (s1 + pow(lambda, dim) * s2);
        }
    }
    free(x_t2);
    free(y_t2);
    res *= pow(ms, nu);
    // apply correction to matrix scaling if nu = d + 2k
    unsigned int k = (unsigned int)fmax(0., nearbyint((nu - (double)dim) / 2));
    if (reg && (nu == (dim + 2 * (double)k))) {
        if (alphaAbs) {
            res -= pow(M_PI, (double)k + ((double)dim / 2)) /
                   tgamma((double)k + ((double)dim / 2)) * negative_one_pow(k + 1) *
                   polynomial_y_der(k, dim, y, alpha, alphaAbs, k) * log(ms * ms) /
                   vol * inverse_imaginary_int_pow(alphaAbs) /
                   real_int_pow(-2 * M_PI, alphaAbs);
        } else {
            double ySquared = 0;
            for (int i = 0; i < dim; i++) {
                ySquared += y[i] * y[i];
            }
            res -= pow(M_PI, (double)(2 * k) + ((double)dim / 2)) /
                   tgamma((double)k + ((double)dim / 2)) * negative_one_pow(k + 1) /
                   tgamma((double)k + 1) * real_int_pow(ySquared, k) * log(ms * ms) /
                   vol;
        }
    }
    return res;
}
#undef G_BOUND
