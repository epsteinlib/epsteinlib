// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024-2026 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file tools.c
 * @brief  Minimal linear algebra for matrix vector operations, binomials and
 * factorials.
 */

#include "tools.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

/*!
 * @brief Minimal distance of two vector elements considered unequal.
 */
#define EPS ldexp(1, -32)

/**
 * @brief Square matrix transpose.
 * @param[in] dim: dimension of the square matrix.
 * @param[in,out] m: square matrix.
 */
void transpose(unsigned int dim, double *m) {
    double swap;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < i; j++) {
            swap = m[(dim * i) + j];
            m[(dim * i) + j] = m[(dim * j) + i];
            m[(dim * j) + i] = swap;
        }
    }
}

/**
 * @brief Check if two vectors are equal.
 * @param[in] dim: dimension of the vectors.
 * @param[in] v1: first vector.
 * @param[in] v2: second vector.
 * @return true if the vectors are equal, false if the vectors are not equal.
 */
bool equals(unsigned int dim, const double *v1, const double *v2) {
    bool eq = true;
    for (int i = 0; i < dim && eq; i++) {
        eq = eq && fabs(v1[i] - v2[i]) < EPS;
    }
    return eq;
}

/**
 * @brief Invert matrix.
 * @param[in] dim: dimension of the vectors.
 * @param[in, out] m: matrix to invert. overwritten bei LU-decomposition.
 * @param[out] p: permutation vector.
 * @param[out] r: where inverse matrix is stored.
 */
// fast inversion where lots of parameters are reused
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void invert(unsigned int dim, double *m, int *p, double *r) {
    // initialize p
    for (int i = 0; i < dim; i++) {
        p[i] = i;
    }
    for (int i = 0; i < dim; i++) {
        // column pivot search
        int pivot = i;
        for (int j = i + 1; j < dim; j++) {
            if (fabs(m[(j * dim) + i]) > fabs(m[(pivot * dim) + i])) {
                pivot = j;
            }
        }
        if (i != pivot) {
            int zw = p[i];
            p[i] = p[pivot];
            p[pivot] = zw;
        }
        // permute accordingly
        for (int k = 0; k < dim; k++) {
            double zw = m[(i * dim) + k];
            m[(i * dim) + k] = m[(pivot * dim) + k];
            m[(pivot * dim) + k] = zw;
        }
        // standard LU-decomposition
        for (int k = i + 1; k < dim; k++) {
            m[(k * dim) + i] /= m[(i * dim) + i]; // l-value
            for (int j = i + 1; j < dim; j++) {
                m[(k * dim) + j] -= m[(k * dim) + i] * m[(i * dim) + j];
            }
        }
    }
    // Compute inverse
    // the analyzer does not know that dim > 0
    // NOLINTNEXTLINE(clang-analyzer-core.VLASize)
    double y[dim];
    for (int i = 0; i < dim; i++) {
        // Solve Ly = P e_i
        for (int k = 0; k < dim; k++) {
            y[k] = (p[k] == i);
            for (int j = 0; j < k; j++) {
                y[k] -= m[(k * dim) + j] * y[j];
            }
        } // Solve Rx=y
        for (int j = (int)dim - 1; j >= 0; j--) {
            r[(j * dim) + i] = y[j]; // every entry of p[i] < dim
            for (int k = j + 1; k < (int)dim; k++) {
                r[(j * dim) + i] -= m[(j * dim) + k] * r[(k * dim) + i];
            }
            r[(j * dim) + i] /= m[(j * dim) + j];
        }
    }
}

/**
 * @brief Compute infinity norm (maximum sum row norm).
 * @param[in] dim: dimension of the vectors.
 * @param[in] m: matrix to compute infinity norm of.
 */
double inf_norm(unsigned int dim, const double *m) { // NOLINT
    double r = 0;
    for (int j = 0; j < dim; j++) {
        r += fabs(m[j]);
    }
    for (int i = 1; i < dim; i++) {
        double w = 0;
        for (int j = 0; j < dim; j++) {
            w += fabs(m[(i * dim) + j]);
        }
        if (w > r) {
            r = w;
        }
    }
    return r;
}

/**
 * @brief Compute absolute value of multi-index, that is the sum of its components.
 * @param[in] dim: dimension of alpha end vec.
 * @param[in] alpha: multi-index.
 * @return absolute values of alpha.
 * @return absolute values of alpha.
 */
unsigned int mult_abs(unsigned int dim, const unsigned int *alpha) {
    unsigned int n = 0;
    for (int i = 0; i < dim; i++) {
        n += alpha[i];
    }
    return n;
}

/**
 * @brief Compute the binomial coefficient bionm(n,k).
 * @param[in] n: non-negative integer greater or equal k.
 * @param[in] k: non-negative integer smaller or equal n.
 * @return binom(n)(k).
 */
unsigned long long binom(unsigned long long n, unsigned long long k) {
    // Calculate binom(n)(n-k) if n - k is smaller than k
    if (k > n - k) {
        k = n - k;
    }
    unsigned long long res = 1;
    for (unsigned int i = 1; i <= k; i++) {
        res = res * (n - k + i) / i;
    }
    return res;
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
#undef EPS
