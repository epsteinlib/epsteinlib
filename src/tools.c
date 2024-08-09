// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan.busse@dlr.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file tools.c
 * @brief  Minimal linear algebra for matrix vector operations.
 */

#include <math.h>
#include <stdbool.h>

/*!
 * @brief minimal distance of two vector elements considered unequal.
 */
#define EPS ldexp(1, -32)

/**
 * @brief euclidean dot product.
 * @param[in] dim: dimension of the input vectors
 * @param[in] v1: first vector.
 * @param[in] v2: second vector.
 * @return dot product of v1 and v2.
 */
double dot(unsigned int dim, const double *v1, const double *v2) {
    double r = 0;
    for (int i = 0; i < dim; i++) {
        r += v1[i] * v2[i];
    }
    return r;
}

/**
 * @brief matrix - (integer) vector multiplication.
 * @param[in] dim: dimension of the square matrix and the integer vector.
 * @param[in] m: square matrix.
 * @param[in] v: integer vector.
 * @param[in,out] res: solution vector of the vector matrix multiplication.
 */
void matrix_intVector(unsigned int dim, const double *m, const int *v, double *res) {
    for (int i = 0; i < dim; i++) {
        res[i] = 0;
        for (int j = 0; j < dim; j++) {
            res[i] += m[i * dim + j] * v[j];
        }
    }
}

/**
 * @brief square matrix transpose.
 * @param[in] dim: dimension of the square matrix.
 * @param[in,out] m: square matrix.
 */
void transpose(unsigned int dim, double *m) {
    double swap;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < i; j++) {
            swap = m[dim * i + j];
            m[dim * i + j] = m[dim * j + i];
            m[dim * j + i] = swap;
        }
    }
}

/**
 * @brief check if two vectors are equal.
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
 * @brief check if vector is zero.
 * @param[in] dim: dimension of the vectors.
 * @param[in] v: vector.
 * @return true if the vector is zero.
 */
bool equalsZero(unsigned int dim, const double *v) {
    bool eq = true;
    for (int i = 0; i < dim && eq; i++) {
        eq = eq && fabs(v[i]) < EPS;
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
void invert(unsigned int dim, double *m, int *p, double *r) { // NOLINT
    // initialize p
    for (int i = 0; i < dim; i++) {
        p[i] = i;
    }
    for (int i = 0; i < dim; i++) {
        // column pivot search
        int r = i;
        for (int j = i + 1; j < dim; j++) {
            if (fabs(m[i * dim + j]) > fabs(m[i * dim + r])) {
                r = j;
            }
        }
        if (i != r) {
            int zw = p[i];
            p[i] = p[r];
            p[r] = zw;
        }
        // permute accordingly
        for (int k = 0; k < dim; k++) {
            double zw = m[i * dim + k];
            m[i * dim + k] = m[r * dim + k];
            m[r * dim + k] = zw;
        }
        // standard LU-decomposition
        for (int k = i + 1; k < dim; k++) {
            m[k * dim + i] /= m[i * dim + i]; // l-value
            for (int j = i + 1; j < dim; j++) {
                m[k * dim + j] -= m[k * dim + i] * m[i * dim + j];
            }
        }
    }
    // Compute inverse
    double y[dim]; // NOLINT user has to provide dim > 0
    for (int i = 0; i < dim; i++) {
        // Solve Ly=e_p[i]
        for (int j = 0; j < p[i]; j++) {
            y[j] = 0;
        }
        y[p[i]] = 1;
        for (int k = p[i] + 1; k < dim; k++) {
            y[k] = 0;
            for (int j = p[i]; j < k; j++) {
                y[k] -= m[k * dim + j] * y[j];
            }
        }
        // Solve Rx=y
        for (int j = (int)dim - 1; j >= 0; j--) {
            r[j * dim + i] = y[j]; // NOLINT every entry of p[i] < dim
            for (int k = j + 1; k < (int)dim; k++) {
                r[j * dim + i] -= m[j * dim + k] * r[k * dim + i];
            }
            r[j * dim + i] /= m[j * dim + j];
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
            w += fabs(m[i * dim + j]);
        }
        if (w > r) {
            r = w;
        }
    }
    return r;
}
#undef EPS
