// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
// SPDX-FileCopyrightText: 2024 Ruben Gutendorf
// <ruben.gutendorf@uni-saarland.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file tools.h
 * @brief  Minimal linear algebra for matrix vector and multi-index operations.
 */

#ifndef EPSTEIN_TOOLS
#define EPSTEIN_TOOLS
#include <complex.h>
#include <stdbool.h>

/** @brief Maximum number of 32-bit limbs in fixed-size representation (256 bits
 * total). */
#define APINT_MAX_LIMBS 8

/** @brief Fixed-size arbitrary-precision integer with separate sign and binary
 * exponent.
 *
 * Representation: value = sign × (mantissa in base 2^32) × 2^exp2
 * Limbs are little-endian: limb[0] is least significant.
 * Invariant: if n > 0, then limb[n-1] != 0 (no leading zeros).
 */
typedef struct {
    signed char sign; // +1 or -1
    signed char n;    // number of used limbs (0 => zero)
    int exp2;         // global power of two
    unsigned int limb[APINT_MAX_LIMBS];
} apint_t;

/** @brief Count trailing zero bits in a 64-bit unsigned integer.
 * @param[in] x: input value
 * @return number of trailing zeros (0..63); returns 64 if x == 0
 */
unsigned ctz64(unsigned long long x);

/** @brief Find index of most significant bit in a 32-bit unsigned integer.
 * @param[in] x: input value
 * @return bit index 0..31; returns 0 for x == 0 (defensive fallback)
 */
int msb32(unsigned int x);
double dot(unsigned int dim, const double *v1, const double *v2);
void matrix_intVector(unsigned int dim, const double *m, const int *v, double *res);
void transpose(unsigned int dim, double *m);
bool equals(unsigned int dim, const double *v1, const double *v2);
bool equalsZero(unsigned int dim, const double *v);
void invert(unsigned int dim, double *m, int *p, double *r);
double inf_norm(unsigned int dim, const double *m);
unsigned int mult_abs(unsigned int dim, const unsigned int *alpha);
unsigned long long binom(unsigned long long n, unsigned long long k);
double complex int_pow(double complex base, unsigned int exp);
#endif
