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

/** @brief Maximum number of 32-bit limbs in fixed-size representation (1024 bits
 * total).
 * */
#define APINT_MAX_LIMBS 32

/** @brief Fixed-size arbitrary-precision integer with separate sign and binary
 * exponent.
 *
 * Representation: value = sign × (mantissa in base 2^32) × 2^exp2
 * Supports both integers (exp2 ≥ 0) and dyadic rationals (exp2 < 0).
 * Limbs are little-endian: limb[0] is least significant.
 * Maximum precision: APINT_MAX_LIMBS × 32 bits.
 * Invariant: if n > 0, then limb[n-1] != 0 (no leading zeros).
 * Invariant: if n > 0, then limb[0] != 0 (no trailing zeros after normalization).
 */
typedef struct {
    signed char sign; // +1 or -1
    unsigned char n;  // number of used limbs (0 => zero)
    int exp2;         // global power of two (negative for dyadic rationals)
    unsigned int limb[APINT_MAX_LIMBS];
} apint_t;

/**
 * @brief Compute the integer power of a double by squaring.
 * Uses switch for small exponents to avoid loop overhead.
 * @param[in] base: the base value.
 * @param[in] exp: non-negative integer exponent.
 * @return base ** exp.
 */
static inline double real_int_pow(double base, unsigned int exp) {
    double b2;
    switch (exp) {
    case 0:
        return 1.0;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        b2 = base * base;
        return b2 * b2;
    case 5:
        b2 = base * base;
        return b2 * b2 * base;
    case 6:
        b2 = base * base;
        return b2 * b2 * b2;
    case 7:
        b2 = base * base;
        return b2 * b2 * b2 * base;
    case 8:
        b2 = base * base;
        b2 = b2 * b2;
        return b2 * b2;
    default:
        break;
    }
    double res = 1.0;
    while (1) {
        if (exp & 1) {
            res *= base;
        }
        exp >>= 1;
        if (!exp) {
            break;
        }
        base *= base;
    }
    return res;
}

/**
 * @brief Compute integer powers of the imaginary unit I.
 * Uses the fact that I^4 = 1, so only exp % 4 matters.
 * @param[in] exp: non-negative integer exponent.
 * @return I ** exp.
 */
static inline double complex imaginary_int_pow(unsigned int exp) {
    static const double complex powers[4] = {1.0, I, -1.0, -I};
    return powers[exp & 3];
}

/**
 * @brief Compute integer powers of -1.
 * Uses the fact that (-1)^exp alternates between 1 and -1.
 * @param[in] exp: non-negative integer exponent.
 * @return (-1) ** exp.
 */
static inline double negative_one_pow(unsigned int exp) {
    if (exp & 1) {
        return -1.;
    }
    return 1.;
}

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

/** @brief Initialize apint from unsigned long long with sign.
 * @param[out] a: destination apint
 * @param[in] x: unsigned integer value
 * @param[in] sign: +1 or -1
 */
void apint_set_ull(apint_t *a, unsigned long long x, signed char sign);

/** @brief Normalize apint: trim leading zeros and left-shift mantissa so MSB of top
 * limb is 1
 * @param[in,out] a: pointer to apint to normalize
 * @return void
 */
void apint_normalize(apint_t *a);

/** @brief Multiply two apints: out = a * b
 * @param[out] out: result (may alias a or b)
 * @param[in] a: first operand
 * @param[in] b: second operand
 */
void apint_mul(apint_t *out, const apint_t *a, const apint_t *b);

/** @brief Right-shift mantissa with sticky bit (value-preserving)
 *
 * Shifts mantissa right by `bits`, increments exp2 by `bits` to preserve value.
 * Any 1-bit shifted out ORs into bit 0 (sticky). Does NOT normalize.
 *
 * @param[out] dst: destination apint
 * @param[in] src: source apint
 * @param[in] bits: number of bits to shift right
 */
void apint_shr_bits_sticky(apint_t *dst, const apint_t *src, unsigned int bits);

/** @brief Left-shift mantissa by specified number of bits (value-preserving)
 * @param[out] dst: destination apint
 * @param[in] src: source apint
 * @param[in] bits: number of bits to shift left
 * @return void
 */
void apint_shl_bits(apint_t *dst, const apint_t *src, int bits);

/** @brief Add two apints: out = a + b
 * @param[out] out: result (supports aliasing)
 * @param[in] a: first operand
 * @param[in] b: second operand
 */
void apint_add(apint_t *out, const apint_t *a, const apint_t *b);

/** @brief Convert apint to double with round-to-nearest-even.
 * @param[in] a: pointer to normalized apint
 * @return closest double representation; ±DBL_MAX on overflow
 */
double apint_to_double(const apint_t *a);

double dot(unsigned int dim, const double *v1, const double *v2);
void matrix_intVector(unsigned int dim, const double *m, const int *v, double *res);
void transpose(unsigned int dim, double *m);
bool equals(unsigned int dim, const double *v1, const double *v2);
bool equalsZero(unsigned int dim, const double *v);
void invert(unsigned int dim, double *m, int *p, double *r);
double inf_norm(unsigned int dim, const double *m);
unsigned int mult_abs(unsigned int dim, const unsigned int *alpha);
unsigned long long binom(unsigned long long n, unsigned long long k);
#endif
