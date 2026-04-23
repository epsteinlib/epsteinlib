// SPDX-FileCopyrightText: 2025-2026 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

/**
 * @file hpdyad.h
 * @brief High-precision dyadic number operations.
 */

#ifndef EPSTEIN_HPDYAD
#define EPSTEIN_HPDYAD

/** @brief Maximum number of 32-bit limbs in fixed-size representation (1024 bits
 * total).
 * */
#define HPDYAD_MAX_LIMBS 32

/** @brief Fixed-size high-precision dyadic number (hpdyad) with separate sign and
 * binary exponent.
 *
 * Representation: value = sign × (mantissa in base 2^32) × 2^exp2
 * Supports both integers (exp2 ≥ 0) and dyadic rationals (exp2 < 0).
 * Limbs are little-endian: limb[0] is least significant.
 * Maximum precision: HPDYAD_MAX_LIMBS × 32 bits.
 * Invariant: if n > 0, then limb[n-1] != 0 (no leading zeros).
 * Invariant: if n > 0, then limb[0] != 0 (no trailing zeros after normalization).
 */
typedef struct {
    signed char sign; // +1 or -1
    unsigned char n;  // number of used limbs (0 => zero)
    int exp2;         // global power of two (negative for dyadic rationals)
    unsigned int limb[HPDYAD_MAX_LIMBS];
} hpdyad_t;
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

/** @brief Initialize hpdyad from unsigned long long with sign.
 * @param[out] a: destination hpdyad
 * @param[in] x: unsigned integer value
 * @param[in] sign: +1 or -1
 */
void hpdyad_set_ull(hpdyad_t *a, unsigned long long x, signed char sign);

/** @brief Normalize hpdyad: trim leading zeros and left-shift mantissa so MSB of top
 * limb is 1
 * @param[in,out] a: pointer to hpdyad to normalize
 * @return void
 */
void hpdyad_normalize(hpdyad_t *a);

/** @brief Multiply two high precision dyadic numbers (hpdyad's): out = a * b
 * @param[out] out: result (may alias a or b)
 * @param[in] a: first operand
 * @param[in] b: second operand
 */
void hpdyad_mul(hpdyad_t *out, const hpdyad_t *a, const hpdyad_t *b);

/** @brief Right-shift mantissa with sticky bit (value-preserving)
 *
 * Shifts mantissa right by `bits`, increments exp2 by `bits` to preserve value.
 * Any 1-bit shifted out ORs into bit 0 (sticky). Does NOT normalize.
 *
 * @param[out] dst: destination hpdyad
 * @param[in] src: source hpdyad
 * @param[in] bits: number of bits to shift right
 */
void hpdyad_shr_bits_sticky(hpdyad_t *dst, const hpdyad_t *src, unsigned int bits);

/** @brief Left-shift mantissa by specified number of bits (value-preserving)
 * @param[out] dst: destination hpdyad
 * @param[in] src: source hpdyad
 * @param[in] bits: number of bits to shift left
 * @return void
 */
void hpdyad_shl_bits(hpdyad_t *dst, const hpdyad_t *src, int bits);

/** @brief Add two high precision dyadic numbers (hpdyad's): out = a + b
 * @param[out] out: result (supports aliasing)
 * @param[in] a: first operand
 * @param[in] b: second operand
 */
void hpdyad_add(hpdyad_t *out, const hpdyad_t *a, const hpdyad_t *b);

/** @brief Convert hpdyad to double with round-to-nearest-even.
 * @param[in] a: pointer to normalized hpdyad
 * @return closest double representation; ±DBL_MAX on overflow
 */
double hpdyad_to_double(const hpdyad_t *a);

#endif
