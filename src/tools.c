// SPDX-FileCopyrightText: 2024 Andreas Buchheit <buchheit@num.uni-sb.de>
// SPDX-FileCopyrightText: 2024 Jonathan Busse <jonathan@jbusse.de>
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
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>

/*!
 * @brief minimal distance of two vector elements considered unequal.
 */
#define EPS ldexp(1, -32)

/** @brief Count trailing zero bits in a 64-bit unsigned integer.
 * @param[in] x: input value
 * @return number of trailing zeros (0..63); returns 64 if x == 0
 */
unsigned ctz64(unsigned long long x) {
    if (x == 0) {
        return 64;
    }
    unsigned n = 0;
    while ((x & 1) == 0) {
        x >>= 1;
        n++;
    }
    return n;
}

/** @brief Find index of most significant bit in a 32-bit unsigned integer.
 * @param[in] x: input value
 * @return bit index 0..31; returns 0 for x == 0 (defensive fallback)
 */
int msb32(unsigned int x) {
    if (x == 0) {
        return 0;
    }
    int n = 0;
    if (x >= 0x10000) {
        n += 16;
        x >>= 16;
    }
    if (x >= 0x100) {
        n += 8;
        x >>= 8;
    }
    if (x >= 0x10) {
        n += 4;
        x >>= 4;
    }
    if (x >= 0x4) {
        n += 2;
        x >>= 2;
    }
    if (x >= 0x2) {
        n += 1;
    }
    return n;
}

/** @brief Initialize apint from unsigned long long with sign.
 * @param[out] a: destination apint
 * @param[in] x: unsigned integer value
 * @param[in] sign: +1 or -1
 */
void hpdyad_set_ull(hpdyad_t *a, unsigned long long x, signed char sign) {
    if (x == 0) {
        a->n = 0;
        a->exp2 = 0;
        a->sign = 1; // canonical zero
        return;
    }

    // Extract and remove trailing zeros
    unsigned tz = ctz64(x);
    x >>= tz;

    // Split into limbs (little-endian)
    a->limb[0] = (unsigned int)(x & 0xFFFFFFFFULL);
    a->limb[1] = (unsigned int)(x >> 32);

    // Set number of used limbs
    a->n = (a->limb[1] != 0) ? 2 : 1;

    // Store exponent and sign
    a->exp2 = (int)tz;
    a->sign = sign;
}

/** @brief Normalize apint: trim leading zeros and left-shift mantissa so MSB of top
 * limb is 1
 * @param[in,out] a: pointer to apint to normalize
 * @return void
 */
void hpdyad_normalize(hpdyad_t *a) {
    int i;
    int s;

    // Trim leading zero limbs
    while (a->n > 0 && a->limb[a->n - 1] == 0) {
        a->n--;
    }

    // Zero canonicalization
    if (a->n == 0) {
        a->sign = 1;
        a->exp2 = 0;
        return;
    }

    // Trim LOW-end zero limbs
    while (a->n > 0 && a->limb[0] == 0) {
        // Shift limbs right by one position
        for (i = 0; i < a->n - 1; i++) {
            a->limb[i] = a->limb[i + 1];
        }
        a->limb[a->n - 1] = 0;
        a->n--;
        a->exp2 += 32;
    }

    // Check for zero after trimming
    if (a->n == 0) {
        a->sign = 1;
        a->exp2 = 0;
        return;
    }

    // Invariant: limb[n-1] != 0 now guaranteed
    // Left-shift to normalize MSB of limb[n-1]
    s = 31 - msb32(a->limb[a->n - 1]);

    if (s == 0) {
        return; // Already normalized
    }

    // Shift entire mantissa left by s bits
    for (i = a->n - 1; i >= 1; i--) {
        a->limb[i] = (a->limb[i] << s) | (a->limb[i - 1] >> (32 - s));
    }
    a->limb[0] <<= s;

    // Adjust exponent
    a->exp2 -= s;

    // Trim LOW-end zeros created by left-shift
    while (a->n > 0 && a->limb[0] == 0) {
        // Shift limbs right by one position
        for (i = 0; i < a->n - 1; i++) {
            a->limb[i] = a->limb[i + 1];
        }
        a->limb[a->n - 1] = 0;
        a->n--;
        a->exp2 += 32;
    }
}

/** @brief Multiply two apints: out = a * b
 * @param[out] out: result (may alias a or b)
 * @param[in] a: first operand
 * @param[in] b: second operand
 */
void hpdyad_mul(hpdyad_t *out, const hpdyad_t *a, const hpdyad_t *b) {
    hpdyad_t result;
    unsigned int i;
    unsigned int j;
    unsigned int k;
    unsigned long long acc;
    unsigned long long prod;
    unsigned int carry;

    // Zero handling
    if (a->n == 0 || b->n == 0) {
        out->n = 0;
        out->sign = 1;
        out->exp2 = 0;
        for (k = 0; k < APINT_MAX_LIMBS; k++) {
            out->limb[k] = 0;
        }
        return;
    }

    // Initialize result limbs to zero
    for (k = 0; k < APINT_MAX_LIMBS; k++) {
        result.limb[k] = 0;
    }

    // Schoolbook multiplication with immediate carry propagation
    for (i = 0; i < a->n; i++) {
        carry = 0;
        for (j = 0; j < b->n; j++) {
            k = i + j;
            if (k >= APINT_MAX_LIMBS) {
                break; // Truncate high limbs
            }

            prod = (unsigned long long)a->limb[i] * b->limb[j];
            acc = (unsigned long long)result.limb[k] + prod + carry;
            result.limb[k] = (unsigned int)acc;
            carry = (unsigned int)(acc >> 32);
        }

        // Propagate final carry
        while (carry && k + 1 < APINT_MAX_LIMBS) {
            k++;
            acc = (unsigned long long)result.limb[k] + carry;
            result.limb[k] = (unsigned int)acc;
            carry = (unsigned int)(acc >> 32);
        }
    }

    // Set result size (min of natural size and max limbs)
    result.n = a->n + b->n;
    if (result.n > APINT_MAX_LIMBS) {
        result.n = APINT_MAX_LIMBS;
    }

    // Set sign and exponent
    result.sign = (a->sign == b->sign) ? 1 : -1;
    result.exp2 = a->exp2 + b->exp2;

    // Normalize
    hpdyad_normalize(&result);

    // Copy to output (handles aliasing)
    *out = result;
}

/** @brief Right-shift mantissa with sticky bit (value-preserving)
 *
 * Shifts mantissa right by `bits`, increments exp2 by `bits` to preserve value.
 * Any 1-bit shifted out ORs into bit 0 (sticky). Does NOT normalize.
 *
 * @param[out] dst: destination apint
 * @param[in] src: source apint
 * @param[in] bits: number of bits to shift right
 */
void hpdyad_shr_bits_sticky(hpdyad_t *dst, const hpdyad_t *src, // NOLINT
                            unsigned int bits) {
    hpdyad_t tmp;
    hpdyad_t *target = (dst == src) ? &tmp : dst;

    // Handle zero source
    if (src->n == 0) {
        target->sign = 1;
        target->n = 0;
        target->exp2 = 0;
        for (unsigned char i = 0; i < APINT_MAX_LIMBS; i++) {
            target->limb[i] = 0;
        }
        if (dst == src) {
            *dst = tmp;
        }
        return;
    }

    // Handle complete shift-out
    if (bits >= 32 * src->n) {
        target->sign = src->sign;
        target->exp2 = src->exp2 + (int)bits;
        target->n = 1;
        target->limb[0] = 1;
        for (unsigned char i = 1; i < APINT_MAX_LIMBS; i++) {
            target->limb[i] = 0;
        }
        if (dst == src) {
            *dst = tmp;
        }
        return;
    }

    unsigned int limb_shift = bits / 32;
    unsigned int bit_shift = bits % 32;
    unsigned char sticky = 0;

    // Collect sticky from dropped limbs
    for (unsigned int i = 0; i < limb_shift; i++) {
        if (src->limb[i] != 0) {
            sticky = 1;
            break;
        }
    }

    // Initialize target limbs to zero
    for (unsigned char i = 0; i < APINT_MAX_LIMBS; i++) {
        target->limb[i] = 0;
    }

    // Shift limbs
    if (bit_shift == 0) {
        // Simple limb shift (no bit alignment needed)
        for (unsigned int i = 0; i < src->n - limb_shift; i++) {
            target->limb[i] = src->limb[i + limb_shift];
        }
    } else {

        // Shift with bit alignment
        for (unsigned int i = 0; i < src->n - limb_shift; i++) {
            unsigned int lo = src->limb[i + limb_shift];
            unsigned int hi = ((unsigned int)(i + limb_shift + 1) < src->n)
                                  ? src->limb[i + limb_shift + 1]
                                  : 0;
            // Collect sticky ONLY from the lowest limb (i==0)
            if (i == 0 && (lo & ((1U << bit_shift) - 1)) != 0) {
                sticky = 1;
            }
            target->limb[i] = (lo >> bit_shift) | (hi << (32 - bit_shift));
        }
    }

    // OR sticky into lowest limb
    if (sticky) {
        target->limb[0] |= 1;
    }

    // Zero upper limbs
    for (unsigned char i = src->n - limb_shift; i < APINT_MAX_LIMBS; i++) {
        target->limb[i] = 0;
    }

    // Set fields (conservatively, value-preserving)
    target->n = src->n - limb_shift;
    target->exp2 = src->exp2 + (int)bits;
    target->sign = src->sign;

    // Handle aliasing
    if (dst == src) {
        *dst = tmp;
    }
}

/** @brief Compare magnitudes of two apints (assumes same exp2 after alignment)
 * @param[in] a: first apint
 * @param[in] b: second apint
 * @return 1 if |a| > |b|, -1 if |a| < |b|, 0 if equal
 */
static int compare_magnitudes(const hpdyad_t *a, const hpdyad_t *b) {
    unsigned char max_n;
    int i;

    // Conservative n handling: scan from highest possible limb
    max_n = (a->n > b->n) ? a->n : b->n;

    for (i = (int)max_n - 1; i >= 0; i--) {
        unsigned int a_limb = (i < a->n) ? a->limb[i] : 0;
        unsigned int b_limb = (i < b->n) ? b->limb[i] : 0;

        if (a_limb > b_limb) {
            return 1;
        }
        if (a_limb < b_limb) {
            return -1;
        }
    }

    return 0; // Equal magnitude
}

/** @brief Add two unsigned mantissas (assumes same exp2)
 * @param[out] result: output apint (only limb[] and n are set)
 * @param[in] a: first addend
 * @param[in] b: second addend
 */
static void add_mantissas_unsigned(hpdyad_t *result, const hpdyad_t *a,
                                   const hpdyad_t *b) {
    unsigned char max_n;
    unsigned long long sum;
    unsigned int carry;
    unsigned char i;

    max_n = (a->n > b->n) ? a->n : b->n;
    carry = 0;

    // Add limbs with carry propagation
    for (i = 0; i < APINT_MAX_LIMBS; i++) {
        // Stop when no more non-zero limbs and no carry
        if (i >= max_n && carry == 0) {
            break;
        }

        unsigned int a_limb = (i < a->n) ? a->limb[i] : 0;
        unsigned int b_limb = (i < b->n) ? b->limb[i] : 0;

        sum = (unsigned long long)a_limb + b_limb + carry;
        result->limb[i] = (unsigned int)sum;
        carry = (unsigned int)(sum >> 32);
    }

    result->n = i;

    // Zero out unused limbs (for cleanliness)
    for (; i < APINT_MAX_LIMBS; i++) {
        result->limb[i] = 0;
    }
}

/** @brief Subtract smaller unsigned mantissa from larger (assumes same exp2)
 * @param[out] result: output apint (only limb[] and n are set)
 * @param[in] larger: larger magnitude operand
 * @param[in] smaller: smaller magnitude operand (assumed |larger| >= |smaller|)
 */
static void subtract_mantissas_unsigned(hpdyad_t *result, const hpdyad_t *larger,
                                        const hpdyad_t *smaller) {
    long long diff;
    unsigned int borrow;
    unsigned char i;

    borrow = 0;

    // Subtract limbs with borrow propagation
    for (i = 0; i < larger->n; i++) {
        unsigned int larger_limb = larger->limb[i];
        unsigned int smaller_limb = (i < smaller->n) ? smaller->limb[i] : 0;

        // Compute difference using signed arithmetic to detect borrow
        diff = (long long)larger_limb - (long long)smaller_limb - (long long)borrow;

        if (diff < 0) {
            // Need to borrow from next limb
            result->limb[i] = (unsigned int)(diff + 0x100000000LL);
            borrow = 1;
        } else {
            result->limb[i] = (unsigned int)diff;
            borrow = 0;
        }
    }

    result->n = larger->n;

    // Zero out unused limbs
    for (i = larger->n; i < APINT_MAX_LIMBS; i++) {
        result->limb[i] = 0;
    }

    // Result may have leading zeros; normalization will clean up
}

/** @brief Left-shift mantissa by specified number of bits (value-preserving)
 * @param[out] dst: destination apint
 * @param[in] src: source apint
 * @param[in] bits: number of bits to shift left
 * @return void
 */
void hpdyad_shl_bits(hpdyad_t *dst, const hpdyad_t *src, int bits) {
    unsigned int limb_shift;
    unsigned int bit_shift;
    unsigned char i;

    // Handle zero input
    if (src->n == 0) {
        dst->sign = 1;
        dst->exp2 = 0;
        dst->n = 0;
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            dst->limb[i] = 0;
        }
        return;
    }
    // Copy sign
    dst->sign = src->sign;
    // Adjust exponent (value-preserving: multiply mantissa, divide by 2^bits)
    dst->exp2 = src->exp2 - bits;
    limb_shift = bits / 32;
    bit_shift = bits % 32;

    // Check if shift exceeds available space
    if (limb_shift >= APINT_MAX_LIMBS) {
        // Complete overflow - saturate to max
        dst->n = APINT_MAX_LIMBS;
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            dst->limb[i] = 0xFFFFFFFF;
        }
        return;
    }
    // Zero out destination
    for (i = 0; i < APINT_MAX_LIMBS; i++) {
        dst->limb[i] = 0;
    }
    // Shift by whole limbs first
    if (bit_shift == 0) {
        // Simple limb shift
        for (i = 0; i < src->n && (i + limb_shift) < APINT_MAX_LIMBS; i++) {
            dst->limb[i + limb_shift] = src->limb[i];
        }
        dst->n = (src->n + limb_shift < APINT_MAX_LIMBS) ? src->n + limb_shift
                                                         : APINT_MAX_LIMBS;
    } else {
        // Shift with bit offset
        for (i = 0; i < src->n && (i + limb_shift) < APINT_MAX_LIMBS; i++) {
            unsigned long long temp = (unsigned long long)src->limb[i] << bit_shift;
            dst->limb[i + limb_shift] |= (unsigned int)temp;
            if (i + limb_shift + 1 < APINT_MAX_LIMBS) {
                dst->limb[i + limb_shift + 1] = (unsigned int)(temp >> 32);
            }
        }
        dst->n = (src->n + limb_shift + 1 < APINT_MAX_LIMBS)
                     ? src->n + limb_shift + 1
                     : APINT_MAX_LIMBS;
    }
}

/** @brief Add two apints: out = a + b
 * @param[out] out: result (supports aliasing)
 * @param[in] a: first operand
 * @param[in] b: second operand
 */
void hpdyad_add(hpdyad_t *out, const hpdyad_t *a, const hpdyad_t *b) { // NOLINT
    hpdyad_t temp_result;
    hpdyad_t *result;
    const hpdyad_t *higher_exp;
    const hpdyad_t *lower_exp;
    hpdyad_t aligned_higher;
    hpdyad_t aligned_lower;
    int d;
    unsigned char i;
    int bits_higher;
    int bits_available;

    // Handle aliasing
    result = (out == a || out == b) ? &temp_result : out;

    // Handle zeros explicitly
    if (a->n == 0) {
        *result = *b;
        if (result == &temp_result) {
            *out = temp_result;
        }
        return;
    }
    if (b->n == 0) {
        *result = *a;
        if (result == &temp_result) {
            *out = temp_result;
        }
        return;
    }

    // Determine which has larger exp2
    if (a->exp2 >= b->exp2) {
        higher_exp = a;
        lower_exp = b;
        d = a->exp2 - b->exp2;
    } else {
        higher_exp = b;
        lower_exp = a;
        d = b->exp2 - a->exp2;
    }

    // Smart alignment strategy
    bits_higher =
        32 * (higher_exp->n - 1) + msb32(higher_exp->limb[higher_exp->n - 1]) + 1;
    bits_available = 32 * APINT_MAX_LIMBS;

    if (d > 0 && d + bits_higher <= bits_available) {
        // Result of shift fits in available space - exact arithmetic
        aligned_lower = *lower_exp;
        hpdyad_shl_bits(&aligned_higher, higher_exp, d);
        result->exp2 = lower_exp->exp2;

    } else {
        // Result wouldn't fit, use sticky bit approach
        aligned_higher = *higher_exp;
        hpdyad_shr_bits_sticky(&aligned_lower, lower_exp, d);
        result->exp2 = higher_exp->exp2;
    }

    // Perform addition or subtraction based on signs
    if (aligned_higher.sign == aligned_lower.sign) {
        // Same sign: add magnitudes
        add_mantissas_unsigned(result, &aligned_higher, &aligned_lower);
        result->sign = aligned_higher.sign;
    } else {
        // Different signs: subtract magnitudes
        int cmp = compare_magnitudes(&aligned_higher, &aligned_lower);

        if (cmp == 0) {
            // Equal magnitude: result is canonical zero
            result->sign = 1;
            result->n = 0;
            result->exp2 = 0;
            for (i = 0; i < APINT_MAX_LIMBS; i++) {
                result->limb[i] = 0;
            }
        } else {
            // Determine which has larger magnitude
            const hpdyad_t *larger_mag =
                (cmp > 0) ? &aligned_higher : &aligned_lower;
            const hpdyad_t *smaller_mag =
                (cmp > 0) ? &aligned_lower : &aligned_higher;

            subtract_mantissas_unsigned(result, larger_mag, smaller_mag);
            result->sign = larger_mag->sign;
        }
    }

    // Normalize result
    hpdyad_normalize(result);

    // Copy back if aliasing
    if (result == &temp_result) {
        *out = temp_result;
    }
}

/** @brief Convert apint to double with round-to-nearest-even.
 * @param[in] a: pointer to normalized apint
 * @return closest double representation; ±DBL_MAX on overflow
 */
double hpdyad_to_double(const hpdyad_t *a) {
    int i;
    int top_bit;
    int adj_exp;
    int guard;
    int sticky;
    unsigned long long top53;
    double result;

    // Zero → +0.0
    if (a->n == 0) {
        return 0.0;
    }

    // MSB position of value (assumes normalized: msb of limb[n-1] is 31)
    top_bit = a->exp2 + 32 * a->n - 1;

    // Overflow
    if (top_bit >= 1024) {
        return (a->sign > 0) ? DBL_MAX : -DBL_MAX;
    }

    if (a->n == 1) {
        // ≤32 bits, always exact
        result = ldexp((double)a->limb[0], a->exp2);
    } else {
        // Extract top 53 bits
        top53 = ((unsigned long long)a->limb[a->n - 1] << 21) |
                (a->limb[a->n - 2] >> 11);

        // Guard bit (bit 10 of limb[n-2]) and sticky bits (bits 9..0)
        guard = ((a->limb[a->n - 2] & 0x400U) != 0);
        sticky = ((a->limb[a->n - 2] & 0x3FFU) != 0);

        // OR in lower limbs for sticky
        for (i = 0; i < a->n - 2; i++) {
            if (a->limb[i] != 0) {
                sticky = 1;
                break;
            }
        }

        // Round to nearest even
        if (guard && (sticky || (top53 & 1))) {
            top53++;
        }

        adj_exp = a->exp2 + 32 * (a->n - 2) + 11;
        result = ldexp((double)top53, adj_exp);

        // Clamp overflow from rounding
        if (result > DBL_MAX) {
            result = DBL_MAX;
        }
    }

    return (a->sign > 0) ? result : -result;
}

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
            res[i] += m[(i * dim) + j] * v[j];
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
            swap = m[(dim * i) + j];
            m[(dim * i) + j] = m[(dim * j) + i];
            m[(dim * j) + i] = swap;
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
    double y[dim]; // NOLINT user has to provide dim > 0
    for (int i = 0; i < dim; i++) {
        // Solve Ly = P e_i
        for (int k = 0; k < dim; k++) {
            y[k] = (p[k] == i);
            for (int j = 0; j < k; j++) {
                y[k] -= m[(k * dim) + j] * y[j];
            }
        } // Solve Rx=y
        for (int j = (int)dim - 1; j >= 0; j--) {
            r[j * dim + i] = y[j]; // NOLINT every entry of p[i] < dim
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
    unsigned long long res = 1;

    // Calculate binom(n)(n-k) if n - k is closer to smaller than k
    if (k > n - k) {
        k = n - k;
    }

    for (unsigned int i = 1; i <= k; i++) {
        res = res * (n - k + i) / i;
    }
    return res;
}

/**
 * @brief Compute the integer power of a double by squaring.
 * @param[in] base: non-negative integer greater or equal k.
 * @param[in] exp: non-negative integer smaller or equal n.
 * @return base ** exp.
 */
double complex int_pow(double complex base, unsigned int exp) {
    double complex res = 1.;
    while (true) {
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

#undef EPS
