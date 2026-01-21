// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "../tools.h"
#include "utils.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

#ifndef BASE_PATH
#define BASE_PATH "csv"
#endif

/*!
 * @brief Benchmarks Count Trailing Zeros (64-bit) and Most Significant Bit (32-bit)
 * functions
 *
 * @return number of failed tests.
 * */
int test_ctzmsb(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/ctzmsb_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    unsigned long long x;
    unsigned int ctz_expected;
    int msb_expected;
    unsigned int ctz_result;

    int scanResult;
    char line[256];
    int testsPassed = 0;
    int totalTests = 0;

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult =
            sscanf(line, "%llu,%u,%d", &x, &ctz_expected, &msb_expected); // NOLINT
        if (scanResult != 3) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 3", scanResult);
            continue;
        }

        ctz_result = ctz64(x);
        int m_is_odd = (x == 0) ? 1 : (int)((x >> ctz_result) & 1);

        // Validate msb32
        int msb_result;
        if (x == 0) {
            msb_result = 0;
        } else {
            unsigned long long m = x >> ctz_result;
            if (m < (1ULL << 32)) {
                msb_result = msb32((unsigned int)m);
            } else {
                msb_result = 32 + msb32((unsigned int)(m >> 32));
            }
        }

        if (ctz_result == ctz_expected && m_is_odd && msb_result == msb_expected) {
            testsPassed++;
        } else {
            printf("\n\nWarning! ");
            printf(
                "x=%llu: ctz=%u (expected %u), msb_m=%d (expected %d), m_odd=%d\n",
                x, ctz_result, ctz_expected, msb_result, msb_expected, m_is_odd);
        }
        totalTests++;
    }

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
    }

    printf("\n\t ... %d out of %d tests passed.", testsPassed, totalTests);
    printf("\n");
    return totalTests - testsPassed;
}

/*!
 * @brief Tests apint_set_ull function against Mathematica reference values
 *
 * @return number of failed tests.
 * */
int test_apint_set_ull(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/apint_set_ull_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    unsigned long long x;
    int sign_input;
    unsigned int limb0_expected;
    unsigned int limb1_expected;
    int n_expected;
    int exp2_expected;
    int sign_expected;
    int scanResult;
    char line[256];
    int testsPassed = 0;
    int totalTests = 0;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult =
            sscanf(line, "%llu,%d,%u,%u,%d,%d,%d", &x, &sign_input, // NOLINT
                   &limb0_expected, &limb1_expected, &n_expected, &exp2_expected,
                   &sign_expected);
        if (scanResult != 7) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 7", scanResult);
            continue;
        }

        apint_t a;
        apint_set_ull(&a, x, (signed char)sign_input);

        int passed;
        if (x == 0) {
            // For zero: only check n, exp2, sign
            passed = (a.n == n_expected && a.exp2 == exp2_expected &&
                      a.sign == sign_expected);
        } else {
            // For non-zero: check all fields
            passed = (a.limb[0] == limb0_expected && a.limb[1] == limb1_expected &&
                      a.n == n_expected && a.exp2 == exp2_expected &&
                      a.sign == sign_expected);
        }

        if (passed) {
            testsPassed++;
        } else {
            printf("\n\nWarning! ");
            printf(
                "x=%llu, sign=%d: limb[0]=%u (expected %u), limb[1]=%u (expected "
                "%u), "
                "n=%d (expected %d), exp2=%d (expected %d), sign=%d (expected %d)\n",
                x, sign_input, a.limb[0], limb0_expected, a.limb[1], limb1_expected,
                a.n, n_expected, a.exp2, exp2_expected, a.sign, sign_expected);
        }
        totalTests++;
    }
    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
    }
    printf("\n\t ... %d out of %d tests passed.", testsPassed, totalTests);
    printf("\n");
    return totalTests - testsPassed;
}

/*!
 * @brief Tests apint_normalize function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_apint_normalize(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/apint_normalize_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    int sign_in;
    int exp2_in;
    unsigned int limbs_in[APINT_MAX_LIMBS];
    int sign_out;
    int exp2_out;
    unsigned int limbs_out[APINT_MAX_LIMBS];
    int scanResult;
    char line[512];
    int testsPassed = 0;
    int totalTests = 0;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf( // NOLINT
            line, "%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d,%u,%u,%u,%u,%u,%u,%u,%u",
            &sign_in, &exp2_in, &limbs_in[0], &limbs_in[1], &limbs_in[2],
            &limbs_in[3], &limbs_in[4], &limbs_in[5], &limbs_in[6], &limbs_in[7],
            &sign_out, &exp2_out, &limbs_out[0], &limbs_out[1], &limbs_out[2],
            &limbs_out[3], &limbs_out[4], &limbs_out[5], &limbs_out[6],
            &limbs_out[7]);
        if (scanResult != 20) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 20", scanResult);
            continue;
        }
        apint_t a;
        a.sign = (signed char)sign_in;
        a.exp2 = exp2_in;
        a.n = APINT_MAX_LIMBS;
        for (int i = 0; i < APINT_MAX_LIMBS; i++) {
            a.limb[i] = limbs_in[i];
        }
        apint_normalize(&a);
        int passed = 1;
        if (a.sign != sign_out || a.exp2 != exp2_out) {
            passed = 0;
        }
        for (int i = 0; i < APINT_MAX_LIMBS; i++) {
            if (a.limb[i] != limbs_out[i]) {
                passed = 0;
                break;
            }
        }
        if (passed) {
            testsPassed++;
        } else {
            printf("\n\nWarning! ");
            printf("Input: sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   sign_in, exp2_in, limbs_in[0], limbs_in[1], limbs_in[2],
                   limbs_in[3], limbs_in[4], limbs_in[5], limbs_in[6], limbs_in[7]);
            printf("Expected: sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   sign_out, exp2_out, limbs_out[0], limbs_out[1], limbs_out[2],
                   limbs_out[3], limbs_out[4], limbs_out[5], limbs_out[6],
                   limbs_out[7]);
            printf("Got:      sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   a.sign, a.exp2, a.limb[0], a.limb[1], a.limb[2], a.limb[3],
                   a.limb[4], a.limb[5], a.limb[6], a.limb[7]);
        }
        totalTests++;
    }
    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
    }
    printf("\n\t ... %d out of %d tests passed.", testsPassed, totalTests);
    printf("\n");
    return totalTests - testsPassed;
}

/*!
 * @brief Tests apint_mul function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_apint_mul(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/apint_mul_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    unsigned long long valA;
    unsigned long long valB;
    int signA;
    int signB;
    int sign_out;
    int exp2_out;
    unsigned int limbs_out[APINT_MAX_LIMBS];
    int scanResult;
    char line[512];
    int testsPassed = 0;
    int totalTests = 0;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf( // NOLINT
            line, "%llu,%d,%llu,%d,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d", &valA, &signA,
            &valB, &signB, &limbs_out[0], &limbs_out[1], &limbs_out[2],
            &limbs_out[3], &limbs_out[4], &limbs_out[5], &limbs_out[6],
            &limbs_out[7], &exp2_out, &sign_out);
        if (scanResult != 14) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 14", scanResult);
            continue;
        }
        apint_t a;
        apint_t b;
        apint_t result;
        apint_set_ull(&a, valA, (signed char)signA);
        apint_set_ull(&b, valB, (signed char)signB);
        apint_mul(&result, &a, &b);

        int passed = 1;
        if (result.sign != sign_out || result.exp2 != exp2_out) {
            passed = 0;
        }
        for (int i = 0; i < APINT_MAX_LIMBS; i++) {
            if (result.limb[i] != limbs_out[i]) {
                passed = 0;
                break;
            }
        }
        if (passed) {
            testsPassed++;
        } else {
            printf("\n\nWarning! ");
            printf("Input: valA=%llu, signA=%d, valB=%llu, signB=%d\n", valA, signA,
                   valB, signB);
            printf("Expected: sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   sign_out, exp2_out, limbs_out[0], limbs_out[1], limbs_out[2],
                   limbs_out[3], limbs_out[4], limbs_out[5], limbs_out[6],
                   limbs_out[7]);
            printf("Got:      sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   result.sign, result.exp2, result.limb[0], result.limb[1],
                   result.limb[2], result.limb[3], result.limb[4], result.limb[5],
                   result.limb[6], result.limb[7]);
        }
        totalTests++;
    }
    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
    }
    printf("\n\t ... %d out of %d tests passed.", testsPassed, totalTests);
    printf("\n");
    return totalTests - testsPassed;
}

/*!
 * @brief Tests apint_shr_bits_sticky function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_apint_shr_bits_sticky(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/apint_shr_sticky_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    int sign_in;
    int exp2_in;
    unsigned int limbs_in[APINT_MAX_LIMBS];
    unsigned int bits;
    int sign_out;
    int exp2_out;
    unsigned int limbs_out[APINT_MAX_LIMBS];
    int scanResult;
    char line[512];
    int testsPassed = 0;
    int totalTests = 0;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf( // NOLINT
            line, "%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d,%u,%u,%u,%u,%u,%u,%u,%u",
            &sign_in, &exp2_in, &limbs_in[0], &limbs_in[1], &limbs_in[2],
            &limbs_in[3], &limbs_in[4], &limbs_in[5], &limbs_in[6], &limbs_in[7],
            &bits, &sign_out, &exp2_out, &limbs_out[0], &limbs_out[1], &limbs_out[2],
            &limbs_out[3], &limbs_out[4], &limbs_out[5], &limbs_out[6],
            &limbs_out[7]);
        if (scanResult != 21) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 21", scanResult);
            continue;
        }
        apint_t src;
        apint_t dst;
        src.sign = (signed char)sign_in;
        src.exp2 = exp2_in;
        src.n = 0;
        for (int i = APINT_MAX_LIMBS - 1; i >= 0; i--) {
            src.limb[i] = limbs_in[i];
            if (limbs_in[i] != 0 && src.n == 0) {
                src.n = (unsigned char)(i + 1);
            }
        }
        apint_shr_bits_sticky(&dst, &src, bits);
        int passed = 1;
        if (dst.sign != sign_out || dst.exp2 != exp2_out) {
            passed = 0;
        }
        for (int i = 0; i < APINT_MAX_LIMBS; i++) {
            if (dst.limb[i] != limbs_out[i]) {
                passed = 0;
                break;
            }
        }
        if (passed) {
            testsPassed++;
        } else {
            printf("\n\nWarning! ");
            printf("Input: sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u], "
                   "bits=%u\n",
                   sign_in, exp2_in, limbs_in[0], limbs_in[1], limbs_in[2],
                   limbs_in[3], limbs_in[4], limbs_in[5], limbs_in[6], limbs_in[7],
                   bits);
            printf("Expected: sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   sign_out, exp2_out, limbs_out[0], limbs_out[1], limbs_out[2],
                   limbs_out[3], limbs_out[4], limbs_out[5], limbs_out[6],
                   limbs_out[7]);
            printf("Got:      sign=%d, exp2=%d, limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   dst.sign, dst.exp2, dst.limb[0], dst.limb[1], dst.limb[2],
                   dst.limb[3], dst.limb[4], dst.limb[5], dst.limb[6], dst.limb[7]);
        }
        totalTests++;
    }
    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
    }
    printf("\n\t ... %d out of %d tests passed.", testsPassed, totalTests);
    printf("\n");
    return totalTests - testsPassed;
}

/*!
 * @brief Main function to run all set zeta derivatives function tests.
 *
 * @return number of failed tests.
 */
int main() {
    printf("start ");
    int failed = 0;
    failed += test_ctzmsb();
    failed += test_apint_set_ull();
    failed += test_apint_normalize();
    failed += test_apint_mul();
    failed += test_apint_shr_bits_sticky();
    return failed;
}
