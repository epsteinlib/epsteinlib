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
    char line[2048];
    int testsPassed = 0;
    int totalTests = 0;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf( // NOLINT
            line,
            "%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%"
            "u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,"
            "%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u",
            &sign_in, &exp2_in, &limbs_in[0], &limbs_in[1], &limbs_in[2],
            &limbs_in[3], &limbs_in[4], &limbs_in[5], &limbs_in[6], &limbs_in[7],
            &limbs_in[8], &limbs_in[9], &limbs_in[10], &limbs_in[11], &limbs_in[12],
            &limbs_in[13], &limbs_in[14], &limbs_in[15], &limbs_in[16],
            &limbs_in[17], &limbs_in[18], &limbs_in[19], &limbs_in[20],
            &limbs_in[21], &limbs_in[22], &limbs_in[23], &limbs_in[24],
            &limbs_in[25], &limbs_in[26], &limbs_in[27], &limbs_in[28],
            &limbs_in[29], &limbs_in[30], &limbs_in[31], &sign_out, &exp2_out,
            &limbs_out[0], &limbs_out[1], &limbs_out[2], &limbs_out[3],
            &limbs_out[4], &limbs_out[5], &limbs_out[6], &limbs_out[7],
            &limbs_out[8], &limbs_out[9], &limbs_out[10], &limbs_out[11],
            &limbs_out[12], &limbs_out[13], &limbs_out[14], &limbs_out[15],
            &limbs_out[16], &limbs_out[17], &limbs_out[18], &limbs_out[19],
            &limbs_out[20], &limbs_out[21], &limbs_out[22], &limbs_out[23],
            &limbs_out[24], &limbs_out[25], &limbs_out[26], &limbs_out[27],
            &limbs_out[28], &limbs_out[29], &limbs_out[30], &limbs_out[31]);
        if (scanResult != 68) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 68", scanResult);
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
            printf("Input: sign=%d, exp2=%d, "
                   "limbs=[%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   sign_in, exp2_in, limbs_in[0], limbs_in[1], limbs_in[2],
                   limbs_in[3], limbs_in[4], limbs_in[5], limbs_in[6], limbs_in[7],
                   limbs_in[8], limbs_in[9], limbs_in[10], limbs_in[11],
                   limbs_in[12], limbs_in[13], limbs_in[14], limbs_in[15]);
            printf("Expected: sign=%d, exp2=%d, "
                   "limbs=[%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   sign_out, exp2_out, limbs_out[0], limbs_out[1], limbs_out[2],
                   limbs_out[3], limbs_out[4], limbs_out[5], limbs_out[6],
                   limbs_out[7], limbs_out[8], limbs_out[9], limbs_out[10],
                   limbs_out[11], limbs_out[12], limbs_out[13], limbs_out[14],
                   limbs_out[15]);
            printf("Got:      sign=%d, exp2=%d, "
                   "limbs=[%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   a.sign, a.exp2, a.limb[0], a.limb[1], a.limb[2], a.limb[3],
                   a.limb[4], a.limb[5], a.limb[6], a.limb[7], a.limb[8], a.limb[9],
                   a.limb[10], a.limb[11], a.limb[12], a.limb[13], a.limb[14],
                   a.limb[15]);
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
    char line[1024];
    int testsPassed = 0;
    int totalTests = 0;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf( // NOLINT
            line,
            "%llu,%d,%llu,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%"
            "u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d",
            &valA, &signA, &valB, &signB, &limbs_out[0], &limbs_out[1],
            &limbs_out[2], &limbs_out[3], &limbs_out[4], &limbs_out[5],
            &limbs_out[6], &limbs_out[7], &limbs_out[8], &limbs_out[9],
            &limbs_out[10], &limbs_out[11], &limbs_out[12], &limbs_out[13],
            &limbs_out[14], &limbs_out[15], &limbs_out[16], &limbs_out[17],
            &limbs_out[18], &limbs_out[19], &limbs_out[20], &limbs_out[21],
            &limbs_out[22], &limbs_out[23], &limbs_out[24], &limbs_out[25],
            &limbs_out[26], &limbs_out[27], &limbs_out[28], &limbs_out[29],
            &limbs_out[30], &limbs_out[31], &exp2_out, &sign_out);
        if (scanResult != 38) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 38", scanResult);
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
    char line[1024];
    int testsPassed = 0;
    int totalTests = 0;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf( // NOLINT
            line,
            "%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%"
            "u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,"
            "%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u",
            &sign_in, &exp2_in, &limbs_in[0], &limbs_in[1], &limbs_in[2],
            &limbs_in[3], &limbs_in[4], &limbs_in[5], &limbs_in[6], &limbs_in[7],
            &limbs_in[8], &limbs_in[9], &limbs_in[10], &limbs_in[11], &limbs_in[12],
            &limbs_in[13], &limbs_in[14], &limbs_in[15], &limbs_in[16],
            &limbs_in[17], &limbs_in[18], &limbs_in[19], &limbs_in[20],
            &limbs_in[21], &limbs_in[22], &limbs_in[23], &limbs_in[24],
            &limbs_in[25], &limbs_in[26], &limbs_in[27], &limbs_in[28],
            &limbs_in[29], &limbs_in[30], &limbs_in[31], &bits, &sign_out, &exp2_out,
            &limbs_out[0], &limbs_out[1], &limbs_out[2], &limbs_out[3],
            &limbs_out[4], &limbs_out[5], &limbs_out[6], &limbs_out[7],
            &limbs_out[8], &limbs_out[9], &limbs_out[10], &limbs_out[11],
            &limbs_out[12], &limbs_out[13], &limbs_out[14], &limbs_out[15],
            &limbs_out[16], &limbs_out[17], &limbs_out[18], &limbs_out[19],
            &limbs_out[20], &limbs_out[21], &limbs_out[22], &limbs_out[23],
            &limbs_out[24], &limbs_out[25], &limbs_out[26], &limbs_out[27],
            &limbs_out[28], &limbs_out[29], &limbs_out[30], &limbs_out[31]);
        if (scanResult != 69) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 69", scanResult);
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
 * @brief Tests apint_add function against Mathematica reference values
 *
 * @return number of failed tests.
 * */
int test_apint_add(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/apint_add_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    int sign_a;
    int exp2_a;
    unsigned int limb_a[APINT_MAX_LIMBS];
    int sign_b;
    int exp2_b;
    unsigned int limb_b[APINT_MAX_LIMBS];
    int sign_expected;
    int exp2_expected;
    unsigned int limb_expected[APINT_MAX_LIMBS];

    int scanResult;
    char line[2048];
    int testsPassed = 0;
    int totalTests = 0;

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        scanResult = sscanf( // NOLINT
            line,
            "%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%"
            "u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,"
            "%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%d,%d,%"
            "u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,"
            "%u,%u,%u,%u,%u,%u,%u,%u",
            &sign_a, &exp2_a, &limb_a[0], &limb_a[1], &limb_a[2], &limb_a[3],
            &limb_a[4], &limb_a[5], &limb_a[6], &limb_a[7], &limb_a[8], &limb_a[9],
            &limb_a[10], &limb_a[11], &limb_a[12], &limb_a[13], &limb_a[14],
            &limb_a[15], &limb_a[16], &limb_a[17], &limb_a[18], &limb_a[19],
            &limb_a[20], &limb_a[21], &limb_a[22], &limb_a[23], &limb_a[24],
            &limb_a[25], &limb_a[26], &limb_a[27], &limb_a[28], &limb_a[29],
            &limb_a[30], &limb_a[31], &sign_b, &exp2_b, &limb_b[0], &limb_b[1],
            &limb_b[2], &limb_b[3], &limb_b[4], &limb_b[5], &limb_b[6], &limb_b[7],
            &limb_b[8], &limb_b[9], &limb_b[10], &limb_b[11], &limb_b[12],
            &limb_b[13], &limb_b[14], &limb_b[15], &limb_b[16], &limb_b[17],
            &limb_b[18], &limb_b[19], &limb_b[20], &limb_b[21], &limb_b[22],
            &limb_b[23], &limb_b[24], &limb_b[25], &limb_b[26], &limb_b[27],
            &limb_b[28], &limb_b[29], &limb_b[30], &limb_b[31], &sign_expected,
            &exp2_expected, &limb_expected[0], &limb_expected[1], &limb_expected[2],
            &limb_expected[3], &limb_expected[4], &limb_expected[5],
            &limb_expected[6], &limb_expected[7], &limb_expected[8],
            &limb_expected[9], &limb_expected[10], &limb_expected[11],
            &limb_expected[12], &limb_expected[13], &limb_expected[14],
            &limb_expected[15], &limb_expected[16], &limb_expected[17],
            &limb_expected[18], &limb_expected[19], &limb_expected[20],
            &limb_expected[21], &limb_expected[22], &limb_expected[23],
            &limb_expected[24], &limb_expected[25], &limb_expected[26],
            &limb_expected[27], &limb_expected[28], &limb_expected[29],
            &limb_expected[30], &limb_expected[31]);
        if (scanResult != 102) {
            printf("\n\t Error reading line: %s", line);
            printf("\t Scanned %d values instead of 102", scanResult);
            continue;
        } // Build apint structures from parsed data
        apint_t a;
        apint_t b;
        apint_t expected;
        apint_t result;
        unsigned char i;

        // Build a
        a.sign = (signed char)sign_a;
        a.exp2 = exp2_a;
        a.n = 0;
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            a.limb[i] = limb_a[i];
            if (limb_a[i] != 0) {
                a.n = i + 1;
            }
        }
        if (a.n == 0) {
            a.sign = 1;
            a.exp2 = 0;
        }

        // Build b
        b.sign = (signed char)sign_b;
        b.exp2 = exp2_b;
        b.n = 0;
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            b.limb[i] = limb_b[i];
            if (limb_b[i] != 0) {
                b.n = i + 1;
            }
        }
        if (b.n == 0) {
            b.sign = 1;
            b.exp2 = 0;
        }

        // Build expected
        expected.sign = (signed char)sign_expected;
        expected.exp2 = exp2_expected;
        expected.n = 0;
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            expected.limb[i] = limb_expected[i];
            if (limb_expected[i] != 0) {
                expected.n = i + 1;
            }
        }
        if (expected.n == 0) {
            expected.sign = 1;
            expected.exp2 = 0;
        }

        // Perform addition
        apint_add(&result, &a, &b);

        // Compare result with expected
        int passed = 1;
        if (result.sign != expected.sign) {
            passed = 0;
        }
        if (result.exp2 != expected.exp2) {
            passed = 0;
        }
        if (result.n != expected.n) {
            passed = 0;
        }
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            if (result.limb[i] != expected.limb[i]) {
                passed = 0;
                break;
            }
        }

        if (passed) {
            testsPassed++;
        } else {
            printf("\n\nWarning! Test %d failed:", totalTests + 1);
            printf("\n  a: sign=%d exp2=%d limbs=[%u,%u,%u,%u,%u,%u,%u,%u]", a.sign,
                   a.exp2, a.limb[0], a.limb[1], a.limb[2], a.limb[3], a.limb[4],
                   a.limb[5], a.limb[6], a.limb[7]);
            printf("\n  b: sign=%d exp2=%d limbs=[%u,%u,%u,%u,%u,%u,%u,%u]", b.sign,
                   b.exp2, b.limb[0], b.limb[1], b.limb[2], b.limb[3], b.limb[4],
                   b.limb[5], b.limb[6], b.limb[7]);
            printf(
                "\n  result:   sign=%d exp2=%d n=%d limbs=[%u,%u,%u,%u,%u,%u,%u,%u]",
                result.sign, result.exp2, result.n, result.limb[0], result.limb[1],
                result.limb[2], result.limb[3], result.limb[4], result.limb[5],
                result.limb[6], result.limb[7]);
            printf("\n  expected: sign=%d exp2=%d n=%d "
                   "limbs=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
                   expected.sign, expected.exp2, expected.n, expected.limb[0],
                   expected.limb[1], expected.limb[2], expected.limb[3],
                   expected.limb[4], expected.limb[5], expected.limb[6],
                   expected.limb[7]);
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
 * @brief Tests apint_to_double function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_apint_to_double(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/apint_to_double_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    int sign_input;
    int exp2_input;
    unsigned int limbs_input[APINT_MAX_LIMBS];
    double ref_double;
    char line[4096];
    int testsPassed = 0;
    int totalTests = 0;
    int i;
    printf("\n\t ... ");
    printf("processing %s ", path);
    while (fgets(line, sizeof(line), data) != NULL) {
        char *ptr = line;
        char *endptr;
        // Parse sign
        sign_input = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;
        // Parse exp2
        exp2_input = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;
        // Parse 32 limbs
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            limbs_input[i] = (unsigned int)strtoul(ptr, &endptr, 10);
            if (*endptr != ',') {
                break;
            }
            ptr = endptr + 1;
        }
        if (i < APINT_MAX_LIMBS) {
            continue;
        }
        // Parse reference double
        ref_double = strtod(ptr, &endptr);
        // Build apint
        apint_t a;
        a.sign = (signed char)sign_input;
        a.exp2 = exp2_input;
        for (i = 0; i < APINT_MAX_LIMBS; i++) {
            a.limb[i] = limbs_input[i];
        }
        // Compute n (number of used limbs)
        a.n = 0;
        for (i = APINT_MAX_LIMBS - 1; i >= 0; i--) {
            if (a.limb[i] != 0) {
                a.n = (unsigned char)(i + 1);
                break;
            }
        }
        // Test
        double result_double = apint_to_double(&a);
        int passed;
        if (ref_double == 0.0) {
            passed = (result_double == 0.0);
        } else {
            passed = (result_double == ref_double);
        }
        if (passed) {
            testsPassed++;
        } else {
            printf("\n\nWarning! ");
            printf("sign=%d, exp2=%d, n=%d: result=%.17g (expected %.17g)\n",
                   sign_input, exp2_input, a.n, result_double, ref_double);
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
    failed += test_apint_add();
    failed += test_apint_to_double();
    return failed;
}
