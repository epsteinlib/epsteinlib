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
 * @brief Main function to run all set zeta derivatives function tests.
 *
 * @return number of failed tests.
 */
int main() {
    printf("start ");
    int failed = 0;
    failed += test_ctzmsb();
    failed += test_apint_set_ull();
    return failed;
}
