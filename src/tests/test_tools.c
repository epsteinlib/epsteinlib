// SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
//
// SPDX-License-Identifier: AGPL-3.0-only

#include "../tools.h"
#include "utils.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 1024
#endif

#ifndef BASE_PATH
#define BASE_PATH "csv"
#endif

#ifndef MAX_DIM
#define MAX_DIM 12
#endif

/*!
 * @brief Tests matrix invert-transpose against Mathematica reference values.
 * CSV format per line: dim, m[0..dim²-1], ref_transpose_inv[0..dim²-1]
 *
 * @return number of failed tests.
 */
int test_matrix_transpose_inverse(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/sdv_Ref.csv", BASE_PATH);
    if (result < 0 || result >= (int)sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    char line[4096];
    int testsPassed = 0;
    int totalTests = 0;
    const double tol = 5e-14;
    double errMin = NAN;
    double errMax = NAN;
    double errSum = 0.;

    printf("\n\t ... ");
    printf("processing %s ", path);

    while (fgets(line, sizeof(line), data) != NULL) {
        char *ptr = line;
        char *endptr = NULL;

        // parse dimension
        unsigned int dim = (unsigned int)strtoul(ptr, &endptr, 10);
        if (endptr == ptr || *endptr != ',' || dim == 0 || dim > MAX_DIM) {
            printf("Error reading line: %s\n", line);
            printf("Invalid dimension: %u\n", dim);
            continue;
        }
        ptr = endptr + 1;

        unsigned int n = dim * dim;

        // parse input matrix m (row-major, dim x dim)
        double m[dim * dim]; // NOLINT
        int ok = 1;
        for (unsigned int i = 0; i < n; i++) {
            m[i] = strtod(ptr, &endptr);
            if (endptr == ptr || *endptr != ',') {
                ok = 0;
                break;
            }
            ptr = endptr + 1;
        }
        if (!ok) {
            printf("Error reading line: %s\n", line);
            printf("Could not parse input matrix (expected %u values)\n", n);
            continue;
        }

        // parse reference: transpose(inverse(m))
        double ref[dim * dim];
        for (unsigned int i = 0; i < n; i++) {
            ref[i] = strtod(ptr, &endptr);
            if (endptr == ptr) {
                ok = 0;
                break;
            }
            if (i < n - 1) {
                if (*endptr != ',') {
                    ok = 0;
                    break;
                }
                ptr = endptr + 1;
            }
        }
        if (!ok) {
            printf("Error reading line: %s\n", line);
            printf("Could not parse reference matrix (expected %u values)\n", n);
            continue;
        }

        // compute transpose(inverse(m)) via invert()
        double m_copy[dim * dim];
        double m_fourier[dim * dim];
        int p[dim];
        for (unsigned int i = 0; i < n; i++) {
            m_copy[i] = m[i];
        }
        invert(dim, m_copy, p, m_fourier);
        transpose(dim, m_fourier);

        // infinity norm of (m_fourier - ref)
        double diff[dim * dim];
        for (unsigned int i = 0; i < n; i++) {
            diff[i] = m_fourier[i] - ref[i];
        }
        double err = inf_norm(dim, diff);

        errMin = (errMin < err) ? errMin : err;
        errMax = (errMax > err) ? errMax : err;
        errSum += err;

        if (err < tol) {
            testsPassed++;
        } else {
            printf("\n");
            printf("Warning! ");
            printf("matrix transpose inverse (dim=%u):\n", dim);
            printf("inf_norm(m_fourier - ref): %E !< %E  (tolerance)\n", err, tol);
            printf("\n");
            printMatrixUnitTest("m:", m, dim);
            printMatrixUnitTest("reci:", m_fourier, dim);
            printMatrixUnitTest("ref:", ref, dim);
            printf("\n");
        }
        totalTests++;
    }

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d\n", errno);
    }
    printf("\n\t ... ");
    printf("%d out of %d tests passed with tolerance %E.", testsPassed, totalTests,
           tol);
    printf("\t    ");
    printf("[ Error →  min: %E | max: %E | avg: %E ]", errMin, errMax,
           errSum / totalTests);
    printf("\n");
    return totalTests - testsPassed;
}

/*!
 * @brief Helper function to test hpdyad_add against a CSV reference file
 *
 * @param filename Name of the CSV file containing test data
 * @return Number of failed tests
 */
static int test_hpdyad_add_from_file(const char *filename) { // NOLINT
    printf("\t ... processing %s ", filename);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/%s", BASE_PATH, filename);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }

    char line[4096];
    int testsPassed = 0;
    int totalTests = 0;
    int i;

    while (fgets(line, sizeof(line), data) != NULL) {
        char *ptr = line;
        char *endptr;

        // Parse a: sign, exp2, 32 limbs
        int sign_a = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        int exp2_a = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        unsigned int limb_a[HPDYAD_MAX_LIMBS];
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            limb_a[i] = (unsigned int)strtoul(ptr, &endptr, 10);
            if (*endptr != ',') {
                break;
            }
            ptr = endptr + 1;
        }
        if (i < HPDYAD_MAX_LIMBS) {
            continue;
        }

        // Parse b: sign, exp2, 32 limbs
        int sign_b = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        int exp2_b = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        unsigned int limb_b[HPDYAD_MAX_LIMBS];
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            limb_b[i] = (unsigned int)strtoul(ptr, &endptr, 10);
            if (*endptr != ',') {
                break;
            }
            ptr = endptr + 1;
        }
        if (i < HPDYAD_MAX_LIMBS) {
            continue;
        }

        // Parse expected: sign, exp2, 32 limbs
        int sign_expected = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        int exp2_expected = (int)strtol(ptr, &endptr, 10);
        if (*endptr != ',') {
            continue;
        }
        ptr = endptr + 1;

        unsigned int limb_expected[HPDYAD_MAX_LIMBS];
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            limb_expected[i] = (unsigned int)strtoul(ptr, &endptr, 10);
            if (i < HPDYAD_MAX_LIMBS - 1 && *endptr != ',') {
                break;
            }
            ptr = endptr + 1;
        }
        if (i < HPDYAD_MAX_LIMBS - 1) {
            continue;
        }

        // Build hpdyad a
        hpdyad_t a;
        a.sign = (signed char)sign_a;
        a.exp2 = exp2_a;
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            a.limb[i] = limb_a[i];
        }
        a.n = 0;
        for (i = HPDYAD_MAX_LIMBS - 1; i >= 0; i--) {
            if (a.limb[i] != 0) {
                a.n = (unsigned char)(i + 1);
                break;
            }
        }
        if (a.n == 0) {
            a.sign = 1;
            a.exp2 = 0;
        }

        // Build hpdyad b
        hpdyad_t b;
        b.sign = (signed char)sign_b;
        b.exp2 = exp2_b;
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            b.limb[i] = limb_b[i];
        }
        b.n = 0;
        for (i = HPDYAD_MAX_LIMBS - 1; i >= 0; i--) {
            if (b.limb[i] != 0) {
                b.n = (unsigned char)(i + 1);
                break;
            }
        }
        if (b.n == 0) {
            b.sign = 1;
            b.exp2 = 0;
        }

        // Build hpdyad expected
        hpdyad_t expected;
        expected.sign = (signed char)sign_expected;
        expected.exp2 = exp2_expected;
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            expected.limb[i] = limb_expected[i];
        }
        expected.n = 0;
        for (i = HPDYAD_MAX_LIMBS - 1; i >= 0; i--) {
            if (expected.limb[i] != 0) {
                expected.n = (unsigned char)(i + 1);
                break;
            }
        }
        if (expected.n == 0) {
            expected.sign = 1;
            expected.exp2 = 0;
        }

        // Perform addition
        hpdyad_t result_hpdyad;
        hpdyad_add(&result_hpdyad, &a, &b);

        // Compare
        int passed =
            (result_hpdyad.sign == expected.sign &&
             result_hpdyad.exp2 == expected.exp2 && result_hpdyad.n == expected.n);
        if (passed) {
            for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
                if (result_hpdyad.limb[i] != expected.limb[i]) {
                    passed = 0;
                    break;
                }
            }
        }

        if (passed) {
            testsPassed++;
        } else {
            printf("\n\nWarning! Test %d failed:", totalTests + 1);
            printf("\n  a: sign=%d exp2=%d n=%d", a.sign, a.exp2, a.n);
            printf("\n  b: sign=%d exp2=%d n=%d", b.sign, b.exp2, b.n);
            printf("\n  result:   sign=%d exp2=%d n=%d", result_hpdyad.sign,
                   result_hpdyad.exp2, result_hpdyad.n);
            printf("\n  expected: sign=%d exp2=%d n=%d\n", expected.sign,
                   expected.exp2, expected.n);
        }
        totalTests++;
    }

    if (fclose(data) != 0) {
        return fprintf(stderr, "Error closing file: %d", errno);
    }

    printf("\n\t ... %d out of %d tests passed.\n", testsPassed, totalTests);
    return totalTests - testsPassed;
}

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
 * @brief Tests hpdyad_set_ull function against Mathematica reference values
 *
 * @return number of failed tests.
 * */
int test_hpdyad_set_ull(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/hpdyad_set_ull_Ref.csv", BASE_PATH);
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

        hpdyad_t a;
        hpdyad_set_ull(&a, x, (signed char)sign_input);

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
 * @brief Tests hpdyad_normalize function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_hpdyad_normalize(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/hpdyad_normalize_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    int sign_in;
    int exp2_in;
    unsigned int limbs_in[HPDYAD_MAX_LIMBS];
    int sign_out;
    int exp2_out;
    unsigned int limbs_out[HPDYAD_MAX_LIMBS];
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
        hpdyad_t a;
        a.sign = (signed char)sign_in;
        a.exp2 = exp2_in;
        a.n = HPDYAD_MAX_LIMBS;
        for (int i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            a.limb[i] = limbs_in[i];
        }
        hpdyad_normalize(&a);
        int passed = 1;
        if (a.sign != sign_out || a.exp2 != exp2_out) {
            passed = 0;
        }
        for (int i = 0; i < HPDYAD_MAX_LIMBS; i++) {
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
 * @brief Tests hpdyad_mul function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_hpdyad_mul(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result = snprintf(path, sizeof(path), "%s/hpdyad_mul_Ref.csv", BASE_PATH);
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
    unsigned int limbs_out[HPDYAD_MAX_LIMBS];
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
        hpdyad_t a;
        hpdyad_t b;
        hpdyad_t result;
        hpdyad_set_ull(&a, valA, (signed char)signA);
        hpdyad_set_ull(&b, valB, (signed char)signB);
        hpdyad_mul(&result, &a, &b);

        int passed = 1;
        if (result.sign != sign_out || result.exp2 != exp2_out) {
            passed = 0;
        }
        for (int i = 0; i < HPDYAD_MAX_LIMBS; i++) {
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
 * @brief Tests hpdyad_shr_bits_sticky function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_hpdyad_shr_bits_sticky(void) {
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/hpdyad_shr_sticky_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    int sign_in;
    int exp2_in;
    unsigned int limbs_in[HPDYAD_MAX_LIMBS];
    unsigned int bits;
    int sign_out;
    int exp2_out;
    unsigned int limbs_out[HPDYAD_MAX_LIMBS];
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
        hpdyad_t src;
        hpdyad_t dst;
        src.sign = (signed char)sign_in;
        src.exp2 = exp2_in;
        src.n = 0;
        for (int i = HPDYAD_MAX_LIMBS - 1; i >= 0; i--) {
            src.limb[i] = limbs_in[i];
            if (limbs_in[i] != 0 && src.n == 0) {
                src.n = (unsigned char)(i + 1);
            }
        }
        hpdyad_shr_bits_sticky(&dst, &src, bits);
        int passed = 1;
        if (dst.sign != sign_out || dst.exp2 != exp2_out) {
            passed = 0;
        }
        for (int i = 0; i < HPDYAD_MAX_LIMBS; i++) {
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
 * @brief Tests hpdyad_add function against Mathematica reference values (integers)
 *
 * @return Number of failed tests
 */
int test_hpdyad_add(void) {
    printf("%s\n", __func__);
    return test_hpdyad_add_from_file("hpdyad_add_Ref.csv");
}

/*!
 * @brief Tests hpdyad_add function against Mathematica reference values (rationals)
 *
 * @return Number of failed tests
 */
int test_hpdyad_add_fractional(void) {
    printf("%s\n", __func__);
    return test_hpdyad_add_from_file("hpdyad_add_fractional_Ref.csv");
}

/*!
 * @brief Tests hpdyad_to_double function against Mathematica reference values
 *
 * @return number of failed tests.
 */
int test_hpdyad_to_double(void) { // NOLINT
    printf("%s ", __func__);
    char path[MAX_PATH_LENGTH];
    int result =
        snprintf(path, sizeof(path), "%s/hpdyad_to_double_Ref.csv", BASE_PATH);
    if (result < 0 || result >= sizeof(path)) {
        return fprintf(stderr, "Error creating file path\n");
    }
    FILE *data = fopen(path, "r");
    if (data == NULL) {
        return fprintf(stderr, "Error opening file: %s\n", path);
    }
    int sign_input;
    int exp2_input;
    unsigned int limbs_input[HPDYAD_MAX_LIMBS];
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
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            limbs_input[i] = (unsigned int)strtoul(ptr, &endptr, 10);
            if (*endptr != ',') {
                break;
            }
            ptr = endptr + 1;
        }
        if (i < HPDYAD_MAX_LIMBS) {
            continue;
        }
        // Parse reference double
        ref_double = strtod(ptr, &endptr);
        // Build hpdyad
        hpdyad_t a;
        a.sign = (signed char)sign_input;
        a.exp2 = exp2_input;
        for (i = 0; i < HPDYAD_MAX_LIMBS; i++) {
            a.limb[i] = limbs_input[i];
        }
        // Compute n (number of used limbs)
        a.n = 0;
        for (i = HPDYAD_MAX_LIMBS - 1; i >= 0; i--) {
            if (a.limb[i] != 0) {
                a.n = (unsigned char)(i + 1);
                break;
            }
        }
        // Test
        double result_double = hpdyad_to_double(&a);
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
    int failed = 0;
    failed += test_ctzmsb();
    failed += test_hpdyad_set_ull();
    failed += test_hpdyad_normalize();
    failed += test_hpdyad_mul();
    failed += test_hpdyad_shr_bits_sticky();
    failed += test_hpdyad_add();
    failed += test_hpdyad_add_fractional();
    failed += test_hpdyad_to_double();
    failed += run_timed_test(test_matrix_transpose_inverse);
    return failed;
}
