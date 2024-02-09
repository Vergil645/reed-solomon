/**
 * @file compare_op.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Compares complexity of "+" and "*" operations in field GF(65536).
 * @date 2024-02-07
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <assert.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include <gf65536.h>

#define NUMBER_OF_OPERATIONS 1000000000
#define USE_AVX512 1

void compare_by_clock(const GF_t* gf, const element_t* a, const element_t* b,
                      element_t* c) {
    clock_t start;
    clock_t end;
    double elapsed_time_add;
    double elapsed_time_mul;

    if (USE_AVX512) {
        __m256* a256 = (__m256*)a;
        __m256* b256 = (__m256*)b;
        __m256* c256 = (__m256*)c;
        size_t iter_cnt =
            NUMBER_OF_OPERATIONS * sizeof(element_t) / sizeof(__m256);

        start = clock();
        for (size_t i = 0; i < iter_cnt; ++i) {
            c256[i] = _mm256_xor_ps(a256[i], b256[i]);
        }
        end = clock();

#ifndef NDEBUG
        for (size_t i = 0; i < NUMBER_OF_OPERATIONS; ++i) {
            assert(c[i] == a[i] ^ b[i]);
        }
#endif
    } else {
        start = clock();
        for (size_t i = 0; i < NUMBER_OF_OPERATIONS; ++i) {
            c[i] = a[i] ^ b[i];
        }
        end = clock();
    }

    elapsed_time_add = (end - start) / (double)CLOCKS_PER_SEC;

    start = clock();
    for (size_t i = 0; i < NUMBER_OF_OPERATIONS; ++i) {
        c[i] = gf_mul_ee(gf, a[i], b[i]);
    }
    end = clock();

    elapsed_time_mul = (end - start) / (double)CLOCKS_PER_SEC;

    printf("===== clock() =====\n");
    printf("    ^     time: %.3f secs\n", elapsed_time_add);
    printf("gf_mul_ee time: %.3f secs\n", elapsed_time_mul);
    printf("ratio \"gf_mul_ee/^\": %.3f\n",
           elapsed_time_mul / elapsed_time_add);
}

int main(void) {
    GF_t _gf;
    GF_t* gf = &_gf;
    element_t* a;
    element_t* b;
    element_t* c;
    int ret = 0;

    ret = gf_alloc(gf);
    if (ret) {
        printf("ERROR: gf_alloc returned %d\n", ret);
        return ret;
    }

    gf_init(gf);

    a = (element_t*)aligned_alloc(32, NUMBER_OF_OPERATIONS * sizeof(element_t));
    if (!a) {
        printf("ERROR: couldn't allocate a\n");
        gf_free(gf);
        return 1;
    }

    b = (element_t*)aligned_alloc(32, NUMBER_OF_OPERATIONS * sizeof(element_t));
    if (!b) {
        printf("ERROR: couldn't allocate b\n");
        free(a);
        gf_free(gf);
        return 1;
    }

    c = (element_t*)aligned_alloc(32, NUMBER_OF_OPERATIONS * sizeof(element_t));
    if (!c) {
        printf("ERROR: couldn't allocate c\n");
        free(b);
        free(a);
        gf_free(gf);
        return 1;
    }

    a[0] = rand();
    b[0] = rand();
    for (size_t i = 1; i < NUMBER_OF_OPERATIONS; ++i) {
        a[i] = a[i - 1] + 1;
        b[i] = b[i - 1] + 1;
        c[i] = 0;
    }

    compare_by_clock(gf, a, b, c);

    free(c);
    free(b);
    free(a);
    gf_free(gf);

    return 0;
}