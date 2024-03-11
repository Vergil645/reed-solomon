/**
 * @file gf65536.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief rs/gf65536.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <rs/gf65536.h>

GF_t* gf_create() {
    GF_t* gf;

    gf = (GF_t*)malloc(sizeof(GF_t));
    if (!gf)
        return NULL;
    memset((void*)gf, 0, sizeof(GF_t));

    gf->normal_repr_by_subfield[0] = gf->_normal_repr_by_subfield_memory;
    for (uint8_t i = 1; i < CC_COSET_SIZES_CNT; ++i)
        gf->normal_repr_by_subfield[i] = gf->normal_repr_by_subfield[i - 1] + N;

    element_t* pow_table = gf->pow_table;
    element_t elem;
    uint16_t** normal_repr_by_subfield = gf->normal_repr_by_subfield;
    uint16_t* log_table = gf->log_table;
    poly_t cur_poly = 1;

    for (uint16_t i = 0; i < N; ++i) {
        pow_table[i] = (element_t)cur_poly;
        log_table[pow_table[i]] = i;

        cur_poly <<= 1;
        if (cur_poly & GF_FIELD_SIZE)
            cur_poly ^= GF_PRIMITIVE_POLY;
    }

    for (uint32_t i = N; i < (N << 1) - 1; ++i)
        pow_table[i] = pow_table[i - N];

    for (uint8_t i = 0; i < CC_COSET_SIZES_CNT; ++i) {
        element_t* normal_basis = g_normal_basis_by_subfield[i];
        uint8_t m = 1 << i; // cyclotomic coset size

        memset((void*)normal_repr_by_subfield[i], 0, N * sizeof(uint16_t));

        for (uint32_t repr = 1; repr != (1 << m); ++repr) {
            elem = 0;
            for (uint8_t j = 0; j < m; ++j) {
                if (repr & (1 << j))
                    elem ^= normal_basis[j];
            }

            assert(elem != 0);
            assert(normal_repr_by_subfield[i][log_table[elem]] == 0);

            normal_repr_by_subfield[i][log_table[elem]] = (uint16_t)repr;
        }
    }

    return gf;
}

void gf_destroy(GF_t* gf) {
    assert(gf != NULL);

    free(gf);
}

element_t gf_mul_ee(GF_t* gf, element_t a, element_t b) {
    assert(gf != NULL);

    if (a == 0 || b == 0)
        return 0;

    uint16_t* log_table = gf->log_table;
    uint32_t a_log = (uint32_t)log_table[a];
    uint32_t b_log = (uint32_t)log_table[b];

    return gf->pow_table[a_log + b_log];
}

element_t gf_div_ee(GF_t* gf, element_t a, element_t b) {
    assert(gf != NULL);
    assert(b != 0);

    if (a == 0)
        return 0;

    uint16_t* log_table = gf->log_table;
    uint32_t a_log = (uint32_t)log_table[a];
    uint32_t b_log = (uint32_t)log_table[b];

    return gf->pow_table[(N + a_log - b_log) % N];
}

#define GF_USE_NATIVE_LIBRARY 0

#if GF_USE_NATIVE_LIBRARY == 1
#include <immintrin.h>
#endif

void gf_add(void* a, const void* b, size_t symbol_size) {
    assert(symbol_size % sizeof(element_t) == 0);

// TODO: implement native gf_add
#if GF_USE_NATIVE_LIBRARY == 1
    // element_t* data_1 = (element_t*)a;
    // element_t* data_2 = (element_t*)b;
    // size_t max_idx = symbol_size / sizeof(element_t);
    // register __m256i in, out;

    // for (element_t* end_1 = data_1 + max_idx; data_1 != end_1; ++data_1,
    // ++data_2) {
    //     in = _mm256_loadu_si256((__m256i_u*)data_2);
    //     out = _mm256_loadu_si256((__m256i_u*)data_1);
    //     out = _mm256_xor_si256(in, out);
    //     _mm256_storeu_si256((__m256i_u*)data_1, out);
    // }
#else
    uint64_t* data64_1 = (uint64_t*)a;
    uint64_t* data64_2 = (uint64_t*)b;
    size_t max_64_idx = symbol_size / sizeof(uint64_t);

    for (const uint64_t* end64_1 = data64_1 + max_64_idx; data64_1 != end64_1;
         ++data64_1, ++data64_2)
        *data64_1 ^= *data64_2;

    element_t* data_1 = (element_t*)data64_1;
    element_t* data_2 = (element_t*)data64_2;
    size_t max_idx = symbol_size / sizeof(element_t);

    for (const element_t* end_1 = (element_t*)a + max_idx; data_1 != end_1; ++data_1, ++data_2)
        *data_1 ^= *data_2;
#endif
}

void gf_mul(GF_t* gf, void* a, element_t coef, size_t symbol_size) {
    assert(symbol_size % sizeof(element_t) == 0);

    if (coef == 0) {
        uint64_t* data64 = (uint64_t*)a;
        size_t max_64_idx = symbol_size / sizeof(uint64_t);

        for (const uint64_t* end64_1 = data64 + max_64_idx; data64 != end64_1; ++data64)
            *data64 = 0;

        element_t* data = (element_t*)data64;
        size_t max_idx = symbol_size / sizeof(element_t);

        for (const element_t* end_1 = (element_t*)a + max_idx; data != end_1; ++data)
            *data = 0;

        return;
    }

    if (coef == 1)
        return;

    element_t* pow_table_shifted;
    element_t* data = (element_t*)a;
    uint16_t* log_table = gf->log_table;
    size_t max_idx = symbol_size / sizeof(element_t);

    pow_table_shifted = gf->pow_table + log_table[coef];

    for (const element_t* end = data + max_idx; data != end; ++data) {
        element_t val = *data;
        if (val != 0)
            *data = pow_table_shifted[log_table[val]];
    }
}

void gf_madd(GF_t* gf, void* a, element_t coef, const void* b, size_t symbol_size) {
    assert(symbol_size % sizeof(element_t) == 0);

    if (coef == 0)
        return;

    if (coef == 1) {
        gf_add(a, b, symbol_size);
        return;
    }

    element_t* pow_table_shifted;
    element_t* data_1 = (element_t*)a;
    element_t* data_2 = (element_t*)b;
    uint16_t* log_table = gf->log_table;
    size_t max_idx = symbol_size / sizeof(element_t);

    pow_table_shifted = gf->pow_table + log_table[coef];

    for (const element_t* end_1 = data_1 + max_idx; data_1 != end_1; ++data_1, ++data_2) {
        element_t val_2 = *data_2;
        if (val_2 != 0)
            *data_1 ^= pow_table_shifted[log_table[val_2]];
    }
}
