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

static inline void _gf_fill_normal_bases(element_t* normal_bases) {
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(1) + 0] = 1;

    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(2) + 0] = 44234;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(2) + 1] = 44235;

    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(4) + 0] = 10800;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(4) + 1] = 47860;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(4) + 2] = 34555;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(4) + 3] = 5694;

    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 0] = 16402;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 1] = 53598;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 2] = 44348;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 3] = 63986;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 4] = 22060;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 5] = 64366;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 6] = 6088;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(8) + 7] = 32521;

    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 0] = 2048;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 1] = 2880;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 2] = 7129;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 3] = 30616;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 4] = 2643;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 5] = 6897;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 6] = 29685;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 7] = 7378;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 8] = 30100;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 9] = 2743;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 10] = 20193;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 11] = 36223;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 12] = 24055;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 13] = 41458;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 14] = 41014;
    normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(16) + 15] = 61451;
}

GF_t* gf_create() {
    GF_t* gf;

    gf = (GF_t*)malloc(sizeof(GF_t));
    if (!gf)
        return NULL;
    memset((void*)gf, 0, sizeof(GF_t));

    gf->normal_repr_by_subfield[1] = gf->_normal_repr_by_subfield_memory;
    for (uint8_t i = 1; i < CC_COSET_SIZES_CNT; ++i)
        gf->normal_repr_by_subfield[1 << i] = gf->normal_repr_by_subfield[1 << (i - 1)] + N;

    _gf_fill_normal_bases(gf->normal_bases);

    element_t* pow_table = gf->pow_table;
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
        uint8_t m = 1 << i; // cyclotomic coset size
        element_t* normal_basis = gf->normal_bases + GF_NORMAL_BASES_FIRST_IDX_BY_M(m);

        memset((void*)normal_repr_by_subfield[m], 0, N * sizeof(uint16_t));

        for (uint32_t repr = 1; repr != (1 << m); ++repr) {
            element_t elem = 0;
            for (uint8_t j = 0; j < m; ++j) {
                if (repr & (1 << j))
                    elem ^= normal_basis[j];
            }

            assert(elem != 0);
            assert(normal_repr_by_subfield[m][log_table[elem]] == 0);

            normal_repr_by_subfield[m][log_table[elem]] = (uint16_t)repr;
        }
    }

    return gf;
}

void gf_destroy(GF_t* gf) {
    assert(gf != NULL);

    free(gf);
}

inline element_t gf_get_normal_basis_element(GF_t* gf, uint8_t m, uint8_t i) {
    assert(gf != NULL);
    assert(i < m);

    return gf->normal_bases[GF_NORMAL_BASES_FIRST_IDX_BY_M(m) + i];
}

inline uint16_t gf_get_normal_repr(GF_t* gf, uint8_t m, uint16_t d) {
    assert(gf != NULL);

    return gf->normal_repr_by_subfield[m][d];
}

element_t gf_mul_ee(GF_t* gf, element_t a, element_t b) {
    assert(gf != NULL);

    if (a == 0 || b == 0)
        return 0;

    uint16_t* log_table = gf->log_table;

    return gf->pow_table[(uint32_t)log_table[a] + (uint32_t)log_table[b]];
}

element_t gf_div_ee(GF_t* gf, element_t a, element_t b) {
    assert(gf != NULL);
    assert(b != 0);

    if (a == 0)
        return 0;

    uint16_t* log_table = gf->log_table;

    return gf->pow_table[(N + (uint32_t)log_table[a] - (uint32_t)log_table[b]) % N];
}

void gf_add(void* a, const void* b, size_t symbol_size) {
    assert(symbol_size % sizeof(element_t) == 0);

    uint64_t* data64_1 = (uint64_t*)a;
    uint64_t* data64_2 = (uint64_t*)b;

    for (const uint64_t* end64_1 = data64_1 + symbol_size / sizeof(uint64_t); data64_1 != end64_1;
         ++data64_1, ++data64_2)
        *data64_1 ^= *data64_2;

    element_t* data_1 = (element_t*)data64_1;
    element_t* data_2 = (element_t*)data64_2;

    for (const element_t* end_1 = (element_t*)a + symbol_size / sizeof(element_t); data_1 != end_1; ++data_1, ++data_2)
        *data_1 ^= *data_2;
}

void gf_mul(GF_t* gf, void* a, element_t coef, size_t symbol_size) {
    assert(symbol_size % sizeof(element_t) == 0);

    if (coef == 0) {
        memset(a, 0, symbol_size);
        return;
    }

    if (coef == 1)
        return;

    element_t* pow_table_shifted;
    element_t* data = (element_t*)a;
    uint16_t* log_table = gf->log_table;

    pow_table_shifted = gf->pow_table + log_table[coef];

    for (const element_t* end = data + symbol_size / sizeof(element_t); data != end; ++data) {
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

    pow_table_shifted = gf->pow_table + log_table[coef];

    for (const element_t* end_1 = data_1 + symbol_size / sizeof(element_t); data_1 != end_1; ++data_1, ++data_2) {
        element_t val_2 = *data_2;
        if (val_2 != 0)
            *data_1 ^= pow_table_shifted[log_table[val_2]];
    }
}
