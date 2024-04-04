/**
 * @file fft.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief rs/fft.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <rs/fft.h>

// cppcheck-suppress unusedFunction
void fft_transform(GF_t* gf, const symbol_seq_t* f, const uint16_t* positions, symbol_seq_t* res) {
    assert(gf != NULL);
    assert(f != NULL);
    assert(positions != NULL);
    assert(res != NULL);
    assert(f->symbol_size == res->symbol_size);

    size_t symbol_size = f->symbol_size;
    element_t* pow_table = gf->pow_table;
    element_t coef;

    for (uint16_t j = 0; j < res->length; ++j) {
        memset((void*)res->symbols[j]->data, 0, symbol_size);

        for (uint16_t i = 0; i < f->length; ++i) {
            coef = pow_table[(positions[i] * j) % N];
            gf_madd(gf, (void*)res->symbols[j]->data, coef, (void*)f->symbols[i]->data,
                    symbol_size);
        }
    }
}

int fft_transform_cycl(GF_t* gf, const symbol_seq_t* f, const uint16_t* positions,
                       symbol_seq_t* res) {
    assert(gf != NULL);
    assert(f != NULL);
    assert(positions != NULL);
    assert(res != NULL);
    assert(f->symbol_size == res->symbol_size);

    symbol_seq_t* u;
    size_t symbol_size = f->symbol_size;

    bool* calculated = (bool*)calloc(res->length, sizeof(bool));
    if (!calculated)
        return 1;

    u = seq_create(CC_MAX_COSET_SIZE, symbol_size);
    if (!u) {
        free(calculated);
        return 1;
    }

    for (uint16_t s = 0; s < res->length; ++s) {
        if (calculated[s])
            continue;

        uint8_t m = cc_get_coset_size(s);

        for (uint8_t t = 0; t < m; ++t)
            memset((void*)u->symbols[t]->data, 0, symbol_size);

        for (uint16_t i = 0; i < f->length; ++i) {
            uint16_t repr = gf_get_normal_repr(gf, m, (s * positions[i]) % N);

            for (uint8_t t = 0; t < m; ++t) {
                if (repr & (1 << t))
                    gf_add((void*)u->symbols[t]->data, (void*)f->symbols[i]->data, symbol_size);
            }
        }

        uint16_t idx = s;
        for (uint8_t j = 0; j < m; ++j) {
            if (idx < res->length) {
                memset((void*)res->symbols[idx]->data, 0, symbol_size);

                for (uint8_t t = 0; t < m; ++t) {
                    element_t coef = gf_get_normal_basis_element(gf, m, (j + t) % m);
                    gf_madd(gf, (void*)res->symbols[idx]->data, coef, (void*)u->symbols[t]->data,
                            symbol_size);
                }

                calculated[idx] = true;
            }

            idx = NEXT_COSET_ELEMENT(idx);
        }

        assert(idx == s);
    }

    seq_destroy(u);
    free(calculated);

    return 0;
}

// cppcheck-suppress unusedFunction
void fft_partial_transform(GF_t* gf, const symbol_seq_t* f, const uint16_t* components,
                           symbol_seq_t* res) {
    assert(gf != NULL);
    assert(f != NULL);
    assert(components != NULL);
    assert(res != NULL);
    assert(f->symbol_size == res->symbol_size);

    size_t symbol_size = f->symbol_size;
    element_t* pow_table = gf->pow_table;
    element_t coef;

    for (uint16_t res_idx = 0; res_idx < res->length; ++res_idx) {
        uint16_t j = (N - components[res_idx]) % N;

        memset((void*)res->symbols[res_idx]->data, 0, symbol_size);

        for (uint16_t i = 0; i < f->length; ++i) {
            coef = pow_table[(i * j) % N];
            gf_madd(gf, (void*)res->symbols[res_idx]->data, coef, (void*)f->symbols[i]->data,
                    symbol_size);
        }
    }
}

int fft_partial_transform_cycl(GF_t* gf, const symbol_seq_t* f, const coset_t* cosets,
                               uint16_t cosets_cnt, symbol_seq_t* res) {
    assert(gf != NULL);
    assert(f != NULL);
    assert(cosets != NULL);
    assert(res != NULL);
    assert(f->symbol_size == res->symbol_size);

    symbol_seq_t* u;
    size_t symbol_size = f->symbol_size;
    uint16_t idx = 0;

    u = seq_create(CC_MAX_COSET_SIZE, symbol_size);
    if (!u)
        return 1;

    for (const coset_t* end = cosets + cosets_cnt; cosets != end; ++cosets) {
        coset_t coset = *cosets;

        uint16_t s = N - coset.leader;
        uint8_t m = coset.size;

        for (uint8_t t = 0; t < m; ++t)
            memset((void*)u->symbols[t]->data, 0, symbol_size);

        for (uint16_t i = 0; i < f->length; ++i) {
            uint16_t repr = gf_get_normal_repr(gf, m, (s * i) % N);

            for (uint8_t t = 0; t < m; ++t) {
                if (repr & (1 << t))
                    gf_add((void*)u->symbols[t]->data, f->symbols[i]->data, symbol_size);
            }
        }

        for (uint8_t j = 0; j < m; ++j, ++idx) {
            assert(idx < res->length);

            memset((void*)res->symbols[idx]->data, 0, symbol_size);

            for (uint8_t t = 0; t < m; ++t) {
                element_t coef = gf_get_normal_basis_element(gf, m, (j + t) % m);
                gf_madd(gf, (void*)res->symbols[idx]->data, coef, (void*)u->symbols[t]->data,
                        symbol_size);
            }
        }
    }

    assert(idx == res->length);

    seq_destroy(u);

    return 0;
}