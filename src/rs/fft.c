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
void __attribute_maybe_unused__ fft_transform(GF_t* gf, const symbol_seq_t* f,
                                              const uint16_t* positions, symbol_seq_t* res) {
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
    element_t coef;
    uint16_t* normal_repr;
    uint16_t idx;

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

        uint8_t m = 1;
        uint8_t m_log = 0;
        while (s != (uint16_t)(((uint32_t)s << m) % N)) {
            m <<= 1;
            ++m_log;
        }
        assert(m <= CC_MAX_COSET_SIZE);

        normal_repr = gf->normal_repr_by_subfield[m_log];

        for (uint8_t t = 0; t < m; ++t)
            memset((void*)u->symbols[t]->data, 0, symbol_size);

        for (uint16_t i = 0; i < f->length; ++i) {
            uint16_t repr = normal_repr[(s * positions[i]) % N];

            for (uint8_t t = 0; t < m; ++t) {
                if (repr & (1 << t))
                    gf_add((void*)u->symbols[t]->data, (void*)f->symbols[i]->data, symbol_size);
            }
        }

        idx = s;
        for (uint8_t j = 0; j < m; ++j) {
            if (idx < res->length) {
                memset((void*)res->symbols[idx]->data, 0, symbol_size);

                for (uint8_t t = 0; t < m; ++t) {
                    coef = g_normal_basis_by_subfield[m_log][(j + t) % m];
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
void __attribute_maybe_unused__ fft_partial_transform(GF_t* gf, const symbol_seq_t* f,
                                                      const uint16_t* components,
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
    uint16_t* normal_repr;
    uint16_t idx = 0;
    element_t coef;

    u = seq_create(CC_MAX_COSET_SIZE, symbol_size);
    if (!u)
        return 1;

    for (const coset_t* end = cosets + cosets_cnt; cosets != end; ++cosets) {
        coset_t coset = *cosets;

        uint16_t s = N - coset.leader;
        uint8_t m = coset.size;
        uint8_t m_log = 0;
        while (m != (1 << m_log))
            ++m_log;
        assert(m_log <= 4);

        normal_repr = gf->normal_repr_by_subfield[m_log];

        for (uint8_t t = 0; t < m; ++t)
            memset((void*)u->symbols[t]->data, 0, symbol_size);

        for (uint16_t i = 0; i < f->length; ++i) {
            uint16_t repr = normal_repr[(s * i) % N];

            for (uint8_t t = 0; t < m; ++t) {
                if (repr & (1 << t))
                    gf_add((void*)u->symbols[t]->data, f->symbols[i]->data, symbol_size);
            }
        }

        for (uint8_t j = 0; j < m; ++j, ++idx) {
            assert(idx < res->length);

            memset((void*)res->symbols[idx]->data, 0, symbol_size);

            for (uint8_t t = 0; t < m; ++t) {
                coef = g_normal_basis_by_subfield[m_log][(j + t) % m];
                gf_madd(gf, (void*)res->symbols[idx]->data, coef, (void*)u->symbols[t]->data,
                        symbol_size);
            }
        }
    }

    assert(idx == res->length);

    seq_destroy(u);

    return 0;
}