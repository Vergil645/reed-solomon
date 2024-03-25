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

FFT_t* fft_create() {
    FFT_t* fft;

    fft = (FFT_t*)malloc(sizeof(FFT_t));
    if (!fft)
        return NULL;
    memset((void*)fft, 0, sizeof(FFT_t));

    uint8_t* precalc_scalar = fft->precalc_scalar;
    for (uint16_t a = 0; a < (1 << 8); ++a) {
        uint8_t* row = precalc_scalar + (a << 8);

        for (uint16_t b = a; b < (1 << 8); ++b) {
            for (uint8_t t = 0; t < 8; ++t)
                row[b] ^= ((uint8_t)(a & b) & (1 << t)) >> t;
            precalc_scalar[(b << 8) + a] = row[b];
        }
    }

    return fft;
}

void fft_destroy(FFT_t* fft) {
    assert(fft != NULL);

    free(fft);
}

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

static void _fast_matrix_mul(FFT_t* fft, uint16_t s, uint8_t m, const uint16_t* normal_repr,
                             symbol_t* tmp, const symbol_seq_t* f, const uint16_t* positions,
                             symbol_seq_t* u) {
    assert(fft != NULL);
    assert(normal_repr != NULL);
    assert(tmp != NULL);
    assert(f != NULL);
    assert(positions != NULL);
    assert(u != NULL);
    assert(f->symbol_size == u->symbol_size);

    size_t symbol_size = f->symbol_size;

    for (uint8_t t = 0; t < m; ++t)
        memset((void*)u->symbols[t]->data, 0, symbol_size);

    uint16_t i;

    for (i = 0; i + 8 <= f->length; i += 8) {
        const uint16_t repr[8] = {
            [0] = normal_repr[(s * positions[i + 0]) % N],
            [1] = normal_repr[(s * positions[i + 1]) % N],
            [2] = normal_repr[(s * positions[i + 2]) % N],
            [3] = normal_repr[(s * positions[i + 3]) % N],
            [4] = normal_repr[(s * positions[i + 4]) % N],
            [5] = normal_repr[(s * positions[i + 5]) % N],
            [6] = normal_repr[(s * positions[i + 6]) % N],
            [7] = normal_repr[(s * positions[i + 7]) % N],
        };

        uint8_t A[CC_MAX_COSET_SIZE] = {0};
        for (uint8_t t = 0; t < m; ++t) {
            for (uint8_t j = 0; j < 8; ++j) // TODO: disloop?
                A[t] |= ((repr[j] & (1 << t)) >> t) << j;
        }

        uint8_t* B = tmp->data; // length == (symbol_size << 3)
        memset((void*)B, 0, symbol_size << 3);
        for (uint8_t j = 0; j < 8; ++j) {
            uint8_t* data = f->symbols[i + j]->data;
            size_t idx = 0;

            for (const uint8_t* end = data + symbol_size; data != end;
                 ++data, idx += 8) { // TODO: improve?
                uint8_t val = *data;

                B[idx + 0] |= ((val & (1 << 0)) >> 0) << j;
                B[idx + 1] |= ((val & (1 << 1)) >> 1) << j;
                B[idx + 2] |= ((val & (1 << 2)) >> 2) << j;
                B[idx + 3] |= ((val & (1 << 3)) >> 3) << j;
                B[idx + 4] |= ((val & (1 << 4)) >> 4) << j;
                B[idx + 5] |= ((val & (1 << 5)) >> 5) << j;
                B[idx + 6] |= ((val & (1 << 6)) >> 6) << j;
                B[idx + 7] |= ((val & (1 << 7)) >> 7) << j;
            }
        }

        for (uint8_t t = 0; t < m; ++t) {
            uint8_t* data = u->symbols[t]->data;
            uint8_t* row_a = fft->precalc_scalar + (A[t] << 8);
            size_t idx = 0;

            for (const uint8_t* end = data + symbol_size; data != end;
                 ++data, idx += 8) { // TODO: improve?
                *data ^= row_a[B[idx + 0]] << 0 | row_a[B[idx + 1]] << 1 | row_a[B[idx + 2]] << 2 |
                         row_a[B[idx + 3]] << 3 | row_a[B[idx + 4]] << 4 | row_a[B[idx + 5]] << 5 |
                         row_a[B[idx + 6]] << 6 | row_a[B[idx + 7]] << 7;
            }
        }
    }

    for (; i < f->length; ++i) {
        uint16_t repr = normal_repr[(s * positions[i]) % N];

        for (uint8_t t = 0; t < m; ++t) {
            if (repr & (1 << t))
                gf_add((void*)u->symbols[t]->data, (void*)f->symbols[i]->data, symbol_size);
        }
    }

    // for (uint16_t i = 0; i < f->length; ++i) {
    //     uint16_t repr = normal_repr[(s * positions[i]) % N];

    //     for (uint8_t t = 0; t < m; ++t) {
    //         if (repr & (1 << t))
    //             gf_add((void*)u->symbols[t]->data, (void*)f->symbols[i]->data, symbol_size);
    //     }
    // }
}

int fft_transform_cycl(FFT_t* fft, GF_t* gf, const symbol_seq_t* f, const uint16_t* positions,
                       symbol_seq_t* res) {
    assert(gf != NULL);
    assert(f != NULL);
    assert(positions != NULL);
    assert(res != NULL);
    assert(f->symbol_size == res->symbol_size);

    symbol_t* tmp;
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

    tmp = symbol_create(symbol_size << 3);
    if (!tmp) {
        seq_destroy(u);
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

        // for (uint8_t t = 0; t < m; ++t)
        //     memset((void*)u->symbols[t]->data, 0, symbol_size);

        // for (uint16_t i = 0; i < f->length; ++i) {
        //     uint16_t repr = normal_repr[(s * positions[i]) % N];

        //     for (uint8_t t = 0; t < m; ++t) {
        //         if (repr & (1 << t))
        //             gf_add((void*)u->symbols[t]->data, (void*)f->symbols[i]->data, symbol_size);
        //     }
        // }

        _fast_matrix_mul(fft, s, m, normal_repr, tmp, f, positions, u);

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

    symbol_destroy(tmp);
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