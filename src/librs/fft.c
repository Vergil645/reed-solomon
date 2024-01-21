/**
 * @file fft.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief fft.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>

#include "fft.h"

int fft_alloc(FFT_t *fft) {
    // Nothing
    return 0;
}

void fft_init(FFT_t *fft, GF_t *gf) { fft->gf = gf; }

void fft_transform(const FFT_t *fft, symbol_seq_t f, const uint16_t *positions,
                   symbol_seq_t res) {
    assert(positions != NULL);
    assert(f.symbol_size == res.symbol_size);

    GF_t *gf = fft->gf;
    element_t *pow_table = gf->pow_table;
    size_t symbol_size = f.symbol_size;
    element_t x;

    for (uint16_t j = 0; j < res.seq_length; ++j) {
        // Initialization
        for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
            res.symbols[j].data[e_idx] = 0;
        }

        for (uint16_t f_idx = 0; f_idx < f.seq_length; ++f_idx) {
            x = pow_table[(positions[f_idx] * j) % N];

            for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                res.symbols[j].data[e_idx] ^=
                    gf_mul_ee(gf, f.symbols[f_idx].data[e_idx], x);
            }
        }
    }
}

void fft_selective_transform(const FFT_t *fft, symbol_seq_t f, symbol_seq_t res,
                             const uint16_t *components) {
    assert(components != NULL);
    assert(f.symbol_size == res.symbol_size);

    GF_t *gf = fft->gf;
    element_t *pow_table = gf->pow_table;
    size_t symbol_size = f.symbol_size;
    uint16_t j;
    element_t x;

    for (uint16_t res_idx = 0; res_idx < res.seq_length; ++res_idx) {
        // Initialization
        for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
            res.symbols[res_idx].data[e_idx] = 0;
        }

        j = (N - components[res_idx]) % N;

        for (uint16_t i = 0; i < f.seq_length; ++i) {
            x = pow_table[(i * j) % N];

            for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                res.symbols[res_idx].data[e_idx] ^=
                    gf_mul_ee(gf, f.symbols[i].data[e_idx], x);
            }
        }
    }
}