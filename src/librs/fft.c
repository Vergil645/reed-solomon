/**
 * @file fft.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief fft.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <string.h>

#include <fft.h>

int fft_alloc(FFT_t *fft) {
    // Nothing
    return 0;
}

void fft_init(FFT_t *fft, GF_t *gf) { fft->gf = gf; }

void fft_free(FFT_t *fft) {
    // Nothing
}

void fft_transform(const FFT_t *fft, symbol_seq_t f, const uint16_t *positions,
                   symbol_seq_t res) {
    assert(positions != NULL);
    assert(f.symbol_size == res.symbol_size);

    GF_t *gf = fft->gf;
    element_t *pow_table = gf->pow_table;
    size_t symbol_size = f.symbol_size;
    element_t c;

    for (uint16_t j = 0; j < res.length; ++j) {
        // Initialization
        memset((void *)res.symbols[j].data, 0, symbol_size * sizeof(element_t));

        for (uint16_t i = 0; i < f.length; ++i) {
            c = pow_table[(positions[i] * j) % N];

            for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                res.symbols[j].data[e_idx] ^=
                    gf_mul_ee(gf, f.symbols[i].data[e_idx], c);
            }
        }
    }
}

void fft_partial_transform(const FFT_t *fft, symbol_seq_t f, symbol_seq_t res,
                           const uint16_t *components) {
    assert(components != NULL);
    assert(f.symbol_size == res.symbol_size);

    GF_t *gf = fft->gf;
    element_t *pow_table = gf->pow_table;
    size_t symbol_size = f.symbol_size;
    uint16_t j;
    element_t x;

    for (uint16_t res_idx = 0; res_idx < res.length; ++res_idx) {
        // Initialization
        memset((void *)res.symbols[res_idx].data, 0,
               symbol_size * sizeof(element_t));

        j = (N - components[res_idx]) % N;

        for (uint16_t i = 0; i < f.length; ++i) {
            x = pow_table[(i * j) % N];

            for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                res.symbols[res_idx].data[e_idx] ^=
                    gf_mul_ee(gf, f.symbols[i].data[e_idx], x);
            }
        }
    }
}