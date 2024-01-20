/**
 * @file fft.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains functions for efficiently calculating Discrete Fourier
 * transform.
 * @date 2024-01-18
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_FFT_H__
#define __REED_SOLOMON_FFT_H__

#include <assert.h>
#include <stdint.h>

#include "constants.h"
#include "gf65536.h"
#include "symbol.h"

/**
 * @brief Galois field data and pre-computed values.
 */
typedef struct FFT {
    GF_t *gf;
} FFT_t;

/**
 * @brief Allocate memory.
 *
 * @param fft where to place allocated data.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int fft_alloc(FFT_t *fft) {
    // Nothing
    return 0;
}

/**
 * @brief Make necessary for Discrete Fourier transform implementation
 * pre-calculations.
 *
 * @param fft context object.
 * @param gf Galois field data.
 */
void fft_init(FFT_t *fft, GF_t *gf) { fft->gf = gf; }

/**
 * @brief Compute a given number of first components of Discrete Fourier
 * transform of a given sequence.
 * @details \f$\mathcal{F}_{r,\Theta}(f)\f$ - computes \f$F_0, \dots, F_{r-1}\f$
 * for any vector \f$f = (f_0, \dots, f_{N-1})\f$, such that \f$f_i \neq 0\f$
 * for \f$i \in \Theta\f$ and \f$f_i = 0\f$ otherwise.
 *
 * @param fft context object.
 * @param f sequence coefficients.
 * @param positions sequence coefficients indices.
 * @param res where to place the result.
 */
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

/**
 * @brief Compute some components of Discrete Fourier transform of a given
 * sequence.
 * @details \f$\tilde{\mathcal{F}}_{\Omega, d}(f)\f$ - computes \f$F_j =
 * f(a^{-j})\f$, \f$j \in \Omega\f$ for any polynomial \f$f(x)\f$ of degree
 * \f$d-1\f$.
 *
 * @param fft context object.
 * @param f sequence coefficients.
 * @param res where to place the result.
 * @param components components of Discrete Fourier transform to be computed.
 */
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

#endif