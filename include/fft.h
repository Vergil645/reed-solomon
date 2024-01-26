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

#include <stdint.h>

#include <gf65536.h>
#include <seq.h>
#include <symbol.h>

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
int fft_alloc(FFT_t *fft);

/**
 * @brief Make necessary for Discrete Fourier transform implementation
 * pre-calculations.
 *
 * @param fft context object.
 * @param gf Galois field data.
 */
void fft_init(FFT_t *fft, GF_t *gf);

/**
 * @brief Deallocate context object.
 *
 * @param fft context object.
 */
void fft_free(FFT_t *fft);

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
                   symbol_seq_t res);

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
 * @param components negative components of the discrete Fourier transform to be
 * calculated.
 */
void fft_partial_transform(const FFT_t *fft, symbol_seq_t f, symbol_seq_t res,
                           const uint16_t *components);

#endif