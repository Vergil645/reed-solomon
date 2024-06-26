/**
 * @file fft.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains functions for efficiently calculating Discrete Fourier transform.
 * @date 2024-01-18
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_FFT_H__
#define __REED_SOLOMON_FFT_H__

#include <stdint.h>

#include "gf65536.h"
#include <memory/seq.h>
#include <memory/symbol.h>

/**
 * @brief Compute a given number of first components of Discrete Fourier transform of a given sequence.
 * @details \f$\mathcal{F}_{r,\Theta}(f)\f$ - computes \f$F_0, \dots, F_{r-1}\f$ for any vector
 * \f$f = (f_0, \dots, f_{N-1})\f$, such that \f$f_i \neq 0\f$ for \f$i \in \Theta\f$ and \f$f_i = 0\f$ otherwise.
 *
 * @param gf Galois field data.
 * @param f sequence coefficients.
 * @param positions sequence coefficients indices.
 * @param res where to place the result.
 */
void fft_transform(GF_t* gf, const symbol_seq_t* f, const uint16_t* positions, symbol_seq_t* res);

/**
 * @brief Compute a given number of first components of Discrete Fourier transform of a given sequence using cyclotomic
 * FFT algorithm.
 *
 * @param gf Galois field data.
 * @param f sequence coefficients.
 * @param positions sequence coefficients indices.
 * @param res where to place the result.
 * @return 0 on success, 1 on memory allocation error.
 */
int fft_transform_cycl(GF_t* gf, const symbol_seq_t* f, const uint16_t* positions, symbol_seq_t* res);

/**
 * @brief Compute some components of Discrete Fourier transform of a given sequence.
 * @details \f$\tilde{\mathcal{F}}_{\Omega, d}(f)\f$ - computes \f$F_j = f(a^{-j})\f$, \f$j \in \Omega\f$ for any
 * polynomial \f$f(x)\f$ of degree \f$d-1\f$.
 *
 * @param gf Galois field data.
 * @param f sequence coefficients.
 * @param components negative components of the discrete Fourier transform to be computed.
 * @param res where to place the result.
 */
void fft_partial_transform(GF_t* gf, const symbol_seq_t* f, const uint16_t* components, symbol_seq_t* res);

/**
 * @brief Compute some components of Discrete Fourier transform of a given sequence using cyclotomic FFT algorithm.
 *
 * @param gf Galois field data.
 * @param f sequence coefficients.
 * @param cosets cyclotomic cosets that forms negative components to be computed.
 * @param cosets_cnt number of cyclotomic cosets.
 * @param res where to place the result.
 * @return 0 on success, 1 on memory allocation error.
 */
int fft_partial_transform_cycl(GF_t* gf, const symbol_seq_t* f, const coset_t* cosets, uint16_t cosets_cnt,
                               symbol_seq_t* res);

#endif