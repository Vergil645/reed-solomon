/**
 * @file rlc_gf256.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains implementation of a Galois field of size 256.
 * @date 2024-02-27
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __RLC_GF256_H__
#define __RLC_GF256_H__

#include <stdint.h>

void gf256_init();

void assign_mul(uint8_t** array);

void assign_inv(uint8_t* array);

void gf256_symbol_add_scaled(void* symbol1, uint8_t coef, const void* symbol2, uint32_t symbol_size,
                             uint8_t** mul);

void gf256_symbol_mul(void* symbol1, uint8_t coef, uint32_t symbol_size, uint8_t** mul);

void gf256_symbol_add(void* symbol1, const void* symbol2, uint32_t symbol_size);

#endif