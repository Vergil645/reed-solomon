/**
 * @file util.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Utility functions for tests implementation.
 * @date 2024-01-23
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __TEST_UTIL_H__
#define __TEST_UTIL_H__

#include <stdint.h>
#include <stdlib.h>

#include <seq.h>
#include <symbol.h>

/**
 * @brief Check that arrays contains same elements and have same length.
 *
 * @param a first array.
 * @param a_len length of the first array.
 * @param b second array.
 * @param b_len length of the second array.
 * @return 0 if arrays are not equal or one of them is NULL,\n
 *         1 otherwise.
 */
int util_u16_array_eq(const uint16_t *a, size_t a_len, const uint16_t *b,
                      size_t b_len);

/**
 * @brief Print array to stdin.
 *
 * @param a array.
 * @param a_len array length.
 */
void util_u16_array_printf(const uint16_t *a, size_t a_len);

/**
 * @brief Check that two symbols are equal
 *
 * @param symbol_size symbol size.
 * @param a first symbol.
 * @param b second symbol.
 * @return 0 if symbols are not equal or one of them has .data == NULL,\n
 *         1 otherwise.
 */
int util_symbol_eq(size_t symbol_size, symbol_t a, symbol_t b);

/**
 * @brief Print symbol to stdin.
 *
 * @param symbol_size symbol size.
 * @param s symbol.
 */
void util_symbol_printf(size_t symbol_size, symbol_t s);

/**
 * @brief Check that two symbol sequences are equal.
 *
 * @param a first symbol sequence.
 * @param b second symbol sequence.
 * @return 0 if sequences are not equal or one of them has .symbols == NULL,\n
 *         1 otherwise.
 */
int util_seq_eq(symbol_seq_t a, symbol_seq_t b);

/**
 * @brief Print symbol sequence to stdin.
 *
 * @param seq symbol sequence.
 */
void util_seq_printf(symbol_seq_t seq);

#endif