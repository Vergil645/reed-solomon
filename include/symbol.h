/**
 * @file symbol.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains symbol_t definition and functions for interaction with
 * symbols.
 * @date 2024-01-20
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_SYMBOL_H__
#define __REED_SOLOMON_SYMBOL_H__

#include <stdlib.h>

#include <gf65536.h>

/**
 * @brief Symbol data type. Can be passed to functions by value.
 */
typedef struct symbol {
    /**
     * @brief Symbol data.
     */
    element_t *data;
} symbol_t;

/**
 * @brief Symbol sequence data type. Can be passed to functions by value.
 */
typedef struct symbol_seq {
    /**
     * @brief Symbol size (similar for each symbol in sequence).
     */
    size_t symbol_size;

    /**
     * @brief Sequence length.
     */
    size_t seq_length;

    /**
     * @brief Sequence symbols.
     */
    symbol_t *symbols;
} symbol_seq_t;

/**
 * @brief Allocate symbol.
 *
 * @param symbol_size symbol size.
 * @param s symbol.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int alloc_symbol(size_t symbol_size, symbol_t *s);

/**
 * @brief Deallocate symbol.
 *
 * @param s symbol.
 */
void free_symbol(symbol_t *s);

/**
 * @brief Allocate sequence.
 *
 * @param symbol_size symbol size.
 * @param seq_length sequence length.
 * @param seq sequence.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int alloc_seq(size_t symbol_size, size_t seq_length, symbol_seq_t *seq);

/**
 * @brief Deallocate sequence.
 *
 * @param seq sequence.
 */
void free_seq(symbol_seq_t *seq);

#endif