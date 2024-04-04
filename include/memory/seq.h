/**
 * @file seq.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains symbol_seq_t definition and functions for interaction with symbol sequences.
 * @date 2024-01-23
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __MEMORY_SEQ_H__
#define __MEMORY_SEQ_H__

#include <stdbool.h>
#include <stddef.h>

#include "symbol.h"

/**
 * @brief Symbol sequence data type.
 */
typedef struct {
    /**
     * @brief Sequence length.
     */
    size_t length;

    /**
     * @brief Symbol size (similar for each symbol in sequence).
     */
    size_t symbol_size;

    /**
     * @brief Sequence symbols.
     */
    symbol_t** symbols;
} symbol_seq_t;

/**
 * @brief Create sequence.
 *
 * @param symbol_size symbol size.
 * @param length sequence length.
 * @return pointer to created sequence on success and NULL otherwise.
 */
symbol_seq_t* seq_create(size_t length, size_t symbol_size);

/**
 * @brief Destroy sequence.
 *
 * @param seq sequence.
 */
void seq_destroy(symbol_seq_t* seq);

/**
 * @brief Check whether sequences are equal.
 *
 * @param a first sequence.
 * @param b second sequence.
 * @return true if sequences are equal, false otherwise.
 */
bool seq_eq(const symbol_seq_t* a, const symbol_seq_t* b);

/**
 * @brief Write sequence to stdout.
 *
 * @param seq sequence.
 */
void seq_printf(const symbol_seq_t* seq);

#endif