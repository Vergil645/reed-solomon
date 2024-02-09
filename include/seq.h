/**
 * @file seq.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains symbol_seq_t definition and functions for interaction with
 * symbol sequences.
 * @date 2024-01-23
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_SEQ_H__
#define __REED_SOLOMON_SEQ_H__

#include <stdlib.h>

#include <symbol.h>

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
    size_t length;

    /**
     * @brief Sequence symbols.
     */
    symbol_t* symbols;
} symbol_seq_t;

/**
 * @brief Allocate sequence.
 *
 * @param symbol_size symbol size.
 * @param length sequence length.
 * @param seq sequence.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int seq_alloc(size_t symbol_size, size_t length, symbol_seq_t* seq);

/**
 * @brief Deallocate sequence.
 *
 * @param seq sequence.
 */
void seq_free(symbol_seq_t* seq);

#endif