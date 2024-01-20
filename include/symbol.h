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

#include "gf65536.h"

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
 * @brief Allocate symbol memory.
 *
 * @param s symbol.
 * @param symbol_size symbol size.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int alloc_symbol(symbol_t *s, size_t symbol_size) {
    assert(s != NULL);

    s->data = (element_t *)malloc(symbol_size * sizeof(element_t));
    return s->data == NULL ? 1 : 0;
}

/**
 * @brief Deallocate symbol memory.
 *
 * @param s symbol.
 */
void free_symbol(symbol_t *s) {
    assert(s != NULL);

    free(s->data);
}

/**
 * @brief Create a symbol_seq_t object.
 *
 * @param seq sequence.
 * @param symbol_size symbol size.
 * @param seq_length sequence length.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int create_seq(symbol_seq_t *seq, size_t symbol_size, size_t seq_length) {
    assert(seq != NULL);

    symbol_t *symbols;
    int ret;

    symbols = (symbol_t *)malloc(seq_length * sizeof(symbol_t));
    if (!symbols) {
        return 1;
    }

    for (size_t i = 0; i < seq_length; ++i) {
        ret = alloc_symbol(symbols + i, symbol_size);
        if (ret) {
            for (size_t j = 0; j < i; ++j) {
                free_symbol(symbols + i);
            }
            free(symbols);
            return ret;
        }
    }

    *seq = {
        .symbol_size = symbol_size,
        .seq_length = seq_length,
        .symbols = symbols,
    };
}

/**
 * @brief Deallocate sequence memory.
 *
 * @param seq sequence.
 */
void free_seq(symbol_seq_t *seq) {
    size_t seq_length = seq->seq_length;
    symbol_t *symbols = seq->symbols;

    for (size_t i = 0; i < seq_length; ++i) {
        free_symbol(symbols + i);
    }
    free(symbols);
}

#endif