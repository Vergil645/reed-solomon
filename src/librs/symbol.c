/**
 * @file symbol.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief symbol.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>

#include "symbol.h"

int alloc_symbol(size_t symbol_size, symbol_t *s) {
    assert(s != NULL);

    element_t *data;

    data = (element_t *)malloc(symbol_size * sizeof(element_t));
    if (!data) {
        return 1;
    }

    s->data = data;

    return 0;
}

void free_symbol(symbol_t *s) {
    assert(s != NULL);

    free(s->data);
}

int alloc_seq(size_t symbol_size, size_t seq_length, symbol_seq_t *seq) {
    assert(seq != NULL);

    symbol_t *symbols;
    int ret;

    symbols = (symbol_t *)malloc(seq_length * sizeof(symbol_t));
    if (!symbols) {
        return 1;
    }

    for (size_t i = 0; i < seq_length; ++i) {
        ret = alloc_symbol(symbol_size, symbols + i);
        if (ret) {
            for (size_t j = 0; j < i; ++j) {
                free_symbol(symbols + j);
            }
            free(symbols);

            return ret;
        }
    }

    seq->symbol_size = symbol_size;
    seq->seq_length = seq_length;
    seq->symbols = symbols;

    return 0;
}

void free_seq(symbol_seq_t *seq) {
    size_t seq_length = seq->seq_length;
    symbol_t *symbols = seq->symbols;

    for (size_t i = 0; i < seq_length; ++i) {
        free_symbol(symbols + i);
    }
    free(symbols);
}