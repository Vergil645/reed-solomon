/**
 * @file seq.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief seq.h implementation.
 * @date 2024-01-23
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>

#include <seq.h>

int seq_alloc(size_t symbol_size, size_t length, symbol_seq_t *seq) {
    assert(seq != NULL);

    symbol_t *symbols;
    int ret;

    symbols = (symbol_t *)malloc(length * sizeof(symbol_t));
    if (!symbols) {
        return 1;
    }

    for (size_t i = 0; i < length; ++i) {
        ret = symbol_alloc(symbol_size, symbols + i);
        if (ret) {
            for (size_t j = 0; j < i; ++j) {
                symbol_free(symbols + j);
            }
            free(symbols);

            return ret;
        }
    }

    seq->symbol_size = symbol_size;
    seq->length = length;
    seq->symbols = symbols;

    return 0;
}

void seq_free(symbol_seq_t *seq) {
    size_t length = seq->length;
    symbol_t *symbols = seq->symbols;

    for (size_t i = 0; i < length; ++i) {
        symbol_free(symbols + i);
    }
    free(symbols);
}