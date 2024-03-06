/**
 * @file seq.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief memory/seq.h implementation.
 * @date 2024-01-23
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <memory/seq.h>

symbol_seq_t* seq_create(size_t length, size_t symbol_size) {
    symbol_seq_t* seq;

    seq = (symbol_seq_t*)malloc(sizeof(symbol_seq_t));
    if (!seq)
        return NULL;
    memset((void*)seq, 0, sizeof(symbol_seq_t));

    seq->length = length;
    seq->symbol_size = symbol_size;

    seq->symbols = (symbol_t**)calloc(length, sizeof(symbol_t*));
    if (!seq->symbols) {
        free(seq);
        return NULL;
    }

    for (size_t i = 0; i < length; ++i) {
        seq->symbols[i] = symbol_create(symbol_size);
        if (!seq->symbols[i]) {
            for (size_t j = 0; j < i; ++j)
                symbol_destroy(seq->symbols[j]);
            free(seq->symbols);
            free(seq);
            return NULL;
        }
    }

    return seq;
}

void seq_destroy(symbol_seq_t* seq) {
    assert(seq != NULL);

    for (size_t i = 0; i < seq->length; ++i)
        symbol_destroy(seq->symbols[i]);
    free(seq->symbols);
    free(seq);
}

bool seq_eq(const symbol_seq_t* a, const symbol_seq_t* b) {
    if (a == NULL || b == NULL || a->symbols == NULL || b->symbols == NULL)
        return false;

    if (a->length != b->length || a->symbol_size != b->symbol_size)
        return false;

    size_t length = a->length;
    size_t symbol_size = a->symbol_size;

    for (size_t i = 0; i < length; ++i) {
        if (!symbol_eq(a->symbols[i], b->symbols[i], symbol_size))
            return false;
    }

    return true;
}

void seq_printf(const symbol_seq_t* seq) {
    if (seq == NULL || seq->symbols == NULL) {
        printf("NULL");
    } else if (seq->length == 0) {
        printf("[]");
    } else {
        size_t length = seq->length;
        size_t symbol_size = seq->symbol_size;
        symbol_t** symbols = seq->symbols;

        printf("[");
        for (symbol_t** end = symbols + length - 1; symbols != end; ++symbols) {
            symbol_printf(*symbols, symbol_size);
            printf(", ");
        }
        symbol_printf(*symbols, symbol_size);
        printf("]");
    }
}