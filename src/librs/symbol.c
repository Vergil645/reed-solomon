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

int symbol_alloc(size_t symbol_size, symbol_t *s) {
    assert(s != NULL);

    element_t *data;

    data = (element_t *)malloc(symbol_size * sizeof(element_t));
    if (!data) {
        return 1;
    }

    s->data = data;

    return 0;
}

void symbol_free(symbol_t *s) {
    assert(s != NULL);

    free(s->data);
}