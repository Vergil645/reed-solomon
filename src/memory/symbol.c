/**
 * @file symbol.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief memory/symbol.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <memory/symbol.h>

symbol_t* symbol_create(size_t symbol_size) {
    symbol_t* s;

    s = (symbol_t*)malloc(sizeof(symbol_t));
    if (!s)
        return NULL;
    memset((void*)s, 0, sizeof(symbol_t));

    s->data = (uint8_t*)calloc(symbol_size, sizeof(uint8_t));
    if (!s->data) {
        free(s);
        return NULL;
    }

    return s;
}

void symbol_destroy(symbol_t* s) {
    assert(s != NULL);

    free(s->data);
    free(s);
}

bool symbol_eq(const symbol_t* a, const symbol_t* b, size_t symbol_size) {
    if (a == NULL || b == NULL || a->data == NULL || b->data == NULL)
        return false;

    uint8_t* data_1 = a->data;
    uint8_t* data_2 = b->data;

    for (uint8_t* end_1 = data_1 + symbol_size; data_1 != end_1; ++data_1, ++data_2) {
        if (*data_1 != *data_2)
            return false;
    }

    return true;
}

void symbol_printf(const symbol_t* s, size_t symbol_size) {
    if (s == NULL || s->data == NULL) {
        printf("NULL");
    } else if (symbol_size == 0) {
        printf("[]");
    } else {
        uint8_t* data = s->data;

        printf("[");
        for (uint8_t* end = data + symbol_size - 1; data != end; ++data)
            printf("%u, ", *data);
        printf("%u]", *data);
    }
}