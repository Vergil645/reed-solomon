/**
 * @file util.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief util.h implementation.
 * @date 2024-01-23
 *
 * @copyright Copyright (c) 2024
 */

#include <stdio.h>

#include <util.h>

int util_u16_array_eq(const uint16_t* a, size_t a_len, const uint16_t* b,
                      size_t b_len) {
    if (a == NULL || b == NULL || a_len != b_len) {
        return 0;
    }

    for (uint16_t i = 0; i < a_len; ++i) {
        if (a[i] != b[i]) {
            return 0;
        }
    }

    return 1;
}

void util_u16_array_printf(const uint16_t* a, size_t a_len) {
    if (a == NULL) {
        printf("NULL");
    } else if (a_len == 0) {
        printf("[]");
    } else {
        printf("[");
        for (size_t i = 0; i < a_len - 1; ++i) {
            printf("%u, ", a[i]);
        }
        printf("%u]", a[a_len - 1]);
    }
}

int util_symbol_eq(size_t symbol_size, symbol_t a, symbol_t b) {
    return util_u16_array_eq(a.data, symbol_size, b.data, symbol_size);
}

void util_symbol_printf(size_t symbol_size, symbol_t s) {
    util_u16_array_printf(s.data, symbol_size);
}

int util_seq_eq(symbol_seq_t a, symbol_seq_t b) {
    size_t symbol_size;

    if (a.symbol_size != b.symbol_size || a.length != b.length ||
        a.symbols == NULL || b.symbols == NULL) {
        return 0;
    }

    symbol_size = a.symbol_size;

    for (size_t i = 0; i < a.length; ++i) {
        if (!util_symbol_eq(symbol_size, a.symbols[i], b.symbols[i])) {
            return 0;
        }
    }

    return 1;
}

void util_seq_printf(symbol_seq_t seq) {
    if (seq.symbols == NULL) {
        printf("NULL");
    } else if (seq.length == 0) {
        printf("[]");
    } else {
        printf("[");
        for (size_t i = 0; i < seq.length - 1; ++i) {
            util_symbol_printf(seq.symbol_size, seq.symbols[i]);
            printf(", ");
        }
        util_symbol_printf(seq.symbol_size, seq.symbols[seq.length - 1]);
        printf("]");
    }
}