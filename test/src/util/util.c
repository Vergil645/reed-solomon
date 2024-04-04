#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <test/util/util.h>

bool util_u16_array_eq(const uint16_t* a, size_t a_len, const uint16_t* b, size_t b_len) {
    if (a == NULL || b == NULL || a_len != b_len)
        return false;

    size_t length = a_len;
    uint16_t* data_1 = (uint16_t*)a;
    uint16_t* data_2 = (uint16_t*)b;

    for (const uint16_t* end_1 = data_1 + length; data_1 != end_1; ++data_1, ++data_2) {
        if (*data_1 != *data_2)
            return false;
    }

    return true;
}

void util_u16_array_printf(const uint16_t* a, size_t a_len) {
    if (a == NULL) {
        printf("NULL");
    } else if (a_len == 0) {
        printf("[]");
    } else {
        printf("[");
        for (const uint16_t* end = a + a_len - 1; a != end; ++a)
            printf("%u, ", *a);
        printf("%u]", *a);
    }
}

void util_generate_inf_symbols(const symbol_seq_t* inf_symbols) {
    size_t symbol_size = inf_symbols->symbol_size;
    size_t length = inf_symbols->length;

    for (size_t i = 0; i < length; ++i) {
        for (size_t j = 0; j < symbol_size; ++j)
            inf_symbols->symbols[i]->data[j] = (uint8_t)rand();
    }
}

void util_init_rcv_symbols(const symbol_seq_t* src_symbols, symbol_seq_t* rcv_symbols) {
    assert(src_symbols->symbol_size == rcv_symbols->symbol_size);
    assert(src_symbols->length == rcv_symbols->length);

    size_t symbol_size = src_symbols->symbol_size;
    size_t length = src_symbols->length;

    for (size_t i = 0; i < length; ++i) {
        memcpy((void*)rcv_symbols->symbols[i]->data, (void*)src_symbols->symbols[i]->data, symbol_size);
    }
}

void util_choose_and_erase_symbols(symbol_seq_t* rcv_symbols, uint16_t t, bool* is_erased) {
    size_t symbol_size = rcv_symbols->symbol_size;
    uint16_t n = (uint16_t)rcv_symbols->length;

    while (t > 0) {
        for (uint16_t i = 0; i < n && t > 0; ++i) {
            if (is_erased[i])
                continue;
            if (rand() % 2 != 0) {
                is_erased[i] = true;
                --t;
            }
        }
    }

    for (uint16_t i = 0; i < n; ++i) {
        if (!is_erased[i])
            continue;
        memset((void*)rcv_symbols->symbols[i]->data, 0, symbol_size);
    }
}