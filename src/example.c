#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rs/reed_solomon.h>

#define SEED 78934
#define TESTS_CNT 100

static void generate_inf_symbols(const symbol_seq_t* inf_symbols) {
    size_t symbol_size = inf_symbols->symbol_size;
    size_t length = inf_symbols->length;

    for (size_t i = 0; i < length; ++i) {
        for (size_t j = 0; j < symbol_size; ++j)
            inf_symbols->symbols[i]->data[j] = (uint8_t)rand();
    }
}

static void init_rcv_symbols(const symbol_seq_t* src_symbols, symbol_seq_t* rcv_symbols) {
    assert(src_symbols->symbol_size == rcv_symbols->symbol_size);
    assert(src_symbols->length == rcv_symbols->length);

    size_t symbol_size = src_symbols->symbol_size;
    size_t length = src_symbols->length;

    for (size_t i = 0; i < length; ++i) {
        memcpy((void*)rcv_symbols->symbols[i]->data, (void*)src_symbols->symbols[i]->data, symbol_size);
    }
}

static void choose_erased(uint16_t n, uint16_t t, bool* is_erased) {
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
}

static void erase_symbols(symbol_seq_t* rcv_symbols, uint16_t t, const bool* is_erased) {
    size_t symbol_size = rcv_symbols->symbol_size;
    uint16_t n = (uint16_t)rcv_symbols->length;

    for (uint16_t i = 0; i < n; ++i) {
        if (!is_erased[i])
            continue;
        memset((void*)rcv_symbols->symbols[i]->data, 0, symbol_size);
    }
}

static int init_data(size_t symbol_size, uint16_t k, uint16_t r, uint16_t t, symbol_seq_t** src_symbols,
                     symbol_seq_t** rcv_symbols, bool** is_erased, symbol_seq_t* inf_symbols,
                     symbol_seq_t* rep_symbols) {
    assert(src_symbols != NULL);
    assert(rcv_symbols != NULL);
    assert(is_erased != NULL);
    assert(inf_symbols != NULL);
    assert(rep_symbols != NULL);

    *src_symbols = seq_create(k + r, symbol_size);
    if (!*src_symbols) {
        printf("ERROR: seq_create returned NULL\n");
        return 1;
    }

    *rcv_symbols = seq_create(k + r, symbol_size);
    if (!*rcv_symbols) {
        printf("ERROR: seq_create returned NULL\n");
        seq_destroy(*src_symbols);
        return 1;
    }

    *is_erased = (bool*)calloc(k + r, sizeof(bool));
    if (!*is_erased) {
        printf("ERROR: couldn't allocate is_erased\n");
        seq_destroy(*rcv_symbols);
        seq_destroy(*src_symbols);
        return 1;
    }

    inf_symbols->symbol_size = symbol_size;
    inf_symbols->length = k;
    inf_symbols->symbols = (*src_symbols)->symbols;

    generate_inf_symbols(inf_symbols);

    rep_symbols->symbol_size = symbol_size;
    rep_symbols->length = r;
    rep_symbols->symbols = (*src_symbols)->symbols + k;

    choose_erased(k + r, t, *is_erased);

    return 0;
}

static void destroy_data(symbol_seq_t* src_symbols, symbol_seq_t* rcv_symbols, bool* is_erased) {
    free(is_erased);
    seq_destroy(rcv_symbols);
    seq_destroy(src_symbols);
}

int main(void) {
    RS_t* rs;

    rs = rs_create();
    if (!rs) {
        printf("ERROR: rs_create returned NULL\n");
        return 1;
    }

    srand(SEED);

    size_t symbol_size = 10; // size of the symbol
    uint16_t k = 100;        // nummber of information (source) symbols
    uint16_t r = 10;         // number of repair (redundant) symbols
    uint16_t t = r;          // number of erased symbols

    symbol_seq_t* src_symbols;
    symbol_seq_t* rcv_symbols;
    symbol_seq_t inf_symbols;
    symbol_seq_t rep_symbols;
    bool* is_erased;
    int err;

    err = init_data(symbol_size, k, r, t, &src_symbols, &rcv_symbols, &is_erased, &inf_symbols, &rep_symbols);
    if (err) {
        rs_destroy(rs);
        return err;
    }

    /* ===== ENCODING =====*/
    err = rs_generate_repair_symbols(rs, &inf_symbols, &rep_symbols);
    if (err) {
        printf("ERROR: rs_generate_repair_symbols returned %d\n", err);
        destroy_data(src_symbols, rcv_symbols, is_erased);
        rs_destroy(rs);
        return err;
    }

    /* ===== ERASING SYMBOLS ===== */
    init_rcv_symbols(src_symbols, rcv_symbols);
    erase_symbols(rcv_symbols, t, is_erased);

    /* ===== RECOVERING ===== */
    err = rs_restore_symbols(rs, k, r, rcv_symbols, is_erased, t);
    if (err) {
        printf("ERROR: rs_restore_symbols returned %d\n", err);
        destroy_data(src_symbols, rcv_symbols, is_erased);
        rs_destroy(rs);
        return err;
    }

    assert(seq_eq(src_symbols, rcv_symbols));

    destroy_data(src_symbols, rcv_symbols, is_erased);
    rs_destroy(rs);

    return 0;
}