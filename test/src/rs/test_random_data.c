#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rs/reed_solomon.h>
#include <test/util/util.h>

#define SEED 234546127
#define TESTS_CNT 100

#define TEST_WRAPPER(_rs, _symbol_size, _k, _r, _t)                                                                    \
    do {                                                                                                               \
        if (test((_rs), (_symbol_size), (_k), (_r), (_t))) {                                                           \
            rs_destroy((_rs));                                                                                         \
            return 1;                                                                                                  \
        }                                                                                                              \
    } while (0)

static int test(RS_t* rs, size_t symbol_size, uint16_t k, uint16_t r, uint16_t t) {
    assert(t <= r);

    symbol_seq_t* src_symbols;
    symbol_seq_t* rcv_symbols;
    symbol_seq_t inf_symbols;
    symbol_seq_t rep_symbols;
    symbol_seq_t rcv_inf_symbols;
    bool* is_erased;
    int err;

    src_symbols = seq_create(k + r, symbol_size);
    if (!src_symbols) {
        printf("ERROR: seq_create returned NULL\n");
        return 1;
    }

    rcv_symbols = seq_create(k + r, symbol_size);
    if (!rcv_symbols) {
        printf("ERROR: seq_create returned NULL\n");
        seq_destroy(src_symbols);
        return 1;
    }

    is_erased = (bool*)calloc(k + r, sizeof(bool));
    if (!is_erased) {
        printf("ERROR: couldn't allocate is_erased\n");
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return 1;
    }

    inf_symbols.symbol_size = symbol_size;
    inf_symbols.length = k;
    inf_symbols.symbols = src_symbols->symbols;

    util_generate_inf_symbols(&inf_symbols);

    rep_symbols.symbol_size = symbol_size;
    rep_symbols.length = r;
    rep_symbols.symbols = src_symbols->symbols + k;

    rcv_inf_symbols.symbol_size = symbol_size;
    rcv_inf_symbols.length = k;
    rcv_inf_symbols.symbols = rcv_symbols->symbols;

    err = rs_generate_repair_symbols(rs, &inf_symbols, &rep_symbols);
    if (err) {
        printf("ERROR: rs_generate_repair_symbols returned %d\n", err);
        free(is_erased);
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return err;
    }

    util_init_rcv_symbols(src_symbols, rcv_symbols);
    util_choose_and_erase_symbols(rcv_symbols, t, is_erased);
    assert(!seq_eq(src_symbols, rcv_symbols));

    err = rs_restore_symbols(rs, k, r, rcv_symbols, is_erased, t);
    if (err) {
        printf("ERROR: rs_restore_symbols returned %d\n", err);
        free(is_erased);
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return err;
    }

    if (!seq_eq(&inf_symbols, &rcv_inf_symbols)) {
        printf("ERROR: inf_symbols != rcv_inf_symbols after restore:\n");

        printf("\tinf_symbols     = ");
        seq_printf(&inf_symbols);
        printf("\n");

        printf("\trcv_inf_symbols = ");
        seq_printf(&rcv_inf_symbols);
        printf("\n");

        err = 1;
    }

    free(is_erased);
    seq_destroy(rcv_symbols);
    seq_destroy(src_symbols);

    return err;
}

int main(void) {
    RS_t* rs;
    size_t symbol_size;
    uint16_t k;
    uint16_t r;
    uint16_t t;

    rs = rs_create();
    if (!rs) {
        printf("ERROR: rs_create returned NULL\n");
        return 1;
    }

    srand(SEED);

    for (int _i = 0; _i < TESTS_CNT / 2; ++_i) {
        symbol_size = 16;
        k = 100 + rand() % 100;
        r = 50 + rand() % 50;
        t = 11 + rand() % (r - 10);

        TEST_WRAPPER(rs, symbol_size, k, r, t);
    }

    for (int _i = 0; _i < (TESTS_CNT + 1) / 2; ++_i) {
        symbol_size = 16;
        k = 100 + rand() % 100;
        r = 50 + rand() % 50;
        t = r;

        TEST_WRAPPER(rs, symbol_size, k, r, t);
    }

    rs_destroy(rs);

    return 0;
}