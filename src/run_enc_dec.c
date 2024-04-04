#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <rlc/rlc.h>
#include <rs/reed_solomon.h>

#define SEED 78934

static void generate_inf_symbols(const symbol_seq_t* inf_symbols) {
    assert(inf_symbols != NULL);

    size_t symbol_size = inf_symbols->symbol_size;
    size_t length = inf_symbols->length;

    for (size_t i = 0; i < length; ++i) {
        for (size_t j = 0; j < symbol_size; ++j)
            inf_symbols->symbols[i]->data[j] = (uint8_t)rand();
    }
}

static void init_rcv_symbols(const symbol_seq_t* src_symbols, symbol_seq_t* rcv_symbols) {
    assert(src_symbols != NULL);
    assert(rcv_symbols != NULL);
    assert(src_symbols->symbol_size == rcv_symbols->symbol_size);
    assert(src_symbols->length == rcv_symbols->length);

    size_t symbol_size = src_symbols->symbol_size;
    size_t length = src_symbols->length;

    for (size_t i = 0; i < length; ++i) {
        memcpy((void*)rcv_symbols->symbols[i]->data, (void*)src_symbols->symbols[i]->data,
               symbol_size);
    }
}

static void choose_erased(uint16_t n, uint16_t t, bool* is_erased) {
    assert(is_erased != NULL);

    memset((void*)is_erased, 0, n * sizeof(bool));

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
    assert(rcv_symbols != NULL);
    assert(is_erased != NULL);

    size_t symbol_size = rcv_symbols->symbol_size;
    uint16_t n = (uint16_t)rcv_symbols->length;

    for (uint16_t i = 0; i < n; ++i) {
        if (!is_erased[i])
            continue;
        memset((void*)rcv_symbols->symbols[i]->data, 0, symbol_size);
    }
}

static int init_data(size_t symbol_size, uint16_t k, uint16_t r, uint16_t t,
                     symbol_seq_t** src_symbols, symbol_seq_t** rcv_symbols, uint32_t** seeds,
                     bool** is_erased, symbol_seq_t* inf_symbols, symbol_seq_t* rep_symbols) {
    assert(src_symbols != NULL);
    assert(rcv_symbols != NULL);
    assert(seeds != NULL);
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

    *seeds = (uint32_t*)calloc(r, sizeof(uint32_t));
    if (!*seeds) {
        printf("ERROR: couldn't allocate seeds\n");
        seq_destroy(*rcv_symbols);
        seq_destroy(*src_symbols);
        return 1;
    }

    *is_erased = (bool*)calloc(k + r, sizeof(bool));
    if (!*is_erased) {
        printf("ERROR: couldn't allocate is_erased\n");
        free(*seeds);
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

static int run_rs_enc(RS_t* rs, const symbol_seq_t* inf_symbols, symbol_seq_t* rep_symbols) {
    int err;

    err = rs_generate_repair_symbols(rs, inf_symbols, rep_symbols);
    if (err) {
        printf("ERROR: rs_generate_repair_symbols returned %d\n", err);
        return err;
    }

    return 0;
}

static int run_rs_dec(RS_t* rs, uint16_t k, uint16_t r, symbol_seq_t* rcv_symbols,
                      const bool is_erased[], uint16_t t) {
    int err;

    err = rs_restore_symbols(rs, k, r, rcv_symbols, is_erased, t);
    if (err) {
        printf("ERROR: rs_restore_symbols returned %d\n", err);
        return err;
    }

    return 0;
}

static int run_rlc_enc(RLC_t* rlc, const symbol_seq_t* inf_symbols, symbol_seq_t* rep_symbols,
                       uint32_t seeds[]) {
    int err;

    err = rlc_generate_repair_symbols(rlc, inf_symbols, rep_symbols, seeds);
    if (err) {
        printf("ERROR: rlc_generate_repair_symbols returned %d\n", err);
        return err;
    }

    return 0;
}

static int run_rlc_dec(RLC_t* rlc, uint16_t k, uint16_t r, symbol_seq_t* rcv_symbols,
                       const uint32_t seeds[], const bool is_erased[], uint16_t t) {
    int err;

    err = rlc_restore_symbols(rlc, k, r, rcv_symbols, seeds, is_erased, t);
    if (err) {
        printf("ERROR: rlc_restore_symbols returned %d\n", err);
        return err;
    }

    return 0;
}

typedef enum algo_type {
    RS,
    RLC,
    NO,
} algo_type_t;

static int parse_args(int argc, char* argv[], algo_type_t* algo_type, uint16_t* k, uint16_t* r,
                      uint16_t* t) {
    if (argc < 4) {
        printf("Expected at least 3 arguments: <RS/RLC/NO> <k> <r> [t]\n");
        return 2;
    } else if (argc > 5) {
        printf("Expected not more than 4 arguments: <RS/RLC/NO> <k> <r> [t]\n");
        return 2;
    }

    if (strcmp(argv[1], "RS") == 0) {
        *algo_type = RS;
    } else if (strcmp(argv[1], "RLC") == 0) {
        *algo_type = RLC;
    } else if (strcmp(argv[1], "NO") == 0) {
        *algo_type = NO;
    } else {
        printf("Unknown algorithm: %s\n", argv[1]);
        return 3;
    }

    char* end = NULL;
    *k = strtoul(argv[2], &end, 10);
    *r = strtoul(argv[3], &end, 10);
    *t = (argc >= 5) ? strtoul(argv[4], &end, 10) : *r;

    return 0;
}

int main(int argc, char* argv[]) {
    algo_type_t algo_type = NO;
    size_t symbol_size = SYMBOL_SIZE;
    uint16_t k = 0;
    uint16_t r = 0;
    uint16_t t = 0;
    int err;

    err = parse_args(argc, argv, &algo_type, &k, &r, &t);
    if (err)
        return err;

    RS_t* rs;
    RLC_t* rlc;

    rs = rs_create();
    if (!rs) {
        printf("ERROR: rs_create returned NULL\n");
        return 1;
    }

    rlc = rlc_create();
    if (!rlc) {
        printf("ERROR: rlc_create returned NULL\n");
        rs_destroy(rs);
        return 1;
    }

    srand(SEED);

    symbol_seq_t* src_symbols = NULL;
    symbol_seq_t* rcv_symbols = NULL;
    symbol_seq_t inf_symbols = {0};
    symbol_seq_t rep_symbols = {0};
    uint32_t* seeds = NULL;
    bool* is_erased = NULL;

    err = init_data(symbol_size, k, r, t, &src_symbols, &rcv_symbols, &seeds, &is_erased,
                    &inf_symbols, &rep_symbols);
    if (err) {
        rlc_destroy(rlc);
        rs_destroy(rs);
        return err;
    }

    const int ITER_COUNT = 100;

    switch (algo_type) {
    case RS:
        for (int i = 0; i < ITER_COUNT; ++i) {
            err = run_rs_enc(rs, &inf_symbols, &rep_symbols);
            if (err) {
                free(is_erased);
                free(seeds);
                seq_destroy(rcv_symbols);
                seq_destroy(src_symbols);
                rlc_destroy(rlc);
                rs_destroy(rs);
                return err;
            }

            init_rcv_symbols(src_symbols, rcv_symbols);
            erase_symbols(rcv_symbols, t, is_erased);

            err = run_rs_dec(rs, k, r, rcv_symbols, is_erased, t);
            if (err) {
                free(is_erased);
                free(seeds);
                seq_destroy(rcv_symbols);
                seq_destroy(src_symbols);
                rlc_destroy(rlc);
                rs_destroy(rs);
                return err;
            }
        }
        break;
    case RLC:
        for (int i = 0; i < ITER_COUNT; ++i) {
            rlc->current_repair_symbol = 0;

            err = run_rlc_enc(rlc, &inf_symbols, &rep_symbols, seeds);
            if (err) {
                free(is_erased);
                free(seeds);
                seq_destroy(rcv_symbols);
                seq_destroy(src_symbols);
                rlc_destroy(rlc);
                rs_destroy(rs);
                return err;
            }

            init_rcv_symbols(src_symbols, rcv_symbols);
            erase_symbols(rcv_symbols, t, is_erased);

            err = run_rlc_dec(rlc, k, r, rcv_symbols, seeds, is_erased, t);
            if (err) {
                free(is_erased);
                free(seeds);
                seq_destroy(rcv_symbols);
                seq_destroy(src_symbols);
                rlc_destroy(rlc);
                rs_destroy(rs);
                return err;
            }
        }
        break;
    case NO:
        for (int i = 0; i < ITER_COUNT; ++i) {
            init_rcv_symbols(src_symbols, rcv_symbols);
            erase_symbols(rcv_symbols, t, is_erased);
        }
        break;
    default:
        printf("ERROR: unreachable\n");
        return 100;
    }

    free(is_erased);
    free(seeds);
    seq_destroy(rcv_symbols);
    seq_destroy(src_symbols);
    rlc_destroy(rlc);
    rs_destroy(rs);

    return 0;
}