#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <rlc/rlc.h>
#include <rs/reed_solomon.h>

#define SEED 78934
#define TESTS_CNT 100
#define Z 1.959963984540054 // q_{1-0.05/2} = q_{0.975} of nomral distribution

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

static int compare(RS_t* rs, RLC_t* rlc, size_t symbol_size, uint16_t k, uint16_t r, uint16_t t,
                   double* enc_rs_rlc_ratio, double* dec_rs_rlc_ratio) {
    clock_t start;
    clock_t end;
    clock_t elapsed_time_enc_rs;
    clock_t elapsed_time_enc_rlc;
    clock_t elapsed_time_dec_rs;
    clock_t elapsed_time_dec_rlc;
    symbol_seq_t* src_symbols;
    symbol_seq_t* rcv_symbols;
    symbol_seq_t inf_symbols;
    symbol_seq_t rep_symbols;
    uint32_t* seeds;
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

    seeds = (uint32_t*)calloc(r, sizeof(uint32_t));
    if (!seeds) {
        printf("ERROR: couldn't allocate seeds\n");
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return 1;
    }

    is_erased = (bool*)calloc(k + r, sizeof(bool));
    if (!is_erased) {
        printf("ERROR: couldn't allocate is_erased\n");
        free(seeds);
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return 1;
    }

    inf_symbols.symbol_size = symbol_size;
    inf_symbols.length = k;
    inf_symbols.symbols = src_symbols->symbols;

    generate_inf_symbols(&inf_symbols);

    rep_symbols.symbol_size = symbol_size;
    rep_symbols.length = r;
    rep_symbols.symbols = src_symbols->symbols + k;

    choose_erased(k + r, t, is_erased);

    /* ===== RLC ===== */

    start = clock();
    err = rlc_generate_repair_symbols(rlc, &inf_symbols, &rep_symbols, seeds);
    end = clock();
    elapsed_time_enc_rlc = end - start;

    if (err) {
        printf("ERROR: rlc_generate_repair_symbols returned %d\n", err);
        free(is_erased);
        free(seeds);
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return err;
    }

    init_rcv_symbols(src_symbols, rcv_symbols);
    erase_symbols(rcv_symbols, t, is_erased);

    start = clock();
    err = rlc_restore_symbols(rlc, k, r, rcv_symbols, seeds, is_erased, t);
    end = clock();
    elapsed_time_dec_rlc = end - start;

    if (err) {
        printf("ERROR: rlc_restore_symbols returned %d\n", err);
        free(is_erased);
        free(seeds);
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return err;
    }

    /* ===== Reed-Solomon ===== */

    start = clock();
    err = rs_generate_repair_symbols(rs, &inf_symbols, &rep_symbols);
    end = clock();
    elapsed_time_enc_rs = end - start;

    if (err) {
        printf("ERROR: rs_generate_repair_symbols returned %d\n", err);
        free(is_erased);
        free(seeds);
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return err;
    }

    init_rcv_symbols(src_symbols, rcv_symbols);
    erase_symbols(rcv_symbols, t, is_erased);

    start = clock();
    err = rs_restore_symbols(rs, k, r, rcv_symbols, is_erased, t);
    end = clock();
    elapsed_time_dec_rs = end - start;

    if (err) {
        printf("ERROR: rs_restore_symbols returned %d\n", err);
        free(is_erased);
        free(seeds);
        seq_destroy(rcv_symbols);
        seq_destroy(src_symbols);
        return 1;
    }

    *enc_rs_rlc_ratio = (double)elapsed_time_enc_rs / (double)elapsed_time_enc_rlc;
    *dec_rs_rlc_ratio = (double)elapsed_time_dec_rs / (double)elapsed_time_dec_rlc;

    free(is_erased);
    free(seeds);
    seq_destroy(rcv_symbols);
    seq_destroy(src_symbols);

    return 0;
}

void calc_mean_with_delta(const double* x, size_t n, double* mean_res, double* delta_res) {
    assert(n >= 30);

    double mean = 0.0;
    double s = 0.0;
    double delta;

    for (size_t i = 0; i < n; ++i)
        mean += x[i];
    mean /= (double)n;

    for (size_t i = 0; i < n; ++i) {
        double tmp = x[i] - mean;
        s += tmp * tmp;
    }
    s = sqrt(s / (double)(n - 1));

    delta = s * Z / sqrt((double)n);

    *mean_res = mean;
    *delta_res = delta;
}

int main(void) {
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

    double enc_rs_rlc_ratio[TESTS_CNT] = {0};
    double dec_rs_rlc_ratio[TESTS_CNT] = {0};

    for (int i = 0; i < TESTS_CNT; ++i) {
        size_t symbol_size = SYMBOL_SIZE;
        // uint16_t k = 1000 + rand() % 1001;
        // uint16_t r = 10 + rand() % 71;
        uint16_t k = 2000;
        uint16_t r = 40;
        uint16_t t = r;
        int err;

        err = compare(rs, rlc, symbol_size, k, r, t, &enc_rs_rlc_ratio[i], &dec_rs_rlc_ratio[i]);
        if (err) {
            printf("ERROR: compare returned %d\n", err);
            rlc_destroy(rlc);
            rs_destroy(rs);
            return err;
        }
    }

    double enc_rs_rlc_ratio_mean = 0.0;
    double enc_rs_rlc_ratio_delta = 0.0;
    calc_mean_with_delta(enc_rs_rlc_ratio, TESTS_CNT, &enc_rs_rlc_ratio_mean, &enc_rs_rlc_ratio_delta);

    double dec_rs_rlc_ratio_mean = 0.0;
    double dec_rs_rlc_ratio_delta = 0.0;
    calc_mean_with_delta(dec_rs_rlc_ratio, TESTS_CNT, &dec_rs_rlc_ratio_mean, &dec_rs_rlc_ratio_delta);

    printf("===== clock() =====\n");
    printf("encode:\n");
    printf("    time ratio \"RS/RLC\": %.3f+-%.3f\n", enc_rs_rlc_ratio_mean, enc_rs_rlc_ratio_delta);
    printf("    time decreasing, %%: %.0f+-%.1f\n", 100.0 * (1 - enc_rs_rlc_ratio_mean),
           100.0 * enc_rs_rlc_ratio_delta);
    printf("decode:\n");
    printf("    time ratio \"RS/RLC\": %.3f+-%.3f\n", dec_rs_rlc_ratio_mean, dec_rs_rlc_ratio_delta);
    printf("    time decreasing, %%: %.0f+-%.1f\n", 100.0 * (1 - dec_rs_rlc_ratio_mean),
           100.0 * dec_rs_rlc_ratio_delta);

    rlc_destroy(rlc);
    rs_destroy(rs);

    return 0;
}