#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <reed_solomon.h>

#include <util.h>

#define SEED 234546127
#define TESTS_CNT 100

#define TEST_WRAPPER(_rs, _symbol_size, _k, _r, _t)                            \
    do {                                                                       \
        if (test((_rs), (_symbol_size), (_k), (_r), (_t))) {                   \
            free_all((_rs));                                                   \
            return 1;                                                          \
        }                                                                      \
    } while (0)

static void generate_src_symbols(symbol_seq_t src_symbols) {
    size_t symbol_size = src_symbols.symbol_size;
    size_t length = src_symbols.length;

    for (size_t i = 0; i < length; ++i) {
        for (size_t j = 0; j < symbol_size; ++j) {
            src_symbols.symbols[i].data[j] = 1 + ((element_t)rand()) % UINT16_MAX;

            assert(src_symbols.symbols[i].data[j] != 0);
        }
    }
}

static void init_rcv_symbols(symbol_seq_t src_symbols,
                             symbol_seq_t rcv_symbols) {
    assert(src_symbols.symbol_size == rcv_symbols.symbol_size);
    assert(src_symbols.length == rcv_symbols.length);

    size_t symbol_size = src_symbols.symbol_size;
    size_t length = src_symbols.length;

    for (size_t i = 0; i < length; ++i) {
        memcpy((void *)rcv_symbols.symbols[i].data,
               (void *)src_symbols.symbols[i].data,
               symbol_size * sizeof(element_t));
    }
}

static void erase_symbols(symbol_seq_t rcv_symbols, uint16_t t,
                          uint16_t *erased_indices) {
    size_t symbol_size = rcv_symbols.symbol_size;
    uint16_t n = (uint16_t)rcv_symbols.length;
    uint16_t s = 0;

    for (uint16_t i = 0; i < t; ++i) {
        assert(n - s - t + i + 1 > 0);

        erased_indices[i] = s + (uint16_t)rand() % (n - s - t + i + 1);
        s = erased_indices[i] + 1;

        assert(erased_indices[i] < n);
    }

    for (uint16_t i = 0; i < t; ++i) {
        memset((void *)rcv_symbols.symbols[erased_indices[i]].data, 0,
               symbol_size * sizeof(element_t));
    }
}

static int test(const RS_t *rs, size_t symbol_size, uint16_t k, uint16_t r,
                uint16_t t) {
    assert(t <= r);

    symbol_seq_t src_symbols;
    symbol_seq_t inf_symbols;
    symbol_seq_t rep_symbols;
    symbol_seq_t rcv_symbols;
    uint16_t *erased_indices;
    int ret;

    ret = seq_alloc(symbol_size, k + r, &src_symbols);
    if (ret) {
        printf("ERROR: seq_alloc returned %d\n", ret);
        return ret;
    }

    ret = seq_alloc(symbol_size, k + r, &rcv_symbols);
    if (ret) {
        printf("ERROR: seq_alloc returned %d\n", ret);
        seq_free(&src_symbols);
        return ret;
    }

    erased_indices = (uint16_t *)malloc(t * sizeof(uint16_t));
    if (!erased_indices) {
        printf("ERROR: couldn't allocate erased_indices\n");
        seq_free(&rcv_symbols);
        seq_free(&src_symbols);
        return 1;
    }

    inf_symbols.symbol_size = symbol_size;
    inf_symbols.length = k;
    inf_symbols.symbols = src_symbols.symbols;

    rep_symbols.symbol_size = symbol_size;
    rep_symbols.length = r;
    rep_symbols.symbols = src_symbols.symbols + k;

    generate_src_symbols(src_symbols);
    rs_generate_repair_symbols(rs, inf_symbols, rep_symbols);

    init_rcv_symbols(src_symbols, rcv_symbols);
    erase_symbols(rcv_symbols, t, erased_indices);
    assert(!util_seq_eq(src_symbols, rcv_symbols));
    rs_restore_symbols(rs, k, r, rcv_symbols, erased_indices, t);

    if (!util_seq_eq(src_symbols, rcv_symbols)) {
        printf("ERROR: src_symbols != rcv_symbols after restore:\n");

        printf("\tsrc_symbols = ");
        util_seq_printf(src_symbols);
        printf("\n");

        printf("\trcv_symbols  = ");
        util_seq_printf(rcv_symbols);
        printf("\n");

        ret = 1;
    }

    free(erased_indices);
    seq_free(&rcv_symbols);
    seq_free(&src_symbols);

    return ret;
}

static void free_all(RS_t *rs) {
    GF_t *gf = rs->gf;
    CC_t *cc = rs->cc;
    FFT_t *fft = rs->fft;

    rs_free(rs);
    fft_free(fft);
    cc_free(cc);
    gf_free(gf);
}

int main(void) {
    GF_t _gf;
    CC_t _cc;
    FFT_t _fft;
    RS_t _rs;
    RS_t *rs = &_rs;
    size_t symbols_size;
    uint16_t k;
    uint16_t r;
    uint16_t t;

    {
        GF_t *gf = &_gf;
        CC_t *cc = &_cc;
        FFT_t *fft = &_fft;

        int ret = 0;

        ret = gf_alloc(gf);
        if (ret) {
            printf("ERROR: gf_alloc returned %d\n", ret);
            return ret;
        }

        gf_init(gf);

        ret = cc_alloc(cc);
        if (ret) {
            printf("ERROR: cc_alloc returned %d\n", ret);
            gf_free(gf);
            return ret;
        }

        cc_init(cc);

        ret = fft_alloc(fft);
        if (ret) {
            printf("ERROR: fft_alloc returned %d\n", ret);
            cc_free(cc);
            gf_free(gf);
            return ret;
        }

        fft_init(fft, gf);

        ret = rs_alloc(rs);
        if (ret) {
            printf("ERROR: rs_alloc returned %d\n", ret);
            fft_free(fft);
            cc_free(cc);
            gf_free(gf);
            return ret;
        }

        rs_init(rs, gf, cc, fft);
    }

    srand(SEED);

    for (int _i = 0; _i < TESTS_CNT; ++_i) {
        symbols_size = 1 + rand() % 5;
        k = 100 + rand() % 100;
        r = 50 + rand() % 50;
        t = 11 + rand() % (r - 10);

        TEST_WRAPPER(rs, symbols_size, k, r, t);
    }

    free_all(rs);

    return 0;
}