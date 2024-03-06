/**
 * @file rlc.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief rlc/rlc.h implementation.
 * @date 2024-02-27
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <rlc/equation.h>
#include <rlc/gf256.h>
#include <rlc/rlc.h>
#include <rlc/system.h>

#include "prng/tinymt32.c"

#define ALIGNMENT 32
static __always_inline size_t align(size_t val) {
    return (val + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;
}

RLC_t* rlc_create() {
    gf256_init();

    RLC_t* rlc = (RLC_t*)malloc(sizeof(RLC_t));
    if (!rlc)
        return NULL;
    memset(rlc, 0, sizeof(RLC_t));

    uint8_t* inv_table = (uint8_t*)malloc(256 * sizeof(uint8_t));
    if (!inv_table) {
        free(rlc);
        return NULL;
    }
    memset(inv_table, 0, 256 * sizeof(uint8_t));
    assign_inv(inv_table);

    uint8_t** mul_table = (uint8_t**)malloc(256 * sizeof(uint8_t*));
    if (!mul_table) {
        free(inv_table);
        free(rlc);
        return NULL;
    }

    for (int i = 0; i < 256; ++i) {
        mul_table[i] = (uint8_t*)malloc(256 * sizeof(uint8_t));
        if (!mul_table[i]) {
            for (int j = 0; j < i; ++j)
                free(mul_table[j]);
            free(mul_table);
            free(inv_table);
            free(rlc);
            return NULL;
        }
        memset(mul_table[i], 0, 256 * sizeof(uint8_t));
    }
    assign_mul(mul_table);

    rlc->inv_table = inv_table;
    rlc->mul_table = mul_table;

    return rlc;
}

void rlc_destroy(RLC_t* rlc) {
    uint8_t* inv_table = rlc->inv_table;
    uint8_t** mul_table = rlc->mul_table;

    for (int i = 0; i < 256; ++i)
        free(mul_table[i]);
    free(mul_table);
    free(inv_table);
    free(rlc);
}

static inline void get_coefs(tinymt32_t* prng, uint32_t seed, int n, uint8_t coefs[n]) {
    tinymt32_init(prng, seed);
    int i;
    for (i = 0; i < n; i++) {
        coefs[i] = (uint8_t)tinymt32_generate_uint32(prng);
        if (coefs[i] == 0)
            coefs[i] = 1;
    }
}

static int get_one_coded_symbol(RLC_t* rlc, const symbol_seq_t* inf_symbols, symbol_t* rep_symbol,
                                uint32_t* seed) {
    tinymt32_t prng;
    prng.mat1 = 0x8f7011ee;
    prng.mat2 = 0xfc78ff1f;
    prng.tmat = 0x3793fdff;
    uint8_t** mul = rlc->mul_table;
    size_t symbol_size = inf_symbols->symbol_size;

    uint8_t* coefs = (uint8_t*)calloc(align(inf_symbols->length), sizeof(uint8_t));
    if (!coefs)
        return 1;

    *seed = rlc->current_repair_symbol++;

    get_coefs(&prng, *seed, inf_symbols->length, coefs);

    memset((void*)rep_symbol->data, 0, symbol_size);

    for (int j = 0; j < inf_symbols->length; j++) {
        gf256_symbol_add_scaled((void*)rep_symbol->data, coefs[j],
                                (void*)inf_symbols->symbols[j]->data, symbol_size, mul);
    }

    free(coefs);

    return 0;
}

int rlc_generate_repair_symbols(RLC_t* rlc, const symbol_seq_t* inf_symbols,
                                const symbol_seq_t* rep_symbols, uint32_t* seeds) {
    int err;

    for (uint16_t i = 0; i < rep_symbols->length; ++i) {
        err = get_one_coded_symbol(rlc, inf_symbols, rep_symbols->symbols[i], &seeds[i]);
        if (err)
            return err;
    }

    return 0;
}

static int receive_repair_symbol(RLC_t* rlc, system_t* system, const symbol_seq_t* inf_symbols,
                                 symbol_t* rep_symbol, uint32_t seed, const bool* is_erased) {
    uint16_t k = (uint16_t)inf_symbols->length;

    equation_t* eq = (equation_t*)malloc(sizeof(equation_t));
    if (!eq)
        return 1;
    eq->pivot = 0;
    eq->last_non_zero_id = k - 1;
    eq->n_coefs = k;
    eq->n_protected_symbols = k;
    eq->symbol_size = inf_symbols->symbol_size;
    eq->constant_term = rep_symbol;

    eq->coefs = (uint8_t*)calloc(align(k), sizeof(uint8_t));
    eq->_coefs_allocated_size = align(k);
    if (!eq->coefs) {
        free(eq);
        return 1;
    }

    tinymt32_t prng;
    prng.mat1 = 0x8f7011ee;
    prng.mat2 = 0xfc78ff1f;
    prng.tmat = 0x3793fdff;
    uint8_t** mul_table = rlc->mul_table;
    uint8_t* inv_table = rlc->inv_table;

    get_coefs(&prng, seed, eq->n_coefs, eq->coefs);

    for (uint16_t i = 0; i < k; ++i) {
        if (is_erased[i])
            continue;
        gf256_symbol_add_scaled((void*)eq->constant_term->data, eq->coefs[i],
                                (void*)inf_symbols->symbols[i]->data, inf_symbols->symbol_size,
                                mul_table);
        eq->coefs[i] = 0;
    }

    equation_adjust_non_zero_bounds(eq);
    if (equation_has_one_id(eq)) {
        equation_multiply(eq, inv_table[equation_get_coef(eq, eq->pivot)], mul_table);
    }
    if (equation_is_zero(eq)) {
        free(eq->coefs);
        free(eq);
        return 0;
    }

    // if (system->first_id_id == ID_NONE && eq->pivot != 0) {
    //     // this can happen when there are unknown source symbols for which no
    //     // repair symbol has been received yet in that case, we add the
    //     // knowledge of symbols from 0 to pivot - 1 to the system
    //     system_set_bounds(system, 0, eq->pivot - 1);
    // }
    if (system->first_id_id == ID_NONE)
        system_set_bounds(system, 0, 0);

    int decoded = 0;
    equation_t* removed = NULL;
    int used_in_system = 0;

    system_add_with_elimination(system, eq, inv_table, mul_table, &decoded, &removed,
                                &used_in_system);
    if (!used_in_system) {
        free(eq->coefs);
        free(eq);
    }

    return 0;
}

int rlc_restore_symbols(RLC_t* rlc, uint16_t k, uint16_t r, symbol_seq_t* rcv_symbols,
                        const uint32_t* seeds, const bool* is_erased, uint16_t t) {
    symbol_seq_t inf_symbols;
    symbol_seq_t rep_symbols;
    int err = 0;

    inf_symbols.length = k;
    inf_symbols.symbol_size = rcv_symbols->symbol_size;
    inf_symbols.symbols = rcv_symbols->symbols;

    rep_symbols.length = r;
    rep_symbols.symbol_size = rcv_symbols->symbol_size;
    rep_symbols.symbols = rcv_symbols->symbols + k;

    equation_t** equations = (equation_t**)calloc(k + r, sizeof(equation_t*));
    if (!equations)
        return 1;

    system_t system = {
        .max_equations = k + r,
        .n_equations = 0,
        .first_id_id = ID_NONE,
        .last_symbol_id = ID_NONE,
        .equations = equations,
    };

    for (uint16_t i = 0; i < r; ++i) {
        if (is_erased[k + i])
            continue;
        err = receive_repair_symbol(rlc, &system, &inf_symbols, rep_symbols.symbols[i], seeds[i],
                                    is_erased);
        if (err) {
            for (uint16_t j = 0; j < system.n_equations; ++j) {
                free(system.equations[j]->coefs);
                free(system.equations[j]);
            }
            return err;
        }
    }

    for (uint16_t i = 0; i < k; ++i) {
        if (!is_erased[i])
            continue;
        assert(system.equations[i] != NULL);

        equation_t* eq = system.equations[i];

        assert(equation_has_one_id(eq));
        assert(equation_get_coef(eq, i) == 1);

        memcpy((void*)rcv_symbols->symbols[i]->data, (void*)eq->constant_term->data,
               rcv_symbols->symbol_size);
    }

    for (uint16_t i = 0; i < system.max_equations; ++i) {
        if (!system.equations[i])
            continue;
        free(system.equations[i]->coefs);
        free(system.equations[i]);
    }
    free(equations);

    return 0;
}