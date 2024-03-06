#include <assert.h>
#include <stdbool.h>

#include <rlc/equation.h>
#include <rlc/gf256.h>

#define MIN(a, b) (((a) <= (b)) ? (a) : (b))

// #define ALIGNMENT 32
// static __always_inline size_t align(size_t val) {
//     return (val + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;
// }

/**
 * @brief get the minimum source index that appears in the symbol
 *        ID_NONE if there is none (e.g. symbol is 0)
 */
uint32_t equation_get_min_symbol_id(equation_t* eq) { return eq->pivot; }

/**
 * @brief get the maximum source index that appears in the symbol
 *        returns ID_NONE if there is none (e.g. symbol is 0)
 */
uint32_t equation_get_max_symbol_id(equation_t* eq) { return eq->last_non_zero_id; }

uint32_t equation_count_allocated_coef(equation_t* equation) {
    if (equation->pivot == ID_NONE) {
        assert(equation->last_non_zero_id == ID_NONE);
        return 0;
    } else {
        return equation->n_protected_symbols;
    }
}

uint8_t equation_get_coef(equation_t* eq, uint16_t i) {
    // if out of bounds, the coef is 0
    if (i > eq->n_protected_symbols) {
        return 0;
    }
    return eq->coefs[i];
}

/**
 * update full symbol max coef
 */
bool full_symbol_adjust_max_coef(equation_t* eq) {
    assert(eq->n_protected_symbols > 0);

    bool result = false;
    eq->last_non_zero_id = ID_NONE;

    for (uint16_t j = 0; j < equation_count_allocated_coef(eq); j++) {
        uint16_t i = eq->n_protected_symbols - 1 - j;
        if (equation_get_coef(eq, i) != 0) {
            eq->last_non_zero_id = i;
            result = true;
            break;
        }
    }
    return result;
}

/**
 * update full symbol min coef
 */
bool full_symbol_adjust_min_coef(equation_t* eq) {
    assert(eq->n_protected_symbols > 0);

    bool result = false;
    eq->pivot = ID_NONE;

    bool adjust_fast = true;
    if (adjust_fast) {
        uint8_t* coefs = eq->coefs;
        uint64_t* coefs_64 = (uint64_t*)eq->coefs;
        int64_t i;
        for (i = 0; result != true && i < eq->n_protected_symbols / sizeof(uint64_t); i++) {
            if (coefs_64[i] != 0) {
                for (int64_t j = i * sizeof(uint64_t); j < (i + 1) * sizeof(uint64_t); j++) {
                    if (coefs[j] != 0) {
                        eq->pivot = j;
                        result = true;
                        break;
                    }
                }
            }
        }
        if (result != true) {
            for (int64_t j = i * sizeof(uint64_t);
                 j < (i + 1) * sizeof(uint64_t) && j < eq->n_protected_symbols; j++) {
                if (coefs[j] != 0) {
                    eq->pivot = j;
                    result = true;
                    break;
                }
            }
        }

    } else {
        for (uint16_t i = 0; i < eq->n_protected_symbols; i++) {
            if (equation_get_coef(eq, i) != 0) {
                eq->pivot = i;
                result = true;
                break;
            }
        }
    }
    return result;
}

/**
 *update full symbol min and max coefs
 *  returns whether there exists non-zero coefs
 */

bool equation_adjust_non_zero_bounds(equation_t* eq) {
    if (eq->n_protected_symbols == 0) {
        eq->pivot = ID_NONE;
        eq->last_non_zero_id = ID_NONE;
        return false;
    }
    bool result1 = full_symbol_adjust_min_coef(eq);
    __attribute_maybe_unused__ bool result2 = full_symbol_adjust_max_coef(eq);
    assert(result1 == result2);
    return result1;
}

/**
 * @brief Returns whether the symbol is an empty symbol
 */
bool equation_is_zero(equation_t* eq) { return equation_get_min_symbol_id(eq) == ID_NONE; }

/**
 * @brief Returns whether the symbol is an empty symbol
 */
bool equation_has_one_id(equation_t* full_symbol) {
    return (!equation_is_zero(full_symbol)) &&
           (equation_get_min_symbol_id(full_symbol) == equation_get_max_symbol_id(full_symbol));
}

void equation_multiply(equation_t* eq, uint8_t coef, uint8_t** mul_table) {
    // multiply the coefficients of the equation (we can do it with one call)
    gf256_symbol_mul(eq->coefs, coef, eq->n_coefs, mul_table);
    // multiply the constant term of the equation
    gf256_symbol_mul(eq->constant_term->data, coef, eq->symbol_size, mul_table);
}

static void add_coefs(equation_t* eq1, equation_t* eq2, uint16_t from, uint16_t to) {
    // from = MAX(from, eq2->constant_term.metadata.first_id);
    // to = MIN(to, repair_symbol_last_id(&eq2->constant_term));
    // uint8_t *eq1_coefs_buffer = &eq1->coefs[from -
    // eq1->constant_term.metadata.first_id]; uint8_t *eq2_coefs_buffer =
    // &eq2->coefs[from - eq2->constant_term.metadata.first_id];

    uint8_t* eq1_coefs_buffer = &eq1->coefs[from];
    uint8_t* eq2_coefs_buffer = &eq2->coefs[from];

    gf256_symbol_add((void*)eq1_coefs_buffer, (void*)eq2_coefs_buffer,
                     MIN(to + 1 - from, eq1->n_protected_symbols));
}

static void full_symbol_add_base(equation_t* eq1, equation_t* eq2) {
    assert(eq1->constant_term->data != NULL && eq2->constant_term->data != NULL);

    //    uint32_t first_coef_index;
    //    uint32_t last_coef_index;

    // XXX: should not be NONE
    if (eq1->pivot == ID_NONE && eq2->pivot == ID_NONE) {
        eq1->pivot = ID_NONE;
        eq1->last_non_zero_id = ID_NONE;
    }
    // if (eq1->last_non_zero_id > eq2->last_non_zero_id){
    //     equation_clear_unused_coefs(eq2);
    // }
    bool add_fast = true;
    if (!add_fast) {
        // source_symbol_id_t first_rs_id =
        // eq1->constant_term.metadata.first_id;
        uint16_t first_rs_id = 0;
        for (uint32_t i = eq2->pivot; i <= eq2->last_non_zero_id; i++) {
            eq1->coefs[i - first_rs_id] = equation_get_coef(eq1, i) ^ equation_get_coef(eq2, i);
        }
    } else {
        add_coefs(eq1, eq2, eq2->pivot, eq2->last_non_zero_id);
    }

    equation_adjust_non_zero_bounds(eq1);
    gf256_symbol_add((void*)eq1->constant_term->data, (void*)eq2->constant_term->data,
                     eq2->symbol_size);
}

int equation_add(equation_t* eq1, equation_t* eq2) {
    // here were a lot of memory magic

    // results stored in eq1
    full_symbol_add_base(eq1, eq2);

    return 0;
}
