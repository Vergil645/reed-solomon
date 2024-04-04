#ifndef __RLC_EQUATION_H__
#define __RLC_EQUATION_H__

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <memory/symbol.h>

#define ID_NONE 0xffu

typedef struct {
    uint16_t pivot; // id of the pivot (i.e. first non-zero coefficient in the
                    // equation)
    uint16_t last_non_zero_id;
    uint16_t n_coefs;

    uint16_t n_protected_symbols;
    size_t symbol_size; // must be multiplied by 2
    symbol_t* constant_term;

    uint8_t* coefs;
    size_t _coefs_allocated_size;
} equation_t;

/**
 * @brief get the minimum source index that appears in the symbol
 *        ID_NONE if there is none (e.g. symbol is 0)
 */
uint32_t equation_get_min_symbol_id(equation_t* eq);

/**
 * @brief get the maximum source index that appears in the symbol
 *        returns ID_NONE if there is none (e.g. symbol is 0)
 */
uint32_t equation_get_max_symbol_id(equation_t* eq);

uint32_t equation_count_allocated_coef(equation_t* equation);

uint8_t equation_get_coef(equation_t* eq, uint16_t i);

/**
 * update full symbol max coef
 */
bool full_symbol_adjust_max_coef(equation_t* eq);

/**
 * update full symbol min coef
 */
bool full_symbol_adjust_min_coef(equation_t* eq);

/**
 *update full symbol min and max coefs
 *  returns whether there exists non-zero coefs
 */

bool equation_adjust_non_zero_bounds(equation_t* eq);

/**
 * @brief Returns whether the symbol is an empty symbol
 */
bool equation_is_zero(equation_t* eq);

/**
 * @brief Returns whether the symbol is an empty symbol
 */
bool equation_has_one_id(equation_t* full_symbol);

void equation_multiply(equation_t* eq, uint8_t coef, uint8_t** mul_table);

int equation_add(equation_t* eq1, equation_t* eq2);

#endif