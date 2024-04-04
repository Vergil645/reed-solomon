/**
 * @file rlc.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains implementaion of Random Linear codes over GF(256).
 * @date 2024-02-27
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __RLC_H__
#define __RLC_H__

#include <stdbool.h>

#include <memory/seq.h>

typedef struct {
    uint32_t current_repair_symbol;
    uint8_t* inv_table;
    uint8_t** mul_table;
} RLC_t;

/**
 * @brief Create context object.
 *
 * @return pointer to context object or NULL if error occured.
 */
RLC_t* rlc_create();

/**
 * @brief Destroy context object.
 *
 * @param rlc context object.
 */
void rlc_destroy(RLC_t* rlc);

/**
 * @brief Generate repair symbols for the given information symbols.
 *
 * @param rlc context object.
 * @param inf_symbols information symbols.
 * @param rep_symbols where to place the result.
 * @param seeds where to place random seeds used to generate coefficients.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int rlc_generate_repair_symbols(RLC_t* rlc, const symbol_seq_t* inf_symbols, const symbol_seq_t* rep_symbols,
                                uint32_t* seeds);

/**
 * @brief Restore erased symbols. Assume that rcv_symbols[i] = 0 for erased i.
 *
 * @param rlc context object.
 * @param k number of information symbols.
 * @param r number of repair symbols.
 * @param rcv_symbols received symbols, restored symbols will be written here.
 * @param seeds random seeds used to generate coefficients.
 * @param is_erased indicates which symbols has been erased.
 * @param t number of erases.
 * @return 0 on success,\n
 *         1 on memory allocation error,\n
 *         or RS_ERR_CANNOT_RESTORE.
 */
int rlc_restore_symbols(RLC_t* rlc, uint16_t k, uint16_t r, symbol_seq_t* rcv_symbols, const uint32_t* seeds,
                        const bool* is_erased, uint16_t t);

#endif