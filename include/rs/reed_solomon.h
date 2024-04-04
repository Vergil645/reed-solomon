/**
 * @file reed_solomon.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains implementaion of Reed-Solomon codes over GF(65536).
 * @date 2024-01-18
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_H__
#define __REED_SOLOMON_H__

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

#include "cyclotomic_coset.h"
#include "gf65536.h"
#include <memory/seq.h>

/**
 * @brief Maximum number of cyclotomic coset locator polynomial coefficients.
 */
#define RS_COSET_LOCATOR_MAX_LEN (CC_MAX_COSET_SIZE + 1)

/**
 * @brief Return code for cases when erases cannot be restored due to code
 * parameters.
 */
#define RS_ERR_CANNOT_RESTORE 100

/**
 * @brief Context data.
 */
typedef struct {
    GF_t* gf;
    CC_t* cc;
} RS_t;

/**
 * @brief Create context object.
 *
 * @return pointer to created context object on success and NULL otherwise.
 */
RS_t* rs_create();

/**
 * @brief Destroy context object.
 *
 * @param rs context object.
 */
void rs_destroy(RS_t* rs);

/**
 * @brief Generate repair symbols for the given information symbols.
 *
 * @param rs context object.
 * @param inf_symbols information symbols.
 * @param rep_symbols where to place the result.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int rs_generate_repair_symbols(RS_t* rs, const symbol_seq_t* inf_symbols,
                               symbol_seq_t* rep_symbols);

/**
 * @brief Restore erased symbols. Assume that rcv_symbols[i] = 0 for erased i.
 *
 * @param rs context object.
 * @param k number of information symbols.
 * @param r number of repair symbols.
 * @param rcv_symbols received symbols, restored symbols will be written here.
 * @param is_erased indicates which symbols has been erased.
 * @param t number of erases.
 * @return 0 on success,\n
 *         1 on memory allocation error,\n
 *         or RS_ERR_CANNOT_RESTORE.
 */
int rs_restore_symbols(RS_t* rs, uint16_t k, uint16_t r, symbol_seq_t* rcv_symbols,
                       const bool* is_erased, uint16_t t);

#endif