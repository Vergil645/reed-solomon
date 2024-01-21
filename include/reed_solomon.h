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
#include <stdint.h>

#include "constants.h"
#include "cyclotomic_coset.h"
#include "fft.h"
#include "symbol.h"

/**
 * @brief Maximum number of cyclotomic coset locator polynomial coefficients.
 */
#define RS_COSET_LOCATOR_MAX_LEN 17

/**
 * @brief Context data.
 */
typedef struct RS {
    GF_t *gf;
    CC_t *cc;
    FFT_t *fft;
} RS_t;

/**
 * @brief Allocate memory.
 *
 * @param rs where to place allocated data.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int rs_alloc(RS_t *rs);

/**
 * @brief Make necessary pre-calculations.
 *
 * @param rs context object.
 * @param gf Galois field data
 * @param cc cyclotomic cosets data.
 * @param fft fft data.
 */
void rs_init(RS_t *rs, GF_t *gf, CC_t *cc, FFT_t *fft);

/**
 * @brief Generate repair symbols for the given information symbols.
 *
 * @param rs context object.
 * @param inf_symbols information symbols.
 * @param rep_symbols where to place the result.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int rs_generate_repair_symbols(const RS_t *rs, symbol_seq_t inf_symbols,
                               symbol_seq_t rep_symbols);

#endif