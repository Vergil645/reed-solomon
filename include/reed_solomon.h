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

#include <cyclotomic_coset.h>
#include <fft.h>
#include <prelude.h>
#include <seq.h>
#include <symbol.h>

/**
 * @brief Maximum number of cyclotomic coset locator polynomial coefficients.
 */
#define RS_COSET_LOCATOR_MAX_LEN 17

/**
 * @brief Return code for cases when erases cannot be restored due to code
 * parameters.
 */
#define RS_ERR_CANNOT_RESTORE 100

/**
 * @brief Context data.
 */
typedef struct RS {
    GF_t* gf;
    CC_t* cc;
    FFT_t* fft;
} RS_t;

/**
 * @brief Allocate memory.
 *
 * @param rs where to place allocated data.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int rs_alloc(RS_t* rs);

/**
 * @brief Make necessary pre-calculations.
 *
 * @param rs context object.
 * @param gf Galois field data
 * @param cc cyclotomic cosets data.
 * @param fft fft data.
 */
void rs_init(RS_t* rs, GF_t* gf, CC_t* cc, FFT_t* fft);

/**
 * @brief Deallocate context object memory.
 *
 * @param rs context object.
 */
void rs_free(RS_t* rs);

/**
 * @brief Generate repair symbols for the given information symbols.
 *
 * @param rs context object.
 * @param inf_symbols information symbols.
 * @param rep_symbols where to place the result.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int rs_generate_repair_symbols(const RS_t* rs, symbol_seq_t inf_symbols,
                               symbol_seq_t rep_symbols);

/**
 * @brief Restore erased symbols. Assume that rcv_symbols[i] = 0 for erased i.
 *
 * @param rs context object.
 * @param k number of information symbols.
 * @param r number of repair symbols.
 * @param rcv_symbols received symbols, restored symbols will be written here.
 * @param erased_indices indices of erased symbols in rcv_symbols.
 * @param t number of erases.
 * @return 0 on success,\n
 *         1 on memory allocation error,\n
 *         or RS_ERR_CANNOT_RESTORE.
 */
int rs_restore_symbols(const RS_t* rs, uint16_t k, uint16_t r,
                       symbol_seq_t rcv_symbols, const uint16_t* erased_indices,
                       uint16_t t);

// Internal functions not for public use.
// Defined here only in Debug mode in order to test them.
#ifndef NDEBUG

/**
 * @brief Compute syndrome polynomial.
 *
 * @param rs context object.
 * @param seq symbol sequence.
 * @param positions symbol positions.
 * @param syndrome_poly where to place the result.
 */
void _rs_get_syndrome_poly(const RS_t* rs, symbol_seq_t seq,
                           const uint16_t* positions,
                           symbol_seq_t syndrome_poly);

/**
 * @brief Compute locator polynomial.
 *
 * @param rs context object.
 * @param positions positions.
 * @param positions_cnt number of positions.
 * @param locator_poly where to place locator polynomial coefficients.
 * @param locator_max_len max number of locator polynomial coefficients.
 */
void _rs_get_locator_poly(const RS_t* rs, const uint16_t* positions,
                          uint16_t positions_cnt, element_t* locator_poly,
                          uint16_t locator_max_len);

/**
 * @brief Compute repair symbols locator polynomial.
 * @details All locator polynomial coefficients will belongs to GF(2) subfield
 * ({0, 1}) of GF(65536). Locator polynomial will have degree equal to number of
 * repair symbols.
 *
 * @param rs context object.
 * @param r number of repair symbols.
 * @param rep_cosets cyclotomic cosets that form repair symbol positions.
 * @param rep_cosets_cnt number of cyclotomic cosets.
 * @param locator_poly where to place locator polynomial coefficients.
 * @param locator_max_len max number of locator polynomial coefficients.
 */
void _rs_get_rep_symbols_locator_poly(const RS_t* rs, uint16_t r,
                                      const coset_t* rep_cosets,
                                      uint16_t rep_cosets_cnt,
                                      element_t* locator_poly,
                                      uint16_t locator_max_len);

/**
 * @brief Compute Forney coefficient for given symbol postion.
 *
 * @param rs context object.
 * @param locator_poly locator polynomial.
 * @param d degree of locator polynomial (number of repair symbols or erasures).
 * @param pos symbol position.
 * @return Forney coefficient.
 */
element_t _rs_get_forney_coef(const RS_t* rs, const element_t* locator_poly,
                              uint16_t d, uint16_t pos);

/**
 * @brief Compute evaluator polynomial modulo x^t (t - number of repair symbols
 * or erasures).
 *
 * @param rs context object.
 * @param syndrome_poly information symbols syndrome polynomial (deg == t - 1).
 * @param locator_poly repair symbols locator polynomial (deg == t).
 * @param evaluator_poly where to place the result.
 */
void _rs_get_evaluator_poly(const RS_t* rs, symbol_seq_t syndrome_poly,
                            const element_t* locator_poly,
                            symbol_seq_t evaluator_poly);

/**
 * @brief Compute repair symbols.
 *
 * @param rs context object.
 * @param locator_poly repair symbols locator polynomial.
 * @param evaluator_poly repair symbols evaluator polynomial.
 * @param rep_symbols where to place the result.
 * @param rep_positions repair symbol positions.
 */
void _rs_get_repair_symbols(const RS_t* rs, const element_t* locator_poly,
                            symbol_seq_t evaluator_poly,
                            symbol_seq_t rep_symbols,
                            const uint16_t* rep_positions);

/**
 * @brief Restore erased symbols if it is possible.
 *
 * @param rs context object.
 * @param locator_poly erased symbols locator polynomial.
 * @param evaluator_poly erased symbols evaluator polynomial.
 * @param positions positions of all symbols.
 * @param rcv_symbols received symbols, restored symbols will be written here.
 * @param erased_indices indices of erased symbols in rcv_symbols.
 */
void _rs_restore_erased(const RS_t* rs, const element_t* locator_poly,
                        symbol_seq_t evaluator_poly, const uint16_t* positions,
                        symbol_seq_t rcv_symbols,
                        const uint16_t* erased_indices);

#endif

#endif