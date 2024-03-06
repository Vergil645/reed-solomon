/**
 * @file cyclotomic_coset.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains functions that uses cyclotomic cosets over GF(2) modulo N
 * (65535).
 * @date 2024-01-17
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_CYCLOTOMIC_COSET_H__
#define __REED_SOLOMON_CYCLOTOMIC_COSET_H__

#include <stdint.h>

/**
 * @brief Number of different sizes of cyclotomic cosets.
 */
#define CC_COSET_SIZES_CNT 5

/**
 * @brief Number of different cyclotomic cosets.
 */
#define CC_COSETS_CNT 4115

/**
 * @brief Maximum cyclotomic coset size.
 */
#define CC_MAX_COSET_SIZE 16

/**
 * @brief Number of leaders of cyclotomic cosets of size 1.
 */
#define CC_LEADERS_1_CNT 1

/**
 * @brief Number of leaders of cyclotomic cosets of size 2.
 */
#define CC_LEADERS_2_CNT 1

/**
 * @brief Number of leaders of cyclotomic cosets of size 4.
 */
#define CC_LEADERS_4_CNT 3

/**
 * @brief Number of leaders of cyclotomic cosets of size 8.
 */
#define CC_LEADERS_8_CNT 30

/**
 * @brief Number of leaders of cyclotomic cosets of size 16.
 */
#define CC_LEADERS_16_CNT 4080

/**
 * @brief If (r > *value*), we have to use cyclotomic cosets of size 1.
 */
#define CC_THRESHOLD_1 0

/**
 * @brief If (r > *value*), we have to use cyclotomic cosets of size 2.
 */
#define CC_THRESHOLD_2 1

/**
 * @brief If (r > *value*), we have to use cyclotomic cosets of size 4.
 */
#define CC_THRESHOLD_4 3

/**
 * @brief If (r > *value*), we have to use cyclotomic cosets of size 8.
 */
#define CC_THRESHOLD_8 15

/**
 * @brief If (r > *value*), we have to use cyclotomic cosets of size 16.
 */
#define CC_THRESHOLD_16 255

/**
 * @brief Next cyclotomic coset element.
 * @details \f$s_{i+1} = \equiv s_i * 2 \mod N\f$
 *
 * @param _s current coset element.
 * @return cyclotomic coset element following _s.
 */
#define NEXT_COSET_ELEMENT(_s) ((uint16_t)(((uint32_t)(_s) << 1) % 65535))

/**
 * @brief Cyclotomic coset over GF(2) modulo N.
 */
typedef struct coset {
    uint16_t leader;
    uint8_t size;
} coset_t;

/**
 * @brief Cyclotomic cosets over GF(2) modulo N pre-computed data.
 */
typedef struct CC {
    /**
     * @brief Leaders of cyclotomic cosets.
     * @details leaders[i] - leaders of cyclotomic cosets of size \f$2^i\f$.\n
     * All members are pointers to different parts of one allocated memory
     * fragment. leaders[0] - pointer to the beginning of this fragment.
     */
    uint16_t* leaders[CC_COSET_SIZES_CNT];

    /**
     * @brief .leaders memory.
     */
    uint16_t _leaders_memory[CC_COSETS_CNT];
} CC_t;

/**
 * @brief Create cyclotomic cosets data structure.
 *
 * @return pointer to created cyclotomic cosets data structure on success and
 * NULL otherwise.
 */
CC_t* cc_create();

/**
 * @brief Destroy cyclotomic cosets data structure.
 *
 * @param cc cyclotomic cosets data.
 */
void cc_destroy(CC_t* cc);

/**
 * @brief Estimate upper limits on the number of cyclotomic cosets that will be
 * selected by cc_select_cosets(...).
 *
 * @param k number of information symbols.
 * @param r number of repair symbols.
 * @param inf_max_cnt upper limit on the number of information symbol cosets.
 * @param rep_max_cnt upper limit on the number of repair symbol cosets.
 */
void cc_estimate_cosets_cnt(uint16_t k, uint16_t r, uint16_t* inf_max_cnt, uint16_t* rep_max_cnt);

/**
 * @brief Select cyclotomic cosets over GF(2) modulo N that form information and
 * repair symbol positions in virtual codeword.
 *
 * @param cc cyclotomic cosets data.
 * @param k number of information symbols.
 * @param r number of repair symbols.
 * @param inf_cosets where to place information symbol cosets.
 * @param inf_max_cnt max number of elements that can be written to inf_leaders.
 * @param inf_cosets_cnt where to place number of written information symbol
 * cosets.
 * @param rep_cosets where to place repair symbol cosets.
 * @param rep_max_cnt max number of elements that can be written to rep_leaders.
 * @param rep_cosets_cnt where to place number of written repair symbol cosets.
 */
void cc_select_cosets(CC_t* cc, uint16_t k, uint16_t r, coset_t* inf_cosets, uint16_t inf_max_cnt,
                      uint16_t* inf_cosets_cnt, coset_t* rep_cosets, uint16_t rep_max_cnt,
                      uint16_t* rep_cosets_cnt);

/**
 * @brief Convert list of cyclotomic cosets to list of symbol positions.
 *
 * @param cosets cyclotomic cosets.
 * @param cosets_cnt number of cyclotomic cosets.
 * @param positions where to place symbol positions.
 * @param positions_cnt required number of positions to be written.
 */
void cc_cosets_to_positions(const coset_t* cosets, uint16_t cosets_cnt, uint16_t* positions,
                            uint16_t positions_cnt);

#endif