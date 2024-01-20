/**
 * @file cyclotomic_coset.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains functions that uses cyclotomic cosets over GF(2) modulo N.
 * @date 2024-01-17
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_CYCLOTOMIC_COSET_H__
#define __REED_SOLOMON_CYCLOTOMIC_COSET_H__

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "constants.h"

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
 * @brief Leaders of cyclotomic cosets of size 1.
 */
#define CC_LEADERS_1                                                           \
    { 0 }

/**
 * @brief Leaders of cyclotomic cosets of size 2.
 */
#define CC_LEADERS_2                                                           \
    { 21845 }

/**
 * @brief Leaders of cyclotomic cosets of size 4.
 */
#define CC_LEADERS_4                                                           \
    { 4369, 13107, 30583 }

/**
 * @brief Leaders of cyclotomic cosets of size 8.
 */
#define CC_LEADERS_8                                                           \
    {                                                                          \
        257, 771, 1285, 1799, 2313, 2827, 3341, 3855, 4883, 5397, 5911, 6425,  \
            6939, 7453, 7967, 9509, 10023, 11051, 11565, 12079, 13621, 14135,  \
            15163, 15677, 16191, 22359, 23387, 24415, 28527, 32639             \
    }

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
 * @brief Minimum of 2 arguments.
 */
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/**
 * @brief Cyclotomic coset over GF(2) modulo N.
 */
typedef struct coset {
    uint16_t leader;
    uint16_t size;
} coset_t;

/**
 * @brief Cyclotomic cosets over GF(2) modulo N pre-computed data.
 */
typedef struct CC {
    /**
     * @brief Leaders of cyclotomic cosets of size 1.
     */
    uint16_t leaders_1[CC_LEADERS_1_CNT];

    /**
     * @brief Leaders of cyclotomic cosets of size 2.
     */
    uint16_t leaders_2[CC_LEADERS_2_CNT];

    /**
     * @brief Leaders of cyclotomic cosets of size 4.
     */
    uint16_t leaders_4[CC_LEADERS_4_CNT];

    /**
     * @brief Leaders of cyclotomic cosets of size 8.
     */
    uint16_t leaders_8[CC_LEADERS_8_CNT];

    /**
     * @brief Leaders of cyclotomic cosets of size 16.
     */
    uint16_t *leaders_16;
} CC_t;

/**
 * @brief Allocate memory for cyclotomic cosets data.
 *
 * @param cc where to place allocated data.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int cc_alloc(CC_t *cc) {
    assert(cc != NULL);

    uint16_t *leaders_16;

    leaders_16 = (uint16_t *)malloc(CC_LEADERS_16_CNT * sizeof(uint16_t));
    if (!leaders_16) {
        return 1;
    }

    *cc = {
        .leaders_1 = CC_LEADERS_1,
        .leaders_2 = CC_LEADERS_2,
        .leaders_4 = CC_LEADERS_4,
        .leaders_8 = CC_LEADERS_8,
        .leaders_16 = leaders_16,
    };
    return 0;
}

/**
 * @brief Make necessary pre-calculations.
 * @details Pre-calculations list:
 * 1. find leaders of cyclotomic cosets of size 8;
 * 2. find leaders of cyclotomic cosets of size 16.
 *
 * @param cc cyclotomic cosets data.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int cc_init(CC_t *cc) {
    assert(cc != NULL);

    bool *processed;
    uint16_t *leaders_16 = cc->leaders_16;
    uint16_t idx_16 = 0;
    uint16_t coset_size;
    uint16_t cur_coset_elem;

    processed = (bool *)malloc(N * sizeof(bool));
    if (!processed) {
        return 1;
    }

    for (uint16_t s = 1; s < N; ++s) {
        if (processed[s]) {
            continue;
        }
        processed[s] = true;

        cur_coset_elem = (s << 1) % N;
        for (coset_size = 1; cur_coset_elem != s; ++coset_size) {
            processed[cur_coset_elem] = true;
            cur_coset_elem = (cur_coset_elem << 1) % N;
        }

        if (coset_size == 16) {
            leaders_16[idx_16++] = s;
        }
    }

    free(processed);

    assert(idx_16 == CC_LEADERS_16_CNT);

    return 0;
}

/**
 * @brief Deallocates cyclotomic cosets data.
 *
 * @param cc cyclotomic cosets data.
 */
void cc_free(CC_t *cc) {
    assert(cc != NULL);

    free(cc->leaders_16);
}

/**
 * @brief Calculate a number of cyclotomic cosets the union of which has a given
 * size.
 *
 * @param r union size.
 * @return number of cyclotomic cosets.
 */
uint16_t cc_get_cosets_cnt(uint16_t r) {
    uint16_t cosets_cnt = 0;
    uint16_t cosets_cnt_inc;

    if (r > CC_THRESHOLD_16) {
        cosets_cnt_inc = (r - CC_THRESHOLD_16 + 15) >> 4;
        cosets_cnt += cosets_cnt_inc;
        r -= cosets_cnt_inc << 4;
    }

    if (r > CC_THRESHOLD_8) {
        cosets_cnt_inc = (r - CC_THRESHOLD_8 + 7) >> 3;
        cosets_cnt += cosets_cnt_inc;
        r -= cosets_cnt_inc << 3;
    }

    if (r > CC_THRESHOLD_4) {
        cosets_cnt_inc = (r - CC_THRESHOLD_4 + 3) >> 2;
        cosets_cnt += cosets_cnt_inc;
        r -= cosets_cnt_inc << 2;
    }

    if (r > CC_THRESHOLD_2) {
        cosets_cnt += 1;
        r -= 2;
    }

    if (r > CC_THRESHOLD_1) {
        cosets_cnt += 1;
        r -= 1;
    }

    assert(r == 0);

    return cosets_cnt;
}

/**
 * @brief Estimate upper limits on the number of cyclotomic cosets that will be
 * selected by cc_select_cosets(...).
 *
 * @param k number of information symbols.
 * @param r number of repair symbols.
 * @param inf_max_cnt upper limit on the number of information symbol cosets.
 * @param rep_max_cnt upper limit on the number of repair symbol cosets.
 */
void cc_estimate_cosets_cnt(uint16_t k, uint16_t r, uint16_t *inf_max_cnt,
                            uint16_t *rep_max_cnt) {
    *inf_max_cnt = cc_get_cosets_cnt(k); // upper limit
    *rep_max_cnt = cc_get_cosets_cnt(r); // accurate estimation
}

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
void cc_select_cosets(const CC_t *cc, uint16_t k, uint16_t r,
                      coset_t *inf_cosets, uint16_t inf_max_cnt,
                      uint16_t *inf_cosets_cnt, coset_t *rep_cosets,
                      uint16_t rep_max_cnt, uint16_t *rep_cosets_cnt) {
    assert(k + r <= N);
    assert(inf_cosets != NULL);
    assert(rep_cosets != NULL);

    uint16_t leaders_idx[5] = {0}; // leaders_idx[i] - index in leaders_<2^i>
    uint16_t threshold[5] = {
        CC_THRESHOLD_1, CC_THRESHOLD_2,  CC_THRESHOLD_4,
        CC_THRESHOLD_8, CC_THRESHOLD_16,
    }; // threshold[i] - threshold for using cyclotomic cosets of size 2^i for
       // information symbols
    uint16_t inf_idx = 0;
    uint16_t rep_idx = 0;

    while (r > CC_THRESHOLD_16 && rep_idx < rep_max_cnt) {
        assert(leaders_idx[4] < CC_LEADERS_16_CNT);
        rep_cosets[rep_idx++] = {
            .leader = cc->leaders_16[leaders_idx[4]++],
            .size = 16,
        };
        r -= 16;
    }

    while (r > CC_THRESHOLD_8 && rep_idx < rep_max_cnt) {
        assert(leaders_idx[3] < CC_LEADERS_8_CNT);
        rep_cosets[rep_idx++] = {
            .leader = cc->leaders_8[leaders_idx[3]++],
            .size = 8,
        };
        threshold[4] -= 8;
        r -= 8;
    }

    while (r > CC_THRESHOLD_4 && rep_idx < rep_max_cnt) {
        assert(leaders_idx[2] < CC_LEADERS_4_CNT);
        rep_cosets[rep_idx++] = {
            .leader = cc->leaders_4[leaders_idx[2]++],
            .size = 4,
        };
        threshold[4] -= 4;
        threshold[3] -= 4;
        r -= 4;
    }

    if (r > CC_THRESHOLD_2 && rep_idx < rep_max_cnt) {
        rep_cosets[rep_idx++] = {
            .leader = cc->leaders_2[leaders_idx[1]++],
            .size = 2,
        };
        threshold[4] -= 2;
        threshold[3] -= 2;
        threshold[2] -= 2;
        r -= 2;
    }

    if (r > CC_THRESHOLD_1 && rep_idx < rep_max_cnt) {
        rep_cosets[rep_idx++] = {
            .leader = cc->leaders_1[leaders_idx[0]++],
            .size = 1,
        };
        threshold[4] -= 1;
        threshold[3] -= 1;
        threshold[2] -= 1;
        threshold[1] -= 1;
        r -= 1;
    }

    assert(r == 0);
    *rep_cosets_cnt = rep_idx;

    while (k > threshold[4] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[4] < CC_LEADERS_16_CNT);
        inf_cosets[inf_idx++] = {
            .leader = cc->leaders_16[leaders_idx[4]++],
            .size = 16,
        };
        k -= MIN(k, 16);
    }

    while (k > threshold[3] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[3] < CC_LEADERS_8_CNT);
        inf_cosets[inf_idx++] = {
            .leader = cc->leaders_8[leaders_idx[3]++],
            .size = 8,
        };
        k -= MIN(k, 8);
    }

    while (k > threshold[2] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[2] < CC_LEADERS_4_CNT);
        inf_cosets[inf_idx++] = {
            .leader = cc->leaders_4[leaders_idx[2]++],
            .size = 4,
        };
        k -= MIN(k, 4);
    }

    if (k > threshold[1] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[1] < CC_LEADERS_2_CNT);
        inf_cosets[inf_idx++] = {
            .leader = cc->leaders_2[leaders_idx[1]++],
            .size = 2,
        };
        k -= MIN(k, 2);
    }

    if (k > threshold[0] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[0] < CC_LEADERS_1_CNT);
        inf_cosets[inf_idx++] = {
            .leader = cc->leaders_1[leaders_idx[0]++],
            .size = 1,
        };
        k -= MIN(k, 1);
    }

    assert(k == 0);
    *inf_cosets_cnt = inf_idx;
}

/**
 * @brief Convert list of cyclotomic cosets to list of symbol positions.
 *
 * @param cosets cyclotomic cosets.
 * @param cosets_cnt number of cyclotomic cosets.
 * @param positions where to place symbol positions.
 * @param positions_cnt required number of positions to be written.
 */
void cc_cosets_to_positions(const coset_t *cosets, uint16_t cosets_cnt,
                            uint16_t *positions, uint16_t positions_cnt) {
    assert(cosets != NULL);
    assert(positions != NULL);

    uint16_t s;
    uint16_t cur_coset_elem;
    uint16_t positions_idx = 0;

    for (uint16_t i = 0; i < cosets_cnt && positions_idx < positions_cnt; ++i) {
        s = cosets[i].leader;
        positions[positions_idx++] = s;
        cur_coset_elem = (s << 1) % N;

        while (cur_coset_elem != s && positions_idx < positions_cnt) {
            positions[positions_idx++] = cur_coset_elem;
            cur_coset_elem = (cur_coset_elem << 1) % N;
        }
    }

    assert(positions_idx == positions_cnt);
}

#endif