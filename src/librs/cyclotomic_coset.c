/**
 * @file cyclotomic_coset.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief cyclotomic_coset.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#include "cyclotomic_coset.h"

/**
 * @brief Initialize cyclotomic coset by given values.
 *
 * @param _coset cyclotomic coset structure to be initialized.
 * @param _leader cyclotomic coset leader.
 * @param _size cyclotomic coset size.
 *
 * @warning Macro defines local variable with name "__coset_tmp__" in the inner
 * scope.
 */
#define INIT_COSET(_coset, _leader, _size)                                     \
    do {                                                                       \
        coset_t __coset_tmp__ = {                                              \
            .leader = (_leader),                                               \
            .size = (_size),                                                   \
        };                                                                     \
        (_coset) = __coset_tmp__;                                              \
    } while (0);

/**
 * @brief Minimum of 2 values.
 *
 * @param _x first value.
 * @param _y second value.
 * @return min(_x, _y).
 *
 * @warning Expression with side-effects multiply evaluated by macro.
 */
#define MIN(_x, _y) (((_x) < (_y)) ? (_x) : (_y))

int cc_alloc(CC_t *cc) {
    assert(cc != NULL);

    CC_t cc_tmp = {
        .leaders_1 = CC_LEADERS_1,
        .leaders_2 = CC_LEADERS_2,
        .leaders_4 = CC_LEADERS_4,
        .leaders_8 = CC_LEADERS_8,
        .leaders_16 = NULL,
    };

    cc_tmp.leaders_16 =
        (uint16_t *)malloc(CC_LEADERS_16_CNT * sizeof(uint16_t));
    if (!cc_tmp.leaders_16) {
        return 1;
    }

    *cc = cc_tmp;

    return 0;
}

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

        cur_coset_elem = NEXT_COSET_ELEMENT(s);
        for (coset_size = 1; cur_coset_elem != s; ++coset_size) {
            processed[cur_coset_elem] = true;
            cur_coset_elem = NEXT_COSET_ELEMENT(cur_coset_elem);
        }

        if (coset_size == 16) {
            leaders_16[idx_16++] = s;
        }
    }

    free(processed);

    assert(idx_16 == CC_LEADERS_16_CNT);

    return 0;
}

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
static uint16_t cc_get_cosets_cnt(uint16_t r) {
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

void cc_estimate_cosets_cnt(uint16_t k, uint16_t r, uint16_t *inf_max_cnt,
                            uint16_t *rep_max_cnt) {
    *inf_max_cnt = cc_get_cosets_cnt(k); // upper limit
    *rep_max_cnt = cc_get_cosets_cnt(r); // accurate estimation
}

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
        INIT_COSET(rep_cosets[rep_idx++], cc->leaders_16[leaders_idx[4]++], 16);
        r -= 16;
    }

    while (r > CC_THRESHOLD_8 && rep_idx < rep_max_cnt) {
        assert(leaders_idx[3] < CC_LEADERS_8_CNT);
        INIT_COSET(rep_cosets[rep_idx++], cc->leaders_8[leaders_idx[3]++], 8);
        threshold[4] -= 8;
        r -= 8;
    }

    while (r > CC_THRESHOLD_4 && rep_idx < rep_max_cnt) {
        assert(leaders_idx[2] < CC_LEADERS_4_CNT);
        INIT_COSET(rep_cosets[rep_idx++], cc->leaders_4[leaders_idx[2]++], 4);
        threshold[4] -= 4;
        threshold[3] -= 4;
        r -= 4;
    }

    if (r > CC_THRESHOLD_2 && rep_idx < rep_max_cnt) {
        INIT_COSET(rep_cosets[rep_idx++], cc->leaders_2[leaders_idx[1]++], 2);
        threshold[4] -= 2;
        threshold[3] -= 2;
        threshold[2] -= 2;
        r -= 2;
    }

    if (r > CC_THRESHOLD_1 && rep_idx < rep_max_cnt) {
        INIT_COSET(rep_cosets[rep_idx++], cc->leaders_1[leaders_idx[0]++], 1);
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
        INIT_COSET(inf_cosets[inf_idx++], cc->leaders_16[leaders_idx[4]++], 16);
        k -= MIN(k, 16);
    }

    while (k > threshold[3] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[3] < CC_LEADERS_8_CNT);
        INIT_COSET(inf_cosets[inf_idx++], cc->leaders_8[leaders_idx[3]++], 8);
        k -= MIN(k, 8);
    }

    while (k > threshold[2] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[2] < CC_LEADERS_4_CNT);
        INIT_COSET(inf_cosets[inf_idx++], cc->leaders_4[leaders_idx[2]++], 4);
        k -= MIN(k, 4);
    }

    if (k > threshold[1] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[1] < CC_LEADERS_2_CNT);
        INIT_COSET(inf_cosets[inf_idx++], cc->leaders_2[leaders_idx[1]++], 2);
        k -= MIN(k, 2);
    }

    if (k > threshold[0] && inf_idx < inf_max_cnt) {
        assert(leaders_idx[0] < CC_LEADERS_1_CNT);
        INIT_COSET(inf_cosets[inf_idx++], cc->leaders_1[leaders_idx[0]++], 1);
        k -= MIN(k, 1);
    }

    assert(k == 0);

    *inf_cosets_cnt = inf_idx;
}

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

        cur_coset_elem = NEXT_COSET_ELEMENT(s);
        while (cur_coset_elem != s && positions_idx < positions_cnt) {
            positions[positions_idx++] = cur_coset_elem;
            cur_coset_elem = NEXT_COSET_ELEMENT(cur_coset_elem);
        }
    }

    assert(positions_idx == positions_cnt);
}