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
#include <string.h>

#include <cyclotomic_coset.h>

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
    } while (0)

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

/**
 * @brief  Number of leaders of cyclotomic cosets of a certain size.
 */
static const uint16_t g_leaders_cnt[CC_COSET_SIZES_CNT] = {
    [0] = CC_LEADERS_1_CNT,
    [1] = CC_LEADERS_2_CNT,
    [2] = CC_LEADERS_4_CNT,
    [3] = CC_LEADERS_8_CNT,
    [4] = CC_LEADERS_16_CNT};

/**
 * @brief Thresholds for using cyclotomic cosets of a certain size.
 */
static const uint16_t g_thresholds[CC_COSET_SIZES_CNT] = {[0] = CC_THRESHOLD_1,
                                                          [1] = CC_THRESHOLD_2,
                                                          [2] = CC_THRESHOLD_4,
                                                          [3] = CC_THRESHOLD_8,
                                                          [4] =
                                                              CC_THRESHOLD_16};

int cc_alloc(CC_t *cc) {
    assert(cc != NULL);

    uint16_t **leaders = cc->leaders;
    uint16_t *leaders_ptr;

    leaders_ptr = (uint16_t *)malloc(CC_COSETS_CNT * sizeof(uint16_t));
    if (!leaders_ptr) {
        return 1;
    }

    leaders[0] = leaders_ptr;
    for (uint8_t i = 1; i < CC_COSET_SIZES_CNT; ++i) {
        leaders[i] = leaders[i - 1] + g_leaders_cnt[i - 1];
    }

    return 0;
}

int cc_init(CC_t *cc) {
    assert(cc != NULL);

    bool *processed;
    uint16_t **leaders = cc->leaders;
    uint16_t idx[CC_COSET_SIZES_CNT] = {0}; // idx[i] - index in leaders[1 << i]
    uint16_t coset_size;
    uint16_t cur_coset_elem;
    uint8_t i;

    processed = (bool *)malloc(N * sizeof(bool));
    if (!processed) {
        return 1;
    }
    memset(processed, false, N * sizeof(bool));

    for (uint16_t s = 0; s < N; ++s) {
        if (processed[s]) {
            continue;
        }
        processed[s] = true;

        cur_coset_elem = NEXT_COSET_ELEMENT(s);
        for (coset_size = 1; coset_size <= CC_MAX_COSET_SIZE; ++coset_size) {
            if (cur_coset_elem == s) {
                break;
            }

            processed[cur_coset_elem] = true;
            cur_coset_elem = NEXT_COSET_ELEMENT(cur_coset_elem);
        }

        i = 0;
        while (i < CC_COSET_SIZES_CNT && (1 << i) != coset_size) {
            ++i;
        }

        assert(i != CC_COSET_SIZES_CNT);
        assert(idx[i] < g_leaders_cnt[i]);

        leaders[i][idx[i]++] = s;
    }

    free(processed);

    assert(idx[0] == CC_LEADERS_1_CNT);
    assert(idx[1] == CC_LEADERS_2_CNT);
    assert(idx[2] == CC_LEADERS_4_CNT);
    assert(idx[3] == CC_LEADERS_8_CNT);
    assert(idx[4] == CC_LEADERS_16_CNT);

    return 0;
}

void cc_free(CC_t *cc) {
    assert(cc != NULL);

    free(cc->leaders[0]); // see cc_init(...)
}

_static uint16_t _cc_get_cosets_cnt(uint16_t r) {
    uint16_t cosets_cnt = 0;
    uint16_t cosets_cnt_inc;

    for (uint8_t i = CC_COSET_SIZES_CNT - 1; r != 0; --i) {
        if (r > g_thresholds[i]) {
            cosets_cnt_inc = (r - g_thresholds[i] + (1 << i) - 1) >> i;
            cosets_cnt += cosets_cnt_inc;
            r -= cosets_cnt_inc << i;
        }

        if (i == 0) {
            break;
        }
    }

    assert(r == 0);

    return cosets_cnt;
}

void cc_estimate_cosets_cnt(uint16_t k, uint16_t r, uint16_t *inf_max_cnt,
                            uint16_t *rep_max_cnt) {
    *inf_max_cnt = _cc_get_cosets_cnt(k); // upper limit
    *rep_max_cnt = _cc_get_cosets_cnt(r); // accurate estimation
}

void cc_select_cosets(const CC_t *cc, uint16_t k, uint16_t r,
                      coset_t *inf_cosets, uint16_t inf_max_cnt,
                      uint16_t *inf_cosets_cnt, coset_t *rep_cosets,
                      uint16_t rep_max_cnt, uint16_t *rep_cosets_cnt) {
    assert(cc != NULL);
    assert(k + r <= N);
    assert(inf_cosets != NULL);
    assert(rep_cosets != NULL);

    uint16_t *const *leaders = cc->leaders;
    uint16_t idx[CC_COSET_SIZES_CNT] = {0}; // idx[i] - index in leaders[1 << i]
    uint16_t inf_thresholds[CC_COSET_SIZES_CNT] = {
        0}; // inf_threshold[i] - threshold for
            // using cyclotomic cosets of size
            // (1 << i) for information symbols
    uint16_t inf_idx = 0;
    uint16_t rep_idx = 0;

    for (uint8_t i = CC_COSET_SIZES_CNT - 1; r != 0; --i) {
        while (r > g_thresholds[i] && rep_idx < rep_max_cnt) {
            assert(idx[i] < g_leaders_cnt[i]);

            INIT_COSET(rep_cosets[rep_idx++], leaders[i][idx[i]++], 1 << i);
            r -= 1 << i;
        }

        if (i == 0) {
            break;
        }
    }

    assert(r == 0);

    *rep_cosets_cnt = rep_idx;

    for (uint8_t i = 0; i < CC_COSET_SIZES_CNT; ++i) {
        inf_thresholds[i] = g_thresholds[i];
    }
    for (uint8_t i = 0; i < CC_COSET_SIZES_CNT - 1; ++i) {
        for (uint8_t j = i + 1; j < CC_COSET_SIZES_CNT; ++j) {
            inf_thresholds[j] -= idx[i] << i;
        }
    }

    for (uint8_t i = CC_COSET_SIZES_CNT - 1; k != 0; --i) {
        while (k > inf_thresholds[i] && inf_idx < inf_max_cnt) {
            assert(idx[i] < g_leaders_cnt[i]);

            INIT_COSET(inf_cosets[inf_idx++], leaders[i][idx[i]++], 1 << i);
            k -= MIN(k, 1 << i);
        }

        if (i == 0) {
            break;
        }
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