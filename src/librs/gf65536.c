/**
 * @file gf65536.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief gf65536.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdlib.h>

#include <gf65536.h>

int gf_alloc(GF_t* gf) {
    assert(gf != NULL);

    element_t* pow_table;
    uint16_t* log_table;

    pow_table = (element_t*)malloc(GF_FIELD_SIZE * sizeof(element_t));
    if (!pow_table) {
        return 1;
    }

    log_table = (element_t*)malloc(GF_FIELD_SIZE * sizeof(element_t));
    if (!log_table) {
        free(pow_table);
        return 1;
    }

    gf->pow_table = pow_table;
    gf->log_table = log_table;

    return 0;
}

void gf_init(GF_t* gf) {
    assert(gf != NULL);

    element_t* pow_table = gf->pow_table;
    uint16_t* log_table = gf->log_table;
    poly_t cur_poly = 1;

    for (uint16_t i = 0; i < N; ++i) {
        pow_table[i] = (element_t)cur_poly;
        log_table[(element_t)cur_poly] = i;

        cur_poly <<= 1;
        if (cur_poly & GF_FIELD_SIZE) {
            cur_poly ^= GF_PRIMITIVE_POLY;
        }
    }
}

void gf_free(GF_t* gf) {
    assert(gf != NULL);

    free(gf->log_table);
    free(gf->pow_table);
}

inline element_t gf_mul_ee(const GF_t* gf, element_t a, element_t b) {
    assert(gf != NULL);

    if (a == 0 || b == 0) {
        return 0;
    }

    uint16_t* log_table = gf->log_table;
    uint32_t a_log = (uint32_t)log_table[a];
    uint32_t b_log = (uint32_t)log_table[b];

    return gf->pow_table[(a_log + b_log) % N];
}

inline element_t gf_div_ee(const GF_t* gf, element_t a, element_t b) {
    assert(gf != NULL);
    assert(b != 0);

    if (a == 0) {
        return 0;
    }

    uint16_t* log_table = gf->log_table;
    uint32_t a_log = (uint32_t)log_table[a];
    uint32_t b_log = (uint32_t)log_table[b];

    return gf->pow_table[(N + a_log - b_log) % N];
}