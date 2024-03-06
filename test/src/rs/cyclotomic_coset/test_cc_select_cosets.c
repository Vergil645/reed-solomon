#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <rs/cyclotomic_coset.h>

#define TEST_WRAPPER(_cc, _k, _r, _inf_cosets, _inf_cosets_cnt, _rep_cosets, _rep_cosets_cnt)      \
    do {                                                                                           \
        if (test((_cc), (_k), (_r), (_inf_cosets), (_inf_cosets_cnt), (_rep_cosets),               \
                 (_rep_cosets_cnt))) {                                                             \
            cc_destroy((_cc));                                                                     \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)

static bool cosets_eq(const coset_t* a, uint16_t a_len, const coset_t* b, uint16_t b_len) {
    if (a == NULL || b == NULL || a_len != b_len) {
        return false;
    }

    for (uint16_t i = 0; i < a_len; ++i) {
        if (a[i].leader != b[i].leader || a[i].size != b[i].size)
            return false;
    }

    return true;
}

static void cosets_printf(const coset_t* cosets, uint16_t cosets_cnt) {
    if (cosets == NULL) {
        printf("NULL");
    } else if (cosets_cnt == 0) {
        printf("[]");
    } else {
        printf("[");
        for (uint16_t i = 0; i < cosets_cnt - 1; ++i) {
            printf("{%u, %u}, ", cosets[i].leader, cosets[i].size);
        }
        printf("{%u, %u}]", cosets[cosets_cnt - 1].leader, cosets[cosets_cnt - 1].size);
    }
}

static int test(CC_t* cc, uint16_t k, uint16_t r, const coset_t* inf_cosets,
                uint16_t inf_cosets_cnt, const coset_t* rep_cosets, uint16_t rep_cosets_cnt) {
    coset_t* _inf_cosets;
    coset_t* _rep_cosets;
    uint16_t _inf_cosets_cnt;
    uint16_t _rep_cosets_cnt;
    int ret = 0;

    _inf_cosets = (coset_t*)calloc(inf_cosets_cnt, sizeof(coset_t));
    if (!_inf_cosets) {
        printf("ERROR: cannot allocate memory for _inf_cosets\n");
        return 1;
    }

    _rep_cosets = (coset_t*)calloc(rep_cosets_cnt, sizeof(coset_t));
    if (!_rep_cosets) {
        printf("ERROR: cannot allocate memory for _rep_cosets\n");
        free(_inf_cosets);
        return 1;
    }

    cc_select_cosets(cc, k, r, _inf_cosets, inf_cosets_cnt, &_inf_cosets_cnt, _rep_cosets,
                     rep_cosets_cnt, &_rep_cosets_cnt);

    if (!cosets_eq(_inf_cosets, _inf_cosets_cnt, inf_cosets, inf_cosets_cnt)) {
        printf("ERROR: cc_select_cosets(*, %u, %u, ...): incorrect inf_cosets:\n", k, r);

        printf("\tcorrect = ");
        cosets_printf(inf_cosets, inf_cosets_cnt);
        printf("\n");

        printf("\tactual  = ");
        cosets_printf(_inf_cosets, _inf_cosets_cnt);
        printf("\n");

        ret = 1;
    } else if (!cosets_eq(_rep_cosets, _rep_cosets_cnt, rep_cosets, rep_cosets_cnt)) {
        printf("ERROR: cc_select_cosets(*, %u, %u, ...): incorrect rep_cosets:\n", k, r);

        printf("\tcorrect = ");
        cosets_printf(rep_cosets, rep_cosets_cnt);
        printf("\n");

        printf("\tactual  = ");
        cosets_printf(_rep_cosets, _rep_cosets_cnt);
        printf("\n");

        ret = 1;
    }

    free(_rep_cosets);
    free(_inf_cosets);

    return ret;
}

int main(void) {
    CC_t* cc;

    cc = cc_create();
    if (!cc) {
        printf("ERROR: cc_create returned NULL\n");
        return 1;
    }

    // Test 1
    {
        uint16_t k = 16;
        uint16_t r = 3;
        uint16_t inf_cosets_cnt = 3;
        uint16_t rep_cosets_cnt = 2;
        coset_t inf_cosets[] = {
            {.leader = 257, .size = 8},
            {.leader = 4369, .size = 4},
            {.leader = 13107, .size = 4},
        };
        coset_t rep_cosets[] = {
            {.leader = 21845, .size = 2},
            {.leader = 0, .size = 1},
        };

        TEST_WRAPPER(cc, k, r, inf_cosets, inf_cosets_cnt, rep_cosets, rep_cosets_cnt);
    }

    // Test 2
    {
        uint16_t k = 11;
        uint16_t r = 11;
        uint16_t inf_cosets_cnt = 2;
        uint16_t rep_cosets_cnt = 4;
        coset_t inf_cosets[] = {
            {.leader = 257, .size = 8},
            {.leader = 30583, .size = 4},
        };
        coset_t rep_cosets[] = {
            {.leader = 4369, .size = 4},
            {.leader = 13107, .size = 4},
            {.leader = 21845, .size = 2},
            {.leader = 0, .size = 1},
        };

        TEST_WRAPPER(cc, k, r, inf_cosets, inf_cosets_cnt, rep_cosets, rep_cosets_cnt);
    }

    // Test 3
    {
        uint16_t k = 19;
        uint16_t r = 18;
        uint16_t inf_cosets_cnt = 3;
        uint16_t rep_cosets_cnt = 4;
        coset_t inf_cosets[] = {
            {.leader = 771, .size = 8},
            {.leader = 1285, .size = 8},
            {.leader = 30583, .size = 4},
        };
        coset_t rep_cosets[] = {
            {.leader = 257, .size = 8},
            {.leader = 4369, .size = 4},
            {.leader = 13107, .size = 4},
            {.leader = 21845, .size = 2},
        };

        TEST_WRAPPER(cc, k, r, inf_cosets, inf_cosets_cnt, rep_cosets, rep_cosets_cnt);
    }

    // Test 4
    {
        uint16_t k = 22;
        uint16_t r = 17;
        uint16_t inf_cosets_cnt = 4;
        uint16_t rep_cosets_cnt = 4;
        coset_t inf_cosets[] = {
            {.leader = 771, .size = 8},
            {.leader = 1285, .size = 8},
            {.leader = 30583, .size = 4},
            {.leader = 21845, .size = 2},
        };
        coset_t rep_cosets[] = {
            {.leader = 257, .size = 8},
            {.leader = 4369, .size = 4},
            {.leader = 13107, .size = 4},
            {.leader = 0, .size = 1},
        };

        TEST_WRAPPER(cc, k, r, inf_cosets, inf_cosets_cnt, rep_cosets, rep_cosets_cnt);
    }

    cc_destroy(cc);

    return 0;
}