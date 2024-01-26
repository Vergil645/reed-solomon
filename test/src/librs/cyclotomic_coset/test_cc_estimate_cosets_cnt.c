#include <stdio.h>

#include <cyclotomic_coset.h>

#define TEST_WRAPPER(_k, _r, _inf_max_cnt_lb, _rep_max_cnt_lb)                 \
    do {                                                                       \
        if (test((_k), (_r), (_inf_max_cnt_lb), (_rep_max_cnt_lb))) {          \
            return 1;                                                          \
        }                                                                      \
    } while (0)

static int test(uint16_t k, uint16_t r, uint16_t inf_max_cnt_lb,
                uint16_t rep_max_cnt_lb) {
    uint16_t inf_max_cnt;
    uint16_t rep_max_cnt;

    cc_estimate_cosets_cnt(k, r, &inf_max_cnt, &rep_max_cnt);

    if (inf_max_cnt < inf_max_cnt_lb) {
        printf("ERROR: cc_estimate_cosets_cnt(%u, %u, ...): inf_max_cnt = %u < "
               "%u\n",
               k, r, inf_max_cnt, inf_max_cnt_lb);
        return 1;
    }

    if (rep_max_cnt < rep_max_cnt_lb) {
        printf("ERROR: cc_estimate_cosets_cnt(%u, %u, ...): rep_max_cnt = %u < "
               "%u\n",
               k, r, rep_max_cnt, rep_max_cnt_lb);
        return 1;
    }

    return 0;
}

int main(void) {
    TEST_WRAPPER(19, 0, 5, 0);
    TEST_WRAPPER(255, 0, 35, 0);
    TEST_WRAPPER(389, 0, 42, 0);
    TEST_WRAPPER(16, 3, 3, 2);
    TEST_WRAPPER(11, 11, 2, 4);
    TEST_WRAPPER(19, 18, 3, 4);

    // Duplicating a test to check that a function is deterministic
    TEST_WRAPPER(1034, 389, 66, 42);
    TEST_WRAPPER(1034, 389, 66, 42);

    return 0;
}