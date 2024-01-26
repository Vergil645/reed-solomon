#include <stdio.h>

#include <cyclotomic_coset.h>

#define TEST_WRAPPER(_r, _res)                                                 \
    do {                                                                       \
        if (test((_r), (_res))) {                                              \
            return 1;                                                          \
        }                                                                      \
    } while (0)

static int test(uint16_t r, uint16_t res) {
    uint16_t _res;

    _res = _cc_get_cosets_cnt(r);

    if (_res != res) {
        printf("ERROR: _cc_get_cosets_cnt(%u) = %u != %u\n", r, _res, res);
        return 1;
    }

    return 0;
}

int main(void) {
    TEST_WRAPPER(0, 0);
    TEST_WRAPPER(1, 1);
    TEST_WRAPPER(2, 1);
    TEST_WRAPPER(3, 2);
    TEST_WRAPPER(8, 2);
    TEST_WRAPPER(11, 4);
    TEST_WRAPPER(12, 3);
    TEST_WRAPPER(16, 3);
    TEST_WRAPPER(18, 4);
    TEST_WRAPPER(19, 5);
    TEST_WRAPPER(255, 35);
    TEST_WRAPPER(256, 32);
    TEST_WRAPPER(389, 42);

    // Duplicating a test to check that a function is deterministic
    TEST_WRAPPER(1034, 82);
    TEST_WRAPPER(1034, 82);

    return 0;
}