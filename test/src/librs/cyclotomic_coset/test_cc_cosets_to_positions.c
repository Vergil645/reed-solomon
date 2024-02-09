#include <stdio.h>
#include <stdlib.h>

#include <cyclotomic_coset.h>

#include <util.h>

#define TEST_WRAPPER(_cosets, _cosets_cnt, _positions, _positions_cnt)         \
    do {                                                                       \
        if (test((_cosets), (_cosets_cnt), (_positions), (_positions_cnt))) {  \
            return 1;                                                          \
        }                                                                      \
    } while (0)

static int test(const coset_t* cosets, uint16_t cosets_cnt, uint16_t* positions,
                uint16_t positions_cnt) {
    uint16_t* _positions;
    int ret = 0;

    _positions = (uint16_t*)malloc(positions_cnt * sizeof(uint16_t));
    if (!_positions) {
        printf("ERROR: cannot allocate memory for _positions\n");
        return 1;
    }

    cc_cosets_to_positions(cosets, cosets_cnt, _positions, positions_cnt);

    if (!util_u16_array_eq(_positions, (size_t)positions_cnt, positions,
                           (size_t)positions_cnt)) {
        printf("ERROR: cc_cosets_to_positions(...): incorrect positions:\n");

        printf("\tcorrect = ");
        util_u16_array_printf(positions, (size_t)positions_cnt);
        printf("\n");

        printf("\tactual  = ");
        util_u16_array_printf(_positions, (size_t)positions_cnt);
        printf("\n");

        ret = 1;
    }

    free(_positions);

    return ret;
}

int main(void) {
    // Test 1
    {
        uint16_t cosets_cnt = 2;
        uint16_t positions_cnt = 3;
        coset_t cosets[] = {
            {.leader = 21845, .size = 2},
            {.leader = 0, .size = 1},
        };
        uint16_t positions[] = {21845, 43690, 0};

        TEST_WRAPPER(cosets, cosets_cnt, positions, positions_cnt);
    }

    // Test 2
    {
        uint16_t cosets_cnt = 4;
        uint16_t positions_cnt = 11;
        coset_t cosets[] = {
            {.leader = 4369, .size = 4},
            {.leader = 13107, .size = 4},
            {.leader = 21845, .size = 2},
            {.leader = 0, .size = 1},
        };
        uint16_t positions[] = {4369,  8738,  17476, 34952, 13107, 26214,
                                52428, 39321, 21845, 43690, 0};

        TEST_WRAPPER(cosets, cosets_cnt, positions, positions_cnt);
    }

    // Test 3
    {
        uint16_t cosets_cnt = 2;
        uint16_t positions_cnt = 11;
        coset_t cosets[] = {
            {.leader = 257, .size = 8},
            {.leader = 30583, .size = 4},
        };
        uint16_t positions[] = {257,   514,   1028,  2056,  4112, 8224,
                                16448, 32896, 30583, 61166, 56797};

        TEST_WRAPPER(cosets, cosets_cnt, positions, positions_cnt);
    }

    // Test 4
    {
        uint16_t cosets_cnt = 3;
        uint16_t positions_cnt = 18;
        coset_t cosets[] = {
            {.leader = 771, .size = 8},
            {.leader = 1285, .size = 8},
            {.leader = 30583, .size = 4},
        };
        uint16_t positions[] = {771,   1542,  3084,  6168,  12336, 24672,
                                49344, 33153, 1285,  2570,  5140,  10280,
                                20560, 41120, 16705, 33410, 30583, 61166};

        TEST_WRAPPER(cosets, cosets_cnt, positions, positions_cnt);
    }

    return 0;
}