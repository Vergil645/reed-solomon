#include <stdio.h>

#include <rs/gf65536.h>

#define TEST_WRAPPER(_gf, _a, _b, _res)                                                                                \
    do {                                                                                                               \
        if (test((_gf), (_a), (_b), (_res))) {                                                                         \
            gf_destroy((_gf));                                                                                         \
            return 1;                                                                                                  \
        }                                                                                                              \
    } while (0)

static int test(GF_t* gf, element_t a, element_t b, element_t res) {
    element_t _res;

    _res = gf_mul_ee(gf, a, b);

    if (_res != res) {
        printf("ERROR: gf_mul_ee(*, %u, %u) = %u != %u\n", a, b, _res, res);
        return 1;
    }

    return 0;
}

int main(void) {
    GF_t* gf;

    gf = gf_create();
    if (!gf) {
        printf("ERROR: gf_create returned NULL\n");
        return 1;
    }

    // Test data obtained by SageMath.
    TEST_WRAPPER(gf, 1, 645, 645);
    TEST_WRAPPER(gf, 46478, 0, 0);
    TEST_WRAPPER(gf, 31981, 38739, 42167);
    TEST_WRAPPER(gf, 2491, 54249, 5290);
    TEST_WRAPPER(gf, 60895, 36296, 21017);
    TEST_WRAPPER(gf, 62824, 46526, 6710);
    TEST_WRAPPER(gf, 58263, 29917, 33120);

    gf_destroy(gf);

    return 0;
}