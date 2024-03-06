#include <stdio.h>

#include <rs/gf65536.h>

#define TEST_WRAPPER(_gf, _a, _b, _res)                                                            \
    do {                                                                                           \
        if (test((_gf), (_a), (_b), (_res))) {                                                     \
            gf_destroy((_gf));                                                                     \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)

static int test(GF_t* gf, element_t a, element_t b, element_t res) {
    element_t _res;

    _res = gf_div_ee(gf, a, b);

    if (_res != res) {
        printf("ERROR: gf_div_ee(*, %u, %u) = %u != %u\n", a, b, _res, res);
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
    TEST_WRAPPER(gf, 0, 45687, 0);
    TEST_WRAPPER(gf, 65512, 65512, 1);
    TEST_WRAPPER(gf, 12320, 29623, 11439);
    TEST_WRAPPER(gf, 31193, 63233, 27486);
    TEST_WRAPPER(gf, 21844, 54054, 49588);
    TEST_WRAPPER(gf, 38756, 35149, 10047);
    TEST_WRAPPER(gf, 5768, 15888, 24163);

    gf_destroy(gf);

    return 0;
}