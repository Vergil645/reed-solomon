#include <stdio.h>

#include <gf65536.h>

#define TEST_GF_DIV_EE(_gf, _a, _b, _r)                                        \
    do {                                                                       \
        if (gf_div_ee((_gf), (_a), (_b)) != (_r)) {                            \
            printf("ERROR: gf_div_ee(*, %u, %u) != %u", (_a), (_b), (_r));     \
            gf_free((_gf));                                                    \
            return 1;                                                          \
        }                                                                      \
    } while (0)

int main(void) {
    GF_t _gf;
    GF_t *gf = &_gf;
    int ret;

    ret = gf_alloc(gf);
    if (ret) {
        printf("ERROR: gf_alloc returned %d", ret);
        return ret;
    }

    gf_init(gf);

    // Test data obtained by SageMath.
    TEST_GF_DIV_EE(gf, 0, 45687, 0);
    TEST_GF_DIV_EE(gf, 65512, 65512, 1);
    TEST_GF_DIV_EE(gf, 12320, 29623, 11439);
    TEST_GF_DIV_EE(gf, 31193, 63233, 27486);
    TEST_GF_DIV_EE(gf, 21844, 54054, 49588);
    TEST_GF_DIV_EE(gf, 38756, 35149, 10047);
    TEST_GF_DIV_EE(gf, 5768, 15888, 24163);

    gf_free(gf);

    return 0;
}