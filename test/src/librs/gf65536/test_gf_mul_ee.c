#include <stdio.h>

#include <gf65536.h>

#define TEST_GF_MUL_EE(_gf, _a, _b, _r)                                        \
    do {                                                                       \
        if (gf_mul_ee((_gf), (_a), (_b)) != (_r)) {                            \
            printf("ERROR: gf_mul_ee(*, %u, %u) != %u", (_a), (_b), (_r));     \
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
    TEST_GF_MUL_EE(gf, 1, 645, 645);
    TEST_GF_MUL_EE(gf, 46478, 0, 0);
    TEST_GF_MUL_EE(gf, 31981, 38739, 42167);
    TEST_GF_MUL_EE(gf, 2491, 54249, 5290);
    TEST_GF_MUL_EE(gf, 60895, 36296, 21017);
    TEST_GF_MUL_EE(gf, 62824, 46526, 6710);
    TEST_GF_MUL_EE(gf, 58263, 29917, 33120);

    gf_free(gf);

    return 0;
}