#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include <rlc/gf256.h>
#include <rlc/rlc.h>

#define ARRAY_SIZE 1300
#define NUMBER_OF_OPERATIONS 100000

static void init_data(uint8_t* a, uint8_t* b, uint8_t* coef) {
    *coef = 2 + (uint8_t)rand() % (UINT8_MAX - 1);
    a[0] = (uint8_t)rand();
    b[0] = (uint8_t)rand();
    for (size_t i = 1; i < ARRAY_SIZE; ++i) {
        a[i] = a[i - 1] + (uint8_t)rand() % 3;
        b[i] = b[i - 1] + (uint8_t)rand() % 4;
    }
}

static void compare_by_clock(uint8_t* a, uint8_t* b, uint8_t** mul) {
    clock_t start;
    clock_t end;
    clock_t elapsed_time_add;
    clock_t elapsed_time_mul;
    uint8_t coef;

    init_data(a, b, &coef);

    start = clock();
    for (int i = 0; i < NUMBER_OF_OPERATIONS; ++i)
        gf256_symbol_add_scaled(a, coef, b, ARRAY_SIZE, mul);
    end = clock();
    elapsed_time_mul = (end - start);

    init_data(a, b, &coef);

    start = clock();
    for (int i = 0; i < NUMBER_OF_OPERATIONS; ++i)
        gf256_symbol_add_scaled(a, 1, b, ARRAY_SIZE, mul);
    end = clock();
    elapsed_time_add = (end - start);

    printf("===== clock() =====\n");
    printf("< += >   (coef=1) time: %lu\n", elapsed_time_add);
    printf("< +=.* > (coef>1) time: %lu\n", elapsed_time_mul);
    printf("time ratio \"< +=.* >/< += >\": %.3f\n", (double)elapsed_time_mul / (double)elapsed_time_add);
}

#define USECS_IN_SEC 1000000
static double get_diff_secs(struct timeval begin, struct timeval end) {
    return (double)(end.tv_usec - begin.tv_usec) / USECS_IN_SEC + (double)(end.tv_sec - begin.tv_sec);
}

static void compare_by_gettimeofday(uint8_t* a, uint8_t* b, uint8_t** mul) {
    struct timeval start;
    struct timeval end;
    double elapsed_time_add;
    double elapsed_time_mul;
    uint8_t coef;

    init_data(a, b, &coef);

    gettimeofday(&start, NULL);
    for (int i = 0; i < NUMBER_OF_OPERATIONS; ++i)
        gf256_symbol_add_scaled(a, coef, b, ARRAY_SIZE, mul);
    gettimeofday(&end, NULL);
    elapsed_time_mul = get_diff_secs(start, end);

    init_data(a, b, &coef);

    gettimeofday(&start, NULL);
    for (int i = 0; i < NUMBER_OF_OPERATIONS; ++i)
        gf256_symbol_add_scaled(a, 1, b, ARRAY_SIZE, mul);
    gettimeofday(&end, NULL);
    elapsed_time_add = get_diff_secs(start, end);

    printf("===== gettimeofday(...) =====\n");
    printf("< += >   (coef=1) time: %.3f secs\n", elapsed_time_add);
    printf("< +=.* > (coef>1) time: %.3f secs\n", elapsed_time_mul);
    printf("time ratio \"< +=.* >/< += >\": %.3f\n", elapsed_time_mul / elapsed_time_add);
}

int main(void) {
    RLC_t* rlc;
    uint8_t* a;
    uint8_t* b;

    rlc = rlc_create();
    if (!rlc) {
        printf("ERROR: rlc_create returned NULL\n");
        return 1;
    }

    srand(time(NULL));

    a = (uint8_t*)calloc(ARRAY_SIZE, sizeof(uint8_t));
    if (!a) {
        printf("ERROR: couldn't allocate a\n");
        return 1;
    }

    b = (uint8_t*)calloc(ARRAY_SIZE, sizeof(uint8_t));
    if (!b) {
        printf("ERROR: couldn't allocate b\n");
        free(a);
        return 1;
    }

    compare_by_clock(a, b, rlc->mul_table);
    compare_by_gettimeofday(a, b, rlc->mul_table);

    free(b);
    free(a);
    rlc_destroy(rlc);

    return 0;
}