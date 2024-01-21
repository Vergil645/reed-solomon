/**
 * @file gf65536.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains implementation of a Galois field of size 65536.
 * @date 2024-01-17
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __REED_SOLOMON_GF65536_H__
#define __REED_SOLOMON_GF65536_H__

#include <stdint.h>

#include "constants.h"

/**
 * @brief Primitive polynomial of a Galois field: x^16 + x^5 + x^3 + x^2 + 1.
 * @details Primitive element: \f$\alpha = x\f$.
 */
#define GF_PRIMITIVE_POLY 65581

/**
 * @brief Galois field element type.
 */
typedef uint16_t element_t;

/**
 * @brief Polynomial type.
 */
typedef uint32_t poly_t;

/**
 * @brief Galois field data.
 * @details Field definition: \f$GF(2)[x] / \left<PRIMITIVE\_POLY\right>\f$.\n
 * Primitive element: \f$\alpha = x\f$.
 */
typedef struct GF_t {
    /**
     * @brief Primitive element powers.
     * @details \f$pow\_table_i = \alpha^i\f$
     */
    element_t *pow_table;

    /**
     * @brief Logarithm to the base of a primitive element.
     * @details \f$log\_table_e = d\f$ s.t. \f$\alpha^d = e\f$
     */
    uint16_t *log_table;
} GF_t;

/**
 * @brief Allocate memory for a Galois field implementation.
 *
 * @param gf where to place allocated data.
 * @return 0 on success,\n
 *         1 on memory allocation error.
 */
int gf_alloc(GF_t *gf);

/**
 * @brief Make necessary for Galois field implementation pre-calculations.
 * @details Pre-calculations list:
 * 1. fill pow_table: \f$pow\_table_i = \alpha^i\f$;
 * 2. fill log_table: \f$log\_table_e = d\f$ such that \f$\alpha^d = e\f$.
 *
 * @param gf Galois field data.
 */
void gf_init(GF_t *gf);

/**
 * @brief Deallocate Galois field data.
 *
 * @param gf Galois field data.
 */
void gf_free(GF_t *gf);

/**
 * @brief Compute multiplication of 2 elements in Galois field.
 * @details Use pre-computed data to optimize computation process:
 * \f$a * b = \alpha^{(\log_{\alpha} a + \log_{\alpha} b) \mod N}\f$
 *
 * @param gf Galois field data.
 * @param a first multiplier.
 * @param b second multiplier.
 * @return multiplication result.
 */
element_t gf_mul_ee(const GF_t *gf, element_t a, element_t b);

/**
 * @brief Compute quotient of 2 elements in Galois field.
 * @details Use pre-computed data to optimize computation process:
 * \f$a / b = \alpha^{(\log_{\alpha} a - \log_{\alpha} b) \mod N}\f$
 *
 * @param gf Galois field data.
 * @param a divisible element.
 * @param b divisor.
 * @return division result.
 */
element_t gf_div_ee(const GF_t *gf, element_t a, element_t b);

#endif