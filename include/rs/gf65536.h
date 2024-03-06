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

#include "cyclotomic_coset.h"
#include "prelude.h"

/**
 * @brief Galois field size. Equal to (N + 1).
 */
#define GF_FIELD_SIZE 65536

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
 * @brief Normal basises of subfields GF(2^m), where m - cyclotomic coset sizes.
 * @details g_normal_basis_by_subfield[i][j] - j-th basis element of
 * \f$GF(2^{2^i})\f$ subfield (i=0,...,4).
 */
static __attribute_maybe_unused__ element_t
    g_normal_basis_by_subfield[CC_COSET_SIZES_CNT][CC_MAX_COSET_SIZE] = {
        [0] = {1},
        [1] = {44234, 44235},
        [2] = {10800, 47860, 34555, 5694},
        [3] = {16402, 53598, 44348, 63986, 22060, 64366, 6088, 32521},
        [4] = {2048, 2880, 7129, 30616, 2643, 6897, 29685, 7378, 30100, 2743, 20193, 36223, 24055,
               41458, 41014, 61451},
};

/**
 * @brief Galois field data.
 * @details Field definition: \f$GF(2)[x] / \left<PRIMITIVE\_POLY\right>\f$.\n
 * Primitive element: \f$\alpha = x\f$.
 */
typedef struct GF_t {
    /**
     * @brief Primitive element powers. Power can belongs to range [0; 2*N-2];
     * @details \f$pow\_table_i = \alpha^i\f$
     */
    element_t pow_table[(N << 1) - 1];

    /**
     * @brief Logarithm to the base of a primitive element.
     * @details \f$log\_table_e = d\f$ s.t. \f$\alpha^d = e\f$
     */
    uint16_t log_table[GF_FIELD_SIZE];

    /**
     * @brief Coefficients in normal basis of subfields GF(2^m).
     * @details j-th bit of normal_repr_by_subfield[i][d] - j-th coefficient in
     * normal basis of \f$GF(2^{2^i})\f$ subfield of element \f$\alpha^d\f$.
     */
    uint16_t* normal_repr_by_subfield[CC_COSET_SIZES_CNT];

    /**
     * @brief .normal_repr_by_subfield memory.
     */
    uint16_t _normal_repr_by_subfield_memory[CC_COSET_SIZES_CNT * N];
} GF_t;

/**
 * @brief Create Galois field data structure.
 *
 * @return pointer to created Galois field data structure on success and NULL
 * otherwise.
 */
GF_t* gf_create();

/**
 * @brief Destroy Galois field data structure.
 *
 * @param gf Galois field data.
 */
void gf_destroy(GF_t* gf);

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
element_t gf_mul_ee(GF_t* gf, element_t a, element_t b);

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
element_t gf_div_ee(GF_t* gf, element_t a, element_t b);

/**
 * @brief Compute the sum of 2 elements in Galois field.
 *
 * @param a first element (result will be placed here).
 * @param b second element.
 * @param symbol_size symbol size (must be divisible by 2).
 */
void gf_add(void* a, const void* b, size_t symbol_size);

/**
 * @brief Compute multiplication of element and coefficient in Galois field.
 *
 * @param gf Galois field data.
 * @param a element (result will be placed here).
 * @param coef coefficient.
 * @param symbol_size symbol size (must be divisible by 2).
 */
void gf_mul(GF_t* gf, void* a, element_t coef, size_t symbol_size);

/**
 * @brief Compute "A += c * B" expression in Galois field.
 *
 * @param gf Galois field data.
 * @param a first element (result will be placed here).
 * @param coef coefficient.
 * @param b second element.
 * @param symbol_size symbol size (must be divisible by 2).
 */
void gf_madd(GF_t* gf, void* a, element_t coef, const void* b, size_t symbol_size);

#endif