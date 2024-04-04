/**
 * @file symbol.h
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief Contains symbol_t definition and functions for interaction with symbols.
 * @date 2024-01-20
 *
 * @copyright Copyright (c) 2024
 */

#ifndef __MEMORY_SYMBOL_H__
#define __MEMORY_SYMBOL_H__

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/**
 * @brief Symbol data type.
 */
typedef struct {
    /**
     * @brief Symbol data.
     */
    uint8_t* data;
} symbol_t;

/**
 * @brief Create symbol.
 *
 * @param symbol_size symbol size.
 * @return pointer to created symbol on success and NULL otherwise.
 */
symbol_t* symbol_create(size_t symbol_size);

/**
 * @brief Destroy symbol.
 *
 * @param s symbol.
 */
void symbol_destroy(symbol_t* s);

/**
 * @brief Check whether symbols are equal.
 *
 * @param a first symbol.
 * @param b second symbol.
 * @param symbol_size symbol size.
 * @return true if symbols are equal, false otherwise.
 */
bool symbol_eq(const symbol_t* a, const symbol_t* b, size_t symbol_size);

/**
 * @brief Write symbol to stdout.
 *
 * @param s symbol.
 * @param symbol_size symbol size.
 */
void symbol_printf(const symbol_t* s, size_t symbol_size);

#endif