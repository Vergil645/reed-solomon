#ifndef __TEST_UTIL_H__
#define __TEST_UTIL_H__

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <memory/seq.h>

bool util_u16_array_eq(const uint16_t* a, size_t a_len, const uint16_t* b, size_t b_len);

void util_u16_array_printf(const uint16_t* a, size_t a_len);

void util_generate_inf_symbols(const symbol_seq_t* inf_symbols);

void util_init_rcv_symbols(const symbol_seq_t* src_symbols, symbol_seq_t* rcv_symbols);

void util_choose_and_erase_symbols(symbol_seq_t* rcv_symbols, uint16_t t, bool* is_erased);

#endif