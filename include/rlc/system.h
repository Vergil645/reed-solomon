#ifndef __RLC_SYSTEM_H__
#define __RLC_SYSTEM_H__

#include "equation.h"

#define ENTRY_INDEX_NONE 0xfffffffful

typedef struct {
    int max_equations;
    int n_equations;
    uint16_t first_id_id;
    uint16_t last_symbol_id;
    equation_t** equations;
} system_t;

bool system_set_bounds(system_t* system, uint16_t first, uint16_t last);

int system_add_with_elimination(system_t* system, equation_t* eq, uint8_t* inv_table, uint8_t** mul_table, int* decoded,
                                equation_t** removed, int* used_in_system);

#endif