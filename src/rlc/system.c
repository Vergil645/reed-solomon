#include <assert.h>
#include <stdbool.h>

#include <rlc/system.h>

#define MIN(a, b) (((a) <= (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

bool system_set_bounds(system_t* system, uint16_t first, uint16_t last) {
    if (last + 1 - first > system->max_equations)
        return false;
    if (system->first_id_id < first)
        return false;
    system->first_id_id = first;
    if (last < system->last_symbol_id && system->last_symbol_id != ID_NONE)
        return false;
    system->last_symbol_id = last;
    return true;
}

equation_t* system_get_pivot_for_id(system_t* system, uint16_t id) {
    if (id >= system->first_id_id && id < (system->max_equations + system->first_id_id) &&
        system->equations[id - system->first_id_id]) {
        return system->equations[id - system->first_id_id];
    }
    return NULL;
}

static int reduce_equation(system_t* system, equation_t* eq, uint8_t** mul_table, uint8_t* inv_table) {

    equation_adjust_non_zero_bounds(eq);
    if (eq->pivot == ID_NONE)
        return 0;

    int err = 0;
    for (uint16_t id = eq->pivot; id <= eq->last_non_zero_id && !equation_is_zero(eq); id++) {
        uint8_t coef = equation_get_coef(eq, id);
        if (coef != 0) {
            equation_t* pivot_equation = system_get_pivot_for_id(system, id);
            if (pivot_equation != NULL) {
                // equation_multiply(pivot_equation, coef, mul_table);
                /* we cancel the coef */

                equation_multiply(eq,
                                  mul_table[equation_get_coef(pivot_equation, pivot_equation->pivot)][inv_table[coef]],
                                  mul_table);

                // we reduce the equation and remove its pivot coefficient by
                // adding the multiplied equation and the system's pivot
                // equation
                err = equation_add(eq, pivot_equation);
                if (err) {
                    break;
                }
            }
        }
    }

    return err;
}

static uint32_t system_add(system_t* system, equation_t* eq, equation_t** removed) {
    assert(system != NULL);
    *removed = NULL;
    if (equation_is_zero(eq)) {
        return ENTRY_INDEX_NONE;
    }

    if (system->first_id_id == ID_NONE) {
        system->first_id_id = eq->pivot;
    }

    // uint32_t old_size = system->max_equations;
    uint16_t new_i0 = eq->pivot;
    uint16_t set_i0 = system->first_id_id;

    // printf("system_add: new_i0=%u set_i0=%u\n", new_i0, set_i0);
    // fflush(stdout);
    assert(new_i0 >= set_i0);
    /* The added symbol has first nonzero index outside the current set */
    // ... was deleted

    uint32_t idx_pos = new_i0 - set_i0;
    if (idx_pos < system->max_equations) {
        if (system->equations[idx_pos] != NULL) {
            *removed = system->equations[idx_pos];
            system->n_equations--;
            //            equation_free(cnx, system->equations[idx_pos]);
        }
        system->equations[idx_pos] = eq;
        system->n_equations++;
        if (system->last_symbol_id == ID_NONE) {
            system->last_symbol_id = eq->last_non_zero_id;
        } else {
            system->last_symbol_id = MAX(eq->last_non_zero_id, system->last_symbol_id);
        }
        return idx_pos;
    }

    assert(idx_pos < system->max_equations);

    return ENTRY_INDEX_NONE;
}

static uint32_t system_add_as_pivot(system_t* system, equation_t* eq, uint8_t* inv_table, uint8_t** mul_table,
                                    int* decoded, equation_t** removed) {
    *decoded = 0;
    equation_adjust_non_zero_bounds(eq);
    if (eq->pivot == ID_NONE) {
        return 0;
    }
    uint16_t first_id = eq->pivot;

    int n_non_null_equations = 0;
    for (uint32_t i = 0; i < system->max_equations && n_non_null_equations < system->n_equations; i++) {
        if (system->equations[i]) {
            n_non_null_equations++;
        }
        // TODO: add one temp equation to store the add and mul results  (it
        // could belong to the system) so that we avoid doing inv each time
        // TODO: instead we could also, instead of doing inv to reset at 1,
        // multiplying  by coef*inv[pivot] instead of just coef and remove the
        // multiplication by inv !
        if (system->equations[i] && equation_get_coef(system->equations[i], first_id) != 0) {
            uint8_t coef = equation_get_coef(system->equations[i], first_id); /* XXX: if coef == 0, nothing to do */
            if (coef != 0) {

                //            equation_multiply(eq,
                //            inv_table[equation_get_coef(eq, eq->pivot)],
                //            mul_table);
                uint8_t pivot_coef = equation_get_coef(eq, eq->pivot);
                assert(mul_table[inv_table[pivot_coef]][coef] != 0);
                equation_multiply(eq, mul_table[inv_table[pivot_coef]][coef], mul_table);

                bool has_one_id_before_add = equation_has_one_id(system->equations[i]);
                // int err = equation_add(system->equations[i], eq);
                equation_add(system->equations[i], eq);
                bool is_decoded = !has_one_id_before_add && equation_has_one_id(system->equations[i]);
                if (is_decoded) {
                    uint16_t si = equation_get_min_symbol_id(system->equations[i]);
                    if (equation_get_coef(system->equations[i], si) != 1) {
                        equation_multiply(system->equations[i], inv_table[equation_get_coef(system->equations[i], si)],
                                          mul_table);
                    }
                    assert(equation_get_coef(system->equations[i], si) == 1);
                    *decoded = 1;
                }
            }
        }
    }
    // TODO: at this point, eq's pivot could be different from 1. I don't think
    // this is a problem, the application must just ask for a normalized source
    // symbol
    return system_add(system, eq, removed);
}

int system_add_with_elimination(system_t* system, equation_t* eq, uint8_t* inv_table, uint8_t** mul_table, int* decoded,
                                equation_t** removed, int* used_in_system) {
    *removed = NULL;
    *decoded = 0;
    *used_in_system = 0;

    int err = reduce_equation(system, eq, mul_table, inv_table);

    if (!err && !equation_is_zero(eq)) {
        uint32_t idx = system_add_as_pivot(system, eq, inv_table, mul_table, decoded, removed);

        if (idx == ENTRY_INDEX_NONE || idx >= system->max_equations) {
            return 0;
        }

        *used_in_system = 1;
        equation_t* stored_symbol = system->equations[idx];

        bool is_decoded = equation_has_one_id(stored_symbol);
        if (is_decoded) {
            *decoded = 1;
            uint16_t si = equation_get_min_symbol_id(stored_symbol);
            if (equation_get_coef(stored_symbol, si) != 1) {
                equation_multiply(stored_symbol, inv_table[equation_get_coef(stored_symbol, si)], mul_table);
            }
            assert(equation_get_coef(stored_symbol, si) == 1);
        }
    }

    return err;
}