/**
 * @file reed_solomon.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief reed_solomon.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <string.h>

#include <reed_solomon.h>

int rs_alloc(RS_t *rs) {
    // Nothing;
    return 0;
}

void rs_init(RS_t *rs, GF_t *gf, CC_t *cc, FFT_t *fft) {
    rs->gf = gf;
    rs->cc = cc;
    rs->fft = fft;
}

void rs_free(RS_t *rs) {
    // Nothing
}

_static void _rs_get_syndrome_poly(const RS_t *rs, symbol_seq_t seq,
                                   const uint16_t *positions,
                                   symbol_seq_t syndrome_poly) {
    assert(rs != NULL);
    assert(positions != NULL);
    assert(seq.symbol_size == syndrome_poly.symbol_size);

    fft_transform(rs->fft, seq, positions, syndrome_poly);
}

_static void _rs_get_locator_poly(const RS_t *rs, const uint16_t *positions,
                                  uint16_t positions_cnt,
                                  element_t *locator_poly,
                                  uint16_t locator_max_len) {
    assert(rs != NULL);
    assert(positions != NULL);
    assert(locator_poly != NULL);
    assert(positions_cnt + 1 <= locator_max_len);

    GF_t *gf = rs->gf;
    element_t *pow_table = gf->pow_table;
    uint16_t pos;
    element_t c;

    locator_poly[0] = 1;
    for (uint16_t d = 0; d < positions_cnt; ++d) {
        locator_poly[d + 1] = 0;
        pos = positions[d];
        c = pow_table[pos];

        for (uint16_t i = d + 1; i > 0; --i) {
            locator_poly[i] ^= gf_mul_ee(gf, locator_poly[i - 1], c);
        }
    }
}

_static void _rs_get_rep_symbols_locator_poly(const RS_t *rs, uint16_t r,
                                              const coset_t *rep_cosets,
                                              uint16_t rep_cosets_cnt,
                                              element_t *locator_poly,
                                              uint16_t locator_max_len) {
    assert(rs != NULL);
    assert(rep_cosets != NULL);
    assert(locator_poly != NULL);
    assert(r + 1 <= locator_max_len);

#ifndef NDEBUG
    uint16_t r1 = 0;
    for (uint16_t i = 0; i < rep_cosets_cnt; ++i) {
        r1 += rep_cosets[i].size;
    }
    assert(r1 == r);
#endif

    uint16_t coset_elements[CC_MAX_COSET_SIZE];
    element_t coset_locator_poly[RS_COSET_LOCATOR_MAX_LEN];
    coset_t coset;
    uint16_t d; // Locator polynomial degree

    // Locator polynomial initialization.
    d = 0;
    memset((void *)locator_poly, 0, (r + 1) * sizeof(element_t));
    locator_poly[0] = 1;

    for (uint16_t coset_idx = 0; coset_idx < rep_cosets_cnt; ++coset_idx) {
        coset = rep_cosets[coset_idx];

        coset_elements[0] = coset.leader;
        for (uint16_t i = 1; i < coset.size; ++i) {
            coset_elements[i] = NEXT_COSET_ELEMENT(coset_elements[i - 1]);
        }

        _rs_get_locator_poly(rs, coset_elements, coset.size, coset_locator_poly,
                             RS_COSET_LOCATOR_MAX_LEN);

#ifndef NDEBUG
        for (uint16_t _i = 0; _i <= coset.size; ++_i) {
            assert(coset_locator_poly[_i] == 0 || coset_locator_poly[_i] == 1);
        }
#endif

        for (uint16_t i = d;; --i) {
            if (locator_poly[i] == 1) {
                for (uint16_t j = 1; j <= coset.size; ++j) {
                    locator_poly[i + j] ^= coset_locator_poly[j];
                }
            }

            if (i == 0) {
                break;
            }
        }

        d += coset.size;
        assert(locator_poly[d] == 1);
    }

    assert(d == r);

#ifndef NDEBUG
    for (uint16_t _i = 0; _i <= r; ++_i) {
        assert(locator_poly[_i] == 0 || locator_poly[_i] == 1);
    }
#endif
}

_static element_t _rs_get_forney_coef(const RS_t *rs,
                                      const element_t *locator_poly, uint16_t d,
                                      uint16_t pos) {
    assert(rs != NULL);
    assert(locator_poly != NULL);

    GF_t *gf = rs->gf;
    element_t *pow_table = gf->pow_table;
    element_t p; // divisible element = alpha^{position}
    element_t q; // divisor = locator_poly'(alpha^{-position})
    element_t c;

    p = pow_table[pos];
    q = 0;
    for (uint16_t j = 0; j < d; j += 2) {
        c = locator_poly[j + 1];

        if (c == 0) {
            continue;
        } else if (c == 1) {
            q ^= pow_table[(j * (N - pos)) % N];
        } else {
            q ^= gf_mul_ee(gf, pow_table[(j * (N - pos)) % N], c);
        }
    }

    return gf_div_ee(gf, p, q);
}

_static void _rs_get_evaluator_poly(const RS_t *rs, symbol_seq_t syndrome_poly,
                                    const element_t *locator_poly,
                                    symbol_seq_t evaluator_poly) {
    assert(rs != NULL);
    assert(syndrome_poly.symbol_size == evaluator_poly.symbol_size);

    GF_t *gf = rs->gf;
    size_t symbol_size = syndrome_poly.symbol_size;
    uint16_t r = syndrome_poly.length;
    element_t c;

    // Evaluator polynomial initialization.
    for (uint16_t i = 0; i < r; ++i) {
        memset((void *)evaluator_poly.symbols[i].data, 0,
               symbol_size * sizeof(element_t));
    }

    for (uint16_t i = 0; i < r; ++i) {
        c = locator_poly[i];

        if (c == 0) {
            continue;
        } else if (c == 1) {
            for (uint16_t j = 0; j < r - i; ++j) {
                for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                    evaluator_poly.symbols[i + j].data[e_idx] ^=
                        syndrome_poly.symbols[j].data[e_idx];
                }
            }
        } else {
            for (uint16_t j = 0; j < r - i; ++j) {
                for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                    evaluator_poly.symbols[i + j].data[e_idx] ^=
                        gf_mul_ee(gf, syndrome_poly.symbols[j].data[e_idx], c);
                }
            }
        }
    }
}

_static void _rs_get_repair_symbols(const RS_t *rs,
                                    const element_t *locator_poly,
                                    symbol_seq_t evaluator_poly,
                                    symbol_seq_t rep_symbols,
                                    const uint16_t *rep_positions) {
    assert(rs != NULL);
    assert(locator_poly != NULL);
    assert(rep_positions != NULL);
    assert(evaluator_poly.symbol_size == rep_symbols.symbol_size);

    GF_t *gf = rs->gf;
    size_t symbol_size = evaluator_poly.symbol_size;
    uint16_t r = rep_symbols.length;
    uint16_t pos;
    element_t coef;

    fft_partial_transform(rs->fft, evaluator_poly, rep_symbols, rep_positions);

    // TODO: implement deffered Forney scalling
    for (uint16_t i = 0; i < r; ++i) {
        pos = rep_positions[i];
        coef = _rs_get_forney_coef(rs, locator_poly, r, pos);

        for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
            rep_symbols.symbols[i].data[e_idx] =
                gf_mul_ee(gf, coef, rep_symbols.symbols[i].data[e_idx]);
        }
    }
}

_static void _rs_restore_erased(const RS_t *rs, const element_t *locator_poly,
                                symbol_seq_t evaluator_poly,
                                const uint16_t *positions,
                                symbol_seq_t rcv_symbols,
                                const uint16_t *erased_indices) {
    assert(rs != NULL);
    assert(locator_poly != NULL);
    assert(positions != NULL);
    assert(erased_indices != NULL);
    assert(evaluator_poly.symbol_size == rcv_symbols.symbol_size);

    GF_t *gf = rs->gf;
    element_t *pow_table = gf->pow_table;
    size_t symbol_size = evaluator_poly.symbol_size;
    uint16_t t = evaluator_poly.length;
    uint16_t erased_id;
    uint16_t pos;
    uint16_t j;
    element_t coef;
    element_t val;

    for (uint16_t erasure_idx = 0; erasure_idx < t; ++erasure_idx) {
        erased_id = erased_indices[erasure_idx];
        pos = positions[erased_id];
        coef = _rs_get_forney_coef(rs, locator_poly, t, pos);

        memset((void *)rcv_symbols.symbols[erased_id].data, 0,
               symbol_size * sizeof(element_t));

        j = (N - pos) % N;

        for (uint16_t i = 0; i < t; ++i) {
            val = pow_table[(i * j) % N];

            for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                rcv_symbols.symbols[erased_id].data[e_idx] ^= gf_mul_ee(
                    gf,
                    gf_mul_ee(gf, coef, evaluator_poly.symbols[i].data[e_idx]),
                    val);
            }
        }
    }
}

int rs_generate_repair_symbols(const RS_t *rs, symbol_seq_t inf_symbols,
                               symbol_seq_t rep_symbols) {
    assert(rs != NULL);
    assert(inf_symbols.length + rep_symbols.length <= N);
    assert(inf_symbols.symbol_size == rep_symbols.symbol_size);

    CC_t *cc = rs->cc;
    int ret = 0;
    size_t symbol_size = inf_symbols.symbol_size;
    uint16_t k = inf_symbols.length;
    uint16_t r = rep_symbols.length;
    uint16_t inf_max_cnt = 0;
    uint16_t rep_max_cnt = 0;
    uint16_t inf_cosets_cnt = 0;
    uint16_t rep_cosets_cnt = 0;
    coset_t *cosets;
    coset_t *inf_cosets;
    coset_t *rep_cosets;
    uint16_t *positions;
    uint16_t *inf_positions;
    uint16_t *rep_positions;
    element_t *locator_poly;
    symbol_seq_t syndrome_poly;
    symbol_seq_t evaluator_poly;

    cc_estimate_cosets_cnt(k, r, &inf_max_cnt, &rep_max_cnt);

    cosets = (coset_t *)malloc((inf_max_cnt + rep_max_cnt) * sizeof(coset_t));
    if (!cosets) {
        return 1;
    }
    inf_cosets = cosets;
    rep_cosets = cosets + inf_max_cnt;

    positions = (uint16_t *)malloc((k + r) * sizeof(uint16_t));
    if (!positions) {
        free(cosets);
        return 1;
    }
    inf_positions = positions;
    rep_positions = positions + k;

    locator_poly = (element_t *)malloc((r + 1) * sizeof(element_t));
    if (!locator_poly) {
        free(positions);
        free(cosets);
        return 1;
    }

    ret = seq_alloc(symbol_size, r, &syndrome_poly);
    if (ret) {
        free(locator_poly);
        free(positions);
        free(cosets);
        return ret;
    }

    ret = seq_alloc(symbol_size, r, &evaluator_poly);
    if (ret) {
        seq_free(&syndrome_poly);
        free(locator_poly);
        free(positions);
        free(cosets);
        return ret;
    }

    cc_select_cosets(cc, k, r, inf_cosets, inf_max_cnt, &inf_cosets_cnt,
                     rep_cosets, rep_max_cnt, &rep_cosets_cnt);

    cc_cosets_to_positions(inf_cosets, inf_cosets_cnt, inf_positions, k);
    cc_cosets_to_positions(rep_cosets, rep_cosets_cnt, rep_positions, r);

    _rs_get_syndrome_poly(rs, inf_symbols, inf_positions, syndrome_poly);

    _rs_get_rep_symbols_locator_poly(rs, r, rep_cosets, rep_cosets_cnt,
                                     locator_poly, r + 1);

    _rs_get_evaluator_poly(rs, syndrome_poly, locator_poly, evaluator_poly);

    _rs_get_repair_symbols(rs, locator_poly, evaluator_poly, rep_symbols,
                           rep_positions);

    seq_free(&evaluator_poly);
    seq_free(&syndrome_poly);
    free(locator_poly);
    free(positions);
    free(cosets);

    return 0;
}

int rs_restore_symbols(const RS_t *rs, uint16_t k, uint16_t r,
                       symbol_seq_t rcv_symbols, const uint16_t *erased_indices,
                       uint16_t t) {
    assert(rs != NULL);
    assert(erased_indices != NULL);
    assert((k + r) == rcv_symbols.length);

    CC_t *cc = rs->cc;
    int ret = 0;
    size_t symbol_size = rcv_symbols.symbol_size;
    uint16_t inf_max_cnt = 0;
    uint16_t rep_max_cnt = 0;
    uint16_t inf_cosets_cnt = 0;
    uint16_t rep_cosets_cnt = 0;
    coset_t *cosets;
    coset_t *inf_cosets;
    coset_t *rep_cosets;
    uint16_t *positions;
    uint16_t *inf_positions;
    uint16_t *rep_positions;
    uint16_t *erased_positions;
    element_t *locator_poly;
    symbol_seq_t syndrome_poly;
    symbol_seq_t evaluator_poly;

    if (r < t) {
        // Too many erases - symbols cannot be restored.
        return RS_ERR_CANNOT_RESTORE;
    }

    cc_estimate_cosets_cnt(k, r, &inf_max_cnt, &rep_max_cnt);

    cosets = (coset_t *)malloc((inf_max_cnt + rep_max_cnt) * sizeof(coset_t));
    if (!cosets) {
        return 1;
    }
    inf_cosets = cosets;
    rep_cosets = cosets + inf_max_cnt;

    positions = (uint16_t *)malloc((k + r) * sizeof(uint16_t));
    if (!positions) {
        free(cosets);
        return 1;
    }
    inf_positions = positions;
    rep_positions = positions + k;

    erased_positions = (element_t *)malloc(t * sizeof(uint16_t));
    if (!erased_positions) {
        free(positions);
        free(cosets);
        return 1;
    }

    locator_poly = (element_t *)malloc((t + 1) * sizeof(element_t));
    if (!locator_poly) {
        free(erased_positions);
        free(positions);
        free(cosets);
        return 1;
    }

    ret = seq_alloc(symbol_size, t, &syndrome_poly);
    if (ret) {
        free(locator_poly);
        free(erased_positions);
        free(positions);
        free(cosets);
        return ret;
    }

    ret = seq_alloc(symbol_size, t, &evaluator_poly);
    if (ret) {
        seq_free(&syndrome_poly);
        free(locator_poly);
        free(erased_positions);
        free(positions);
        free(cosets);
        return ret;
    }

    cc_select_cosets(cc, k, r, inf_cosets, inf_max_cnt, &inf_cosets_cnt,
                     rep_cosets, rep_max_cnt, &rep_cosets_cnt);

    cc_cosets_to_positions(inf_cosets, inf_cosets_cnt, inf_positions, k);
    cc_cosets_to_positions(rep_cosets, rep_cosets_cnt, rep_positions, r);

    _rs_get_syndrome_poly(rs, rcv_symbols, positions, syndrome_poly);

    for (uint16_t i = 0; i < t; ++i) {
        erased_positions[i] = positions[erased_indices[i]];
    }

    _rs_get_locator_poly(rs, erased_positions, t, locator_poly, t + 1);

    _rs_get_evaluator_poly(rs, syndrome_poly, locator_poly, evaluator_poly);

    _rs_restore_erased(rs, locator_poly, evaluator_poly, positions, rcv_symbols,
                       erased_indices);

    seq_free(&evaluator_poly);
    seq_free(&syndrome_poly);
    free(locator_poly);
    free(erased_positions);
    free(positions);
    free(cosets);

    return 0;
}