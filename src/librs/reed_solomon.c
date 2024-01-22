/**
 * @file reed_solomon.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief reed_solomon.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>

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

_static void _rs_get_coset_locator_poly(const RS_t *rs, coset_t coset,
                                        element_t *coset_locator_poly,
                                        uint16_t coset_locator_max_len) {
    assert(rs != NULL);
    assert(coset_locator_poly != NULL);
    assert(coset.size + 1 <= coset_locator_max_len);

    GF_t *gf = rs->gf;
    element_t *pow_table = gf->pow_table;
    uint16_t cur_coset_elem;

    coset_locator_poly[0] = 1;
    cur_coset_elem = coset.leader;
    for (uint16_t d = 0; d < coset.size; ++d) {
        coset_locator_poly[d + 1] = 0;
        for (uint16_t i = d + 1; i > 0; --i) {
            coset_locator_poly[i] ^= gf_mul_ee(gf, coset_locator_poly[i - 1],
                                               pow_table[cur_coset_elem]);
        }
        cur_coset_elem = (cur_coset_elem << 1) % N;
    }

    assert(cur_coset_elem == coset.leader);
    assert(coset_locator_poly[coset.size] == 1);

#ifndef NDEBUG
    for (uint16_t i = 0; i <= coset.size; ++i) {
        assert(coset_locator_poly[i] == 0 || coset_locator_poly[i] == 1);
    }
#endif
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

    element_t coset_locator_poly[RS_COSET_LOCATOR_MAX_LEN];
    uint16_t d; // Locator polynomial degree
    uint16_t i = 0;

    // Locator polynomial initialization.
    d = 0;
    locator_poly[0] = 1;
    for (i = 1; i <= r; ++i) {
        locator_poly[i] = 0;
    }

    for (uint16_t coset_idx = 0; coset_idx < rep_cosets_cnt; ++coset_idx) {
        _rs_get_coset_locator_poly(rs, rep_cosets[coset_idx],
                                   coset_locator_poly,
                                   RS_COSET_LOCATOR_MAX_LEN);

        i = d;
        while (1) {
            if (locator_poly[i] == 1) {
                for (uint16_t j = 0; j <= rep_cosets[coset_idx].size; ++j) {
                    locator_poly[i + j] ^= coset_locator_poly[j];
                }
            }
            assert(locator_poly[i] == 0);

            if (i == 0) {
                break;
            } else {
                --i;
            }
        }

        d += rep_cosets[coset_idx].size;
        assert(locator_poly[d] == 1);
    }

    assert(d == r);

#ifndef NDEBUG
    for (uint16_t _i = 0; _i <= r; ++_i) {
        assert(locator_poly[_i] == 0 || locator_poly[_i] == 1);
    }
#endif
}

_static void _rs_get_enc_forney_coefs(const RS_t *rs, uint16_t r,
                                      const element_t *locator_poly,
                                      const uint16_t *rep_positions,
                                      element_t *forney_coefs) {
    assert(rs != NULL);
    assert(locator_poly != NULL);
    assert(rep_positions != NULL);

    GF_t *gf = rs->gf;
    element_t *pow_table = gf->pow_table;
    uint16_t pos;
    element_t p; // divisible element = alpha^{position}
    element_t q; // divisor = locator_poly'(alpha^{-position})

    for (uint16_t i = 0; i < r; ++i) {
        pos = rep_positions[i];

        p = pow_table[pos];
        q = 0;
        for (uint16_t j = 0; j < r; j += 2) {
            if (locator_poly[j + 1] == 0) {
                continue;
            }
            assert(locator_poly[j + 1] == 1);

            q ^= pow_table[(j * (N - pos)) % N];
        }

        forney_coefs[i] = gf_div_ee(gf, p, q);
    }
}

_static void _rs_get_inf_symbols_syndrome_poly(const RS_t *rs,
                                               symbol_seq_t inf_symbols,
                                               const uint16_t *inf_positions,
                                               symbol_seq_t syndrome_poly) {
    assert(rs != NULL);
    assert(inf_positions != NULL);
    assert(inf_symbols.symbol_size == syndrome_poly.symbol_size);

    fft_transform(rs->fft, inf_symbols, inf_positions, syndrome_poly);
}

_static void _rs_get_repair_symbols_evaluator_poly(
    const RS_t *rs, symbol_seq_t syndrome_poly, const element_t *locator_poly,
    symbol_seq_t evaluator_poly) {
    assert(rs != NULL);
    assert(syndrome_poly.symbol_size == evaluator_poly.symbol_size);

    size_t symbol_size = syndrome_poly.symbol_size;
    uint16_t r = syndrome_poly.seq_length;

    // Evaluator polynomial initialization.
    for (uint16_t i = 0; i < r; ++i) {
        for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
            evaluator_poly.symbols[i].data[e_idx] = 0;
        }
    }
    for (uint16_t i = 0; i < r; ++i) {
        if (locator_poly[i] == 0) {
            continue;
        }
        assert(locator_poly[i] == 1);

        for (uint16_t j = 0; j < r - i; ++j) {
            for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
                evaluator_poly.symbols[i + j].data[e_idx] ^=
                    syndrome_poly.symbols[j].data[e_idx];
            }
        }
    }
}

_static void _rs_get_repair_symbols(const RS_t *rs,
                                    const element_t *forney_coefs,
                                    symbol_seq_t evaluator_poly,
                                    symbol_seq_t rep_symbols,
                                    const uint16_t *rep_positions) {
    assert(rs != NULL);
    assert(forney_coefs != NULL);
    assert(rep_positions != NULL);
    assert(evaluator_poly.symbol_size == rep_symbols.symbol_size);

    GF_t *gf = rs->gf;
    size_t symbol_size = evaluator_poly.symbol_size;
    element_t coef;

    fft_selective_transform(rs->fft, evaluator_poly, rep_symbols,
                            rep_positions);

    // TODO: implement deffered Forney scalling
    for (uint16_t i = 0; i < rep_symbols.seq_length; ++i) {
        coef = forney_coefs[i];

        for (size_t e_idx = 0; e_idx < symbol_size; ++e_idx) {
            rep_symbols.symbols[i].data[e_idx] =
                gf_mul_ee(gf, coef, rep_symbols.symbols[i].data[e_idx]);
        }
    }
}

int rs_generate_repair_symbols(const RS_t *rs, symbol_seq_t inf_symbols,
                               symbol_seq_t rep_symbols) {
    assert(rs != NULL);
    assert(inf_symbols.seq_length + rep_symbols.seq_length <= N);

    CC_t *cc = rs->cc;
    int ret;
    uint16_t k = inf_symbols.seq_length;
    uint16_t r = rep_symbols.seq_length;
    uint16_t inf_max_cnt;
    uint16_t rep_max_cnt;
    uint16_t inf_cosets_cnt;
    uint16_t rep_cosets_cnt;
    coset_t *inf_cosets;
    coset_t *rep_cosets;
    uint16_t *inf_positions;
    uint16_t *rep_positions;
    element_t *locator_poly;
    element_t *forney_coefs;
    symbol_seq_t syndrome_poly;
    symbol_seq_t evaluator_poly;

    cc_estimate_cosets_cnt(k, r, &inf_max_cnt, &rep_max_cnt);

    inf_cosets = (coset_t *)malloc(inf_max_cnt * sizeof(coset_t));
    if (!inf_cosets) {
        return 1;
    }

    rep_cosets = (coset_t *)malloc(rep_cosets_cnt * sizeof(coset_t));
    if (!rep_cosets) {
        free(inf_cosets);
        return 1;
    }

    inf_positions = (uint16_t *)malloc(k * sizeof(uint16_t));
    if (!inf_positions) {
        free(rep_cosets);
        free(inf_cosets);
        return 1;
    }

    rep_positions = (uint16_t *)malloc(r * sizeof(uint16_t));
    if (!rep_positions) {
        free(inf_positions);
        free(rep_cosets);
        free(inf_cosets);
        return 1;
    }

    locator_poly = (element_t *)malloc((r + 1) * sizeof(element_t));
    if (!locator_poly) {
        free(rep_positions);
        free(inf_positions);
        free(rep_cosets);
        free(inf_cosets);
        return 1;
    }

    forney_coefs = (element_t *)malloc(r * sizeof(element_t));
    if (!forney_coefs) {
        free(locator_poly);
        free(rep_positions);
        free(inf_positions);
        free(rep_cosets);
        free(inf_cosets);
        return 1;
    }

    ret = alloc_seq(SYMBOL_SIZE, r, &syndrome_poly);
    if (ret) {
        free(forney_coefs);
        free(locator_poly);
        free(rep_positions);
        free(inf_positions);
        free(rep_cosets);
        free(inf_cosets);
        return ret;
    }

    ret = alloc_seq(SYMBOL_SIZE, r, &evaluator_poly);
    if (ret) {
        free_seq(&syndrome_poly);
        free(forney_coefs);
        free(locator_poly);
        free(rep_positions);
        free(inf_positions);
        free(rep_cosets);
        free(inf_cosets);
        return ret;
    }

    cc_select_cosets(cc, k, r, inf_cosets, inf_max_cnt, &inf_cosets_cnt,
                     rep_cosets, rep_max_cnt, &rep_cosets_cnt);

    cc_cosets_to_positions(inf_cosets, inf_cosets_cnt, inf_positions, k);
    cc_cosets_to_positions(rep_cosets, rep_cosets_cnt, rep_positions, r);

    _rs_get_rep_symbols_locator_poly(rs, r, rep_cosets, rep_cosets_cnt,
                                     locator_poly, r + 1);

    _rs_get_enc_forney_coefs(rs, r, locator_poly, rep_positions, forney_coefs);

    _rs_get_inf_symbols_syndrome_poly(rs, inf_symbols, inf_positions,
                                      syndrome_poly);

    _rs_get_repair_symbols_evaluator_poly(rs, syndrome_poly, locator_poly,
                                          evaluator_poly);

    _rs_get_repair_symbols(rs, forney_coefs, evaluator_poly, rep_symbols,
                           rep_positions);

    free_seq(&evaluator_poly);
    free_seq(&syndrome_poly);
    free(forney_coefs);
    free(locator_poly);
    free(rep_positions);
    free(inf_positions);
    free(rep_cosets);
    free(inf_cosets);

    return 0;
}