/**
 * @file reed_solomon.c
 * @author Matvey Kolesov (kolesov645@gmail.com)
 * @brief rs/reed_solomon.h implementation.
 * @date 2024-01-21
 *
 * @copyright Copyright (c) 2024
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <rs/fft.h>
#include <rs/reed_solomon.h>

RS_t* rs_create() {
    RS_t* rs;

    rs = (RS_t*)malloc(sizeof(RS_t));
    if (!rs)
        return NULL;
    memset((void*)rs, 0, sizeof(RS_t));

    rs->gf = gf_create();
    if (!rs->gf) {
        free(rs);
        return NULL;
    }

    rs->cc = cc_create();
    if (!rs->cc) {
        gf_destroy(rs->gf);
        free(rs);
        return NULL;
    }

    rs->fft = fft_create();
    if (!rs->fft) {
        cc_destroy(rs->cc);
        gf_destroy(rs->gf);
        free(rs);
        return NULL;
    }

    return rs;
}

void rs_destroy(RS_t* rs) {
    assert(rs != NULL);

    fft_destroy(rs->fft);
    cc_destroy(rs->cc);
    gf_destroy(rs->gf);
    free(rs);
}

/**
 * @brief Compute syndrome polynomial.
 *
 * @param rs context object.
 * @param seq symbol sequence.
 * @param positions symbol positions.
 * @param syndrome_poly where to place the result.
 * @return 0 on success,\n
 *         !0 on error.
 */
static int _rs_get_syndrome_poly(RS_t* rs, const symbol_seq_t* seq, const uint16_t* positions,
                                 symbol_seq_t* syndrome_poly) {
    assert(rs != NULL);
    assert(seq != NULL);
    assert(positions != NULL);
    assert(syndrome_poly != NULL);

    int err;

    // fft_transform(rs->gf, seq, positions, syndrome_poly);
    err = fft_transform_cycl(rs->fft, rs->gf, seq, positions, syndrome_poly);
    if (err)
        return err;

    return 0;
}

/**
 * @brief Compute locator polynomial.
 *
 * @param rs context object.
 * @param positions positions.
 * @param positions_cnt number of positions.
 * @param locator_poly where to place locator polynomial coefficients.
 * @param locator_max_len max number of locator polynomial coefficients.
 */
static void _rs_get_locator_poly(const RS_t* rs, const uint16_t* positions, uint16_t positions_cnt,
                                 element_t* locator_poly, uint16_t locator_max_len) {
    assert(rs != NULL);
    assert(positions != NULL);
    assert(locator_poly != NULL);
    assert(positions_cnt + 1 <= locator_max_len);

    GF_t* gf = rs->gf;
    element_t* pow_table = gf->pow_table;

    locator_poly[0] = 1;
    for (uint16_t d = 0; d < positions_cnt; ++d) {
        locator_poly[d + 1] = 0;
        uint16_t pos = positions[d];
        element_t coef = pow_table[pos];

        for (uint16_t i = d + 1; i > 0; --i)
            locator_poly[i] ^= gf_mul_ee(gf, locator_poly[i - 1], coef);
    }
}

/**
 * @brief Compute repair symbols locator polynomial.
 * @details All locator polynomial coefficients will belongs to GF(2) subfield
 * ({0, 1}) of GF(65536). Locator polynomial will have degree equal to number of
 * repair symbols.
 *
 * @param rs context object.
 * @param r number of repair symbols.
 * @param rep_cosets cyclotomic cosets that form repair symbol positions.
 * @param rep_cosets_cnt number of cyclotomic cosets.
 * @param locator_poly where to place locator polynomial coefficients.
 * @param locator_max_len max number of locator polynomial coefficients.
 */
static void _rs_get_rep_symbols_locator_poly(const RS_t* rs, uint16_t r, const coset_t* rep_cosets,
                                             uint16_t rep_cosets_cnt, element_t* locator_poly,
                                             uint16_t locator_max_len) {
    assert(rs != NULL);
    assert(rep_cosets != NULL);
    assert(locator_poly != NULL);
    assert(r + 1 <= locator_max_len);

#ifndef NDEBUG
    uint16_t _r1 = 0;
    for (uint16_t _i = 0; _i < rep_cosets_cnt; ++_i) {
        _r1 += rep_cosets[_i].size;
    }
    assert(_r1 == r);
#endif

    element_t coset_locator_poly[RS_COSET_LOCATOR_MAX_LEN] = {0};
    coset_t coset;
    uint16_t coset_elements[CC_MAX_COSET_SIZE];
    uint16_t d; // locator polynomial degree

    d = 0;
    memset((void*)locator_poly, 0, (r + 1) * sizeof(element_t));
    locator_poly[0] = 1;

    for (uint16_t coset_idx = 0; coset_idx < rep_cosets_cnt; ++coset_idx) {
        coset = rep_cosets[coset_idx];

        coset_elements[0] = coset.leader;
        for (uint16_t i = 1; i < coset.size; ++i)
            coset_elements[i] = NEXT_COSET_ELEMENT(coset_elements[i - 1]);

        _rs_get_locator_poly(rs, coset_elements, coset.size, coset_locator_poly,
                             RS_COSET_LOCATOR_MAX_LEN);

#ifndef NDEBUG
        for (uint16_t _i = 0; _i <= coset.size; ++_i)
            assert(coset_locator_poly[_i] == 0 || coset_locator_poly[_i] == 1);
#endif

        for (uint16_t i = d;; --i) {
            if (locator_poly[i] == 1) {
                for (uint16_t j = 1; j <= coset.size; ++j)
                    locator_poly[i + j] ^= coset_locator_poly[j];
            }
            if (i == 0)
                break;
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

/**
 * @brief Compute Forney coefficient for given symbol postion.
 *
 * @param rs context object.
 * @param locator_poly locator polynomial.
 * @param d degree of locator polynomial (number of repair symbols or erasures).
 * @param pos symbol position.
 * @return Forney coefficient.
 */
static element_t _rs_get_forney_coef(const RS_t* rs, const element_t* locator_poly, uint16_t d,
                                     uint16_t pos) {
    assert(rs != NULL);
    assert(locator_poly != NULL);

    GF_t* gf = rs->gf;
    element_t* pow_table = gf->pow_table;
    element_t p; // divisible element = alpha^{position}
    element_t q; // divisor = locator_poly'(alpha^{-position})

    p = pow_table[pos];
    q = 0;
    for (uint16_t j = 0; j < d; j += 2) {
        element_t coef = locator_poly[j + 1];

        if (coef == 0) {
            continue;
        } else if (coef == 1) {
            q ^= pow_table[(j * (N - pos)) % N];
        } else {
            q ^= gf_mul_ee(gf, pow_table[(j * (N - pos)) % N], coef);
        }
    }

    return gf_div_ee(gf, p, q);
}

/**
 * @brief Compute evaluator polynomial modulo x^t (t - number of repair symbols
 * or erasures).
 *
 * @param rs context object.
 * @param syndrome_poly information symbols syndrome polynomial (deg == t - 1).
 * @param locator_poly repair symbols locator polynomial (deg == t).
 * @param evaluator_poly where to place the result.
 */
static void _rs_get_evaluator_poly(const RS_t* rs, const symbol_seq_t* syndrome_poly,
                                   const element_t* locator_poly, symbol_seq_t* evaluator_poly) {
    assert(rs != NULL);
    assert(syndrome_poly != NULL);
    assert(locator_poly != NULL);
    assert(evaluator_poly != NULL);
    assert(syndrome_poly->symbol_size == evaluator_poly->symbol_size);

    GF_t* gf = rs->gf;
    size_t symbol_size = syndrome_poly->symbol_size;
    uint16_t r = syndrome_poly->length;

    // Evaluator polynomial initialization.
    for (uint16_t i = 0; i < r; ++i)
        memset((void*)evaluator_poly->symbols[i]->data, 0, symbol_size);

    for (uint16_t i = 0; i < r; ++i) {
        element_t coef = locator_poly[i];

        if (coef == 0)
            continue;

        for (uint16_t j = 0; j < r - i; ++j)
            gf_madd(gf, (void*)evaluator_poly->symbols[i + j]->data, coef,
                    (void*)syndrome_poly->symbols[j]->data, symbol_size);
    }
}

/**
 * @brief Compute repair symbols.
 *
 * @param rs context object.
 * @param locator_poly repair symbols locator polynomial.
 * @param evaluator_poly repair symbols evaluator polynomial.
 * @param rep_positions repair symbol positions.
 * @param rep_cosets repair symbol positions in form of cyclotomic cosets union.
 * @param rep_cosets_cnt number of repair symbol cyclotomic cosets.
 * @param rep_symbols where to place the result.
 * @return 0 on success,\n
 *         !0 on error.
 */
static int _rs_get_repair_symbols(const RS_t* rs, const element_t* locator_poly,
                                  const symbol_seq_t* evaluator_poly, const uint16_t* rep_positions,
                                  const coset_t* rep_cosets, uint16_t rep_cosets_cnt,
                                  symbol_seq_t* rep_symbols) {
    assert(rs != NULL);
    assert(locator_poly != NULL);
    assert(evaluator_poly != NULL);
    assert(rep_positions != NULL);
    assert(rep_cosets != NULL);
    assert(rep_symbols != NULL);
    assert(evaluator_poly->symbol_size == rep_symbols->symbol_size);

    GF_t* gf = rs->gf;
    size_t symbol_size = evaluator_poly->symbol_size;
    uint16_t r = rep_symbols->length;
    int err;

    // fft_partial_transform(gf, evaluator_poly, rep_symbols,
    // rep_positions);
    err = fft_partial_transform_cycl(gf, evaluator_poly, rep_cosets, rep_cosets_cnt, rep_symbols);
    if (err)
        return err;

    // TODO: implement deffered Forney scalling
    for (uint16_t i = 0; i < r; ++i) {
        element_t coef = _rs_get_forney_coef(rs, locator_poly, r, rep_positions[i]);
        gf_mul(gf, (void*)rep_symbols->symbols[i]->data, coef, symbol_size);
    }

    return 0;
}

/**
 * @brief Restore erased symbols if it is possible.
 *
 * @param rs context object.
 * @param k number of information symbols.
 * @param locator_poly erased symbols locator polynomial.
 * @param evaluator_poly erased symbols evaluator polynomial.
 * @param positions positions of all symbols.
 * @param is_erased indicates which symbols has been erased.
 * @param rcv_symbols received symbols, restored symbols will be written here.
 */
static void _rs_restore_erased(const RS_t* rs, uint16_t k, const element_t* locator_poly,
                               const symbol_seq_t* evaluator_poly, const uint16_t* positions,
                               const bool* is_erased, symbol_seq_t* rcv_symbols) {
    assert(rs != NULL);
    assert(locator_poly != NULL);
    assert(evaluator_poly != NULL);
    assert(positions != NULL);
    assert(is_erased != NULL);
    assert(rcv_symbols != NULL);
    assert(evaluator_poly->symbol_size == rcv_symbols->symbol_size);

    GF_t* gf = rs->gf;
    size_t symbol_size = evaluator_poly->symbol_size;
    element_t* pow_table = gf->pow_table;
    element_t forney_coef;
    element_t coef;
    uint16_t t = evaluator_poly->length;
    uint16_t pos;
    uint16_t j;

    for (uint16_t id = 0; id < k; ++id) {
        if (!is_erased[id])
            continue;

        pos = positions[id];
        forney_coef = _rs_get_forney_coef(rs, locator_poly, t, pos);

        memset((void*)rcv_symbols->symbols[id]->data, 0, symbol_size);

        j = (N - pos) % N;

        for (uint16_t i = 0; i < t; ++i) {
            coef = gf_mul_ee(gf, forney_coef, pow_table[(i * j) % N]);
            gf_madd(gf, (void*)rcv_symbols->symbols[id]->data, coef,
                    (void*)evaluator_poly->symbols[i]->data, symbol_size);
        }
    }
}

int rs_generate_repair_symbols(RS_t* rs, const symbol_seq_t* inf_symbols,
                               symbol_seq_t* rep_symbols) {
    assert(rs != NULL);
    assert(inf_symbols != NULL);
    assert(rep_symbols != NULL);
    assert(inf_symbols->length + rep_symbols->length <= N);
    assert(inf_symbols->symbol_size == rep_symbols->symbol_size);

    CC_t* cc = rs->cc;
    size_t symbol_size = inf_symbols->symbol_size;
    uint16_t k = inf_symbols->length;
    uint16_t r = rep_symbols->length;
    uint16_t inf_max_cnt = 0;
    uint16_t rep_max_cnt = 0;
    uint16_t inf_cosets_cnt = 0;
    uint16_t rep_cosets_cnt = 0;
    coset_t* _cosets;
    coset_t* inf_cosets;
    coset_t* rep_cosets;
    uint16_t* positions;
    uint16_t* inf_positions;
    uint16_t* rep_positions;
    element_t* locator_poly;
    symbol_seq_t* syndrome_poly;
    symbol_seq_t* evaluator_poly;
    int err;

    cc_estimate_cosets_cnt(k, r, &inf_max_cnt, &rep_max_cnt);

    _cosets = (coset_t*)calloc(inf_max_cnt + rep_max_cnt, sizeof(coset_t));
    if (!_cosets)
        return 1;
    inf_cosets = _cosets;
    rep_cosets = _cosets + inf_max_cnt;

    positions = (uint16_t*)calloc(k + r, sizeof(uint16_t));
    if (!positions) {
        free(_cosets);
        return 1;
    }
    inf_positions = positions;
    rep_positions = positions + k;

    locator_poly = (element_t*)calloc(r + 1, sizeof(element_t));
    if (!locator_poly) {
        free(positions);
        free(_cosets);
        return 1;
    }

    syndrome_poly = seq_create(r, symbol_size);
    if (!syndrome_poly) {
        free(locator_poly);
        free(positions);
        free(_cosets);
        return 1;
    }

    evaluator_poly = seq_create(r, symbol_size);
    if (!evaluator_poly) {
        seq_destroy(syndrome_poly);
        free(locator_poly);
        free(positions);
        free(_cosets);
        return 1;
    }

    cc_select_cosets(cc, k, r, inf_cosets, inf_max_cnt, &inf_cosets_cnt, rep_cosets, rep_max_cnt,
                     &rep_cosets_cnt);

    cc_cosets_to_positions(inf_cosets, inf_cosets_cnt, inf_positions, k);
    cc_cosets_to_positions(rep_cosets, rep_cosets_cnt, rep_positions, r);

    err = _rs_get_syndrome_poly(rs, inf_symbols, inf_positions, syndrome_poly);
    if (err) {
        seq_destroy(evaluator_poly);
        seq_destroy(syndrome_poly);
        free(locator_poly);
        free(positions);
        free(_cosets);
        return err;
    }

    _rs_get_rep_symbols_locator_poly(rs, r, rep_cosets, rep_cosets_cnt, locator_poly, r + 1);

    _rs_get_evaluator_poly(rs, syndrome_poly, locator_poly, evaluator_poly);

    err = _rs_get_repair_symbols(rs, locator_poly, evaluator_poly, rep_positions, rep_cosets,
                                 rep_cosets_cnt, rep_symbols);
    if (err) {
        seq_destroy(evaluator_poly);
        seq_destroy(syndrome_poly);
        free(locator_poly);
        free(positions);
        free(_cosets);
        return err;
    }

    seq_destroy(evaluator_poly);
    seq_destroy(syndrome_poly);
    free(locator_poly);
    free(positions);
    free(_cosets);

    return 0;
}

int rs_restore_symbols(RS_t* rs, uint16_t k, uint16_t r, symbol_seq_t* rcv_symbols,
                       const bool* is_erased, uint16_t t) {
    assert(rs != NULL);
    assert(rcv_symbols != NULL);
    assert(is_erased != NULL);
    assert((k + r) == rcv_symbols->length);

    CC_t* cc = rs->cc;
    size_t symbol_size = rcv_symbols->symbol_size;
    uint16_t inf_max_cnt = 0;
    uint16_t rep_max_cnt = 0;
    uint16_t inf_cosets_cnt = 0;
    uint16_t rep_cosets_cnt = 0;
    coset_t* _cosets;
    coset_t* inf_cosets;
    coset_t* rep_cosets;
    uint16_t* positions;
    uint16_t* inf_positions;
    uint16_t* rep_positions;
    uint16_t* erased_positions;
    element_t* locator_poly;
    symbol_seq_t* syndrome_poly;
    symbol_seq_t* evaluator_poly;
    int err;

    if (r < t) {
        // Too many erases - symbols cannot be restored.
        return RS_ERR_CANNOT_RESTORE;
    }

    cc_estimate_cosets_cnt(k, r, &inf_max_cnt, &rep_max_cnt);

    _cosets = (coset_t*)calloc(inf_max_cnt + rep_max_cnt, sizeof(coset_t));
    if (!_cosets)
        return 1;
    inf_cosets = _cosets;
    rep_cosets = _cosets + inf_max_cnt;

    positions = (uint16_t*)calloc(k + r, sizeof(uint16_t));
    if (!positions) {
        free(_cosets);
        return 1;
    }
    inf_positions = positions;
    rep_positions = positions + k;

    erased_positions = (element_t*)calloc(t, sizeof(uint16_t));
    if (!erased_positions) {
        free(positions);
        free(_cosets);
        return 1;
    }

    locator_poly = (element_t*)calloc(t + 1, sizeof(element_t));
    if (!locator_poly) {
        free(erased_positions);
        free(positions);
        free(_cosets);
        return 1;
    }

    syndrome_poly = seq_create(t, symbol_size);
    if (!syndrome_poly) {
        free(locator_poly);
        free(erased_positions);
        free(positions);
        free(_cosets);
        return 1;
    }

    evaluator_poly = seq_create(t, symbol_size);
    if (!evaluator_poly) {
        seq_destroy(syndrome_poly);
        free(locator_poly);
        free(erased_positions);
        free(positions);
        free(_cosets);
        return 1;
    }

    cc_select_cosets(cc, k, r, inf_cosets, inf_max_cnt, &inf_cosets_cnt, rep_cosets, rep_max_cnt,
                     &rep_cosets_cnt);

    cc_cosets_to_positions(inf_cosets, inf_cosets_cnt, inf_positions, k);
    cc_cosets_to_positions(rep_cosets, rep_cosets_cnt, rep_positions, r);

    err = _rs_get_syndrome_poly(rs, rcv_symbols, positions, syndrome_poly);
    if (err) {
        seq_destroy(evaluator_poly);
        seq_destroy(syndrome_poly);
        free(locator_poly);
        free(erased_positions);
        free(positions);
        free(_cosets);
        return err;
    }

    uint16_t idx = 0;
    for (uint16_t i = 0; i < k + r; ++i) {
        if (!is_erased[i])
            continue;
        erased_positions[idx++] = positions[i];
    }

    _rs_get_locator_poly(rs, erased_positions, t, locator_poly, t + 1);

    _rs_get_evaluator_poly(rs, syndrome_poly, locator_poly, evaluator_poly);

    _rs_restore_erased(rs, k, locator_poly, evaluator_poly, positions, is_erased, rcv_symbols);

    seq_destroy(evaluator_poly);
    seq_destroy(syndrome_poly);
    free(locator_poly);
    free(erased_positions);
    free(positions);
    free(_cosets);

    return 0;
}