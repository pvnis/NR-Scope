/*
 * SPDX-License-Identifier: MIT
 *
 * GoldenSniffer port: soft-information error correction of the GF(2) DCI decoder,
 * ported from the LGS ("lineares Gleichungssystem") error correction of 5GDescrambler
 * (https://github.com/kit-ps/5GDescrambler, pdcch/src/lgs_error_correction.rs and the
 * error budgets of pdcch_macros/src/decode_mat.rs).
 *
 * Algorithm (5GDescrambler, mapped onto the GoldenSniffer system x G = y of gs_gf2.c):
 *   1. hard decision; if the GF(2) system is consistent, done (exactly gs_dci_try).
 *   2. otherwise the t least reliable bits (smallest |LLR|) become unknowns next to
 *      the payload and the RNTI. Each parity check of the code (a non-pivot bit of the
 *      reduced generator) is one linear equation in them; Gaussian elimination
 *      (unknowns in ascending reliability) solves for the error pattern on those bits.
 *      Unknowns the equations do not determine ("free") keep their hard decision.
 *   3. the corrected word must satisfy every parity check; payload and RNTI are then
 *      read from it as in the hard decoder.
 *   4. confidence: sum over the t bits of +|LLR| (kept), -|LLR| (flipped or free),
 *      divided by 2 R (5GDescrambler's "simplified linear SNR", info.metric_raw); the
 *      returned metric divides the sum by t mean|LLR| instead, so it does not depend on
 *      the LLR scale (~0 for a random solution, ~0.3-0.6 for typical corrections).
 *
 * Mapping and adaptations (see gs_ec.c for detail):
 *   - 5GDescrambler does not know N_ID and solves for c_init too (K + 32 unknowns);
 *     here N_ID comes from the PDCCH DM-RS, so the unknowns are payload + RNTI
 *     (N + 16) and every one of the remaining R - N - 16 bits is a parity check. The
 *     RNTI MSB is therefore never "free" (5GDescrambler's NUM_FREE = 1) and the
 *     syndrome-pattern second solve of its NUM_FREE > 1 path has no counterpart.
 *   - Repetition (E > N, AL8/16): like 5GDescrambler, only the first R = N bits (one
 *     copy of every coded bit) enter the soft solve ("rate matching limit").
 *   - Error budget t: 5GDescrambler's 10% of E (n_RNTI = RNTI) / 5% of E (n_RNTI = 0),
 *     capped so that at least min_parity checks stay unused by the unknowns (its AL1
 *     cap "errors <= 97 - K" plays the same role for its 32 extra unknowns).
 *   - Acceptance: 5GDescrambler accepts any consistent solution at AL >= 2 and at AL1
 *     only "errors + free < 2" (MAX_ERRORS_AND_FREE_AL_1). Here a solution is accepted
 *     when the effective parity  Q - log2(sum_{i<=w} C(r, i))  is >= min_parity
 *     (Q checks, r = rank of the unknowns, w = bits flipped): a random input passes
 *     with probability <= 2^-min_parity, which bounds false alarms the same way at
 *     every AL and size (a bound for inputs whose syndrome is random: noise, other
 *     sizes and windows). With n_RNTI = RNTI the RNTI must be in 1..0xFFF2 (the
 *     5GDescrambler plausibility test).
 *   - Added (not in 5GDescrambler): RNTI bits that do not enter the scrambling (all 16
 *     with n_RNTI = 0, the MSB with n_RNTI = RNTI; the latter is 5GDescrambler's free
 *     "RNTI MSB") have low-weight generator rows (4..32 bits), so the corrected word
 *     plus such a row is a nearby codeword with another RNTI. The decode is rejected
 *     when that neighbour agrees with the LLRs at least as well (min_nb_margin): this
 *     removed the large-budget miscorrections with n_RNTI = RNTI, at ~1e-4 of the
 *     decodes. 5GDescrambler instead leaves the MSB free and sets it from the syndrome.
 *
 * LLR convention: llr[i] > 0 means bit i is more likely 0 (hard bit = llr < 0), bit
 * order as gs_dci_try (PDCCH REs in mapping order, Re then Im). Any positive scaling
 * works (only signs and the order of |llr| are used; metric is scale free). From
 * ZF-equalised QPSK with noise variance sigma^2 per RE: llr = 2 sqrt(2) Re/Im(y) / sigma^2;
 * where the channel gain varies across the candidate, weighting by |h|^2 (i.e.
 * llr = 2 sqrt(2) Re/Im(y conj(h)) / sigma_n^2 before equalisation) orders |llr| better.
 */
#ifndef NRSCOPE_GS_EC_H
#define NRSCOPE_GS_EC_H

#include "nrscope/hdr/golden_sniffer/gs_core.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GS_EC_EFFORT_AUTO UINT32_MAX /* 5GDescrambler error budget */
/* Largest budget the acceptance allows, t = Q - min_parity (<= GS_EC_T_MAX), with
 * n_RNTI = RNTI only. With n_RNTI = 0 the RNTI enters the CRC alone: its generator
 * rows are low-weight polar codewords, and large budgets then correct into the
 * neighbouring codeword (another RNTI) at low SNR (~0.1-0.2% of the corrected DCIs at
 * t >= 20% of E, none seen at 5GDescrambler's 5%), so that case keeps the AUTO budget. */
#define GS_EC_EFFORT_MAX (UINT32_MAX - 1)
#define GS_EC_T_MAX 192              /* largest error budget (AL16: 10% of 1728 = 172) */

typedef struct {
  float    frac_ue;    /* error budget, fraction of E, n_RNTI = RNTI (5GDescrambler: 0.10) */
  float    frac_css;   /* error budget, fraction of E, n_RNTI = 0 (5GDescrambler: 0.05) */
  uint32_t min_parity; /* least effective parity bits (hard path: the caller's gate N + 16 + min_parity <= E) */
  uint32_t max_flips;  /* most bits corrected, 0 = no cap (5GDescrambler AL1: 1) */
  float    min_metric; /* least metric, -INFINITY = off */
  float    min_nb_margin; /* least neighbour margin (info.nb_margin), default 0, -INFINITY = off */
  bool     hard_first; /* try gs_dci_try first (default true; false only for self-tests) */
} gs_ec_args_t;

typedef struct {
  bool     hard;       /* decoded by the hard decision */
  uint32_t R;          /* bits used by the soft solve (E, or N under repetition) */
  uint32_t Q;          /* parity checks: R - nof_bits - 16 */
  uint32_t t;          /* least reliable bits taken as unknowns */
  uint32_t rank;       /* of their parity-check rows */
  uint32_t flips;      /* bits corrected */
  float    eff_parity; /* Q - log2(sum_{i<=flips} C(rank, i)) */
  float    flip_llr;   /* sum of |llr| of the corrected bits */
  float    metric_raw; /* 5GDescrambler's "simplified linear SNR": sum over the t bits of
                          +|llr| (kept) / -|llr| (flipped or free), over 2 R */
  float    metric;     /* the same sum over t * mean|llr| of the R bits: scale free, ~0 for
                          a random solution, larger for few corrections of unreliable bits */
  float    nb_margin;  /* how much likelier the decode is than its nearest RNTI-only
                          neighbour (corrected word + one RNTI row), in mean |llr| */
} gs_ec_info_t;

void gs_ec_args_default(gs_ec_args_t* a);

/* Soft-decision GF(2) DCI decoder: same contract as gs_dci_try (nof_bits payload bits
 * in E = 108 L LLRs, n_rnti selects the scrambling), but from LLRs, correcting bit
 * errors when the hard decision is inconsistent. effort = number of least reliable
 * bits taken as unknowns (0: hard decision only, GS_EC_EFFORT_AUTO: 5GDescrambler's
 * budget, GS_EC_EFFORT_MAX: see above; any t is capped to Q - min_parity). Uses q->args.min_parity (q may be NULL: 16). Returns true with payload, rnti
 * and metric (may be NULL) on an accepted decode. Thread safe; tables are built on
 * first use of each (size, AL). */
bool gs_dci_try_soft(gs_t* q, const float* llr, uint32_t E, uint32_t nof_bits, uint32_t n_id, bool n_rnti,
                     uint32_t effort, uint8_t* payload, uint16_t* rnti, float* metric);

/* Same with explicit acceptance parameters; info may be NULL. */
bool gs_dci_try_soft_ex(const gs_ec_args_t* a, const float* llr, uint32_t E, uint32_t nof_bits, uint32_t n_id,
                        bool n_rnti, uint32_t effort, uint8_t* payload, uint16_t* rnti, gs_ec_info_t* info);

/* Builds the tables of one size for all ALs ahead of time (optional). */
void gs_ec_prepare(uint32_t nof_bits);

#ifdef __cplusplus
}
#endif

#endif /* NRSCOPE_GS_EC_H */
