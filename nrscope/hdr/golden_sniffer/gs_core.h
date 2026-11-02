/*
 * SPDX-License-Identifier: MIT
 *
 * GoldenSniffer blind PDCCH detection, ported to C for NR-Scope.
 *
 * Port of the MATLAB reference https://github.com/5G-GoldenSniffer/5G-GoldenSniffer
 * (Copyright 2024-2026 the authors, MIT License): Palama, Mangione, Dino, Focarelli,
 * Bianchi, "Golden Sniffer: Low-Complexity 5G New Radio Signal Interception Exploiting
 * Gold Sequences", MobiCom '26.
 *
 * What it does, per downlink slot, without knowing any RNTI or search space:
 *   1. power-clusters the busy subcarriers of the first PDCCH symbols into runs
 *      (power_detect_REGs.m) and enumerates PDCCH candidates from them
 *      (enum_candidates.m);
 *   2. removes the QPSK modulation by phase unwrapping (unwrap_QPSK.m) and reads the
 *      PDCCH DM-RS c_init, hence the scrambling ID N_ID, by inverting the Gold sequence
 *      (DMRS_process.m; the encrypted PDCCH_DMRS_solve_underdet.p is replaced by
 *      gs_dmrs_solve_underdet());
 *   3. equalizes with the recovered DM-RS (PDCCH_equalize.m) and recovers the DCI
 *      payload and the RNTI with one GF(2) linear solve per DCI size: with the
 *      scrambling known, polar code + CRC + RNTI mask are an affine map of
 *      [payload | RNTI] (make_G_dci.m, make_E_dci.m, DCIdecode_try.m).
 *
 * The GF(2) model is built from srsRAN's own PDCCH encoder chain, so it matches what
 * srsran_pdcch_nr_encode() transmits bit for bit.
 *
 * Grid convention: symbol l, subcarrier k at grid[l * nof_sc + k]. Two reference
 * points are given in grid subcarriers: point A (CRB0 subcarrier 0, the DM-RS
 * reference of every CORESET but CORESET#0) and the first subcarrier of CORESET#0
 * (its own DM-RS reference). Either may lie outside the grid (negative index).
 */
#ifndef NRSCOPE_GS_CORE_H
#define NRSCOPE_GS_CORE_H

#include "srsran/config.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GS_MAX_DURATION 3
#define GS_MAX_L_LOG2 4      /* aggregation levels 1..16 */
#define GS_MAX_DCI_BITS 140  /* polar K <= 164 for DCI */
#define GS_MAX_IV 3          /* runs (intervals) per candidate: interleaved triples */
#define GS_MAX_DCI_PER_SLOT 64
#define GS_MAX_SC (275 * 12)
#define GS_NO_REF INT32_MIN  /* reference point not available */

typedef struct {
  int32_t first; /* grid subcarrier, inclusive */
  int32_t last;  /* grid subcarrier, inclusive */
} gs_iv_t;

typedef struct {
  uint32_t nof_iv;
  gs_iv_t  iv[GS_MAX_IV];
} gs_cand_t;

typedef struct {
  uint32_t slot_idx;  /* slot within the frame */
  uint32_t duration;  /* detected CORESET duration (symbols) */
  uint32_t L;         /* aggregation level */
  gs_cand_t cand;     /* where it was found (grid subcarriers) */
  bool     coreset0;  /* DM-RS referenced to the CORESET#0 start (else point A) */
  bool     ref_ambiguous; /* both references coincide: coreset0 is a guess */
  uint32_t n_id;      /* PDCCH DM-RS scrambling ID read from the DM-RS */
  uint16_t rnti;      /* RNTI recovered from the CRC mask */
  bool     n_rnti;    /* scrambling used n_RNTI = RNTI (UE search space, pdcch-DMRS-ScramblingID) */
  uint32_t nof_bits;  /* DCI payload size */
  uint8_t  payload[GS_MAX_DCI_BITS]; /* a_0 .. a_{nof_bits-1}, one bit per byte */
  float    evm;       /* rms distance of the equalized REs from QPSK, relative */
  float    epre;      /* mean |Y|^2 of the candidate REs */
} gs_dci_t;

typedef struct {
  uint32_t max_duration;   /* PDCCH symbols scanned (1..3), default 2 */
  uint32_t size_min;       /* DCI sizes tried, default 20..84 */
  uint32_t size_max;
  const uint32_t* sizes;   /* if nof_sizes > 0, only these sizes */
  uint32_t nof_sizes;
  /* Sizes accepted when the scrambling is n_RNTI = 0 (common search spaces, or a UE
   * search space without pdcch-DMRS-ScramblingID). Such a decode is weak: polar codes
   * of different sizes are nested, so a size-N codeword also solves every larger
   * size with only 8 CRC bits left to check (1/256 per size). Sizes are tried in
   * ascending order, which leaves ~0.4% of these DCIs at a wrong size; listing the
   * sizes actually used (e.g. DCI 1_0 in CORESET#0) removes that. Applied to the
   * candidates that can be in CORESET#0: DM-RS referenced to it (or to point A when
   * the two coincide) and, when its width is known, lying inside it. 0 = no restriction. */
  const uint32_t* css_sizes;
  uint32_t nof_css_sizes;
  uint32_t max_L_log2;     /* default 4 (AL16) */
  /* Least number of observed bits beyond the N + 16 unknowns (payload, RNTI) for a
   * decode to count: a wrong candidate solves with probability 2^-(parity bits). The
   * MATLAB allows 8, which at AL1 and large sizes gave false DCIs; default 16. */
  uint32_t min_parity;
} gs_args_t;

typedef struct {
  const cf_t* grid;     /* symbol l at grid + l * nof_sc, l < max_duration */
  uint32_t    nof_sc;
  uint32_t    slot_idx; /* n_s,f^mu: slot index within the frame */
  int32_t     ref_crb0; /* grid subcarrier index of point A (usually <= 0), GS_NO_REF if unknown */
  int32_t     ref_cs0;  /* grid subcarrier index of the first CORESET#0 subcarrier, GS_NO_REF if none */
  uint32_t    cs0_nof_sc; /* CORESET#0 width in subcarriers; 0 if unknown */
  uint32_t    pci;
} gs_slot_t;

typedef struct gs_s gs_t;

/* Builds the shared GF(2) tables on the first call (thread safe), then a per-thread
 * searcher. Returns NULL on failure. */
gs_t* gs_create(const gs_args_t* args);
void  gs_destroy(gs_t* q);

/* Blind search of one slot. Returns the number of DCIs written to out (<= max_out). */
int gs_search_slot(gs_t* q, const gs_slot_t* s, gs_dci_t* out, uint32_t max_out);

/* Default arguments */
void gs_args_default(gs_args_t* args);

/* ---- building blocks, exposed for the self-test ---- */

/* Gold sequence c(n), n in [offset, offset + len), as one bit per byte. */
void gs_gold_bits(uint32_t c_init, uint32_t offset, uint32_t len, uint8_t* out);

/* PDCCH DM-RS c_init (TS 38.211 7.4.1.3.1) */
uint32_t gs_pdcch_dmrs_cinit(uint32_t slot_idx, uint32_t l, uint32_t n_id);

/* DMRS_process.m: c_est are the DM-RS REs up to a common QPSK ambiguity (on the
 * axes), D the bit offset of the first one from the reference. Returns the number
 * of rotations (0..4) whose bits satisfy the x2 recursion; c_init[i] and rot[i]
 * (the rotation index, exp(j(pi/4 + rot pi/2))) for each. Needs n >= 16. */
int gs_dmrs_process(const cf_t* c_est, uint32_t n, uint32_t D, uint32_t c_init[4], uint32_t rot[4]);

/* c_init_PDCCH_DMRS_decode.m */
bool gs_pdcch_cinit_check(uint32_t c_init, uint32_t slot_idx, uint32_t l, uint32_t* n_id);

/* Replacement of the encrypted PDCCH_DMRS_solve_underdet.p: the DM-RS of `dur`
 * symbols (n per symbol, column l at c_est + l * n, sharing one QPSK ambiguity),
 * each too short to invert alone, jointly determine N_ID. Returns true if exactly
 * one (N_ID, rotation) explains all symbols. n * 2 must be <= 64. */
bool gs_dmrs_solve_underdet(const cf_t* c_est, uint32_t n, uint32_t dur, uint32_t D, uint32_t slot_idx,
                            uint32_t* n_id, uint32_t* rot);

/* unwrap_QPSK.m */
void gs_unwrap_qpsk(const cf_t* y, uint32_t n, cf_t* y_corr);

/* GF(2) DCI decoder for nof_bits payload bits in E = 108 L hard bits (one per byte).
 * n_rnti selects the scrambling: false -> c_init = n_id (common search space),
 * true -> c_init = (RNTI << 16) + n_id. Returns true on a consistent solution. */
bool gs_dci_try(const uint8_t* bits, uint32_t E, uint32_t nof_bits, uint32_t n_id, bool n_rnti, uint8_t* payload,
                uint16_t* rnti);

#ifdef __cplusplus
}
#endif

#endif /* NRSCOPE_GS_CORE_H */
