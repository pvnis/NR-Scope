/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
/*
 * \brief One snapshot per stream per slot from every grant of the slot (sensing.join_grants).
 *
 * The PDSCH DM-RS of a symbol is one sequence over the whole carrier, indexed by
 * absolute CRB (38.211 7.4.1.1), so the despread estimates of several UEs' grants in
 * the same symbol sit on one pilot lattice and can be laid side by side. Without this,
 * every grant was pushed on its own and dropped when narrower than half the carrier:
 * four UEs of 68 PRBs each in one slot gave nothing.
 *
 * Per slot and worker: begin, then every grant's estimates are added (from any of the
 * worker's decoder threads), then the worker flushes once all of them are done.
 *
 * What may share a buffer (the key): receive chain, DM-RS port bitmap, layer, and the
 * DM-RS power offset (beta). Ports fix the comb (k_step, k_offset). Rows are per OFDM
 * symbol, so a grant with an extra DM-RS symbol only adds to that row. On an overlap
 * the first grant's estimates are kept.
 *
 * At flush, per row: the grants are walked in frequency. Where two grants of different
 * RNTIs meet, the channel has to be continuous across the junction (same effective
 * channel: an unprecoded cell, or the same precoder): a line through one PRB of pilots
 * on each side, the row's phase slope removed, extrapolated to the junction, must agree
 * within the lines' own noise plus join_max_step of |H|. A gap of at most
 * join_max_gap_prb (0 = none) between grants is filled from those lines, under the
 * same check. At low SNR one junction cannot tell two precoders apart, so the cell is
 * judged too (join_learn): when the clean junctions between UEs mostly fail, the cell
 * beamforms each UE and UEs are no longer joined; that verdict is printed and revised. The widest
 * run of joined pilots is the row's snapshot, kept when it covers half the carrier, as
 * before. A single full-band grant therefore gives the very same snapshot as before.
 *
 * Every row, junction and gap is counted (nrscope_sensing_join_report()) and, with a
 * map dump, written to <dump>.join.csv: what a commercial scheduler does with its
 * PRBs (gaps, adjacency, precoding steps) is read from there.
 */
#ifndef NRSCOPE_SENSING_JOIN_H
#define NRSCOPE_SENSING_JOIN_H

#include <stdbool.h>
#include <stdint.h>

#include "srsran/config.h"

#ifdef __cplusplus
extern "C" {
#endif

struct nrscope_sensing_s;
typedef struct nrscope_sensing_join_s nrscope_sensing_join_t;

/* Knobs, the sensing: block of the yaml (load_config.cc). */
typedef struct {
  bool  enable;        // join_grants, default true
  int   max_gap_prb;   // join_max_gap_prb, default 0: grants must touch
  bool  check;         // join_check, default true: continuity test at every junction
  float max_step;      // join_max_step, default 0.3: step allowed past the noise, share of |H|
  bool  learn;         // join_learn, default true: stop joining UEs in a cell that precodes each
} nrscope_sensing_join_args_t;

extern nrscope_sensing_join_args_t nrscope_sensing_join_args;

/* One per worker. grid_crb0 / n_sc_grid: the frame every row is kept in, the carrier
   (offsetToCarrier, carrierBandwidth * 12), which is the DCI decoders' sensing grid. */
nrscope_sensing_join_t* nrscope_sensing_join_alloc(uint32_t grid_crb0, uint32_t n_sc_grid);
void                    nrscope_sensing_join_free(nrscope_sensing_join_t* j);

/* Start a slot: forget the previous one's rows. */
void nrscope_sensing_join_begin(nrscope_sensing_join_t* j, uint32_t sfn, uint32_t slot_idx);

/* One grant's estimates of one layer on one chain, rows of dmrs_mask, in the caller's
   grid (n_sc_grid subcarriers from CRB grid_crb0). Thread safe. Returns 1 when taken,
   0 when it cannot be (no slot open, a grid outside the frame, no free key): the caller
   then pushes it on its own as before. */
int nrscope_sensing_join_add(nrscope_sensing_join_t* j,
                             uint32_t                aarx,
                             uint16_t                ports,
                             int                     layer,
                             int                     n_layers,
                             float                   beta,
                             uint16_t                rnti,
                             uint16_t                dmrs_mask,
                             const cf_t*             H,
                             const bool*             valid,
                             uint32_t                n_sc_grid,
                             uint32_t                grid_crb0,
                             uint32_t                sfn,
                             uint32_t                slot_idx,
                             uint64_t                slot_abs);

/* The sensing context the grants go to (nrscope_sensing_process_grant_join sets it). */
void nrscope_sensing_join_set_context(nrscope_sensing_join_t* j, struct nrscope_sensing_s* s);

/* End the slot: build every row and push them, chain 0 first. s NULL: the context the
   grants came with. Returns snapshots pushed. */
int nrscope_sensing_join_flush(nrscope_sensing_join_t* j, struct nrscope_sensing_s* s);

/* Totals over every worker, printed at exit. */
void nrscope_sensing_join_report(void);

/* The cell's verdict (join_learn): 0 undecided, 1 coherent across UEs, 2 precoded per UE.
   And, for the tests, forget every total and the verdict. */
int  nrscope_sensing_join_cell(void);
void nrscope_sensing_join_reset(void);

/* For the tests: the row a flush builds, without the sensing context. Runs the same
   code as nrscope_sensing_join_flush() but hands each key's rows to the callback. */
typedef void (*nrscope_sensing_join_sink_t)(void*       ctx,
                                            uint32_t    aarx,
                                            uint16_t    ports,
                                            int         layer,
                                            uint16_t    mask,
                                            const cf_t* H,
                                            const bool* valid,
                                            uint32_t    n_sc);
int nrscope_sensing_join_flush_to(nrscope_sensing_join_t* j, nrscope_sensing_join_sink_t sink, void* ctx);

/* Continuity across a junction of a row (exposed for the tests): the pilots ending at kl
   and the ones starting at kr (subcarriers, same comb), m of them on each side, the row's
   phase slope ph (radians per lattice step) removed. Each side's line is extrapolated to
   the middle; the step between them passes when within the fits' noise (chi-square 99%)
   plus join_max_step of |H|. *db, *deg: the step; *informative: the fits are clean
   enough to see a precoder change. Returns 1 pass, 0 cut, -1 no pilots. */
int nrscope_sensing_join_test(const cf_t* H,
                              const bool* valid,
                              int         n_sc,
                              int         kl,
                              int         kr,
                              int         k_step,
                              int         m,
                              float       ph,
                              float*      db,
                              float*      deg,
                              bool*       informative);

#ifdef __cplusplus
}
#endif

#endif
