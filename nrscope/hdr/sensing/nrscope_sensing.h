/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
/*
 * \brief The sensing pipeline behind one context, fed per decoded downlink grant.
 *
 * Everything in the sensing directory is the algorithm moved from
 * OpenAirInterface5G (tag isac-reference-for-nrscope-port). This file is the glue
 * that feeds it from NR-Scope's receiver, and the only sensing interface the rest
 * of NR-Scope sees.
 *
 * Per grant, on the worker that decoded it:
 *   1. the grid is the slot's full-carrier resource grid the DCI decoder already
 *      demodulates for the DM-RS check (dmrs_check.h), so nothing is transformed
 *      twice;
 *   2. each layer's channel is estimated on the DM-RS symbols by the despreader
 *      (nr_ue_dmrs_despread.h), on the pilots the DM-RS check verifies on air;
 *   3. nr_ue_sensing_slot_profile() turns each DM-RS symbol into a delay response
 *      and pushes it, stamped with its absolute sample time, into the history the
 *      range-Doppler maps are built from.
 *
 *   4. when a stream's newest snapshot is one observation window past the end of its
 *      last map (nr_ue_sensing_window_s(), the window the transform gathers over),
 *      its history is copied and handed to a map thread, which runs OAI's map task
 *      (nr_ue_sensing_map_task in OAI's phy_procedures_nr_ue.c): range-Doppler map,
 *      TDD detector, AoA, MUSIC, and the map dump, with OAI's prints.
 *
 * Scratch buffers belong to the caller (one per DCI decoder) so workers estimate in
 * parallel; only the push into the shared history takes its lock.
 */
#ifndef NRSCOPE_SENSING_H
#define NRSCOPE_SENSING_H

#include <stdbool.h>
#include <stdint.h>

#include "nrscope/hdr/sensing/nr_ue_dmrs_despread.h"
#include "srsran/phy/phch/phch_cfg_nr.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Runtime configuration, the sensing: block of the yaml. One field per sensing
   option of OAI's nr-uesoftmodem (--sensing-*), same meaning and same defaults. */
typedef struct {
  bool enable;
  /* --sensing-symbols: the most reference symbols one map may use, which caps the
  cost of a transform (O(n_freq * n_bins * symbols)) and nothing else. The map is
  triggered by time, not by this count, so a value under what the window holds only
  thins a map; NR_SENSING_MIN_SNAPSHOTS..NR_SENSING_HISTORY_DEPTH. */
  int symbols;
  /// --sensing-max-speed: largest speed the map window is sized for, m/s
  double max_speed_ms;
  /// --sensing-dump: map CSV (plus <dump>.snap.csv once, .music.csv, .L1.csv); empty for none
  char dump[256];
  /// --sensing-clutter-removal: remove the static scene (slow-time mean)
  bool clutter_removal;
  /// --sensing-clutter-kernel: remove it with each snapshot's point spread function instead
  bool clutter_kernel;
  /// --sensing-antenna-avg: one map averaged over the Rx antennas, and the AoA
  bool antenna_avg;
  /// --sensing-layer-avg: one map averaged over the layers
  bool layer_avg;
  /// --sensing-random-drop: break up the periodicity of the slow-time sampling
  bool random_drop;
  /// --sensing-tdd-detect: run the TDD detector on each map
  bool tdd_detect;
  /// --sensing-music: also build each map with Doppler MUSIC
  bool music;
  /* NR-Scope only, on: set aside as uncertain (verdict 3, not drawn, not localised) a
  detection whose mirror at -v on the same range is within NR_TDD_MIRROR_DB of it. Off
  keeps every detection the TDD detector accepts. */
  bool mirror_reject;
  /* test_record_samples: record-only mode. Write one file per window holding every
  (chain, layer) of its slow-time samples to <dump>.rec.<NNNNN>.csv, up to
  record_max_files, and skip the whole live map pipeline (range_doppler, detector, AoA,
  map dump). A Python reprocessor (scripts/sensing/replay_samples.py)
  rebuilds the maps offline and sweeps any parameter without a rebuild or recapture. */
  bool   test_record_samples;
  int    record_max_files;
  /* avg_maps: sliding non-coherent average of the last avg_maps maps of the same stream
  (1 = off, at most NR_SENSING_AVG_MAX). Each map is still built from its own window;
  the dumped map, and the AoA's own CFAR when tdd_detect is off, see the average. The
  TDD detector works on the window's slow-time samples and is not affected. Maps of a
  stream are all built on the first one's Doppler grid while this is on.
  avg_vcomp: shift every Doppler column of an older map by its path-length rate
  (lambda * f * dt) before averaging, so a mover stays in its cell. */
  int    avg_maps;
  bool   avg_vcomp;
  /* map_max_gap_ms: a stream whose snapshots jump by more than this lost samples upstream
  (capture outage, overflow). The window that was filling is discarded and the next one
  starts after the gap, so a map is only built from a window without a hole and short
  maps from a few snapshots cannot happen. 0 disables it. */
  double map_max_gap_ms;
} nrscope_sensing_args_t;

/// Most maps avg_maps may average over
#define NR_SENSING_AVG_MAX 16

void nrscope_sensing_default_args(nrscope_sensing_args_t* args);

/* The args the yaml gave, set once by load_config before the decoders start. */
extern nrscope_sensing_args_t nrscope_sensing_args;

typedef struct nrscope_sensing_s nrscope_sensing_t;

/* Per-caller scratch: one slot of channel estimates and their valid mask, over
   the whole carrier. About 400 kB, so one per DCI decoder, allocated once. */
typedef struct nrscope_sensing_scratch_s nrscope_sensing_scratch_t;
nrscope_sensing_scratch_t* nrscope_sensing_scratch_alloc(uint32_t n_sc_grid);
void                       nrscope_sensing_scratch_free(nrscope_sensing_scratch_t* sc);

/* The shared context, created on first use by whichever decoder gets there first
   and then returned to every caller. NULL when sensing is disabled.
   carrier_hz  : the downlink carrier centre, which sets the Doppler-to-speed scale
   srate_hz    : sample rate of the capture, which sets the sample clock of t_sample
   ofdm_size   : FFT size of the grid (4096 at 122.88 Msps and 30 kHz)
   nof_antennas: receive chains captured, 1 .. NR_AOA_MAX_ANT; each gets its own rings */
nrscope_sensing_t* nrscope_sensing_get(uint64_t carrier_hz, double srate_hz, uint32_t ofdm_size, uint32_t scs_hz,
                                       uint32_t nof_antennas);

/* One decoded downlink grant, as received on one chain.
   aarx        : the receive chain grid came from. Call once per chain for a slot,
                 in chain order, chain 0 first: a map averaging the chains is
                 triggered by the last one, and the alignment gives the other chains
                 the correction chain 0's symbols got, which keeps the phase between
                 the chains that the AoA reads
   grid        : the slot's full-carrier grid on that chain, n_sc_grid subcarriers
                 per symbol, symbol-major, starting at CRB place->grid_crb0
   cfg         : the PDSCH configuration derived from the DCI (grant + DM-RS)
   dci_ports   : the DCI's antenna-ports field (TS 38.212 7.3.1.2.2), or -1 for a
                 DCI 1_0, which has none and means port 1000
   pci         : physical cell identity, the DM-RS scrambling fallback
   sfn, slot_idx : where the slot sits in the frame
   Returns snapshots pushed, 0 when the grant is not used (too narrow, a DM-RS
   configuration the estimator does not model), negative on error. */
int nrscope_sensing_process_grant(nrscope_sensing_t*          s,
                                  nrscope_sensing_scratch_t*  sc,
                                  uint32_t                    aarx,
                                  const cf_t*                 grid,
                                  uint32_t                    n_sc_grid,
                                  const nr_dmrs_placement_t*  place,
                                  const srsran_sch_cfg_nr_t*  cfg,
                                  int                         dci_ports,
                                  uint32_t                    pci,
                                  uint32_t                    sfn,
                                  uint32_t                    slot_idx);

/* Block until every queued map is built and dumped. For tests. */
void nrscope_sensing_wait_maps(nrscope_sensing_t* s);

/* DM-RS ports behind the DCI's antenna-ports value, TS 38.212 table 7.3.1.2.2-1
   (configuration type 1, maxLength 1, one codeword): a bitmap, bit i meaning
   antenna port 1000 + i. Sets *cdm_groups to the CDM groups without data the same
   row gives. Returns 0 for a reserved value. */
uint16_t nrscope_dmrs_ports_type1_len1(uint32_t value, uint32_t* cdm_groups);

#ifdef __cplusplus
}
#endif

#endif
