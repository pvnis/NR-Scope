#ifndef NRSCOPE_RUN_RECORDER_H
#define NRSCOPE_RUN_RECORDER_H

#include <stdint.h>
#include <string>

#include "srsran/phy/phch/dci_nr.h"
#include "srsran/phy/phch/phch_cfg_nr.h"

#include "nrscope/hdr/dmrs_check.h"

/* Per-run CSV recordings for offline statistics, switched on by
  log_config.recording_mode in the yaml. Each run writes, under the project
  root rather than the build tree:

    msg4/msg4_<date>_<time>_pci<N>.csv  one row per RRCSetup decoded
    DCIs/dci_<date>_<time>_pci<N>.csv   one row per DCI found for a known UE

  <date>_<time> is when the run started, so the two files of a run pair up. A
  file is only created with its first row: a run that never sees an RRCSetup,
  or never a UE's DCI, leaves no file of that kind.

  While recording, the terminal keeps only the "Found DCI" and "rrc_setup,
  hooray" lines; the grant dumps and the JSON go to the files instead. */
namespace RunRecorder {

void init(bool enable);
bool enabled();

/* One-shot snapshot of the cell configuration recovered from the air, used to
  draw the static cell map (SSB / CORESET#0 / SIB1 / initial DL BWP against
  Point A). Every field comes from what NR-Scope itself decoded (MIB, SIB1,
  CORESET#0 zero tables) or computed from it, never from a gNB-side log. All
  frequencies are in Hz; RB counts use the resource block of the subcarrier
  spacing named alongside them, except offset_to_point_a_rb which, per TS 38.211,
  counts RBs of 15 kHz. */
struct CellSummary {
  uint32_t pci                       = 0;
  double   dl_center_freq_hz         = 0.0; // capture centre
  double   ssb_center_freq_hz        = 0.0; // SS/PBCH block centre
  double   coreset0_lower_freq_hz    = 0.0; // lowest subcarrier of CORESET#0
  double   coreset0_center_freq_hz   = 0.0;
  uint32_t ssb_scs_khz               = 0;
  uint32_t common_scs_khz            = 0;   // scs of the initial DL BWP / carrier
  uint32_t k_ssb                     = 0;   // subcarrier offset SSB -> CRB grid
  uint32_t ssb_idx                   = 0;
  uint32_t coreset0_idx              = 0;   // pdcch-ConfigSIB1 4 MSBs
  uint32_t ss0_idx                   = 0;   // pdcch-ConfigSIB1 4 LSBs
  uint32_t coreset0_offset_rb        = 0;   // from Point A, in common-scs RBs
  uint32_t coreset0_bw_rb            = 0;
  uint32_t coreset0_duration_symbols = 0;
  uint32_t coreset0_first_symbol     = 0;   // first OFDM symbol in its slot (table 13-11)
  uint32_t coreset0_slot_n0          = 0;   // monitoring slot n_0 within a frame
  uint32_t coreset0_sfn_c            = 0;   // 0 = even SFN, 1 = odd SFN
  char     ssb_pattern[8]            = {0}; // "A".."E", sets the SSB symbol positions
  uint32_t offset_to_point_a_rb      = 0;   // SIB1, in 15 kHz RBs
  uint32_t carrier_bw_rb             = 0;   // SIB1 scs-SpecificCarrier, common-scs RBs
  uint32_t carrier_offset_to_carrier = 0;   // SIB1 scs-SpecificCarrier
  /* Where SIB1's own PDSCH landed, from the SI-RNTI DCI in CORESET#0. RBs are
    on the CORESET#0 grid (RB 0 = CORESET#0's lowest RB); 0 if not captured. */
  uint32_t sib1_prb_start           = 0;
  uint32_t sib1_nof_prb             = 0;
  uint32_t sib1_symbol_start        = 0;
  uint32_t sib1_nof_symbols         = 0;
  uint32_t sib1_slot_idx            = 0;
  /* SIB1's PDCCH within CORESET#0, and the CORESET#0 CCE-to-REG mapping needed
    to place those CCEs on the grid (TS 38.211 7.3.2.2). */
  uint32_t sib1_pdcch_agg_level     = 0;    // number of CCEs (1<<L)
  uint32_t sib1_pdcch_ncce          = 0;    // first CCE index
  uint32_t coreset0_reg_bundle_size = 0;    // L, in REGs (2/3/6)
  uint32_t coreset0_interleaver_size = 0;   // R (2/3/6); 0 if non-interleaved
  uint32_t coreset0_shift_index     = 0;
  uint32_t coreset0_interleaved     = 0;    // 1 if interleaved CCE-to-REG mapping
};

/* Written once, the first time SIB1 is decoded for a cell, to
  cell/cell_<date>_<time>_pci<N>.json. No-op unless recording is enabled. */
void record_cell_summary(const CellSummary& c);

/* Flush and close the files. Rows are also flushed periodically, since runs
  usually end by being killed. */
void close();

void record_rrc_setup(uint32_t           pci,
                      uint32_t           sfn,
                      uint32_t           slot_idx,
                      uint16_t           tc_rnti,
                      uint16_t           c_rnti,
                      uint32_t           rrc_transaction_id,
                      const char*        dci_str,
                      uint32_t           rrc_offset,
                      const uint8_t*     msg4_bytes,
                      uint32_t           nof_bytes,
                      const std::string& master_cell_group_json);

/* sch_cfg is the PDSCH/PUSCH configuration derived from the DCI, or NULL when
  the decoder did not derive one (downlink formats other than 1_1). ca_variant
  says the DCI matched the carrier aggregation DCI sizes rather than the normal
  ones, and dci_bits is its raw payload as '0'/'1', so other field layouts can
  be tested offline. dci_str is the line printed as "Found DCI", kept whole. */
void record_dci(uint32_t                   pci,
                uint32_t                   sfn,
                uint32_t                   slot_idx,
                bool                       downlink,
                const srsran_dci_ctx_t&    ctx,
                uint32_t                   freq_alloc,
                uint32_t                   time_alloc,
                uint32_t                   mcs,
                uint32_t                   ndi,
                uint32_t                   rv,
                uint32_t                   harq_id,
                uint32_t                   tpc,
                uint32_t                   ports,
                uint32_t                   dmrs_id,
                uint32_t                   srs_request,
                const srsran_sch_cfg_nr_t* sch_cfg,
                bool                       ca_variant,
                const char*                dci_bits,
                const char*                dci_str,
                const DmrsCheckResult*     dmrs_check = nullptr);

/* Debug: every PDCCH candidate the DCI decoder evaluated for a known RNTI, one
  row per candidate and per DCI size set tried, to PDCCH/pdcch_<run>_pci<N>.csv.
  Switched on by log_config.record_pdcch_candidates, and only with
  recording_mode. At 2000 slots/s this is tens of thousands of rows a second, so
  it is meant for short diagnostic runs.

  stage says how far the candidate got: "no_measure" (the DM-RS measure was not
  a number), "epre" or "corr" (dropped by that DM-RS threshold before
  decoding), "crc_fail" or "crc_ok" (decoded; CRC against the RNTI). The PDCCH
  DM-RS does not depend on the RNTI, so epre / norm_corr say whether a PDCCH was
  physically where NR-Scope looked; the CRC says whether it was this UE's.
  al is the aggregation level in CCEs. */
void enable_pdcch_candidates(bool enable);
bool pdcch_candidates_enabled();
void record_pdcch_candidate(uint32_t    pci,
                            uint32_t    sfn,
                            uint32_t    slot_idx,
                            uint16_t    rnti,
                            bool        ca_variant,
                            const char* ss_type,
                            uint32_t    coreset_id,
                            uint32_t    al,
                            uint32_t    cce,
                            uint32_t    nof_bits,
                            float       epre_dBfs,
                            float       rsrp_dBfs,
                            float       norm_corr,
                            float       sync_error_us,
                            const char* stage);

} // namespace RunRecorder

#endif
