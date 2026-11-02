/*
 * GoldenSniffer in NR-Scope: blind DCI detection without known RNTIs.
 *
 * The C core (golden_sniffer/gs_core.h) finds DCIs in a slot's first PDCCH symbols
 * and reads their RNTI, PDCCH DM-RS scrambling ID, size and payload. This file runs
 * it in every worker from the moment SIB1 is known, on its own full-carrier grid
 * (the part of the carrier inside the capture), without waiting for a RACH, and
 * reports what it finds:
 *   - every DCI to a CSV (golden_sniffer.csv);
 *   - a per-RNTI count, on stdout every report_s seconds and at exit;
 *   - with feed_rntis, the C-RNTIs it finds (seen min_hits times) to the regular
 *     DCI decoders, with the scrambling ID it read, once an RRCSetup gave them
 *     CORESET 1.
 *
 * Config, top-level section:
 *   golden_sniffer:
 *     enable: true
 *     max_duration: 2        # PDCCH symbols scanned
 *     dci_sizes: []          # sizes tried (empty: 20..84)
 *     css_sizes: auto        # sizes allowed for n_RNTI = 0 DCIs whose DM-RS may be CORESET#0's:
 *                            # auto (DCI 1_0 over CORESET#0 and over the initial BWP), a list,
 *                            # or [] for any
 *     max_al: 16
 *     min_parity: 16         # observed bits beyond the unknowns a decode needs
 *     feed_rntis: false
 *     min_hits: 3
 *     csv: <path>            # default <root>/golden_sniffer.csv
 *     report_s: 10
 */
#ifndef NRSCOPE_GOLDEN_SNIFFER_H
#define NRSCOPE_GOLDEN_SNIFFER_H

#include "nrscope/hdr/golden_sniffer/gs_core.h"
#include "nrscope/hdr/nrscope_def.h"
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct GoldenSnifferArgs {
  bool                  enable       = false;
  uint32_t              max_duration = 2;
  std::vector<uint32_t> dci_sizes;
  bool                  css_auto = true;
  std::vector<uint32_t> css_sizes;
  uint32_t              max_al     = 16;
  uint32_t              min_parity = 16;
  bool                  feed_rntis = false;
  uint32_t              min_hits   = 3;
  std::string           csv;
  double                report_s = 10.0;
};

extern GoldenSnifferArgs nrscope_gs_args;

void load_golden_sniffer_config(const YAML::Node& config_yaml);

/* Per-RNTI summary on stdout and CSV flush; called by the exit path */
void golden_sniffer_report();

/* Under the scheduler lock: logs the slot's DCIs, and with feed_rntis adds C-RNTIs */
void golden_sniffer_update(WorkState& state, const SlotResult& result, double now);

class GoldenSnifferDecoder
{
public:
  GoldenSnifferDecoder() = default;
  ~GoldenSnifferDecoder();

  /* After SIB1: places the grid, builds the searcher */
  int Init(WorkState* state, cf_t* input0);

  /* One slot (a worker thread entry) */
  void DecodeSlot(srsran_slot_cfg_t* slot, WorkState* state, SlotResult* result);

  bool ready = false;

private:
  srsran_ue_dl_nr_t             ue_dl = {};
  srsran_carrier_nr_t           grid_carrier = {};
  srsran_ue_dl_nr_sratescs_info arg_scs_grid = {};
  gs_t*                         gs       = nullptr;
  int32_t                       ref_crb0 = GS_NO_REF;
  int32_t                       ref_cs0  = GS_NO_REF;
  uint32_t                      pci      = 0;
  uint32_t                      slots_per_frame = 20;
  uint32_t                      first_crb = 0; /* CRB of grid subcarrier 0 */
};

#endif
