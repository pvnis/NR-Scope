/*
 * GoldenSniffer in NR-Scope; see nrscope/hdr/golden_sniffer.h.
 */
#include "nrscope/hdr/golden_sniffer.h"
#include <cmath>
#include <cstdio>
#include <map>
#include <mutex>

GoldenSnifferArgs nrscope_gs_args;

namespace {

struct RntiStats {
  uint64_t                     dcis     = 0;
  uint32_t                     n_id     = 0;
  bool                         n_rnti   = false;
  bool                         coreset0 = false;
  bool                         fed      = false;
  std::map<uint32_t, uint64_t> sizes;
  uint64_t                     al[5] = {};
};

std::mutex                    gs_log_mutex;
FILE*                         gs_csv = nullptr;
std::map<uint16_t, RntiStats> gs_stats;
double                        gs_last_report = 0;
uint64_t                      gs_total       = 0;

const char* rnti_kind(uint16_t rnti, bool coreset0, bool n_rnti)
{
  if (rnti == 0xffff) {
    return "SI";
  }
  if (rnti == 0xfffe) {
    return "P";
  }
  /* RA-RNTI = 1 + s_id + 14 t_id + 14 80 f_id + 14 80 8 ul_carrier_id <= 0x4600 (TS 38.321 5.1.3) */
  if (coreset0 && !n_rnti && rnti >= 1 && rnti <= 0x4600) {
    return "RA/TC";
  }
  return "C";
}

void print_report_locked()
{
  printf("GoldenSniffer: %llu DCIs, %zu RNTIs\n", (unsigned long long)gs_total, gs_stats.size());
  for (const auto& e : gs_stats) {
    const RntiStats& s = e.second;
    printf("  RNTI 0x%04x %-5s N_ID %5u %s%s DCIs %8llu  sizes", e.first, rnti_kind(e.first, s.coreset0, s.n_rnti),
           s.n_id, s.n_rnti ? "n_RNTI " : "", s.coreset0 ? "CORESET#0" : "         ", (unsigned long long)s.dcis);
    for (const auto& z : s.sizes) {
      printf(" %u:%llu", z.first, (unsigned long long)z.second);
    }
    printf("  AL");
    for (int l = 0; l < 5; l++) {
      if (s.al[l]) {
        printf(" %d:%llu", 1 << l, (unsigned long long)s.al[l]);
      }
    }
    printf("%s\n", s.fed ? "  (fed to the DCI decoders)" : "");
  }
  fflush(stdout);
}

} // namespace

void load_golden_sniffer_config(const YAML::Node& config_yaml)
{
  GoldenSnifferArgs& a = nrscope_gs_args;
  a.csv                = std::string(NRSCOPE_ROOT_DIR) + "/golden_sniffer.csv";
  if (!config_yaml["golden_sniffer"]) {
    return;
  }
  const YAML::Node g = config_yaml["golden_sniffer"];
  if (g["enable"]) {
    a.enable = g["enable"].as<bool>();
  }
  if (g["max_duration"]) {
    a.max_duration = g["max_duration"].as<uint32_t>();
    if (a.max_duration < 1 || a.max_duration > GS_MAX_DURATION) {
      std::cerr << "golden_sniffer: max_duration must be 1.." << GS_MAX_DURATION << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  if (g["dci_sizes"]) {
    a.dci_sizes = g["dci_sizes"].as<std::vector<uint32_t> >();
  }
  if (g["css_sizes"]) {
    if (g["css_sizes"].IsScalar() && g["css_sizes"].as<std::string>() == "auto") {
      a.css_auto = true;
    } else {
      a.css_auto  = false;
      a.css_sizes = g["css_sizes"].as<std::vector<uint32_t> >();
    }
  }
  if (g["max_al"]) {
    a.max_al = g["max_al"].as<uint32_t>();
  }
  if (g["feed_rntis"]) {
    a.feed_rntis = g["feed_rntis"].as<bool>();
  }
  if (g["min_parity"]) {
    a.min_parity = g["min_parity"].as<uint32_t>();
  }
  if (g["min_hits"]) {
    a.min_hits = g["min_hits"].as<uint32_t>();
  }
  if (g["csv"]) {
    a.csv = g["csv"].as<std::string>();
  }
  if (g["report_s"]) {
    a.report_s = g["report_s"].as<double>();
  }
  if (a.enable) {
    gs_csv = fopen(a.csv.c_str(), "w");
    if (gs_csv == nullptr) {
      std::cerr << "golden_sniffer: cannot write " << a.csv << std::endl;
    } else {
      fprintf(gs_csv, "timestamp,sfn,slot,rnti,kind,n_id,n_rnti,coreset0,al,duration,first_crb,interleaved,size,"
                      "payload_hex,evm,epre_db\n");
    }
    std::cout << "GoldenSniffer enabled, DCIs to " << a.csv << std::endl;
  }
}

void golden_sniffer_report()
{
  if (!nrscope_gs_args.enable) {
    return;
  }
  std::lock_guard<std::mutex> lock(gs_log_mutex);
  print_report_locked();
  if (gs_csv != nullptr) {
    fflush(gs_csv);
  }
}

void golden_sniffer_update(WorkState& state, const SlotResult& result, double now)
{
  if (!nrscope_gs_args.enable) {
    return;
  }
  std::lock_guard<std::mutex> lock(gs_log_mutex);
  for (const GsDciResult& d : result.gs_dcis) {
    gs_total++;
    RntiStats& s = gs_stats[d.rnti];
    s.dcis++;
    s.n_id     = d.n_id;
    s.n_rnti   = d.n_rnti;
    s.coreset0 = d.coreset0;
    s.sizes[d.nof_bits]++;
    for (int l = 0; l < 5; l++) {
      if (d.L == (1u << l)) {
        s.al[l]++;
      }
    }
    if (gs_csv != nullptr) {
      char     hex[GS_MAX_DCI_BITS / 4 + 2] = {};
      uint32_t nh                           = (d.nof_bits + 3) / 4;
      for (uint32_t h = 0; h < nh; h++) {
        uint32_t v = 0;
        for (uint32_t b = 0; b < 4; b++) {
          uint32_t i = 4 * h + b;
          v          = (v << 1) | (i < d.nof_bits ? d.payload[i] : 0);
        }
        hex[h] = "0123456789abcdef"[v];
      }
      fprintf(gs_csv, "%.6f,%u,%u,0x%04x,%s,%u,%d,%d,%u,%u,%u,%d,%u,%s,%.3f,%.1f\n", now, result.outcome.sfn,
              d.slot_idx, d.rnti, rnti_kind(d.rnti, d.coreset0, d.n_rnti), d.n_id, d.n_rnti, d.coreset0, d.L,
              d.duration, d.first_crb, d.nof_iv > 1, d.nof_bits, hex, d.evm,
              d.epre > 0 ? 10.0 * std::log10(d.epre) : -999.0);
    }

    /* C-RNTIs to the regular DCI decoders, with the scrambling ID read off the DM-RS */
    if (nrscope_gs_args.feed_rntis && !s.fed && s.dcis >= nrscope_gs_args.min_hits && state.rach_found &&
        std::string(rnti_kind(d.rnti, d.coreset0, d.n_rnti)) == "C") {
      bool known = false;
      for (uint16_t r : state.known_rntis) {
        known |= r == d.rnti;
      }
      if (!known) {
        state.known_rntis.push_back(d.rnti);
        state.last_seen.push_back(now);
        state.nof_known_rntis++;
        const auto& sp_ded = state.master_cell_group.sp_cell_cfg.sp_cell_cfg_ded;
        if (d.n_rnti && sp_ded.init_dl_bwp_present && sp_ded.init_dl_bwp.pdcch_cfg_present &&
            sp_ded.init_dl_bwp.pdcch_cfg.type().value == asn1::setup_release_opts::setup) {
          for (const auto& cs : sp_ded.init_dl_bwp.pdcch_cfg.setup().ctrl_res_set_to_add_mod_list) {
            if (cs.pdcch_dmrs_scrambling_id_present) {
              state.pdcch_dmrs_ids_by_rnti[d.rnti][cs.ctrl_res_set_id] = d.n_id;
            }
          }
        }
        printf("GoldenSniffer: RNTI 0x%04x (N_ID %u) added to the DCI decoders\n", d.rnti, d.n_id);
      }
      s.fed = true;
    }
  }
  if (nrscope_gs_args.report_s > 0 && now - gs_last_report >= nrscope_gs_args.report_s) {
    if (gs_last_report > 0) {
      print_report_locked();
    }
    gs_last_report = now;
    if (gs_csv != nullptr) {
      fflush(gs_csv);
    }
  }
}

GoldenSnifferDecoder::~GoldenSnifferDecoder()
{
  if (gs != nullptr) {
    gs_destroy(gs);
  }
  if (ready) {
    srsran_ue_dl_nr_free(&ue_dl);
  }
}

int GoldenSnifferDecoder::Init(WorkState* state, cf_t* input0)
{
  const auto&  cell    = state->cell;
  const auto&  sib1    = state->sib1;
  const auto&  freq_dl = sib1.serving_cell_cfg_common.dl_cfg_common.freq_info_dl;
  const auto&  carrier = freq_dl.scs_specific_carrier_list[0];
  const double scs     = cell.abs_pdcch_scs;

  /* Point A, as DCIDecoder places it (offsetToPointA in 15 kHz RBs in FR1) */
  const double pointA = state->srsran_searcher_cfg_t.ssb_freq_hz - (SRSRAN_SSB_BW_SUBC / 2) * cell.abs_ssb_scs -
                        cell.k_ssb * SRSRAN_SUBC_SPACING_NR(srsran_subcarrier_spacing_15kHz) -
                        freq_dl.offset_to_point_a * SRSRAN_SUBC_SPACING_NR(srsran_subcarrier_spacing_15kHz) *
                            NRSCOPE_NSC_PER_RB_NR;
  const double f_lo  = pointA + carrier.offset_to_carrier * NRSCOPE_NSC_PER_RB_NR * scs; // carrier PRB 0
  const double fc    = state->args_t.base_carrier.dl_center_frequency_hz;
  const double srate = state->arg_scs.srate;
  /* Up to the FFT's edge, not just the flat part: a CORESET at the carrier edge (CRB 0)
    is common, and a run clipped by a few RBs is lost whole. The roll-off costs SNR only. */
  const double half  = srate / 2.0 - NRSCOPE_NSC_PER_RB_NR * scs;
  const double rb_hz = NRSCOPE_NSC_PER_RB_NR * scs;

  /* The carrier, cut to what the radio passes flat */
  int first = (int)std::ceil((fc - half - f_lo) / rb_hz);
  int end   = (int)std::floor((fc + half - f_lo) / rb_hz);
  first     = std::max(first, 0);
  end       = std::min(end, (int)carrier.carrier_bw);
  if (end - first < 6) {
    ERROR("GoldenSniffer: the capture holds no part of the carrier");
    return SRSRAN_ERROR;
  }
  if (first > 0 || end < (int)carrier.carrier_bw) {
    printf("GoldenSniffer: the capture holds carrier PRBs %d..%d of %u; PDCCHs outside are not searched\n", first,
           end - 1, (unsigned)carrier.carrier_bw);
  }
  first_crb            = carrier.offset_to_carrier + (uint32_t)first;
  grid_carrier         = state->args_t.base_carrier;
  grid_carrier.nof_prb = (uint32_t)(end - first);
  const double grid_lo = f_lo + first * rb_hz;
  const double center  = grid_lo + grid_carrier.nof_prb / 2.0 * rb_hz;
  arg_scs_grid         = state->arg_scs;
  arg_scs_grid.coreset_offset_scs = (int)std::lround((fc - center) / scs);
  nrscope_check_grid_in_capture("GoldenSniffer grid", arg_scs_grid.coreset_offset_scs, grid_carrier.nof_prb, srate,
                                (uint32_t)scs);

  srsran_ue_dl_nr_args_t args = {};
  args.nof_rx_antennas        = 1;
  args.nof_max_prb            = 275;
  args.pdsch.sch.disable_simd = false;
  args.pdcch.disable_simd     = false;
  args.pdcch_dmrs_corr_thr    = 0.05;
  cf_t* input[SRSRAN_MAX_PORTS] = {input0};
  SRSRAN_MEM_ZERO(&ue_dl, srsran_ue_dl_nr_t, 1); // no CORESET: only demodulates
  if (srsran_ue_dl_nr_init_nrscope(&ue_dl, input, &args, arg_scs_grid) < SRSRAN_SUCCESS ||
      srsran_ue_dl_nr_set_carrier_nrscope(&ue_dl, &grid_carrier, arg_scs_grid) < SRSRAN_SUCCESS) {
    ERROR("GoldenSniffer: could not set up the carrier grid");
    return SRSRAN_ERROR;
  }

  ref_crb0 = -(int32_t)(first_crb * NRSCOPE_NSC_PER_RB_NR);
  ref_cs0  = (int32_t)std::lround((state->coreset0_args_t.coreset0_lower_freq_hz - grid_lo) / scs);
  pci      = state->cs_ret.ssb_res.N_id;
  slots_per_frame = SRSRAN_NSLOTS_PER_FRAME_NR(state->arg_scs.scs);

  gs_args_t ga;
  gs_args_default(&ga);
  ga.max_duration = nrscope_gs_args.max_duration;
  ga.min_parity   = nrscope_gs_args.min_parity;
  ga.max_L_log2   = 0;
  while ((2u << ga.max_L_log2) <= nrscope_gs_args.max_al && ga.max_L_log2 < GS_MAX_L_LOG2) {
    ga.max_L_log2++;
  }
  ga.sizes     = nrscope_gs_args.dci_sizes.data();
  ga.nof_sizes = (uint32_t)nrscope_gs_args.dci_sizes.size();
  std::vector<uint32_t> css = nrscope_gs_args.css_sizes;
  if (nrscope_gs_args.css_auto) {
    /* DCI 0_0 / 1_0: 28 bits + the frequency assignment (TS 38.212 7.3.1.2.1; 0_0 is
      padded or truncated to it), over the CORESET#0 bandwidth in the common search
      spaces, and over the initial BWP in a UE search space of a UE scrambled with the
      PCI. Both, since the DM-RS cannot tell the CORESETs apart when CORESET#0 starts at
      point A. */
    auto fallback = [](uint32_t n) { return 28u + (uint32_t)std::ceil(std::log2(n * (n + 1) / 2.0)); };
    css           = {fallback(srsran_coreset_get_bw(&state->coreset0_t))};
    const uint32_t riv = sib1.serving_cell_cfg_common.dl_cfg_common.init_dl_bwp.generic_params.location_and_bw;
    const uint32_t N = 275, a = riv / N, b = riv % N;
    const uint32_t bwp = (a + b < N) ? a + 1 : N - a + 1;
    if (fallback(bwp) != css[0]) {
      css.push_back(fallback(bwp));
    }
  }
  ga.css_sizes     = css.data();
  ga.nof_css_sizes = (uint32_t)css.size();
  gs               = gs_create(&ga);
  if (gs == nullptr) {
    ERROR("GoldenSniffer: could not create the searcher");
    srsran_ue_dl_nr_free(&ue_dl);
    return SRSRAN_ERROR;
  }
  ready = true;
  printf("GoldenSniffer ready: grid of %u PRBs from CRB %u, point A at grid subcarrier %d, CORESET#0 at %d, PCI %u, "
         "CSS sizes",
         grid_carrier.nof_prb, first_crb, ref_crb0, ref_cs0, pci);
  for (uint32_t c : css) {
    printf(" %u", c);
  }
  printf("\n");
  return SRSRAN_SUCCESS;
}

void GoldenSnifferDecoder::DecodeSlot(srsran_slot_cfg_t* slot, WorkState* state, SlotResult* result)
{
  if (!ready) {
    return;
  }
  /* one symbol more than searched: the duration detector compares it */
  uint32_t mask[SRSRAN_MAX_PORTS] = {(1u << (nrscope_gs_args.max_duration + 1)) - 1u};
  uint32_t done[SRSRAN_MAX_PORTS] = {};
  srsran_ue_dl_nr_estimate_fft_nrscope_symbols(&ue_dl, slot, arg_scs_grid, mask, done);

  gs_slot_t s = {};
  s.grid      = ue_dl.sf_symbols[0];
  s.nof_sc    = grid_carrier.nof_prb * NRSCOPE_NSC_PER_RB_NR;
  s.slot_idx  = slot->idx % slots_per_frame;
  s.ref_crb0  = ref_crb0;
  s.ref_cs0   = ref_cs0;
  s.pci       = pci;

  /* Debug: GS_DUMP=<file> appends the searched symbols of every slot, as
    [sfn u32, slot u32, nof_sc u32, nof_symb u32, cf32 grid], for offline analysis */
  static FILE*      dump = nullptr;
  static std::mutex dump_mutex;
  if (getenv("GS_DUMP") != nullptr) {
    std::lock_guard<std::mutex> lock(dump_mutex);
    if (dump == nullptr) {
      dump = fopen(getenv("GS_DUMP"), "wb");
    }
    if (dump != nullptr) {
      uint32_t hdr[4] = {state->sfn, s.slot_idx, s.nof_sc, nrscope_gs_args.max_duration + 1};
      fwrite(hdr, sizeof(hdr), 1, dump);
      fwrite(s.grid, sizeof(cf_t), (size_t)s.nof_sc * hdr[3], dump);
    }
  }

  gs_dci_t out[GS_MAX_DCI_PER_SLOT];
  int      n = gs_search_slot(gs, &s, out, GS_MAX_DCI_PER_SLOT);
  for (int i = 0; i < n; i++) {
    GsDciResult r;
    r.slot_idx  = s.slot_idx;
    r.rnti      = out[i].rnti;
    r.n_id      = out[i].n_id;
    r.n_rnti    = out[i].n_rnti;
    r.coreset0  = out[i].coreset0;
    r.L         = out[i].L;
    r.duration  = out[i].duration;
    r.first_crb = first_crb + (uint32_t)out[i].cand.iv[0].first / NRSCOPE_NSC_PER_RB_NR;
    r.nof_iv    = out[i].cand.nof_iv;
    r.nof_bits  = out[i].nof_bits;
    r.payload.assign(out[i].payload, out[i].payload + out[i].nof_bits);
    r.evm  = out[i].evm;
    r.epre = out[i].epre;
    result->gs_dcis.push_back(std::move(r));
  }
}
