/*
 * GoldenSniffer in NR-Scope; see nrscope/hdr/golden_sniffer.h.
 */
#include "nrscope/hdr/golden_sniffer.h"
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <deque>
#include <map>
#include <mutex>
#include <thread>

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

/* What the workers hand over: the searched symbols of one slot, copied out of the
  worker's grid so the worker can go on with the next slot at once. */
struct GsJob {
  uint32_t          sfn = 0, slot_idx = 0, nof_sc = 0, nof_symb = 0, first_crb = 0, pci = 0;
  int32_t           ref_crb0 = GS_NO_REF, ref_cs0 = GS_NO_REF;
  uint32_t          cs0_nof_sc = 0;
  std::vector<cf_t> grid;
};

/* The searcher pool: its own threads and its own bounded queue. GoldenSniffer is never
  on the worker's critical path; when the pool falls behind, it drops its own slots. */
struct GsPool {
  std::once_flag          started;
  std::mutex              m;
  std::condition_variable cv;
  std::deque<GsJob>       q;
  gs_args_t               args = {};
  std::vector<uint32_t>   sizes, css;
  /* statistics, under m */
  uint64_t submitted = 0, searched = 0, dropped = 0;
  double   t_sum_ms = 0, t_max_ms = 0;
};
GsPool gs_pool;

/* C-RNTIs ready to go to the DCI decoders (feed_rntis), under gs_log_mutex */
struct GsFeed {
  uint16_t rnti;
  uint32_t n_id;
  bool     n_rnti;
};
std::vector<GsFeed> gs_feeds;

double wall_now()
{
  return std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
}

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
  uint64_t searched, dropped, submitted;
  double   t_sum, t_max;
  {
    std::lock_guard<std::mutex> lock(gs_pool.m);
    searched  = gs_pool.searched;
    dropped   = gs_pool.dropped;
    submitted = gs_pool.submitted;
    t_sum     = gs_pool.t_sum_ms;
    t_max     = gs_pool.t_max_ms;
  }
  printf("GoldenSniffer: %llu DCIs, %zu RNTIs; %llu slots searched (%.3f ms mean, %.3f ms max), %llu of %llu dropped "
         "(searchers behind)\n",
         (unsigned long long)gs_total, gs_stats.size(), (unsigned long long)searched,
         searched ? t_sum / (double)searched : 0.0, t_max, (unsigned long long)dropped,
         (unsigned long long)submitted);
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
  if (g["threads"]) {
    a.threads = g["threads"].as<uint32_t>();
  }
  if (g["queue"]) {
    a.queue = g["queue"].as<uint32_t>();
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

/* Logs one searched slot's DCIs; called by the searcher threads */
static void gs_log_dcis(const GsJob& job, const gs_dci_t* out, int n)
{
  if (n <= 0) {
    return;
  }
  const double                now = wall_now();
  std::lock_guard<std::mutex> lock(gs_log_mutex);
  for (int i = 0; i < n; i++) {
    const gs_dci_t& d = out[i];
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
          uint32_t k = 4 * h + b;
          v          = (v << 1) | (k < d.nof_bits ? d.payload[k] : 0);
        }
        hex[h] = "0123456789abcdef"[v];
      }
      fprintf(gs_csv, "%.6f,%u,%u,0x%04x,%s,%u,%d,%d,%u,%u,%u,%d,%u,%s,%.3f,%.1f\n", now, job.sfn, job.slot_idx, d.rnti,
              rnti_kind(d.rnti, d.coreset0, d.n_rnti), d.n_id, d.n_rnti, d.coreset0, d.L, d.duration,
              job.first_crb + (uint32_t)d.cand.iv[0].first / NRSCOPE_NSC_PER_RB_NR, d.cand.nof_iv > 1, d.nof_bits, hex,
              d.evm, d.epre > 0 ? 10.0 * std::log10(d.epre) : -999.0);
    }
    if (nrscope_gs_args.feed_rntis && !s.fed && s.dcis >= nrscope_gs_args.min_hits &&
        std::string(rnti_kind(d.rnti, d.coreset0, d.n_rnti)) == "C") {
      gs_feeds.push_back({d.rnti, d.n_id, d.n_rnti});
      s.fed = true;
    }
  }
}

static void gs_searcher_thread(uint32_t id)
{
  char name[16];
  snprintf(name, sizeof(name), "nrs-gs%u", id);
  pthread_setname_np(pthread_self(), name);
  gs_t* gs = gs_create(&gs_pool.args);
  if (gs == nullptr) {
    ERROR("GoldenSniffer: could not create a searcher");
    return;
  }
  gs_dci_t out[GS_MAX_DCI_PER_SLOT];
  while (true) {
    GsJob job;
    {
      std::unique_lock<std::mutex> lock(gs_pool.m);
      gs_pool.cv.wait(lock, [] { return !gs_pool.q.empty(); });
      job = std::move(gs_pool.q.front());
      gs_pool.q.pop_front();
    }
    gs_slot_t sl = {};
    sl.grid      = job.grid.data();
    sl.nof_sc    = job.nof_sc;
    sl.slot_idx  = job.slot_idx;
    sl.ref_crb0  = job.ref_crb0;
    sl.ref_cs0   = job.ref_cs0;
    sl.cs0_nof_sc = job.cs0_nof_sc;
    sl.pci       = job.pci;
    auto t0      = std::chrono::steady_clock::now();
    int  n       = gs_search_slot(gs, &sl, out, GS_MAX_DCI_PER_SLOT);
    double ms    = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    {
      std::lock_guard<std::mutex> lock(gs_pool.m);
      gs_pool.searched++;
      gs_pool.t_sum_ms += ms;
      gs_pool.t_max_ms = std::max(gs_pool.t_max_ms, ms);
    }
    gs_log_dcis(job, out, n);
  }
}

void golden_sniffer_update(WorkState& state, double now)
{
  if (!nrscope_gs_args.enable) {
    return;
  }
  std::lock_guard<std::mutex> lock(gs_log_mutex);
  /* C-RNTIs to the regular DCI decoders, with the scrambling ID read off the DM-RS. They
    need the CORESET from an RRCSetup; until one was seen the RNTIs wait. */
  if (state.rach_found && !gs_feeds.empty()) {
    for (const GsFeed& f : gs_feeds) {
      bool known = false;
      for (uint16_t r : state.known_rntis) {
        known |= r == f.rnti;
      }
      if (known) {
        continue;
      }
      state.known_rntis.push_back(f.rnti);
      state.last_seen.push_back(now);
      state.nof_known_rntis++;
      const auto& sp_ded = state.master_cell_group.sp_cell_cfg.sp_cell_cfg_ded;
      if (f.n_rnti && sp_ded.init_dl_bwp_present && sp_ded.init_dl_bwp.pdcch_cfg_present &&
          sp_ded.init_dl_bwp.pdcch_cfg.type().value == asn1::setup_release_opts::setup) {
        for (const auto& cs : sp_ded.init_dl_bwp.pdcch_cfg.setup().ctrl_res_set_to_add_mod_list) {
          if (cs.pdcch_dmrs_scrambling_id_present) {
            state.pdcch_dmrs_ids_by_rnti[f.rnti][cs.ctrl_res_set_id] = f.n_id;
          }
        }
      }
      printf("GoldenSniffer: RNTI 0x%04x (N_ID %u) added to the DCI decoders\n", f.rnti, f.n_id);
    }
    gs_feeds.clear();
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
  cs0_nof_sc = srsran_coreset_get_bw(&state->coreset0_t) * NRSCOPE_NSC_PER_RB_NR;
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
    /* When CORESET#0 starts at point A the DM-RS cannot tell its DCIs from those of a
      CORESET over the same RBs, and a UE's DCIs of other sizes there (0_1, 1_1) would be
      rejected. Then no list: the most frequent of the sizes that solve settles it. */
    if (ref_cs0 == ref_crb0) {
      css.clear();
      printf("GoldenSniffer: CORESET#0 starts at point A, no CSS size list\n");
    }
  }
  /* The searcher pool, once for all workers (they all see the same cell) */
  std::call_once(gs_pool.started, [&]() {
    gs_pool.sizes          = nrscope_gs_args.dci_sizes;
    gs_pool.css            = css;
    gs_pool.args           = ga;
    gs_pool.args.sizes     = gs_pool.sizes.data();
    gs_pool.args.nof_sizes = (uint32_t)gs_pool.sizes.size();
    gs_pool.args.css_sizes = gs_pool.css.data();
    gs_pool.args.nof_css_sizes = (uint32_t)gs_pool.css.size();
    gs_t* warm = gs_create(&gs_pool.args); // builds the shared GF(2) tables before the threads start
    gs_destroy(warm);
    for (uint32_t i = 0; i < std::max(1u, nrscope_gs_args.threads); i++) {
      std::thread(gs_searcher_thread, i).detach();
    }
  });
  ready = true;
  printf("GoldenSniffer ready: grid of %u PRBs from CRB %u, point A at grid subcarrier %d, CORESET#0 at %d, PCI %u, "
         "CSS sizes",
         grid_carrier.nof_prb, first_crb, ref_crb0, ref_cs0, pci);
  for (uint32_t c : css) {
    printf(" %u", c);
  }
  printf("%s; %u searcher threads, queue %u slots\n", css.empty() ? " none" : "", nrscope_gs_args.threads,
         nrscope_gs_args.queue);
  return SRSRAN_SUCCESS;
}

void GoldenSnifferDecoder::Submit(srsran_slot_cfg_t* slot, WorkState* state)
{
  if (!ready) {
    return;
  }
  /* Demodulate the PDCCH symbols here (cheap: a few FFTs) since the worker's buffer is
    only valid during its slot, then hand a copy to the searchers. One symbol more than
    searched: the duration detector compares it. */
  const uint32_t nsym                    = nrscope_gs_args.max_duration + 1;
  uint32_t       mask[SRSRAN_MAX_PORTS] = {(1u << nsym) - 1u};
  uint32_t       done[SRSRAN_MAX_PORTS] = {};
  srsran_ue_dl_nr_estimate_fft_nrscope_symbols(&ue_dl, slot, arg_scs_grid, mask, done);

  const uint32_t nof_sc = grid_carrier.nof_prb * NRSCOPE_NSC_PER_RB_NR;
  GsJob          job;
  job.sfn       = state->sfn;
  job.slot_idx  = slot->idx % slots_per_frame;
  job.nof_sc    = nof_sc;
  job.nof_symb  = nsym;
  job.first_crb = first_crb;
  job.pci       = pci;
  job.ref_crb0  = ref_crb0;
  job.ref_cs0   = ref_cs0;
  job.cs0_nof_sc = cs0_nof_sc;

  /* Debug: GS_DUMP=<file> appends the searched symbols of every slot, as
    [sfn u32, slot u32, nof_sc u32, nof_symb u32, cf32 grid], for offline analysis */
  static FILE*      dump = nullptr;
  static std::mutex dump_mutex;
  if (getenv("GS_DUMP") != nullptr) {
    std::lock_guard<std::mutex> lock(dump_mutex);
    if (dump == nullptr) {
      dump = fopen(getenv("GS_DUMP"), "wb");
    }
    static uint64_t nof_dumped = 0;
    const uint64_t  max_dump   = getenv("GS_DUMP_SLOTS") ? strtoull(getenv("GS_DUMP_SLOTS"), nullptr, 10) : UINT64_MAX;
    if (dump != nullptr && nof_dumped++ < max_dump) {
      uint32_t hdr[4] = {job.sfn, job.slot_idx, nof_sc, nsym};
      fwrite(hdr, sizeof(hdr), 1, dump);
      fwrite(ue_dl.sf_symbols[0], sizeof(cf_t), (size_t)nof_sc * nsym, dump);
    }
  }

  /* Replaying a file nothing is dropped: wait for room. Live, never wait: drop. */
  while (nrscope_lossless.load(std::memory_order_relaxed) && !nrscope_stop.load(std::memory_order_relaxed)) {
    {
      std::lock_guard<std::mutex> lock(gs_pool.m);
      if (gs_pool.q.size() < std::max(1u, nrscope_gs_args.queue)) {
        break;
      }
    }
    std::this_thread::sleep_for(std::chrono::microseconds(100));
  }
  {
    std::lock_guard<std::mutex> lock(gs_pool.m);
    gs_pool.submitted++;
    if (gs_pool.q.size() >= std::max(1u, nrscope_gs_args.queue)) {
      gs_pool.dropped++;
      return;
    }
  }
  job.grid.assign(ue_dl.sf_symbols[0], ue_dl.sf_symbols[0] + (size_t)nof_sc * nsym);
  {
    std::lock_guard<std::mutex> lock(gs_pool.m);
    gs_pool.q.push_back(std::move(job));
  }
  gs_pool.cv.notify_one();
}
