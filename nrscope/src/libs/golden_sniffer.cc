/*
 * GoldenSniffer in NR-Scope; see nrscope/hdr/golden_sniffer.h.
 */
#include "nrscope/hdr/golden_sniffer.h"
#include "nrscope/hdr/gs_fdra.h"
#include "nrscope/hdr/sensing/nrscope_sensing.h"
#include <algorithm>
#include <atomic>
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
/* Execution time of one GoldenSniffer stage: count, mean, max, and a log histogram
  (10 bins per decade from 1 us) for the median and the 99th percentile */
struct GsTimer {
  /* lock-free: the real-time workers add to it */
  const char*           name;
  std::atomic<uint64_t> n{0}, sum_ns{0}, max_ns{0};
  std::atomic<uint32_t> hist[70];
  explicit GsTimer(const char* nm) : name(nm)
  {
    for (auto& h : hist) {
      h.store(0, std::memory_order_relaxed);
    }
  }
  void add(double ms)
  {
    const uint64_t ns = ms > 0 ? (uint64_t)(ms * 1e6) : 0;
    int            b  = ms > 0 ? (int)std::floor(10.0 * std::log10(ms * 1e3)) : 0;
    b                 = std::max(0, std::min(69, b));
    n.fetch_add(1, std::memory_order_relaxed);
    sum_ns.fetch_add(ns, std::memory_order_relaxed);
    hist[b].fetch_add(1, std::memory_order_relaxed);
    uint64_t m = max_ns.load(std::memory_order_relaxed);
    while (ns > m && !max_ns.compare_exchange_weak(m, ns, std::memory_order_relaxed)) {
    }
  }
  std::string str()
  {
    const uint64_t cnt = n.load(std::memory_order_relaxed);
    if (cnt == 0) {
      return std::string("  ") + name + ": -\n";
    }
    const double mx = 1e-6 * (double)max_ns.load(std::memory_order_relaxed);
    auto         q  = [&](double p) {
      uint64_t want = (uint64_t)std::ceil(p * (double)cnt), acc = 0;
      for (int b = 0; b < 70; b++) {
        acc += hist[b].load(std::memory_order_relaxed);
        if (acc >= want) {
          return std::pow(10.0, (b + 1) / 10.0) * 1e-3; /* bin upper edge, ms */
        }
      }
      return mx;
    };
    char buf[200];
    snprintf(buf, sizeof(buf), "  %-26s n %9llu  mean %8.3f  p50 <%7.3f  p99 <%7.3f  max %8.3f ms\n", name,
             (unsigned long long)cnt, 1e-6 * (double)sum_ns.load(std::memory_order_relaxed) / (double)cnt, q(0.5),
             q(0.99), mx);
    return buf;
  }
};
static GsTimer gs_t_fft("worker: PDCCH FFT"), gs_t_track("worker: tracked decode"), gs_t_enq("worker: copy + queue"),
    gs_t_wait("searcher: queue wait"), gs_t_search0("searcher: slot, no DCI"), gs_t_search1("searcher: slot with DCIs"),
    gs_t_lat_search("latency: slot -> search DCI"), gs_t_lat_track("latency: slot -> tracked DCI"),
    gs_t_update("scheduler: GS update"), gs_t_submit("worker: GS total");

static double gs_steady_s()
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

/* A mutex with priority inheritance: the real-time workers and the searchers (normal
  priority) share the job queue; a searcher preempted while holding it is boosted rather
  than leaving a worker stuck behind it (priority inversion). */
class GsPiMutex
{
public:
  GsPiMutex()
  {
    pthread_mutexattr_t a;
    pthread_mutexattr_init(&a);
    pthread_mutexattr_setprotocol(&a, PTHREAD_PRIO_INHERIT);
    pthread_mutex_init(&m_, &a);
    pthread_mutexattr_destroy(&a);
  }
  void lock() { pthread_mutex_lock(&m_); }
  void unlock() { pthread_mutex_unlock(&m_); }
  bool try_lock() { return pthread_mutex_trylock(&m_) == 0; }

private:
  pthread_mutex_t m_;
};

/* what the searchers need of a converged UE to leave its DCIs to the tracked decode */
struct GsJobTrack {
  uint16_t rnti;
  uint8_t  nof_sizes;
  uint8_t  sizes[8];
};
#define GS_JOB_TRACKS 512
struct GsJob {
  uint32_t sfn = 0, slot_idx = 0, nof_sc = 0, nof_symb = 0, first_crb = 0, pci = 0;
  int32_t  ref_crb0 = GS_NO_REF, ref_cs0 = GS_NO_REF;
  uint32_t cs0_nof_sc = 0;
  cf_t*    grid       = nullptr; /* one of the pool's preallocated buffers */
  int      buf        = -1;
  double   t_submit = 0, t_enq = 0; /* steady clock, s */
  /* converged UEs in this slot: their DCIs of a learned size come from the tracked decode */
  GsJobTrack tracked[GS_JOB_TRACKS];
  uint32_t   nof_tracked = 0;
};
static bool gs_job_tracks(const GsJob& job, const gs_dci_t& d)
{
  for (uint32_t k = 0; k < job.nof_tracked; k++) {
    const GsJobTrack& t = job.tracked[k];
    if (t.rnti == d.rnti && !d.coreset0) {
      for (uint32_t z = 0; z < t.nof_sizes; z++) {
        if (t.sizes[z] == d.nof_bits) {
          return true;
        }
      }
    }
  }
  return false;
}

/* The searcher pool: its own threads and its own bounded queue. GoldenSniffer is never
  on the worker's critical path; when the pool falls behind, it drops its own slots. */
struct GsPool {
  std::once_flag              started;
  GsPiMutex                   m;
  std::condition_variable_any cv;
  /* jobs: preallocated (no allocation on the workers' path); free and queued ones are
    indices into them, in fixed rings */
  std::vector<GsJob>             jobs;
  std::vector<std::vector<cf_t>> bufs;
  std::vector<int>               free_idx; /* stack */
  std::vector<int>               ring;     /* queued, FIFO */
  size_t                         r_head = 0, r_len = 0;
  std::atomic<uint64_t>          busy{0};  /* slots not queued: the lock was taken */
  size_t qsize() const { return r_len; }
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

/* Fallback (golden_sniffer.h) */
std::atomic<bool>     gs_pool_active{true}; // whole-slot search wanted (windows, armed/learning UEs)
std::atomic<bool>     gs_any_converged{false};
std::atomic<double>   gs_ready_time{0};     // when the first worker's GoldenSniffer was ready
std::atomic<uint32_t> gs_spf{0};            // slots per frame
/* the cell's FDRA hypotheses (gs_fdra.h), set once by the first worker */
static GsFdraHyp      gs_fdra_h[GS_FDRA_MAX_HYP];
static std::atomic<uint32_t> gs_fdra_nh{0};
/* the cell's FDRA hypothesis per DL DCI size, once UEs of that size agree (-1: none yet):
   a UE whose config changes its DCI size gets the layout of its new size */
static std::atomic<int16_t> gs_fdra_size_hyp[256];
static struct GsFdraSizeInit {
  GsFdraSizeInit()
  {
    for (auto& v : gs_fdra_size_hyp) {
      v.store(-1);
    }
  }
} gs_fdra_size_init;
static std::once_flag gs_fdra_once;
static FILE*          gs_grants_csv = nullptr;

/* GNU complex helpers for C++ (the C99 conjf/cargf/crealf are not declared here) */
static inline cf_t gs_cj(cf_t a)
{
  cf_t r;
  __real__ r = __real__ a;
  __imag__ r = -__imag__ a;
  return r;
}
static inline float gs_arg(cf_t a)
{
  return atan2f(__imag__ a, __real__ a);
}
static inline float gs_re(cf_t a)
{
  return __real__ a;
}
static FILE*          gs_survey_csv = nullptr; // golden_sniffer.survey_csv
static std::mutex     gs_survey_mutex;
/* grants of tracked UEs: DM-RS ports (bit i = port 1000 + i), CDM groups 2, extra symbols */
static uint64_t gs_grant_ports[4] = {}, gs_grant_cdm2 = 0, gs_grant_add[4] = {};
/* Replay (lossless file input) runs at any speed: GoldenSniffer's timers (startup, sweeps,
  arming, quiet, keep, convergence) then run on air time, from the SFN and slot */
static double  gs_air_ready = -1;
static int64_t gs_air_abs   = -1;
static double gs_air_now(const SlotResult& r)
{
  const uint32_t spf = gs_spf.load();
  if (spf == 0) {
    return 0;
  }
  const int64_t wrap = 1024LL * spf, a = (int64_t)r.outcome.sfn * spf + r.slot.idx % spf;
  if (gs_air_abs < 0) {
    gs_air_abs = a;
  } else {
    int64_t base = gs_air_abs - gs_air_abs % wrap + a, best = base;
    for (int64_t d : {-wrap, wrap}) {
      if (std::llabs(base + d - gs_air_abs) < std::llabs(best - gs_air_abs)) {
        best = base + d;
      }
    }
    if (best > gs_air_abs) { /* results can come out of order: the clock only moves forward */
      gs_air_abs = best;
    }
  }
  return (double)gs_air_abs / (100.0 * spf);
}
std::atomic<uint32_t> gs_first_crb{0};
std::vector<gs_dci_t> gs_sightings;         // searcher C-RNTI DCIs for the scheduler, under gs_log_mutex
std::string           gs_track_summary;     // per-UE fallback states for the report, under gs_log_mutex
double                gs_sweep_until = 0;   // an extra sweep, started when a tracked UE went quiet

const char* state_name(BlindUe::State s)
{
  static const char* n[] = {"RRC", "ARMED", "STALE", "CONVERGED", "DISCOVERED"};
  return n[s];
}

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

void print_report_locked(bool final = false)
{
  uint64_t searched, dropped, submitted;
  double   t_sum, t_max;
  {
    std::lock_guard<GsPiMutex> lock(gs_pool.m);
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
  if (final)
  printf("GoldenSniffer execution time:\n%s%s%s%s%s%s%s%s%s%s", gs_t_submit.str().c_str(), gs_t_fft.str().c_str(),
         gs_t_track.str().c_str(), gs_t_enq.str().c_str(), gs_t_wait.str().c_str(), gs_t_search0.str().c_str(),
         gs_t_search1.str().c_str(), gs_t_update.str().c_str(), gs_t_lat_search.str().c_str(),
         gs_t_lat_track.str().c_str());
  if (!gs_track_summary.empty()) {
    printf("%s", gs_track_summary.c_str());
  }
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
  if (g["size_min"]) {
    a.size_min = g["size_min"].as<uint32_t>();
  }
  if (g["size_max"]) {
    a.size_max = g["size_max"].as<uint32_t>();
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
  if (g["mode"]) {
    const std::string m = g["mode"].as<std::string>();
    if (m != "fallback" && m != "always") {
      std::cerr << "golden_sniffer: mode must be fallback or always" << std::endl;
      exit(EXIT_FAILURE);
    }
    a.always = m == "always";
  }
  if (g["startup_s"]) {
    a.startup_s = g["startup_s"].as<double>();
  }
  if (g["scan_period_s"]) {
    a.scan_period_s = g["scan_period_s"].as<double>();
  }
  if (g["scan_window_s"]) {
    a.scan_window_s = g["scan_window_s"].as<double>();
  }
  if (g["trigger_on_srb"]) {
    a.trigger_on_srb = g["trigger_on_srb"].as<bool>();
  }
  if (g["arm_s"]) {
    a.arm_s = g["arm_s"].as<double>();
  }
  if (g["stale_hits"]) {
    a.stale_hits = g["stale_hits"].as<uint32_t>();
  }
  if (g["converge_hits"]) {
    a.converge_hits = g["converge_hits"].as<uint32_t>();
  }
  if (g["converge_s"]) {
    a.converge_s = g["converge_s"].as<double>();
  }
  if (g["keep_s"]) {
    a.keep_s = g["keep_s"].as<double>();
  }
  if (g["track_full_band"]) {
    a.track_full_band = g["track_full_band"].as<bool>();
  }
  if (g["discover_in_track"]) {
    a.discover_in_track = g["discover_in_track"].as<bool>();
  }
  if (g["searcher_idle"]) {
    a.searcher_idle = g["searcher_idle"].as<bool>();
  }
  if (g["soft"]) {
    a.soft = g["soft"].as<bool>();
  }
  if (g["fdra"]) {
    a.fdra = g["fdra"].as<bool>();
  }
  if (g["grants_csv"]) {
    a.grants_csv = g["grants_csv"].as<std::string>();
  }
  if (g["survey_csv"]) {
    a.survey_csv = g["survey_csv"].as<std::string>();
  }
  if (g["quiet_s"]) {
    a.quiet_s = g["quiet_s"].as<double>();
  }
  if (a.enable) {
    /* the GF(2) tables, built in the background while the radio starts and the cell is
      acquired: no real-time thread builds or waits for one */
    {
      std::vector<uint32_t> first = a.dci_sizes;
      first.insert(first.end(), a.css_sizes.begin(), a.css_sizes.end());
      gs_gf2_prebuild_start(a.size_min, a.size_max, first.data(), (uint32_t)first.size());
    }
    gs_csv = fopen(a.csv.c_str(), "w");
    if (!a.survey_csv.empty()) {
      gs_survey_csv = fopen(a.survey_csv.c_str(), "w");
      if (gs_survey_csv != nullptr) {
        setvbuf(gs_survey_csv, nullptr, _IOFBF, 1 << 22);
        fprintf(gs_survey_csv, "sfn,slot,n_dl,n_ul,rbs_comb0_sym1..13,rbs_comb1_sym1..13,l0_bitmap,l0_occ_median,"
                               "l0_comb1_power_median,l0_bitmap08,l0_type2_rbs,step_inside,step_crb4,step_crb2,"
                               "step_odd (n:deg:dB),csi_rs,ue_nid_tests\n");
      }
    }
    if (a.fdra) {
      std::string gp = a.grants_csv;
      if (gp.empty()) {
        gp            = a.csv;
        const size_t dot = gp.rfind(".csv");
        gp            = (dot == std::string::npos ? gp : gp.substr(0, dot)) + "_grants.csv";
      }
      gs_grants_csv = fopen(gp.c_str(), "w");
      if (gs_grants_csv != nullptr) {
        fprintf(gs_grants_csv,
                "timestamp,sfn,slot,rnti,bwp_start_crb,rb_start,nof_rb,dmrs_coherence,hyp,sensed,ports,cdm,dmrs_add,"
                "occ,comb1_power,comb1_coherence,dmrs_symbols,prb_runs\n");
      }
    }
    if (gs_csv == nullptr) {
      std::cerr << "golden_sniffer: cannot write " << a.csv << std::endl;
    } else {
      fprintf(gs_csv, "timestamp,sfn,slot,rnti,kind,n_id,n_rnti,coreset0,al,duration,first_crb,interleaved,size,"
                      "payload_hex,evm,epre_db,source\n");
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
  print_report_locked(true);
  if (gs_survey_csv != nullptr) {
    std::lock_guard<std::mutex> lk(gs_survey_mutex);
    fflush(gs_survey_csv);
  }
  const uint64_t ng = gs_grant_ports[1] + gs_grant_ports[2] + gs_grant_ports[3]; // grants with DM-RS seen
  if (ng > 0) {
    printf("GoldenSniffer grants: %lu; DM-RS ports on air: 1000 %lu, 1001 %lu, 1000+1001 %lu; CDM groups without data 2: "
           "%lu; extra DM-RS symbols pos1 %lu, pos2 %lu, pos3 %lu\n",
           (unsigned long)ng, (unsigned long)gs_grant_ports[1], (unsigned long)gs_grant_ports[2],
           (unsigned long)gs_grant_ports[3], (unsigned long)gs_grant_cdm2, (unsigned long)gs_grant_add[1],
           (unsigned long)gs_grant_add[2], (unsigned long)gs_grant_add[3]);
  }
  if (gs_csv != nullptr) {
    fflush(gs_csv);
  }
}

/* Logs one DCI (CSV, counts); caller holds gs_log_mutex. source: "search" or "tracked". */
static void gs_log_one_locked(double now, uint32_t sfn, uint32_t slot_idx, uint32_t first_crb, const gs_dci_t& d,
                              const char* source)
{
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
    fprintf(gs_csv, "%.6f,%u,%u,0x%04x,%s,%u,%d,%d,%u,%u,%u,%d,%u,%s,%.3f,%.1f,%s\n", now, sfn, slot_idx, d.rnti,
            rnti_kind(d.rnti, d.coreset0, d.n_rnti), d.n_id, d.n_rnti, d.coreset0, d.L, d.duration,
            first_crb + (uint32_t)d.cand.iv[0].first / NRSCOPE_NSC_PER_RB_NR, d.cand.nof_iv > 1, d.nof_bits, hex, d.evm,
            d.epre > 0 ? 10.0 * std::log10(d.epre) : -999.0, source);
  }
}

/* DCIs found by the searchers, for the scheduler to log and learn from. The searchers run
  at idle priority: they must not hold a lock the pipeline waits for, so they only append
  here (a few instructions, priority inheritance); the scheduler drains it with try_lock. */
struct GsFound {
  double   t;
  uint32_t sfn, slot_idx, first_crb;
  gs_dci_t d;
};
static GsPiMutex            gs_found_m;
static std::vector<GsFound> gs_found, gs_found_drain;
static bool                 gs_found_reserved = [] {
  gs_found.reserve(8192);
  gs_found_drain.reserve(8192);
  return true;
}();

static void gs_log_dcis(const GsJob& job, const gs_dci_t* out, int n)
{
  if (n <= 0) {
    return;
  }
  const double               now = wall_now();
  std::lock_guard<GsPiMutex> lock(gs_found_m);
  for (int i = 0; i < n; i++) {
    const gs_dci_t& d = out[i];
    if (gs_job_tracks(job, d)) {
      continue; // the worker's tracked decode reports it
    }
    if (gs_found.size() < 8192) {
      gs_found.push_back({now, job.sfn, job.slot_idx, job.first_crb, d});
    }
  }
}

/* The scheduler's side, under gs_log_mutex: logs the searchers' DCIs, hands the C-RNTI
  ones to the fallback state machine (gs_sightings), feeds the regular decoders */
static void gs_drain_found_locked()
{
  {
    std::unique_lock<GsPiMutex> l(gs_found_m, std::try_to_lock);
    if (!l) {
      return; /* next slot */
    }
    gs_found_drain.swap(gs_found);
  }
  for (const GsFound& f : gs_found_drain) {
    const gs_dci_t& d = f.d;
    gs_log_one_locked(f.t, f.sfn, f.slot_idx, f.first_crb, d, "search");
    const bool c_rnti = std::string(rnti_kind(d.rnti, d.coreset0, d.n_rnti)) == "C";
    if (c_rnti && gs_sightings.size() < 4096) {
      gs_sightings.push_back(d);
    }
    RntiStats& s = gs_stats[d.rnti];
    if (nrscope_gs_args.feed_rntis && !s.fed && s.dcis >= nrscope_gs_args.min_hits && c_rnti) {
      gs_feeds.push_back({d.rnti, d.n_id, d.n_rnti});
      s.fed = true;
    }
  }
  gs_found_drain.clear(); /* keeps its capacity */
}

static void gs_searcher_thread(uint32_t id)
{
  char name[16];
  snprintf(name, sizeof(name), "nrs-gs%u", id);
  pthread_setname_np(pthread_self(), name);
  /* idle priority: the blind search only takes CPU time nobody else wants. Unpinned and
    at normal priority it landed on the workers' CPUs and slowed them (dropped slots in
    every sweep): NR-Scope's workers run at normal priority too. */
  if (nrscope_gs_args.searcher_idle) {
    struct sched_param sp = {};
    pthread_setschedparam(pthread_self(), SCHED_IDLE, &sp);
  }
  gs_t* gs = gs_create(&gs_pool.args);
  if (gs == nullptr) {
    ERROR("GoldenSniffer: could not create a searcher");
    return;
  }
  gs_dci_t out[GS_MAX_DCI_PER_SLOT];
  while (true) {
    int idx;
    {
      std::unique_lock<GsPiMutex> lock(gs_pool.m);
      gs_pool.cv.wait(lock, [] { return gs_pool.r_len > 0; });
      idx            = gs_pool.ring[gs_pool.r_head];
      gs_pool.r_head = (gs_pool.r_head + 1) % gs_pool.ring.size();
      gs_pool.r_len--;
    }
    const GsJob& job = gs_pool.jobs[idx];
    gs_slot_t sl = {};
    sl.grid      = job.grid;
    sl.nof_sc    = job.nof_sc;
    sl.slot_idx  = job.slot_idx;
    sl.ref_crb0  = job.ref_crb0;
    sl.ref_cs0   = job.ref_cs0;
    sl.cs0_nof_sc = job.cs0_nof_sc;
    sl.pci       = job.pci;
    auto t0      = std::chrono::steady_clock::now();
    gs_t_wait.add(1e3 * (gs_steady_s() - job.t_enq));
    int  n       = gs_search_slot(gs, &sl, out, GS_MAX_DCI_PER_SLOT);
    double ms    = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    (n > 0 ? gs_t_search1 : gs_t_search0).add(ms);
    if (n > 0) {
      gs_t_lat_search.add(1e3 * (gs_steady_s() - job.t_submit));
    }
    gs_log_dcis(job, out, n);
    {
      std::lock_guard<GsPiMutex> lock(gs_pool.m);
      gs_pool.searched++;
      gs_pool.t_sum_ms += ms;
      gs_pool.t_max_ms = std::max(gs_pool.t_max_ms, ms);
      gs_pool.free_idx.push_back(idx); /* capacity reserved: no allocation */
    }
  }
}

/* Learns one DCI of a UE: its PDCCH N_ID and scrambling, the CORESET extent its DCIs were
  seen in, its ALs and its DCI sizes (the most frequent four). */
static void gs_learn_into(BlindUe& e, const gs_dci_t& d, double now)
{
  /* The tracked decode targets the UE-specific CORESET. DCIs in CORESET#0 (Msg4, the
    RRCSetup, fallback DCIs while the UE is still on the common config) are interleaved
    and have their own size: learning them would mark the UE interleaved for good and
    keep it from converging, and add their size and extent to the track. */
  if (d.coreset0) {
    return;
  }
  gs_track_t& t   = e.track;
  int32_t     lo  = d.cand.iv[0].first, hi = d.cand.iv[0].last;
  for (uint32_t i = 1; i < d.cand.nof_iv; i++) {
    lo = std::min(lo, d.cand.iv[i].first);
    hi = std::max(hi, d.cand.iv[i].last);
  }
  if (e.hits == 0) {
    t.sc_lo = lo;
    t.sc_hi = hi;
  } else {
    t.sc_lo = std::min(t.sc_lo, lo);
    t.sc_hi = std::max(t.sc_hi, hi);
  }
  t.rnti     = d.rnti;
  t.n_id     = d.n_id;
  t.n_rnti   = d.n_rnti;
  t.coreset0 = d.coreset0;
  t.duration = d.duration;
  for (uint32_t l = 0; l <= GS_MAX_L_LOG2; l++) {
    if (d.L == (1u << l)) {
      t.al_mask |= 1u << l;
    }
  }
  e.size_count[d.nof_bits]++;
  std::vector<std::pair<uint32_t, uint32_t> > by(e.size_count.begin(), e.size_count.end());
  std::sort(by.begin(), by.end(), [](const auto& x, const auto& y) { return x.second > y.second; });
  t.nof_sizes = 0;
  for (const auto& z : by) {
    if (t.nof_sizes < 8) {
      t.sizes[t.nof_sizes++] = z.first;
    }
  }
  e.interleaved = e.interleaved || d.cand.nof_iv > 1;
  /* The first bit of DCI 0_x/1_x is the format identifier (TS 38.212 7.3.1): 0 UL, 1 DL.
    UL and DL DCIs differ in size, so a UE that converged on one direction alone would
    have the other's size missing from its track and lose those DCIs. */
  if (d.nof_bits > 0) {
    e.dirs |= d.payload[0] ? 2 : 1;
  }
  if (e.hits == 0) {
    e.learn_since = now;
  }
  e.hits++;
}

/* Whether a DCI fits what was learned of the UE: the tracked decode would find it */
static bool gs_agrees(const gs_track_t& t, const gs_dci_t& d)
{
  if (d.n_id != t.n_id || d.n_rnti != t.n_rnti || d.duration != t.duration) {
    return false;
  }
  for (uint32_t z = 0; z < t.nof_sizes; z++) {
    if (t.sizes[z] == d.nof_bits) {
      return true;
    }
  }
  return false;
}

static std::string gs_sizes_str(const gs_track_t& t)
{
  std::string s;
  for (uint32_t z = 0; z < t.nof_sizes; z++) {
    s += (z ? "/" : "") + std::to_string(t.sizes[z]);
  }
  return s.empty() ? "-" : s;
}

static void gs_start_shadow(BlindUe& e, uint16_t rnti, double now, const char* why)
{
  if (e.shadow) {
    return;
  }
  e.shadow       = std::make_shared<BlindUe>();
  e.shadow_since = now;
  printf("GoldenSniffer: UE 0x%04x CONVERGED, re-learning (%s)\n", rnti, why);
}

/* Learns one DCI of a UE. A converged UE learns only what agrees with its track; a DCI
  that does not (a reconfiguration changed its sizes, its CORESET or its scrambling)
  starts a fresh learning in its shadow, which replaces the track once it converges. */
static void gs_set_state(BlindUe& e, uint16_t rnti, BlindUe::State to, const char* why);

static void gs_learn(BlindUe& e, const gs_dci_t& d, double now)
{
  if (d.coreset0) {
    return;
  }
  if (d.nof_bits > 0) {
    e.last_dir[d.payload[0] ? 1 : 0] = now;
  }
  if (e.state != BlindUe::CONVERGED) {
    gs_learn_into(e, d, now);
    return;
  }
  if (d.cand.nof_iv > 1) {
    /* Its CORESET is interleaved: the first DCIs it converged on were AL1 (one REG bundle
      looks contiguous), and the tracked decode (non-interleaved) misses the others. Back
      to learning: the blind search, which handles interleaving, decodes it. */
    e.interleaved = true;
    e.shadow.reset();
    gs_set_state(e, d.rnti, e.from_rach ? BlindUe::STALE : BlindUe::DISCOVERED,
                 "an interleaved DCI: its CORESET is interleaved, left to the blind search");
    gs_learn_into(e, d, now);
    return;
  }
  const bool agrees = gs_agrees(e.track, d);
  if (!agrees) {
    char why[96];
    snprintf(why, sizeof(why), "a DCI of size %u, N_ID %u, %u symbol(s), not in its track %s", d.nof_bits, d.n_id,
             d.duration, gs_sizes_str(e.track).c_str());
    gs_start_shadow(e, d.rnti, now, why);
  }
  if (e.shadow) {
    gs_learn_into(*e.shadow, d, now);
  }
  if (agrees) {
    gs_learn_into(e, d, now);
  }
}

static bool gs_is_known(const WorkState& state, uint16_t rnti)
{
  return std::find(state.known_rntis.begin(), state.known_rntis.end(), rnti) != state.known_rntis.end();
}

static void gs_set_state(BlindUe& e, uint16_t rnti, BlindUe::State to, const char* why)
{
  if (e.state != to) {
    printf("GoldenSniffer: UE 0x%04x %s -> %s (%s)\n", rnti, state_name(e.state), state_name(to), why);
    e.state = to;
  }
}

/* Discovery windows: at start and periodic sweeps, all RNTIs */
static bool gs_discovery_window(double now)
{
  const GoldenSnifferArgs& a     = nrscope_gs_args;
  const double             ready = nrscope_lossless.load() ? gs_air_ready : gs_ready_time.load();
  if (a.always) {
    return true;
  }
  if (ready <= 0) {
    return false;
  }
  if (now < gs_sweep_until) {
    return true;
  }
  const double t = now - ready;
  if (t < a.startup_s) {
    return true;
  }
  return a.scan_period_s > 0 && std::fmod(t - a.startup_s, a.scan_period_s) < a.scan_window_s;
}

void golden_sniffer_update(WorkState&                   state,
                           const SlotResult&            result,
                           const std::vector<uint16_t>& regular_hits,
                           double                       now)
{
  if (!nrscope_gs_args.enable) {
    return;
  }
  struct UpdTimer {
    double t0 = gs_steady_s();
    ~UpdTimer() { gs_t_update.add(1e3 * (gs_steady_s() - t0)); }
  } upd_timer;
  const double now_wall = now; /* the regular RNTI list keeps wall-clock times */
  if (nrscope_lossless.load() && gs_spf.load() > 0) {
    now = gs_air_now(result);
    if (gs_air_ready < 0 && gs_ready_time.load() > 0) {
      gs_air_ready = now;
    }
  }
  const GoldenSnifferArgs&    a = nrscope_gs_args;
  std::lock_guard<std::mutex> lock(gs_log_mutex);
  auto&                       ues = state.blind_ues;
  gs_drain_found_locked();

  /* A fresh RRCSetup: the UE is decoded from it again, whatever it was before. Armed at
    once: a reconfiguration follows within a second, its SRB1 PDU may go with a DCI the
    RACH decoder does not search, and the regular decoders may lose the UE right after. If
    they keep finding it, it falls back to RRC after arm_s. */
  if (result.rach_result && result.found_rach) {
    for (uint16_t r : result.new_rntis_found) {
      auto it = ues.find(r);
      if (it != ues.end() && it->second.state != BlindUe::RRC && it->second.state != BlindUe::ARMED) {
        printf("GoldenSniffer: UE 0x%04x %s -> RRC (new RRCSetup)\n", r, state_name(it->second.state));
      }
      BlindUe e;
      e.from_rach = true;
      e.state     = BlindUe::ARMED;
      e.armed_at  = now;
      ues[r]      = e;
      printf("GoldenSniffer: UE 0x%04x ARMED (new RRCSetup: a reconfiguration follows)\n", r);
    }
  }

  /* The regular decoders found these */
  for (uint16_t r : regular_hits) {
    auto it = ues.find(r);
    if (it != ues.end()) {
      it->second.last_regular    = now;
      it->second.blind_since_reg = 0;
    }
  }

  /* An RRC message to the UE on SRB1, after its RRCSetup: a reconfiguration may follow */
  if (a.trigger_on_srb) {
    for (const SrbEvent& ev : result.srb_events) {
      if (ev.lcid != 1 || !gs_is_known(state, ev.rnti)) {
        continue;
      }
      BlindUe& e = ues[ev.rnti];
      e.from_rach = true;
      e.armed_at  = now;
      if (e.state == BlindUe::RRC) {
        char why[96];
        snprintf(why, sizeof(why), "RRC PDU on SRB1, %u bytes, PDCP SN %d", ev.len, ev.pdcp_sn);
        gs_set_state(e, ev.rnti, BlindUe::ARMED, why);
      } else if (e.state == BlindUe::CONVERGED) {
        char why[96];
        snprintf(why, sizeof(why), "RRC PDU on SRB1, %u bytes, PDCP SN %d", ev.len, ev.pdcp_sn);
        gs_start_shadow(e, ev.rnti, now, why);
      }
    }
  }

  /* DCIs of converged UEs, decoded in the workers */
  const uint32_t first_crb = gs_first_crb.load();
  if (!result.tracked_dcis.empty() && result.gs_t_submit > 0) {
    gs_t_lat_track.add(1e3 * (gs_steady_s() - result.gs_t_submit));
  }
  for (const gs_dci_t& d : result.region_dcis) {
    gs_log_one_locked(now, result.outcome.sfn, d.slot_idx, first_crb, d, "region");
  }
  /* FDRA learning: per UE, the mean lift of every hypothesis over its DL DCIs (only the UEs
    whose PDSCH beam reaches us can show it); one that stands out is kept. Locked hypotheses
    are remembered per DL DCI size, and a UE whose own DCIs never show it takes the one two
    other UEs of its size were locked on. */
  static std::map<uint16_t, std::map<int16_t, uint32_t>> fdra_by_size;
  auto describe = [](const GsFdraHyp& hy) {
    char b[160];
    snprintf(b, sizeof(b), "%s%s, BWP %u PRBs at CRB %u, field at bit %u (%u bits), PDSCH DM-RS N_ID %s",
             hy.type == 1 ? "type 1" : hy.type == 0 ? "type 0" : "dynamic",
             hy.type == 1 ? "" : (hy.rbg_cfg == 1 ? " rbg config 1" : " rbg config 2"), hy.nbwp, hy.bwp_start,
             hy.off, hy.width, hy.nid_ue ? "= its PDCCH DM-RS N_ID" : "= PCI");
    return std::string(b);
  };
  for (const GsFdraObs& o : result.gs_fdra_obs) {
    auto it = ues.find(o.rnti);
    if (it == ues.end() || it->second.fdra_lock >= 0) {
      continue;
    }
    BlindUe& e = it->second;
    if (e.fdra_learn_size == 0) {
      e.fdra_learn_size = o.size;
    } else if (e.fdra_learn_size != o.size) {
      continue; // another DCI size: another layout, not mixed into this one
    }
    for (uint32_t k = 0; k < o.nof && k < GS_FDRA_MAX_HYP; k++) {
      if (o.score[k] > -1.5f) {
        e.fdra_sum[k] += o.score[k];
        e.fdra_cnt[k]++;
      }
    }
    e.fdra_n++;
    if (e.fdra_n < 30 || e.fdra_n % 5 != 0) {
      continue;
    }
    int   best = -1;
    float b1 = -9, b2 = -9;
    for (uint32_t k = 0; k < o.nof && k < GS_FDRA_MAX_HYP; k++) {
      if (e.fdra_cnt[k] < 20) {
        continue;
      }
      const float m = e.fdra_sum[k] / (float)e.fdra_cnt[k];
      if (m > b1) {
        b2   = b1;
        b1   = m;
        best = (int)k;
      } else if (m > b2) {
        b2 = m;
      }
    }
    if (best >= 0 && b1 >= 0.5f && b1 - b2 >= 0.25f) {
      e.fdra_lock      = (int16_t)best;
      e.fdra_lock_size = o.size;
      fdra_by_size[o.size][(int16_t)best]++;
      {
        int16_t  h = -1;
        uint32_t c = 0, tot = 0;
        for (const auto& kv : fdra_by_size[o.size]) {
          tot += kv.second;
          if (kv.second > c) {
            c = kv.second;
            h = kv.first;
          }
        }
        if (o.size < 256 && c >= 2 && c * 4 >= tot * 3) {
          gs_fdra_size_hyp[o.size].store(h);
        }
      }
      printf("GoldenSniffer: UE 0x%04x FDRA: %s; lift %.2f over %u DL DCIs (next hypothesis %.2f), DL size %u\n",
             o.rnti, describe(gs_fdra_h[best]).c_str(), b1, e.fdra_cnt[best], b2, o.size);
      continue;
    }
    /* not shown by its own DCIs (yet): the cell's, for its DCI size */
    auto bs = fdra_by_size.find(o.size);
    if (bs != fdra_by_size.end() && e.fdra_n >= 50) {
      int16_t  h = -1;
      uint32_t c = 0, tot = 0;
      for (const auto& kv : bs->second) {
        tot += kv.second;
        if (kv.second > c) {
          c = kv.second;
          h = kv.first;
        }
      }
      if (h >= 0 && c >= 2 && c * 4 >= tot * 3) {
        e.fdra_lock      = h;
        e.fdra_lock_size = o.size;
        printf("GoldenSniffer: UE 0x%04x FDRA: %s, taken from %u UE(s) with DL size %u (its own DCIs: best lift %.2f)\n",
               o.rnti, describe(gs_fdra_h[h]).c_str(), c, o.size, b1);
        continue;
      }
    }
    if (e.fdra_n >= 600 && !e.fdra_gave_up) {
      e.fdra_gave_up = true;
      printf("GoldenSniffer: UE 0x%04x FDRA not found over %u DL DCIs (best lift %.2f, next %.2f), DL size %u\n",
             o.rnti, e.fdra_n, b1, b2, o.size);
    }
  }
  /* PDSCH DM-RS symbols of tracked UEs: the most symbols a fair share of their grants
    showed (fewer in a special slot or a short PDSCH, never more by chance) */
  for (const GsDmrsObs& o : result.gs_dmrs_obs) {
    auto it = ues.find(o.rnti);
    if (it == ues.end() || it->second.dmrs_n >= 24) {
      continue;
    }
    BlindUe& e = it->second;
    e.dmrs_votes[o.pattern & 3]++;
    if (++e.dmrs_n == 24) {
      int p = 0;
      for (int q = 3; q > 0 && p == 0; q--) {
        p = e.dmrs_votes[q] >= 6 ? q : 0;
      }
      e.dmrs_add = (int8_t)p;
      static const char* sym[4] = {"first only", "+11", "+7,11", "+5,8,11"};
      printf("GoldenSniffer: UE 0x%04x PDSCH DM-RS symbols: %s (dmrs-AdditionalPosition pos%d; over 24 grants %u/%u/%u/%u "
             "for pos0..3)\n",
             o.rnti, sym[p], p, e.dmrs_votes[0], e.dmrs_votes[1], e.dmrs_votes[2], e.dmrs_votes[3]);
    }
  }
  for (const GsGrant& g : result.gs_grants) {
    gs_grant_ports[g.ports & 3]++; // 0: no DM-RS seen on its PRBs
    gs_grant_cdm2 += g.cdm == 2;
    gs_grant_add[g.dmrs_add & 3]++;
    if (gs_grants_csv != nullptr) {
      std::string rs;
      for (uint32_t r = 0; r < g.nruns; r++) {
        rs += (r ? ";" : "") + std::to_string(g.runs[2 * r]) + "-" + std::to_string(g.runs[2 * r + 1]);
      }
      fprintf(gs_grants_csv, "%.6f,%u,%u,0x%04x,%u,%u,%u,%.3f,%d,%d,%u,%u,%u,%.3f,%.3f,%.3f,%u,%s\n", now_wall, g.sfn,
              g.slot, g.rnti, g.bwp_start, g.rb_start, g.nof_rb, g.coherence, g.hyp, g.sensed ? 1 : 0, g.ports, g.cdm,
              g.dmrs_add, g.occ, g.e1, g.coh1, g.sym_seen, rs.c_str());
    }
  }
  for (const gs_dci_t& d : result.tracked_dcis) {
    gs_log_one_locked(now, result.outcome.sfn, d.slot_idx, first_crb, d, "tracked");
    auto it = ues.find(d.rnti);
    if (it != ues.end()) {
      gs_learn(it->second, d, now);
      it->second.last_seen = now;
    }
  }

  /* Blind sightings: from the searchers, and the tracked decode's other C-RNTIs, which
    may discover a UE outside the discovery windows (discover_in_track) */
  std::vector<std::pair<const gs_dci_t*, bool> > sightings;
  for (const gs_dci_t& d : gs_sightings) {
    sightings.push_back({&d, false});
  }
  for (const gs_dci_t& d : result.region_dcis) {
    if (std::string(rnti_kind(d.rnti, d.coreset0, d.n_rnti)) == "C") {
      sightings.push_back({&d, true});
    }
  }
  for (const auto& sp : sightings) {
    const gs_dci_t& d           = *sp.first;
    const bool      from_region = sp.second;
    auto it = ues.find(d.rnti);
    if (gs_is_known(state, d.rnti)) {
      BlindUe& e = ues[d.rnti];
      e.from_rach = true;
      gs_learn(e, d, now);
      e.last_seen = now;
      /* DCIs in CORESET#0 (Msg4, RRCSetup, fallback) are the RACH decoder's, not the
        regular per-UE decoders': they say nothing about a stale UE-specific config */
      if (!d.coreset0) {
        e.blind_since_reg++;
      }
      const uint32_t need = e.state == BlindUe::ARMED ? 2 : a.stale_hits;
      if ((e.state == BlindUe::RRC || e.state == BlindUe::ARMED) && e.blind_since_reg >= need) {
        gs_set_state(e, d.rnti, BlindUe::STALE, "seen blind, not by the regular decoders");
      }
    } else if (it != ues.end()) {
      gs_learn(it->second, d, now);
      it->second.last_seen = now;
    } else if (gs_discovery_window(now) || (from_region && a.discover_in_track)) {
      BlindUe e;
      e.state     = BlindUe::DISCOVERED;
      e.last_seen = now;
      gs_learn(e, d, now);
      ues[d.rnti] = e;
      printf("GoldenSniffer: UE 0x%04x discovered blind (N_ID %u, %s)\n", d.rnti, d.n_id,
             d.n_rnti ? "n_RNTI" : "PCI scrambling");
    }
  }
  gs_sightings.clear();

  /* Convergence, arming timeout, expiry, keep-alive of the regular list */
  bool any_learning = false, any_converged = false;
  for (auto it = ues.begin(); it != ues.end();) {
    BlindUe&       e    = it->second;
    const uint16_t rnti = it->first;
    const bool     known = gs_is_known(state, rnti);
    if ((e.state == BlindUe::STALE || e.state == BlindUe::DISCOVERED) && e.hits >= a.converge_hits &&
        !e.interleaved && e.track.nof_sizes > 0 && (e.dirs == 3 || now - e.learn_since >= a.converge_s)) {
      char why[160];
      char sz[32];
      if (e.track.nof_sizes > 1) {
        snprintf(sz, sizeof(sz), "%u/%u", e.track.sizes[0], e.track.sizes[1]);
      } else {
        snprintf(sz, sizeof(sz), "%u", e.track.sizes[0]);
      }
      snprintf(why, sizeof(why), "learned: N_ID %u, sizes %s, CRB %d-%d, %u symbol(s)", e.track.n_id, sz,
               (int)first_crb + e.track.sc_lo / NRSCOPE_NSC_PER_RB_NR,
               (int)first_crb + e.track.sc_hi / NRSCOPE_NSC_PER_RB_NR, e.track.duration);
      gs_set_state(e, rnti, BlindUe::CONVERGED, why);
      /* test: a reconfiguration that changed one DCI size, as the track would see it */
      static bool test_dropped = false;
      if (!test_dropped && getenv("GS_TEST_DROP_SIZE") != nullptr) {
        const uint32_t drop = (uint32_t)atoi(getenv("GS_TEST_DROP_SIZE"));
        uint32_t       k    = 0;
        for (uint32_t z = 0; z < e.track.nof_sizes; z++) {
          if (e.track.sizes[z] != drop) {
            e.track.sizes[k++] = e.track.sizes[z];
          }
        }
        if (k < e.track.nof_sizes) {
          e.track.nof_sizes = k;
          test_dropped      = true;
          printf("GoldenSniffer: [test] size %u removed from the track of 0x%04x: %s\n", drop, rnti,
                 gs_sizes_str(e.track).c_str());
        }
      }
    }
    /* Re-learning of a converged UE: the shadow replaces the track once it converged */
    if (e.state == BlindUe::CONVERGED && e.shadow) {
      const BlindUe& s = *e.shadow;
      if (!s.interleaved && s.hits >= a.converge_hits && s.track.nof_sizes > 0 &&
          (s.dirs == 3 || now - s.learn_since >= a.converge_s)) {
        printf("GoldenSniffer: UE 0x%04x re-converged: sizes %s -> %s, N_ID %u -> %u, CRB %d-%d, %u symbol(s)\n", rnti,
               gs_sizes_str(e.track).c_str(), gs_sizes_str(s.track).c_str(), e.track.n_id, s.track.n_id,
               (int)first_crb + s.track.sc_lo / NRSCOPE_NSC_PER_RB_NR,
               (int)first_crb + s.track.sc_hi / NRSCOPE_NSC_PER_RB_NR, s.track.duration);
        /* the sizes of both, by count: a re-learning started by a fallback DCI (1_0/0_0 in
          the UE search space) saw few 0_1 and lost that size, which then looked like UL
          stopping and started it again */
        std::map<uint32_t, uint32_t> merged = e.size_count;
        for (const auto& z : s.size_count) {
          merged[z.first] += z.second;
        }
        e.track       = s.track;
        e.size_count  = merged;
        {
          std::vector<std::pair<uint32_t, uint32_t> > by(merged.begin(), merged.end());
          std::sort(by.begin(), by.end(), [](const auto& x, const auto& y) { return x.second > y.second; });
          e.track.nof_sizes = 0;
          for (const auto& z : by) {
            if (e.track.nof_sizes < 8) {
              e.track.sizes[e.track.nof_sizes++] = z.first;
            }
          }
        }
        e.dirs        = s.dirs;
        e.hits        = s.hits;
        e.learn_since = s.learn_since;
        e.shadow.reset();
      } else if (now - e.shadow_since > a.arm_s) {
        printf("GoldenSniffer: UE 0x%04x re-learning ended, %u DCI(s) learned: track kept (%s)\n", rnti, s.hits,
               gs_sizes_str(e.track).c_str());
        e.shadow.reset();
      }
    }
    /* One direction stopped while the other goes on: a reconfiguration may have changed
      that format's size, which the tracked decode cannot see. A sweep finds it. */
    if (e.state == BlindUe::CONVERGED && e.dirs == 3 && e.last_dir[0] > 0 && e.last_dir[1] > 0) {
      const bool ul_q = now - e.last_dir[0] > a.quiet_s, dl_q = now - e.last_dir[1] > a.quiet_s;
      if (ul_q != dl_q && !e.dir_quiet) {
        e.dir_quiet    = true;
        gs_sweep_until = std::max(gs_sweep_until, now + a.scan_window_s);
        printf("GoldenSniffer: UE 0x%04x %s DCIs stopped %.1f s ago, %s go on: sweeping for %.0f s\n", rnti,
               ul_q ? "UL" : "DL", now - e.last_dir[ul_q ? 0 : 1], ul_q ? "DL" : "UL", a.scan_window_s);
      } else if (ul_q == dl_q) {
        e.dir_quiet = false;
      }
    }
    /* A blind-tracked UE that was active and went silent: idle, gone, or back under a
      new C-RNTI (re-attach or reconfiguration with sync) whose RACH was missed. A sweep
      costs little, so start one now rather than at the next period. */
    if (e.state != BlindUe::RRC && e.last_seen > 0) {
      if (!e.quiet && now - e.last_seen > a.quiet_s) {
        e.quiet        = true;
        gs_sweep_until = std::max(gs_sweep_until, now + a.scan_window_s);
        printf("GoldenSniffer: UE 0x%04x quiet for %.1f s, sweeping for %.0f s\n", rnti, now - e.last_seen,
               a.scan_window_s);
      } else if (e.quiet && now - e.last_seen <= a.quiet_s) {
        e.quiet = false;
      }
    }
    if (e.state == BlindUe::ARMED && now - e.armed_at > a.arm_s) {
      gs_set_state(e, rnti, BlindUe::RRC, "no sign of a changed config");
    }
    /* still seen blind: not expired from the regular list either. A converged UE is kept
      there until silent for keep_s: an idle UE sends nothing for a while, and if it was
      dropped, its next DCIs would be lost until a sweep found it again. */
    const bool keep = e.state == BlindUe::CONVERGED && e.last_seen > 0 && now - e.last_seen < a.keep_s;
    if (known && e.last_seen > 0) {
      const double alive = (keep ? now : e.last_seen) + (now_wall - now);
      for (size_t i = 0; i < state.known_rntis.size(); i++) {
        if (state.known_rntis[i] == rnti && state.last_seen[i] < alive) {
          state.last_seen[i] = alive;
        }
      }
    }
    const double last   = std::max(e.last_seen, e.last_regular);
    const double expiry = e.state == BlindUe::CONVERGED ? std::max(a.keep_s, state.rnti_expiry_s) : state.rnti_expiry_s;
    if ((!known && e.from_rach) || (!known && last > 0 && now - last > expiry)) {
      if (e.state != BlindUe::RRC) {
        printf("GoldenSniffer: UE 0x%04x %s dropped (not seen for %.1f s)\n", rnti, state_name(e.state), now - last);
      }
      gs_set_known_rnti(rnti, false);
      it = ues.erase(it);
      continue;
    }
    /* known to the searchers once seen a few times (one-off false RNTIs are not) */
    gs_set_known_rnti(rnti, known || e.state == BlindUe::CONVERGED || e.hits >= 3);
    any_learning  = any_learning || e.state == BlindUe::STALE || e.state == BlindUe::DISCOVERED || e.shadow ||
                   (e.state != BlindUe::RRC && now - e.armed_at < a.arm_s);
    any_converged = any_converged || e.state == BlindUe::CONVERGED;
    ++it;
  }
  gs_pool_active.store(a.always || gs_discovery_window(now) || any_learning);
  gs_any_converged.store(any_converged);

  /* C-RNTIs to the regular DCI decoders (feed_rntis), as before */
  if (state.rach_found && !gs_feeds.empty()) {
    for (const GsFeed& f : gs_feeds) {
      if (gs_is_known(state, f.rnti)) {
        continue;
      }
      state.known_rntis.push_back(f.rnti);
      state.last_seen.push_back(now_wall);
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

  if (a.report_s > 0 && now - gs_last_report >= a.report_s) {
    gs_track_summary.clear();
    for (const auto& u : ues) {
      if (u.second.state == BlindUe::RRC) {
        continue;
      }
      char line[160];
      std::string sz;
      for (uint32_t z = 0; z < u.second.track.nof_sizes; z++) {
        sz += (z ? "/" : "") + std::to_string(u.second.track.sizes[z]);
      }
      snprintf(line, sizeof(line), "  UE 0x%04x %-10s N_ID %5u, %u DCIs learned, sizes %s\n", u.first,
               state_name(u.second.state), u.second.track.n_id, u.second.hits, sz.c_str());
      gs_track_summary += line;
    }
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
  if (track_gs != nullptr) {
    gs_destroy(track_gs);
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
  {
    /* where the carrier is (SIB1), for setting rx_center_freq on an unknown cell */
    static std::atomic<bool> told{false};
    if (!told.exchange(true)) {
      const double c_lo = f_lo, c_hi = c_lo + carrier.carrier_bw * rb_hz; // f_lo: carrier PRB 0
      printf("GoldenSniffer: carrier %u PRBs, %.3f..%.3f MHz (centre %.3f MHz); the capture passes %.3f..%.3f MHz%s",
             (unsigned)carrier.carrier_bw, c_lo / 1e6, c_hi / 1e6, (c_lo + c_hi) / 2e6, (fc - half) / 1e6,
             (fc + half) / 1e6, first > 0 || end < (int)carrier.carrier_bw ? "" : "\n");
      if (first > 0 || end < (int)carrier.carrier_bw) {
        printf(": set rx_center_freq: %.0f\n", std::round((c_lo + c_hi) / 2.0 / 1e3) * 1e3);
      }
    }
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

  dmrs_l0 = cell.mib.dmrs_typeA_pos == srsran_dmrs_sch_typeA_pos_3 ? 3 : 2;
  if (nrscope_gs_args.fdra) {
    std::call_once(gs_fdra_once, [&]() {
      const uint32_t riv = sib1.serving_cell_cfg_common.dl_cfg_common.init_dl_bwp.generic_params.location_and_bw;
      const uint32_t N = 275, ra = riv / N, rb = riv % N;
      const uint32_t bw_start = (ra + rb < N) ? rb : N - 1 - rb, bw_size = (ra + rb < N) ? ra + 1 : N - ra + 1;
      const uint32_t nh = gs_fdra_hypotheses(carrier.offset_to_carrier, carrier.carrier_bw,
                                             carrier.offset_to_carrier + bw_start, bw_size, gs_fdra_h, GS_FDRA_MAX_HYP);
      gs_fdra_nh.store(nh);
      printf("GoldenSniffer: FDRA of DCI 1_1 learned on air: %u hypotheses (carrier %u PRBs at CRB %u, initial BWP %u "
             "at CRB %u), PDSCH DM-RS symbol %u\n",
             nh, (unsigned)carrier.carrier_bw, (unsigned)carrier.offset_to_carrier, bw_size,
             (unsigned)(carrier.offset_to_carrier + bw_start), dmrs_l0);
    });
  }
  if (nrscope_gs_args.fdra && nrscope_sensing_args.enable) {
    const double grid_center = f_lo + (first + grid_carrier.nof_prb / 2.0) * rb_hz;
    sensing = nrscope_sensing_get((uint64_t)grid_center, srate, (uint32_t)std::lround(srate / scs), (uint32_t)scs, 1);
    if (sensing != nullptr) {
      sensing_scratch = nrscope_sensing_scratch_alloc(grid_carrier.nof_prb * NRSCOPE_NSC_PER_RB_NR);
    }
  }
  ref_crb0 = -(int32_t)(first_crb * NRSCOPE_NSC_PER_RB_NR);
  ref_cs0  = (int32_t)std::lround((state->coreset0_args_t.coreset0_lower_freq_hz - grid_lo) / scs);
  cs0_nof_sc = srsran_coreset_get_bw(&state->coreset0_t) * NRSCOPE_NSC_PER_RB_NR;
  pci      = state->cs_ret.ssb_res.N_id;
  slots_per_frame = SRSRAN_NSLOTS_PER_FRAME_NR(state->arg_scs.scs);
  gs_spf.store(slots_per_frame);

  gs_args_t ga;
  gs_args_default(&ga);
  ga.size_min = nrscope_gs_args.size_min;
  ga.soft     = nrscope_gs_args.soft ? 1 : 0;
  ga.size_max = nrscope_gs_args.size_max;
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
    {
      const size_t nq = std::max(1u, nrscope_gs_args.queue);
      gs_pool.jobs.resize(nq);
      gs_pool.bufs.resize(nq);
      gs_pool.free_idx.reserve(nq);
      gs_pool.ring.assign(nq, -1);
      for (size_t k = 0; k < nq; k++) {
        gs_pool.bufs[k].assign((size_t)GS_MAX_SC * (GS_MAX_DURATION + 1), cf_t{});
        gs_pool.jobs[k].grid = gs_pool.bufs[k].data();
        gs_pool.jobs[k].buf  = (int)k;
        gs_pool.free_idx.push_back((int)k);
      }
    }
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
  track_gs = gs_create(&gs_pool.args);
  if (track_gs == nullptr) {
    ERROR("GoldenSniffer: could not create the tracked-UE decoder");
    return SRSRAN_ERROR;
  }
  gs_first_crb.store(first_crb);
  double zero = 0;
  gs_ready_time.compare_exchange_strong(zero, wall_now());
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

void GoldenSnifferDecoder::Submit(srsran_slot_cfg_t* slot, WorkState* state, SlotResult* result)
{
  if (!ready) {
    return;
  }
  const bool search  = gs_pool_active.load(std::memory_order_relaxed);
  const bool tracked = gs_any_converged.load(std::memory_order_relaxed);
  if (!search && !tracked && getenv("GS_DUMP") == nullptr) {
    return; // nothing to do in this slot: not even the FFT
  }
  struct SubTimer {
    double t0 = gs_steady_s();
    ~SubTimer() { gs_t_submit.add(1e3 * (gs_steady_s() - t0)); }
  } sub_timer;
  /* Demodulate the PDCCH symbols here (cheap: a few FFTs) since the worker's buffer is
    only valid during its slot. One symbol more than searched: the duration detector
    compares it. */
  const uint32_t nsym                    = nrscope_gs_args.max_duration + 1;
  uint32_t       mask[SRSRAN_MAX_PORTS] = {(1u << nsym) - 1u};
  uint32_t       done[SRSRAN_MAX_PORTS] = {};
  const double t_fft0 = gs_steady_s();
  result->gs_t_submit = t_fft0;
  srsran_ue_dl_nr_estimate_fft_nrscope_symbols(&ue_dl, slot, arg_scs_grid, mask, done);
  gs_t_fft.add(1e3 * (gs_steady_s() - t_fft0));

  const uint32_t nof_sc = grid_carrier.nof_prb * NRSCOPE_NSC_PER_RB_NR;

  /* Converged UEs: decoded right here, from what was learned of their PDCCH */
  static thread_local std::vector<gs_track_t> tracked_now; /* keeps its capacity */
  tracked_now.clear();
  if (tracked && track_gs != nullptr) {
    gs_slot_t s  = {};
    s.grid       = ue_dl.sf_symbols[0];
    s.nof_sc     = nof_sc;
    s.slot_idx   = slot->idx % slots_per_frame;
    s.ref_crb0   = ref_crb0;
    s.ref_cs0    = ref_cs0;
    s.cs0_nof_sc = cs0_nof_sc;
    s.pci        = pci;
    for (const auto& u : state->blind_ues) {
      if (u.second.state == BlindUe::CONVERGED) {
        tracked_now.push_back(u.second.track);
      }
    }
    if (!tracked_now.empty()) {
      /* all converged UEs at once; while the search is off, the other C-RNTIs found in
        their CORESETs are sightings (the search reports them while it is on) */
      const double t_tr0 = gs_steady_s();
      gs_dci_t     out[GS_MAX_DCI_PER_SLOT];
      uint8_t      trk[GS_MAX_DCI_PER_SLOT];
      const bool   want_other = nrscope_gs_args.discover_in_track && !search;
      int          n = gs_decode_tracked_multi(track_gs, &s, tracked_now.data(), (uint32_t)tracked_now.size(),
                                      nrscope_gs_args.track_full_band, want_other, out, trk, GS_MAX_DCI_PER_SLOT);
      for (int i = 0; i < n; i++) {
        (trk[i] ? result->tracked_dcis : result->region_dcis).push_back(out[i]);
      }
      gs_t_track.add(1e3 * (gs_steady_s() - t_tr0));
      /* survey: the PDSCH DM-RS on every symbol of this slot, PCI-scrambled, and each DL
        DCI's own N_ID on the first DM-RS symbol (cell statistics, offline) */
      if (gs_survey_csv != nullptr) {
        uint32_t all[SRSRAN_MAX_PORTS] = {0x3ffe};
        srsran_ue_dl_nr_estimate_fft_nrscope_symbols(&ue_dl, slot, arg_scs_grid, all, done);
        const uint32_t           np = grid_carrier.nof_prb, sl = slot->idx % slots_per_frame;
        static thread_local float c0[SRSRAN_MAX_PRB_NR], c1[SRSRAN_MAX_PRB_NR], oc[SRSRAN_MAX_PRB_NR],
            q0[SRSRAN_MAX_PRB_NR], q1[SRSRAN_MAX_PRB_NR];
        int n_dl = 0, n_ul = 0;
        for (int i = 0; i < n; i++) {
          (out[i].nof_bits > 0 && out[i].payload[0] == 1 ? n_dl : n_ul)++;
        }
        std::string line = std::to_string(state->sfn) + "," + std::to_string(sl) + "," + std::to_string(n_dl) + "," +
                           std::to_string(n_ul);
        uint32_t cnt1[14] = {};
        std::string bitmap, l0_extra;
        for (uint32_t l = 1; l < 14; l++) {
          gs_dmrs_rb_map(ue_dl.sf_symbols[0] + l * nof_sc, first_crb, np, sl, l, pci, c0, c1, oc, q0, q1);
          uint32_t k0 = 0;
          for (uint32_t c = 0; c < np; c++) {
            k0 += c0[c] > 0.5f;
            cnt1[l] += c1[c] > 0.5f;
          }
          line += "," + std::to_string(k0);
          if (l == dmrs_l0) {
            std::vector<float> o, p;
            char               hx[8];
            for (uint32_t c = 0; c < np; c += 4) {
              uint32_t v = 0;
              for (uint32_t b = 0; b < 4 && c + b < np; b++) {
                v |= (c0[c + b] > 0.5f ? 1u : 0u) << (3 - b);
              }
              snprintf(hx, sizeof(hx), "%x", v);
              bitmap += hx;
            }
            for (uint32_t c = 0; c < np; c++) {
              if (c0[c] > 0.5f) {
                o.push_back(oc[c]);
                p.push_back(q0[c] > 0 ? q1[c] / q0[c] : 0);
              }
            }
            auto med = [](std::vector<float>& v) {
              if (v.empty()) {
                return -1.0f;
              }
              std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
              return v[v.size() / 2];
            };
            char m[64];
            snprintf(m, sizeof(m), ",%.3f,%.3f", med(o), med(p));
            l0_extra = m;
            /* strict occupancy (coherence 0.8: an empty PRB passes 0.5 a third of the time),
              type 2 test, and the precoding: phase/amplitude step between neighbouring PRBs
              that both carry DM-RS, against the step inside a PRB, by CRB class */
            static thread_local cf_t  pr[3 * SRSRAN_MAX_PRB_NR];
            static thread_local float t2[SRSRAN_MAX_PRB_NR];
            gs_dmrs_rb_pairs(ue_dl.sf_symbols[0] + l * nof_sc, first_crb, np, sl, l, pci, pr, t2);
            std::string bm8;
            uint32_t    n_t2 = 0;
            for (uint32_t c = 0; c < np; c += 4) {
              uint32_t v = 0;
              for (uint32_t b = 0; b < 4 && c + b < np; b++) {
                v |= (c0[c + b] > 0.8f ? 1u : 0u) << (3 - b);
              }
              snprintf(hx, sizeof(hx), "%x", v);
              bm8 += hx;
            }
            cf_t rho = 0;
            for (uint32_t c = 0; c < np; c++) {
              n_t2 += t2[c] > 0.8f;
              if (c0[c] > 0.8f) {
                rho += pr[3 * c + 1] * gs_cj(pr[3 * c]) + pr[3 * c + 2] * gs_cj(pr[3 * c + 1]);
              }
            }
            const float phi = gs_arg(rho);
            auto        dev = [phi](cf_t a, cf_t b) -> float { // phase of b against a, the slope removed, |deg|
              float x = gs_arg(b * gs_cj(a)) - phi;
              x       = x > (float)M_PI ? x - 2 * (float)M_PI : x < -(float)M_PI ? x + 2 * (float)M_PI : x;
              return fabsf(x) * 57.29578f;
            };
            auto adb = [](cf_t a, cf_t b) { return fabsf(10.0f * log10f((gs_re(b * gs_cj(b)) + 1e-20f) / (gs_re(a * gs_cj(a)) + 1e-20f))); };
            double   sd[4] = {}, sa[4] = {}; // 0: inside a PRB, 1: CRB%4==0, 2: CRB%4==2, 3: odd CRB
            uint32_t nd[4] = {};
            for (uint32_t c = 0; c < np; c++) {
              if (c0[c] <= 0.8f) {
                continue;
              }
              for (uint32_t n = 0; n < 2; n++) {
                sd[0] += dev(pr[3 * c + n], pr[3 * c + n + 1]);
                sa[0] += adb(pr[3 * c + n], pr[3 * c + n + 1]);
                nd[0]++;
              }
              if (c + 1 < np && c0[c + 1] > 0.8f) {
                const uint32_t crb = first_crb + c + 1, k = crb % 4 == 0 ? 1 : crb % 2 == 0 ? 2 : 3;
                sd[k] += dev(pr[3 * c + 2], pr[3 * c + 3]);
                sa[k] += adb(pr[3 * c + 2], pr[3 * c + 3]);
                nd[k]++;
              }
            }
            char e[256];
            snprintf(e, sizeof(e), ",%s,%u", bm8.c_str(), n_t2);
            l0_extra += e;
            for (int k = 0; k < 4; k++) {
              snprintf(e, sizeof(e), ",%u:%.1f:%.2f", nd[k], nd[k] ? sd[k] / nd[k] : 0, nd[k] ? sa[k] / nd[k] : 0);
              l0_extra += e;
            }
          }
        }
        for (uint32_t l = 1; l < 14; l++) {
          line += "," + std::to_string(cnt1[l]);
        }
        /* CSI-RS: every symbol and subcarrier, the mappings over 0.3 coherence (l:c:mapping:coh),
          and the REs left empty inside the PDSCH (z:l:c:power share) */
        std::string csi;
        for (uint32_t l = 1; l < 14; l++) {
          float co[12][4], pz[12];
          gs_csirs_scan(ue_dl.sf_symbols[0] + l * nof_sc, first_crb, np, sl, l, pci, co, pz);
          char e[48];
          for (uint32_t c = 0; c < 12; c++) {
            for (uint32_t h = 0; h < 4; h++) {
              if (co[c][h] > 0.3f) {
                snprintf(e, sizeof(e), "%s%u:%u:%u:%.2f", csi.empty() ? "" : ";", l, c, h, co[c][h]);
                csi += e;
              }
            }
            if (pz[c] < 0.25f) {
              snprintf(e, sizeof(e), "%sz:%u:%u:%.2f", csi.empty() ? "" : ";", l, c, pz[c]);
              csi += e;
            }
          }
        }
        line += "," + bitmap + l0_extra + "," + csi + ",";
        for (int i = 0; i < n; i++) { // every DL DCI's own N_ID on the first DM-RS symbol
          if (out[i].nof_bits == 0 || out[i].payload[0] != 1 || out[i].n_id == pci) {
            continue;
          }
          gs_dmrs_rb_map(ue_dl.sf_symbols[0] + dmrs_l0 * nof_sc, first_crb, np, sl, dmrs_l0, out[i].n_id, c0, c1, oc,
                         q0, q1);
          int k = 0, lo = -1, hi = -1;
          for (uint32_t c = 0; c < np; c++) {
            if (c0[c] > 0.5f) {
              k++;
              lo = lo < 0 ? (int)c : lo;
              hi = (int)c;
            }
          }
          char e[64];
          snprintf(e, sizeof(e), "%x:%u:%d:%d-%d;", out[i].rnti, out[i].n_id, k, lo, hi);
          line += e;
        }
        std::lock_guard<std::mutex> lk(gs_survey_mutex);
        fprintf(gs_survey_csv, "%s\n", line.c_str());
      }
      /* FDRA and grants of the tracked UEs' DL DCIs (format identifier 1) */
      const uint32_t nh = gs_fdra_nh.load(std::memory_order_relaxed);
      if (nrscope_gs_args.fdra && nh > 0) {
        bool                            have_map = false;
        static thread_local cf_t        rb_s[SRSRAN_MAX_PRB_NR];
        static thread_local float       rb_d[SRSRAN_MAX_PRB_NR];
        static thread_local cf_t        ue_s[SRSRAN_MAX_PRB_NR]; // the same with a UE's own N_ID
        static thread_local float       ue_d[SRSRAN_MAX_PRB_NR];
        static thread_local bool        occ_p[SRSRAN_MAX_PRB_NR], occ_u[SRSRAN_MAX_PRB_NR]; // strict maps
        float                           f_p = 0, f_u = 0;
        uint32_t                        ue_map_nid = UINT32_MAX;
        const uint32_t                  nprb = grid_carrier.nof_prb;
        for (int i = 0; i < n; i++) {
          const gs_dci_t& dd = out[i];
          if (!trk[i] || dd.nof_bits < 2 || dd.payload[0] != 1) {
            continue;
          }
          auto ue = state->blind_ues.find(dd.rnti);
          if (ue == state->blind_ues.end()) {
            continue;
          }
          const BlindUe& b = ue->second;
          if (b.fdra_lock < 0 && b.fdra_n >= 2000) { // gave up long ago: not even a cell hypothesis
            continue;
          }
          if (!have_map) {
            /* the first PDSCH DM-RS symbol of the slot, demodulated if it was not */
            uint32_t mask2[SRSRAN_MAX_PORTS] = {1u << dmrs_l0};
            srsran_ue_dl_nr_estimate_fft_nrscope_symbols(&ue_dl, slot, arg_scs_grid, mask2, done);
            gs_fdra_rb_coherence(ue_dl.sf_symbols[0] + dmrs_l0 * nof_sc, first_crb, nprb,
                                 slot->idx % slots_per_frame, dmrs_l0, pci, rb_s, rb_d);
            f_p      = gs_fdra_occ_map(rb_s, rb_d, nprb, occ_p);
            have_map = true;
          }
          /* and the map under this UE's own PDCCH DM-RS N_ID, when it has one */
          const bool own_nid = dd.n_id != pci;
          if (own_nid && ue_map_nid != dd.n_id) {
            gs_fdra_rb_coherence(ue_dl.sf_symbols[0] + dmrs_l0 * nof_sc, first_crb, nprb,
                                 slot->idx % slots_per_frame, dmrs_l0, dd.n_id, ue_s, ue_d);
            f_u        = gs_fdra_occ_map(ue_s, ue_d, nprb, occ_u);
            ue_map_nid = dd.n_id;
          }
          /* the layout for this DCI's size: the cell's for that size, else the UE's own lock
            when it was learned at this size; another size with neither is learned (or, for a
            UE already locked at another size, skipped) */
          int16_t hsel = dd.nof_bits < 256 ? gs_fdra_size_hyp[dd.nof_bits].load(std::memory_order_relaxed) : -1;
          if (hsel < 0 && b.fdra_lock >= 0 && b.fdra_lock_size == dd.nof_bits) {
            hsel = b.fdra_lock;
          }
          if (hsel < 0 && b.fdra_lock >= 0) {
            continue;
          }
          if (hsel < 0) {
            GsFdraObs o;
            o.rnti = dd.rnti;
            o.nof  = (uint8_t)nh;
            o.size = (uint16_t)dd.nof_bits;
            for (uint32_t k = 0; k < nh; k++) {
              const bool u = gs_fdra_h[k].nid_ue != 0;
              o.score[k]   = u && !own_nid ? -2.0f // the same as its PCI twin
                                           : gs_fdra_lift(gs_fdra_h[k], dd.payload, dd.nof_bits, first_crb, nprb,
                                                          u ? occ_u : occ_p, u ? f_u : f_p);
            }
            result->gs_fdra_obs.push_back(o);
            continue;
          }
          const GsFdraHyp& hy = gs_fdra_h[hsel];
          bool             prb[SRSRAN_MAX_PRB_NR];
          uint32_t         nof = 0;
          if (!gs_fdra_parse(hy, dd.payload, dd.nof_bits, prb, &nof) || nof == 0) {
            continue;
          }
          float          coh    = 0;
          const bool     use_ue = hy.nid_ue != 0 && own_nid;
          const uint32_t dm_nid = use_ue ? dd.n_id : pci; // the PDSCH DM-RS N_ID of this UE
          gs_fdra_score(hy, dd.payload, dd.nof_bits, first_crb, nprb, use_ue ? ue_s : rb_s, use_ue ? ue_d : rb_d, &coh);
          GsGrant g = {};
          g.rnti      = dd.rnti;
          g.sfn       = state->sfn;
          g.slot      = slot->idx % slots_per_frame;
          g.bwp_start = hy.bwp_start;
          g.hyp       = hsel;
          g.coherence = coh;
          uint32_t first_rb = 0;
          while (first_rb < hy.nbwp && !prb[first_rb]) {
            first_rb++;
          }
          g.rb_start = (uint16_t)first_rb;
          g.nof_rb   = (uint16_t)nof;
          for (uint32_t k = 0; k < hy.nbwp && k < SRSRAN_MAX_PRB_NR && g.nruns < 18; k++) {
            if (prb[k] && (k == 0 || !prb[k - 1])) {
              g.runs[2 * g.nruns] = (uint16_t)k;
            }
            if (prb[k] && (k + 1 >= hy.nbwp || !prb[k + 1])) {
              g.runs[2 * g.nruns + 1] = (uint16_t)k;
              g.nruns++;
            }
          }
          /* sensing: the grant as a single-layer PDSCH on port 1000 (its antenna-ports field
            is not known), the DM-RS symbol from the MIB, scrambled with the PCI */
          /* The DM-RS read on air for every grant with DM-RS, sensing or not (cell statistics) */
          if (coh > 0.5f) {
            /* What its antenna-ports field and its RRC would have said, read on air: the
              ports on comb 0, whether comb 1 is left empty (CDM groups without data 2,
              DM-RS 3 dB up), and the DM-RS symbols past the first one */
            static thread_local bool gprb[SRSRAN_MAX_PRB_NR];
            memset(gprb, 0, sizeof(gprb));
            for (uint32_t k = 0; k < hy.nbwp && k < SRSRAN_MAX_PRB_NR; k++) {
              const int32_t c = (int32_t)(hy.bwp_start + k) - (int32_t)first_crb;
              if (prb[k] && c >= 0 && c < (int32_t)nprb) {
                gprb[c] = true;
              }
            }
            GsDmrsShape sh;
            gs_dmrs_shape(ue_dl.sf_symbols[0] + dmrs_l0 * nof_sc, first_crb, gprb, nprb, g.slot, dmrs_l0, dm_nid, &sh);
            uint16_t ports     = 0x1;
            bool     cdm2      = false;
            int      dci_ports = 0;
            gs_dmrs_decide(sh, &ports, &cdm2, &dci_ports);
            int        add      = 0;
            const bool     learning  = b.dmrs_n < 24;
            if (learning || b.dmrs_add > 0) {
              // while learning every candidate symbol, then only the UE's own pattern
              const uint32_t cand = learning ? (1u << 5) | (1u << 7) | (1u << 8) | (1u << 11)
                                    : b.dmrs_add == 3 ? (1u << 5) | (1u << 8) | (1u << 11)
                                    : b.dmrs_add == 2 ? (1u << 7) | (1u << 11)
                                                      : (1u << 11);
              uint32_t mask3[SRSRAN_MAX_PORTS] = {cand};
              srsran_ue_dl_nr_estimate_fft_nrscope_symbols(&ue_dl, slot, arg_scs_grid, mask3, done);
              bool on[14] = {};
              for (uint32_t l : {5u, 7u, 8u, 11u}) {
                if (!(cand & (1u << l))) {
                  continue;
                }
                GsDmrsShape sl;
                gs_dmrs_shape(ue_dl.sf_symbols[0] + l * nof_sc, first_crb, gprb, nprb, g.slot, l, dm_nid, &sl);
                on[l] = sl.coh0 > 0.5f;
              }
              const int seen = (on[5] && on[8] && on[11]) ? 3 : (on[7] && on[11]) ? 2 : on[11] ? 1 : 0;
              if (learning) {
                result->gs_dmrs_obs.push_back(GsDmrsObs{dd.rnti, (uint8_t)seen});
                add = seen;
              } else { // the UE's pattern, when this slot's PDSCH is long enough to hold it
                const bool ok = b.dmrs_add == 3 ? (on[5] && on[8] && on[11]) : b.dmrs_add == 2 ? (on[7] && on[11]) : on[11];
                add           = ok ? b.dmrs_add : 0;
              }
            }
            g.ports    = (uint8_t)ports;
            g.cdm      = cdm2 ? 2 : 1;
            g.dmrs_add = (uint8_t)add;
            g.occ      = sh.occ;
            g.e1       = sh.e1;
            g.coh1     = sh.coh1;
            g.sym_seen = (uint16_t)(1u << dmrs_l0);
            if (learning || b.dmrs_add > 0) {
              for (uint32_t l : {5u, 7u, 8u, 11u}) {
                GsDmrsShape sl;
                if (done[0] & (1u << l)) {
                  gs_dmrs_shape(ue_dl.sf_symbols[0] + l * nof_sc, first_crb, gprb, nprb, g.slot, l, dm_nid, &sl);
                  g.sym_seen |= sl.coh0 > 0.5f ? (uint16_t)(1u << l) : 0;
                }
              }
            }
            srsran_sch_cfg_nr_t cfg = {};
            for (uint32_t k = 0; k < hy.nbwp && k < SRSRAN_MAX_PRB_NR; k++) {
              cfg.grant.prb_idx[k] = prb[k];
            }
            cfg.grant.nof_prb                          = nof;
            cfg.grant.S                                = std::max(dd.duration, 1u);
            cfg.grant.L                                = 14 - cfg.grant.S;
            cfg.grant.mapping                          = srsran_sch_mapping_type_A;
            cfg.grant.k                                = 0;
            cfg.grant.rnti                             = dd.rnti;
            cfg.grant.rnti_type                        = srsran_rnti_type_c;
            cfg.grant.dci_format                       = srsran_dci_format_nr_1_1;
            cfg.grant.nof_layers                       = 2; // despread as 1000+1001 (gs_dmrs_decide)
            cfg.grant.nof_dmrs_cdm_groups_without_data = cdm2 ? 2 : 1;
            cfg.grant.beta_dmrs                        = cdm2 ? (float)M_SQRT2 : 1.0f;
            cfg.dmrs.type                              = srsran_dmrs_sch_type_1;
            cfg.dmrs.typeA_pos      = dmrs_l0 == 3 ? srsran_dmrs_sch_typeA_pos_3 : srsran_dmrs_sch_typeA_pos_2;
            cfg.dmrs.additional_pos = add == 3   ? srsran_dmrs_sch_add_pos_3
                                      : add == 2 ? srsran_dmrs_sch_add_pos_2
                                      : add == 1 ? srsran_dmrs_sch_add_pos_1
                                                 : srsran_dmrs_sch_add_pos_0;
            cfg.dmrs.length         = srsran_dmrs_sch_len_1;
            cfg.dmrs.scrambling_id0_present = use_ue;
            cfg.dmrs.scrambling_id0         = dm_nid;
            const nr_dmrs_placement_t place = {first_crb, hy.bwp_start, 0};
            if (sensing != nullptr && sensing_scratch != nullptr) {
              g.sensed = nrscope_sensing_process_grant_layers(sensing, sensing_scratch, sensing_join, ports, 0,
                                                         ue_dl.sf_symbols[0], nprb * NRSCOPE_NSC_PER_RB_NR, &place,
                                                         &cfg, dci_ports, pci, state->sfn, g.slot) > 0;
            }
          }
          result->gs_grants.push_back(g);
        }
      }
    }
  }

  GsJob job;
  job.sfn        = state->sfn;
  job.slot_idx   = slot->idx % slots_per_frame;
  job.nof_sc     = nof_sc;
  job.nof_symb   = nsym;
  job.first_crb  = first_crb;
  job.pci        = pci;
  job.ref_crb0   = ref_crb0;
  job.ref_cs0    = ref_cs0;
  job.cs0_nof_sc = cs0_nof_sc;
  job.nof_tracked = 0;
  for (const gs_track_t& t : tracked_now) {
    if (job.nof_tracked < GS_JOB_TRACKS) {
      GsJobTrack& jt = job.tracked[job.nof_tracked++];
      jt.rnti        = t.rnti;
      jt.nof_sizes   = (uint8_t)std::min<uint32_t>(t.nof_sizes, 8);
      for (uint32_t z = 0; z < jt.nof_sizes; z++) {
        jt.sizes[z] = (uint8_t)t.sizes[z];
      }
    }
  }

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
  if (!search) {
    return;
  }

  /* A free buffer: replaying a file nothing is dropped (wait for one); live, never wait
    nor block: the lock taken or no buffer free, this slot is not searched. */
  int idx = -1;
  while (idx < 0) {
    std::unique_lock<GsPiMutex> lock(gs_pool.m, std::defer_lock);
    const bool                  lossless = nrscope_lossless.load(std::memory_order_relaxed);
    if (lossless) {
      lock.lock();
    } else if (!lock.try_lock()) {
      gs_pool.busy.fetch_add(1, std::memory_order_relaxed);
      return;
    }
    if (!gs_pool.free_idx.empty()) {
      idx = gs_pool.free_idx.back();
      gs_pool.free_idx.pop_back();
      gs_pool.submitted++;
      break;
    }
    if (!lossless || nrscope_stop.load(std::memory_order_relaxed)) {
      gs_pool.submitted++;
      gs_pool.dropped++;
      return;
    }
    lock.unlock();
    std::this_thread::sleep_for(std::chrono::microseconds(100));
  }
  const double t_enq0 = gs_steady_s();
  GsJob&       dst    = gs_pool.jobs[idx];
  cf_t* const  buf    = dst.grid;
  dst                 = job;
  dst.grid            = buf;
  dst.buf             = idx;
  dst.t_submit        = t_fft0;
  memcpy(dst.grid, ue_dl.sf_symbols[0], sizeof(cf_t) * (size_t)nof_sc * nsym);
  dst.t_enq = gs_steady_s();
  {
    /* the buffer is ours: wait for the lock rather than lose it (held only for a few
      instructions, with priority inheritance) */
    std::lock_guard<GsPiMutex> lock(gs_pool.m);
    gs_pool.ring[(gs_pool.r_head + gs_pool.r_len) % gs_pool.ring.size()] = idx;
    gs_pool.r_len++;
  }
  gs_t_enq.add(1e3 * (gs_steady_s() - t_enq0));
  gs_pool.cv.notify_one();
}
