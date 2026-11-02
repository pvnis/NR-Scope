/*
 * gs_cellsim: synthetic 5G NR downlink with ground truth, for testing GoldenSniffer end to end.
 *
 * Writes the complex float32 IQ (fc32) that NR-Scope's file radio reads (device_name: file),
 * a ground-truth CSV of every DCI placed, an event log, and a matching NR-Scope YAML config.
 * The cell is a real-looking SA cell: SSB (PSS/SSS/PBCH with a MIB), CORESET#0 + Type0 CSS
 * with the SIB1 DCI (SI-RNTI) and its PDSCH, a SIB1 built like srsRAN's gNB builds it (the
 * carrier, the initial BWP, RACH, TDD pattern), and UE DCIs in a CORESET 1:
 *   - non-interleaved or interleaved (REG bundle L, interleaver R, shift index),
 *     duration 1-3, any 6-RB aligned frequency range;
 *   - DM-RS (and PDCCH data) scrambled with the PCI, or with pdcch-DMRS-ScramblingID and
 *     n_RNTI = C-RNTI (UE search space with the scrambling ID configured);
 *   - 1..N C-RNTIs, each with its own DCI 1_1 / 0_1 sizes, fallback 1_0 / 0_0 (in the USS on
 *     CORESET 1, or in a common search space on CORESET#0), an AL distribution, an activity
 *     pattern with idle periods; CCE blocking respected (search space hashing, or free
 *     placement); the first payload bit is the DL/UL format identifier, the rest random.
 *   - events: DCI size changes (RRC reconfiguration), RNTI churn (a UE leaves, a new RNTI
 *     appears after a gap), at given times or at random rates.
 * Channel: AWGN (SNR per RE), CFO, integer-sample multipath, a timing offset (samples
 * dropped at the start), a capture centre offset (offset tuning or a partial capture: only
 * the subcarriers inside the capture's pass band are generated, an ideal filter).
 *
 * The PHY is srsRAN's own transmitter chain (PDCCH polar encoder and mapper, PDCCH DM-RS,
 * PDSCH with its DM-RS, the SSB/PBCH encoder, the DCI packer for SI-RNTI), so what NR-Scope
 * decodes is bit-exact with what a srsRAN gNB transmits. The OFDM modulator is written here,
 * with the TS 38.211 CP layout and per-symbol phase (the long CP every 0.5 ms; srsRAN's own
 * modulator uses the LTE layout, which is wrong at 30 kHz).
 *
 * Usage: gs_cellsim [-c file.ini] [--key value | --key=value] ...
 *   gs_cellsim --help          lists the keys and their defaults
 * The INI file is "key = value" lines (# comments, [sections] ignored); command-line keys
 * override it. Output: out (fc32 file, or a FIFO NR-Scope reads while it is generated),
 * truth (CSV), truth + ".events.csv", yaml (NR-Scope config for this cell).
 */
#include "srsran/asn1/rrc_nr.h"
#include "srsran/common/band_helper.h"
#include "srsran/srsran.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace asn1::rrc_nr;
using cpx = std::complex<float>; // our buffers; srsRAN takes them as cf_t*

/* ------------------------------------------------------------------------------------------ */
/* configuration                                                                               */
/* ------------------------------------------------------------------------------------------ */

struct KeyDoc {
  const char* key;
  const char* def;
  const char* doc;
};

static const KeyDoc keys_doc[] = {
    {"duration_s", "2", "length of the signal (s)"},
    {"sfn0", "0", "SFN of the first frame"},
    {"seed", "1", "random seed (traffic, payloads, channel)"},
    {"out", "cellsim.fc32", "IQ output, interleaved float32 (a file or a FIFO)"},
    {"truth", "truth.csv", "ground-truth CSV of every DCI placed"},
    {"yaml", "", "NR-Scope config written for this cell (empty: <out>.yaml)"},
    {"gs_csv", "golden_sniffer.csv", "golden_sniffer.csv path written into the YAML"},
    {"gs_mode", "always", "golden_sniffer mode written into the YAML (always|fallback)"},
    {"nof_workers", "2", "NR-Scope workers written into the YAML"},
    {"scs_khz", "15", "subcarrier spacing of SSB, CORESET#0 and the carrier (15|30)"},
    {"nof_prb", "106", "carrier bandwidth in PRB (106 = 20 MHz @15 kHz, 51 = 20 MHz @30 kHz, 106 = 40 MHz @30, 273 = 100 MHz @30)"},
    {"srate_hz", "0", "sampling rate (0: smallest standard FFT size that holds the carrier)"},
    {"freq_mhz", "0", "approximate carrier centre (0: 1980 for 15 kHz, 3500 for 30 kHz); moved so that the SSB is on the sync raster"},
    {"rx_offset_hz", "0", "capture centre minus carrier centre (multiple of the SCS); a capture narrower than the carrier keeps the subcarriers in its pass band"},
    {"passband", "0.92", "fraction of the sampling rate kept when the capture is narrower than the carrier"},
    {"band", "0", "NR band (0: from the frequency)"},
    {"pci", "1", "physical cell id"},
    {"tdd_period_ms", "5", "TDD pattern period (ms) when the band is TDD"},
    {"tdd_dl_slots", "7", "full DL slots per TDD period"},
    {"tdd_dl_symbols", "6", "DL symbols of the special slot"},
    {"tdd_ul_slots", "2", "full UL slots per TDD period"},
    {"tdd_ul_symbols", "4", "UL symbols of the special slot"},
    {"coreset0_idx", "6", "CORESET#0 index (MIB controlResourceSetZero, TS 38.213 tables 13-1..13-4)"},
    {"ss0_idx", "0", "search space #0 index (MIB searchSpaceZero, table 13-11, 0..9)"},
    {"ssb_rb", "-1", "CRB (carrier SCS) where the SSB starts (-1: centred)"},
    {"k_ssb", "0", "SSB subcarrier offset k_SSB (15 kHz units; even at 30 kHz)"},
    {"ssb_period_ms", "20", "SSB periodicity (10, 20, 40, 80, 160)"},
    {"dmrs_typeA_pos", "0", "MIB dmrs-TypeA-Position (2|3; 0: 3 if CORESET#0 has 3 symbols, else 2)"},
    {"sib1_period_frames", "2", "SIB1 sent every this many frames (in the Type0 occasions)"},
    {"sib1_al", "4", "aggregation level of the SIB1 DCI"},
    {"sib1_mcs", "2", "SIB1 PDSCH MCS (raised if the SIB1 does not fit)"},
    {"cs1_rb_start", "0", "CORESET 1 first CRB (multiple of 6)"},
    {"cs1_nof_rb", "0", "CORESET 1 width in RB (multiple of 6; 0: the whole carrier)"},
    {"cs1_duration", "1", "CORESET 1 symbols (1..3)"},
    {"cs1_interleaved", "0", "CORESET 1 CCE-to-REG mapping interleaved (1) or not (0)"},
    {"cs1_reg_bundle", "6", "interleaved REG bundle size L (2|3|6)"},
    {"cs1_interleaver_size", "2", "interleaver size R (2|3|6)"},
    {"cs1_shift", "-1", "interleaver shift index (-1: the PCI)"},
    {"cs1_scrambling_id", "-1", "pdcch-DMRS-ScramblingID (-1: not configured, DM-RS on the PCI, n_RNTI = 0)"},
    {"ss_candidates", "2,2,2,2,1", "UE search space candidates at AL 1,2,4,8,16"},
    {"css_candidates", "0,0,2,2,1", "common search space (CORESET#0) candidates at AL 1,2,4,8,16"},
    {"placement", "hash", "hash: TS 38.213 search space hashing; free: any free CCE-aligned position"},
    {"nof_ues", "1", "number of C-RNTIs connected at the start"},
    {"rnti_min", "0x1000", "C-RNTIs are drawn in [rnti_min, rnti_max]"},
    {"rnti_max", "0xffef", ""},
    {"rnti_list", "", "explicit initial C-RNTIs, comma separated (overrides the draw)"},
    {"dl_size_min", "40", "DCI 1_1 size range (bits) per UE"},
    {"dl_size_max", "75", ""},
    {"ul_size_min", "40", "DCI 0_1 size range (bits) per UE"},
    {"ul_size_max", "75", ""},
    {"fallback_prob", "0.05", "fraction of DCIs sent as 1_0 / 0_0"},
    {"css_fallback_prob", "0.5", "fraction of the fallback DCIs sent in the CORESET#0 common search space (else the USS on CORESET 1)"},
    {"al_weights", "1:1,2:2,4:2,8:1,16:0.5", "AL distribution"},
    {"al_sticky", "0.7", "probability a DCI uses the UE's own AL rather than a fresh draw"},
    {"p_dl_min", "0.2", "per-UE probability of a DL DCI in a DL slot, drawn in [min, max]"},
    {"p_dl_max", "0.5", ""},
    {"p_ul_min", "0.1", "per-UE probability of an UL DCI in a DL slot"},
    {"p_ul_max", "0.3", ""},
    {"idle_rate", "0", "rate (1/s) at which an active UE goes idle"},
    {"idle_min_s", "1", "idle period range (s)"},
    {"idle_max_s", "3", ""},
    {"max_dl_per_slot", "8", "DL DCIs per slot at most (the rest wait)"},
    {"max_ul_per_slot", "8", "UL DCIs per slot at most"},
    {"resize_at", "", "size changes: t_s:ue_index,... (the UE gets new DCI 1_1/0_1 sizes)"},
    {"resize_rate", "0", "random size changes per second (cell-wide)"},
    {"churn_at", "", "re-attaches: t_s:ue_index,... (the RNTI leaves, a new one appears after reattach_gap_s)"},
    {"churn_rate", "0", "random re-attaches per second (cell-wide)"},
    {"reattach_gap_s", "0.2", "gap between the old RNTI's last DCI and the new RNTI's first"},
    {"churn_new_sizes", "0", "a re-attached UE gets new DCI sizes (1) or keeps them (0)"},
    {"snr_db", "30", "SNR per RE (signal RE power over noise power per subcarrier); inf: no noise"},
    {"cfo_hz", "0", "carrier frequency offset"},
    {"taps", "", "multipath: delay_samples:gain_db,... (normalised to unit power; random phases)"},
    {"drop_samples", "0", "samples dropped at the start (timing offset: the file starts mid-slot)"},
    {"gain_db", "0", "output scaling"},
    {nullptr, nullptr, nullptr}};

static std::map<std::string, std::string> kv;

static std::string trim(const std::string& s)
{
  size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
  return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

static bool known_key(const std::string& k)
{
  for (const KeyDoc* d = keys_doc; d->key; d++) {
    if (k == d->key) {
      return true;
    }
  }
  return false;
}

static void set_key(const std::string& k, const std::string& v)
{
  if (!known_key(k)) {
    fprintf(stderr, "unknown key '%s' (--help lists them)\n", k.c_str());
    exit(1);
  }
  kv[k] = v;
}

static void load_ini(const char* path)
{
  std::ifstream f(path);
  if (!f) {
    fprintf(stderr, "cannot read %s\n", path);
    exit(1);
  }
  std::string line;
  while (std::getline(f, line)) {
    size_t h = line.find('#');
    if (h != std::string::npos) {
      line = line.substr(0, h);
    }
    line = trim(line);
    if (line.empty() || line[0] == '[') {
      continue;
    }
    size_t eq = line.find('=');
    if (eq == std::string::npos) {
      fprintf(stderr, "%s: bad line '%s'\n", path, line.c_str());
      exit(1);
    }
    set_key(trim(line.substr(0, eq)), trim(line.substr(eq + 1)));
  }
}

static std::string sget(const char* k)
{
  auto it = kv.find(k);
  if (it != kv.end()) {
    return it->second;
  }
  for (const KeyDoc* d = keys_doc; d->key; d++) {
    if (!strcmp(d->key, k)) {
      return d->def;
    }
  }
  fprintf(stderr, "internal: no key %s\n", k);
  exit(1);
}
static double dget(const char* k)
{
  std::string s = sget(k);
  if (s == "inf") {
    return INFINITY;
  }
  return strtod(s.c_str(), nullptr);
}
static long iget(const char* k)
{
  return strtol(sget(k).c_str(), nullptr, 0);
}

static std::vector<std::string> split(const std::string& s, char c)
{
  std::vector<std::string> out;
  std::stringstream        ss(s);
  std::string              item;
  while (std::getline(ss, item, c)) {
    item = trim(item);
    if (!item.empty()) {
      out.push_back(item);
    }
  }
  return out;
}

[[noreturn]] static void die(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  fprintf(stderr, "gs_cellsim: ");
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, "\n");
  va_end(ap);
  exit(1);
}

/* ------------------------------------------------------------------------------------------ */
/* cell                                                                                        */
/* ------------------------------------------------------------------------------------------ */

struct Cell {
  uint32_t                    scs_hz, mu, spf, nof_prb, nfft, pci, band;
  srsran_subcarrier_spacing_t scs;
  double                      srate, fc, fc_rx, pointA, ssb_freq;
  int32_t                     rx_shift_sc; // (fc - fc_rx) / scs
  srsran_duplex_mode_t        duplex;
  srsran_ssb_pattern_t        ssb_pattern;
  uint32_t                    offset_to_point_a, k_ssb; // 15 kHz RBs, 15 kHz subcarriers
  uint32_t                    ssb_sc0;                  // carrier subcarrier of the SSB's first subcarrier
  uint32_t                    ssb_rb_first, ssb_rb_last, ssb_l0;
  uint32_t                    cs0_idx, ss0_idx, cs0_crb, cs0_bw, cs0_dur, n0, sfn_c;
  srsran_dmrs_sch_typeA_pos_t typeA_pos;
  srsran_coreset_t            cs0; // offset_rb = 0: on its own grid
  srsran_coreset_t            cs1;
  uint32_t                    cs1_crb[SRSRAN_MAX_PRB_NR]; // CORESET-relative RB -> CRB
  uint32_t                    cs1_bw;
  // TDD
  bool     tdd;
  uint32_t tdd_period_slots, tdd_dl_slots, tdd_dl_symbols, tdd_ul_slots, tdd_ul_symbols;
  double   tdd_period_ms;
};

/* PDCCH symbols available in slot n (frame-relative): 14 DL, k special, 0 UL */
static uint32_t dl_symbols(const Cell& c, uint32_t slot)
{
  if (!c.tdd) {
    return 14;
  }
  uint32_t s = slot % c.tdd_period_slots;
  if (s < c.tdd_dl_slots) {
    return 14;
  }
  if (s == c.tdd_dl_slots) {
    return c.tdd_dl_symbols;
  }
  return 0;
}

/* Table 13-11 (FR1, multiplexing pattern 1): O, M for indexes 0..9 (first symbol 0 for i = 0) */
static void ss0_table(uint32_t idx, double* O, double* M)
{
  static const double t[10][2] = {{0, 1}, {0, 0.5}, {2, 1}, {2, 0.5}, {5, 1}, {5, 0.5}, {7, 1}, {7, 0.5}, {0, 2}, {5, 2}};
  if (idx > 9) {
    die("ss0_idx %u: only 0..9 (first symbol 0) are supported", idx);
  }
  *O = t[idx][0];
  *M = t[idx][1];
}

/* the global synchronisation raster (TS 38.104 5.4.3.1) nearest to f */
static double raster_nearest(double f)
{
  double best = 0;
  if (f < 3000e6) {
    for (int N = (int)(f / 1.2e6) - 2; N <= (int)(f / 1.2e6) + 2; N++) {
      for (int M : {1, 3, 5}) {
        double g = N * 1.2e6 + M * 50e3;
        if (best == 0 || std::fabs(g - f) < std::fabs(best - f)) {
          best = g;
        }
      }
    }
  } else {
    double N = std::round((f - 3000e6) / 1.44e6);
    best     = 3000e6 + N * 1.44e6;
  }
  return best;
}

static void make_cell(Cell& c)
{
  srsran::srsran_band_helper bh;
  c.scs_hz = (uint32_t)iget("scs_khz") * 1000;
  if (c.scs_hz != 15000 && c.scs_hz != 30000) {
    die("scs_khz must be 15 or 30");
  }
  c.scs     = c.scs_hz == 15000 ? srsran_subcarrier_spacing_15kHz : srsran_subcarrier_spacing_30kHz;
  c.mu      = c.scs_hz == 15000 ? 0 : 1;
  c.spf     = 10u << c.mu;
  c.nof_prb = (uint32_t)iget("nof_prb");
  c.pci     = (uint32_t)iget("pci");
  if (c.nof_prb < 24 || c.nof_prb > 275 || c.pci > 1007) {
    die("nof_prb 24..275, pci 0..1007");
  }
  c.srate = dget("srate_hz");
  if (c.srate <= 0) {
    static const uint32_t sizes[] = {256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096};
    for (uint32_t s : sizes) {
      if (s * 0.85 >= c.nof_prb * 12.0) {
        c.nfft = s;
        break;
      }
    }
    c.srate = (double)c.nfft * c.scs_hz;
  } else {
    c.nfft = (uint32_t)std::lround(c.srate / c.scs_hz);
    if (std::fabs(c.nfft * (double)c.scs_hz - c.srate) > 1e-3 || c.nfft % 128 != 0) {
      die("srate_hz must be a multiple of 128 x SCS");
    }
  }

  double f0 = dget("freq_mhz") * 1e6;
  if (f0 <= 0) {
    f0 = c.scs_hz == 15000 ? 1980e6 : 3500e6;
  }

  // CORESET#0 and its offset to the SSB
  c.cs0_idx = (uint32_t)iget("coreset0_idx");
  c.ss0_idx = (uint32_t)iget("ss0_idx");
  srsran_coreset_t cs0 = {};
  if (srsran_coreset_zero(c.pci, 0, c.scs, c.scs, c.cs0_idx, &cs0) < SRSRAN_SUCCESS) {
    die("coreset0_idx %u is not valid for %u kHz", c.cs0_idx, c.scs_hz / 1000);
  }
  int cs0_off = srsran_coreset0_ssb_offset(c.cs0_idx, c.scs, c.scs);
  c.cs0       = cs0;
  c.cs0.offset_rb = 0;
  c.cs0_bw    = srsran_coreset_get_bw(&cs0);
  c.cs0_dur   = cs0.duration;

  // SSB position: CRB ssb_rb of the carrier, k_ssb; then the nearest raster frequency, and the
  // carrier moved under it (point A from the SSB, the centre from point A)
  c.k_ssb    = (uint32_t)iget("k_ssb");
  long ssb_rb = iget("ssb_rb");
  if (ssb_rb < 0) {
    ssb_rb = ((long)c.nof_prb - 20) / 2;
    ssb_rb = std::max(ssb_rb, (long)cs0_off);
    ssb_rb = std::min(ssb_rb, (long)(c.nof_prb - c.cs0_bw + cs0_off));
  }
  if (c.mu == 1 && c.k_ssb % 2) {
    die("k_ssb must be even at 30 kHz (the SSB on the 30 kHz grid)");
  }
  if (c.k_ssb > 23 || (c.mu == 0 && c.k_ssb > 11)) {
    die("k_ssb 0..11 at 15 kHz, 0..23 at 30 kHz");
  }
  const double rb15   = 12 * 15e3;
  double       pA_apx = f0 - 6.0 * c.nof_prb * c.scs_hz;
  double       sc0    = pA_apx + ssb_rb * 12.0 * c.scs_hz + c.k_ssb * 15e3;
  c.ssb_freq          = raster_nearest(sc0 + 120.0 * c.scs_hz);
  sc0                 = c.ssb_freq - 120.0 * c.scs_hz;
  long O              = std::lround((sc0 - c.k_ssb * 15e3 - pA_apx) / rb15);
  if (c.mu == 1 && O % 2) {
    O += ((sc0 - c.k_ssb * 15e3 - pA_apx) / rb15 > O) ? 1 : -1;
  }
  c.offset_to_point_a = (uint32_t)O;
  c.pointA            = sc0 - c.k_ssb * 15e3 - O * rb15;
  c.fc                = c.pointA + 6.0 * c.nof_prb * c.scs_hz;
  c.ssb_sc0           = (uint32_t)std::lround((sc0 - c.pointA) / c.scs_hz);
  c.ssb_rb_first      = c.ssb_sc0 / 12;
  c.ssb_rb_last       = (c.ssb_sc0 + 239) / 12;
  if (c.ssb_rb_last >= c.nof_prb) {
    die("the SSB (CRB %u..%u) is outside the carrier", c.ssb_rb_first, c.ssb_rb_last);
  }
  // CORESET#0: offset RBs below the CRB overlapping the SSB's first RB (TS 38.213 13)
  long ssb_crb = (long)O * 15000 / (long)c.scs_hz;
  long cs0_crb = ssb_crb - cs0_off;
  if (cs0_crb < 0 || cs0_crb + c.cs0_bw > c.nof_prb) {
    die("CORESET#0 (CRB %ld..%ld) does not fit in the carrier: move ssb_rb", cs0_crb, cs0_crb + c.cs0_bw - 1);
  }
  c.cs0_crb = (uint32_t)cs0_crb;

  // capture
  double rx_off = dget("rx_offset_hz");
  c.fc_rx       = c.fc + rx_off;
  c.rx_shift_sc = -(int32_t)std::lround(rx_off / c.scs_hz);
  if (std::fabs(c.rx_shift_sc * (double)c.scs_hz + rx_off) > 1e-3) {
    die("rx_offset_hz must be a multiple of the SCS");
  }
  if (std::fmod(std::fabs(c.ssb_freq - c.fc_rx), (double)c.scs_hz) > 1e-3 &&
      std::fabs(std::fmod(std::fabs(c.ssb_freq - c.fc_rx), (double)c.scs_hz) - c.scs_hz) > 1e-3) {
    die("internal: SSB not on the capture's subcarrier grid");
  }

  // band, duplex, SSB pattern
  c.band = (uint32_t)iget("band");
  if (c.band == 0) {
    c.band = bh.get_band_from_dl_freq_Hz_and_scs(c.fc, c.scs);
    if (c.band == UINT16_MAX) {
      c.band = bh.get_band_from_dl_freq_Hz(c.fc);
    }
  }
  c.duplex      = bh.get_duplex_mode(c.band);
  c.ssb_pattern = srsran::srsran_band_helper::get_ssb_pattern(c.band, c.scs);
  if (c.ssb_pattern == SRSRAN_SSB_PATTERN_INVALID || c.duplex == SRSRAN_DUPLEX_MODE_INVALID) {
    die("band n%u has no %u kHz SSB", c.band, c.scs_hz / 1000);
  }
  c.ssb_l0 = c.ssb_pattern == SRSRAN_SSB_PATTERN_B ? 4 : 2;
  c.tdd    = c.duplex == SRSRAN_DUPLEX_MODE_TDD;
  if (c.tdd) {
    c.tdd_period_ms  = dget("tdd_period_ms");
    c.tdd_period_slots = (uint32_t)std::lround(c.tdd_period_ms * (1u << c.mu));
    c.tdd_dl_slots   = (uint32_t)iget("tdd_dl_slots");
    c.tdd_dl_symbols = (uint32_t)iget("tdd_dl_symbols");
    c.tdd_ul_slots   = (uint32_t)iget("tdd_ul_slots");
    c.tdd_ul_symbols = (uint32_t)iget("tdd_ul_symbols");
    if (c.tdd_dl_slots + c.tdd_ul_slots + 1 > c.tdd_period_slots || c.tdd_dl_slots < 1 ||
        c.tdd_dl_symbols + c.tdd_ul_symbols > 14 || c.spf % c.tdd_period_slots) {
      die("TDD pattern does not fit its period (%u slots)", c.tdd_period_slots);
    }
  }

  // Type0 occasion
  double O0, M0;
  ss0_table(c.ss0_idx, &O0, &M0);
  c.n0    = (uint32_t)(O0 * (1u << c.mu)) % c.spf;
  c.sfn_c = ((uint32_t)(O0 * (1u << c.mu)) / c.spf) % 2;
  if (dl_symbols(c, c.n0) < 14 || dl_symbols(c, (c.n0 + 1) % c.spf) < 14 || dl_symbols(c, 0) < 14) {
    die("slot 0 and the Type0 slots %u, %u must be DL slots in the TDD pattern", c.n0, c.n0 + 1);
  }
  long tp = iget("dmrs_typeA_pos");
  c.typeA_pos = (tp == 3 || (tp == 0 && c.cs0_dur == 3)) ? srsran_dmrs_sch_typeA_pos_3 : srsran_dmrs_sch_typeA_pos_2;

  // CORESET 1
  srsran_coreset_t& cs1 = c.cs1;
  cs1                   = {};
  cs1.id                = 1;
  cs1.duration          = (uint32_t)iget("cs1_duration");
  long rb0 = iget("cs1_rb_start"), nrb = iget("cs1_nof_rb");
  if (nrb == 0) {
    nrb = ((long)c.nof_prb - rb0) / 6 * 6;
  }
  if (rb0 % 6 || nrb % 6 || nrb <= 0 || rb0 + nrb > (long)c.nof_prb || cs1.duration < 1 || cs1.duration > 3) {
    die("CORESET 1: rb_start and nof_rb multiples of 6 inside the carrier, duration 1..3");
  }
  for (long g = rb0 / 6; g < (rb0 + nrb) / 6; g++) {
    cs1.freq_resources[g] = true;
  }
  c.cs1_bw = (uint32_t)nrb;
  for (uint32_t i = 0; i < c.cs1_bw; i++) {
    c.cs1_crb[i] = (uint32_t)rb0 + i;
  }
  cs1.precoder_granularity = srsran_coreset_precoder_granularity_reg_bundle;
  if (iget("cs1_interleaved")) {
    cs1.mapping_type = srsran_coreset_mapping_type_interleaved;
    auto bsz         = [](long v) {
      switch (v) {
        case 2:
          return srsran_coreset_bundle_size_n2;
        case 3:
          return srsran_coreset_bundle_size_n3;
        case 6:
          return srsran_coreset_bundle_size_n6;
      }
      die("bundle/interleaver sizes are 2, 3 or 6");
    };
    long L = iget("cs1_reg_bundle"), R = iget("cs1_interleaver_size");
    cs1.reg_bundle_size  = bsz(L);
    cs1.interleaver_size = bsz(R);
    long sh              = iget("cs1_shift");
    cs1.shift_index      = (uint32_t)(sh < 0 ? c.pci : sh);
    long nreg            = nrb * cs1.duration;
    if (nreg % (L * R) || L % cs1.duration) {
      die("interleaved CORESET 1: %ld REGs not a multiple of L*R = %ld, or L %% duration != 0", nreg, L * R);
    }
  } else {
    cs1.mapping_type = srsran_coreset_mapping_type_non_interleaved;
  }
  long sid = iget("cs1_scrambling_id");
  if (sid >= 0) {
    cs1.dmrs_scrambling_id_present = true;
    cs1.dmrs_scrambling_id         = (uint32_t)sid;
  }
}

/* ------------------------------------------------------------------------------------------ */
/* SIB1, as srsRAN's gNB fills it (srsgnb/src/stack/rrc/cell_asn1_config.cc)                  */
/* ------------------------------------------------------------------------------------------ */

static uint32_t riv_type1(uint32_t N, uint32_t start, uint32_t len)
{
  return (len - 1) <= N / 2 ? N * (len - 1) + start : N * (N - len + 1) + (N - 1 - start);
}

static std::vector<uint8_t> make_sib1(const Cell& c)
{
  bcch_dl_sch_msg_s msg;
  sib1_s&           sib1 = msg.msg.set_c1().set_sib_type1();

  sib1.cell_sel_info_present            = true;
  sib1.cell_sel_info.q_rx_lev_min       = -70;
  sib1.cell_sel_info.q_qual_min_present = true;
  sib1.cell_sel_info.q_qual_min         = -20;

  sib1.cell_access_related_info.plmn_id_list.resize(1);
  auto& pl = sib1.cell_access_related_info.plmn_id_list[0];
  pl.plmn_id_list.resize(1);
  pl.plmn_id_list[0].mcc_present = true;
  pl.plmn_id_list[0].mcc         = {0, 0, 1};
  pl.plmn_id_list[0].mnc.resize(2);
  pl.plmn_id_list[0].mnc[0] = 0;
  pl.plmn_id_list[0].mnc[1] = 1;
  pl.tac_present            = true;
  pl.tac.from_number(7);
  pl.cell_id.from_number((0x19bULL << 8) + c.pci % 256);
  pl.cell_reserved_for_oper.value = plmn_id_info_s::cell_reserved_for_oper_opts::not_reserved;

  sib1.conn_est_fail_ctrl_present                             = true;
  sib1.conn_est_fail_ctrl.conn_est_fail_count.value           = conn_est_fail_ctrl_s::conn_est_fail_count_opts::n1;
  sib1.conn_est_fail_ctrl.conn_est_fail_offset_validity.value = conn_est_fail_ctrl_s::conn_est_fail_offset_validity_opts::s30;
  sib1.conn_est_fail_ctrl.conn_est_fail_offset_present        = true;
  sib1.conn_est_fail_ctrl.conn_est_fail_offset                = 1;

  sib1.serving_cell_cfg_common_present = true;
  serving_cell_cfg_common_sib_s& sc   = sib1.serving_cell_cfg_common;
  const auto                     scs  = (subcarrier_spacing_opts::options)c.scs;

  // DL
  auto& fdl = sc.dl_cfg_common.freq_info_dl;
  fdl.freq_band_list.resize(1);
  fdl.freq_band_list[0].freq_band_ind_nr_present = true;
  fdl.freq_band_list[0].freq_band_ind_nr         = (uint16_t)c.band;
  fdl.offset_to_point_a                          = (uint16_t)c.offset_to_point_a;
  fdl.scs_specific_carrier_list.resize(1);
  fdl.scs_specific_carrier_list[0].offset_to_carrier  = 0;
  fdl.scs_specific_carrier_list[0].subcarrier_spacing = scs;
  fdl.scs_specific_carrier_list[0].carrier_bw         = (uint16_t)c.nof_prb;
  auto& bwp                                           = sc.dl_cfg_common.init_dl_bwp;
  bwp.generic_params.location_and_bw                  = (uint16_t)riv_type1(275, 0, c.nof_prb);
  bwp.generic_params.subcarrier_spacing               = scs;
  bwp.pdcch_cfg_common_present                        = true;
  auto& pc                                            = bwp.pdcch_cfg_common.set_setup();
  pc.common_search_space_list.resize(1);
  search_space_s& ss1                             = pc.common_search_space_list[0];
  ss1.search_space_id                             = 1;
  ss1.ctrl_res_set_id_present                     = true;
  ss1.ctrl_res_set_id                             = 0;
  ss1.monitoring_slot_periodicity_and_offset_present = true;
  ss1.monitoring_slot_periodicity_and_offset.set_sl1();
  ss1.monitoring_symbols_within_slot_present = true;
  ss1.monitoring_symbols_within_slot.from_number(0b10000000000000);
  ss1.search_space_type_present                                                    = true;
  ss1.search_space_type.set_common().dci_format0_minus0_and_format1_minus0_present = true;
  ss1.nrof_candidates_present                                                      = true;
  {
    std::vector<std::string> cc = split(sget("css_candidates"), ',');
    uint32_t                 n[5] = {};
    for (size_t i = 0; i < 5 && i < cc.size(); i++) {
      n[i] = (uint32_t)atoi(cc[i].c_str());
    }
    asn1::number_to_enum(ss1.nrof_candidates.aggregation_level1, n[0]);
    asn1::number_to_enum(ss1.nrof_candidates.aggregation_level2, n[1]);
    asn1::number_to_enum(ss1.nrof_candidates.aggregation_level4, n[2]);
    asn1::number_to_enum(ss1.nrof_candidates.aggregation_level8, n[3]);
    asn1::number_to_enum(ss1.nrof_candidates.aggregation_level16, n[4]);
  }
  pc.search_space_sib1_present           = true;
  pc.search_space_sib1                   = 0;
  pc.search_space_other_sys_info_present = true;
  pc.search_space_other_sys_info         = 1;
  pc.paging_search_space_present         = true;
  pc.paging_search_space                 = 1;
  pc.ra_search_space_present             = true;
  pc.ra_search_space                     = 1;
  bwp.pdsch_cfg_common_present           = true;
  auto& pdsch                            = bwp.pdsch_cfg_common.set_setup();
  pdsch.pdsch_time_domain_alloc_list.resize(1);
  pdsch.pdsch_time_domain_alloc_list[0].map_type.value       = pdsch_time_domain_res_alloc_s::map_type_opts::type_a;
  pdsch.pdsch_time_domain_alloc_list[0].start_symbol_and_len = 40;
  sc.dl_cfg_common.bcch_cfg.mod_period_coeff.value           = bcch_cfg_s::mod_period_coeff_opts::n4;
  sc.dl_cfg_common.pcch_cfg.default_paging_cycle.value       = paging_cycle_opts::rf128;
  sc.dl_cfg_common.pcch_cfg.nand_paging_frame_offset.set_one_t();
  sc.dl_cfg_common.pcch_cfg.ns.value = pcch_cfg_s::ns_opts::one;

  // UL
  sc.ul_cfg_common_present = true;
  auto& ul                 = sc.ul_cfg_common;
  ul.freq_info_ul.freq_band_list.resize(1);
  ul.freq_info_ul.freq_band_list[0].freq_band_ind_nr_present = true;
  ul.freq_info_ul.freq_band_list[0].freq_band_ind_nr         = (uint16_t)c.band;
  ul.freq_info_ul.scs_specific_carrier_list.resize(1);
  ul.freq_info_ul.scs_specific_carrier_list[0].offset_to_carrier  = 0;
  ul.freq_info_ul.scs_specific_carrier_list[0].subcarrier_spacing = scs;
  ul.freq_info_ul.scs_specific_carrier_list[0].carrier_bw         = (uint16_t)c.nof_prb;
  ul.freq_info_ul.p_max_present                                   = true;
  ul.freq_info_ul.p_max                                           = 23;
  ul.init_ul_bwp.generic_params.location_and_bw                   = (uint16_t)riv_type1(275, 0, c.nof_prb);
  ul.init_ul_bwp.generic_params.subcarrier_spacing                = scs;
  ul.init_ul_bwp.rach_cfg_common_present                          = true;
  rach_cfg_common_s& rach                                         = ul.init_ul_bwp.rach_cfg_common.set_setup();
  rach.rach_cfg_generic.prach_cfg_idx                             = c.tdd ? 8 : 0;
  rach.rach_cfg_generic.msg1_fdm.value                            = rach_cfg_generic_s::msg1_fdm_opts::one;
  rach.rach_cfg_generic.msg1_freq_start                           = 1;
  rach.rach_cfg_generic.zero_correlation_zone_cfg                 = 0;
  rach.rach_cfg_generic.preamb_rx_target_pwr                      = -110;
  rach.rach_cfg_generic.preamb_trans_max.value                    = rach_cfg_generic_s::preamb_trans_max_opts::n7;
  rach.rach_cfg_generic.pwr_ramp_step.value                       = rach_cfg_generic_s::pwr_ramp_step_opts::db4;
  rach.rach_cfg_generic.ra_resp_win.value                         = rach_cfg_generic_s::ra_resp_win_opts::sl10;
  rach.ssb_per_rach_occasion_and_cb_preambs_per_ssb_present       = true;
  asn1::number_to_enum(rach.ssb_per_rach_occasion_and_cb_preambs_per_ssb.set_one(), 64);
  rach.ra_contention_resolution_timer.value = rach_cfg_common_s::ra_contention_resolution_timer_opts::sf64;
  rach.prach_root_seq_idx.set_l839()        = 1;
  rach.restricted_set_cfg.value             = rach_cfg_common_s::restricted_set_cfg_opts::unrestricted_set;
  ul.init_ul_bwp.pusch_cfg_common_present   = true;
  pusch_cfg_common_s& pusch                 = ul.init_ul_bwp.pusch_cfg_common.set_setup();
  pusch.pusch_time_domain_alloc_list.resize(1);
  pusch.pusch_time_domain_alloc_list[0].k2_present           = true;
  pusch.pusch_time_domain_alloc_list[0].k2                   = 4;
  pusch.pusch_time_domain_alloc_list[0].map_type.value       = pusch_time_domain_res_alloc_s::map_type_opts::type_a;
  pusch.pusch_time_domain_alloc_list[0].start_symbol_and_len = 27;
  pusch.p0_nominal_with_grant_present                        = true;
  pusch.p0_nominal_with_grant                                = -76;
  ul.init_ul_bwp.pucch_cfg_common_present                    = true;
  pucch_cfg_common_s& pucch                                  = ul.init_ul_bwp.pucch_cfg_common.set_setup();
  pucch.pucch_res_common_present                             = true;
  pucch.pucch_res_common                                     = 11;
  pucch.pucch_group_hop.value                                = pucch_cfg_common_s::pucch_group_hop_opts::neither;
  pucch.p0_nominal_present                                   = true;
  pucch.p0_nominal                                           = -90;
  ul.time_align_timer_common.value                           = time_align_timer_opts::infinity;

  sc.ssb_positions_in_burst.in_one_group.from_number(0x80);
  uint32_t sp = (uint32_t)iget("ssb_period_ms");
  using sp_e  = serving_cell_cfg_common_sib_s::ssb_periodicity_serving_cell_opts;
  sc.ssb_periodicity_serving_cell.value = sp == 5    ? sp_e::ms5
                                          : sp == 10 ? sp_e::ms10
                                          : sp == 20 ? sp_e::ms20
                                          : sp == 40 ? sp_e::ms40
                                          : sp == 80 ? sp_e::ms80
                                                     : sp_e::ms160;
  sc.n_timing_advance_offset_present    = true;
  sc.n_timing_advance_offset            = serving_cell_cfg_common_sib_s::n_timing_advance_offset_opts::n0;
  if (c.tdd) {
    sc.tdd_ul_dl_cfg_common_present          = true;
    auto& tdd                                = sc.tdd_ul_dl_cfg_common;
    tdd.ref_subcarrier_spacing.value         = scs;
    using per                                = tdd_ul_dl_pattern_s::dl_ul_tx_periodicity_opts;
    const double                     p       = c.tdd_period_ms;
    tdd.pattern1.dl_ul_tx_periodicity.value  = p == 0.5    ? per::ms0p5
                                               : p == 0.625 ? per::ms0p625
                                               : p == 1     ? per::ms1
                                               : p == 1.25  ? per::ms1p25
                                               : p == 2     ? per::ms2
                                               : p == 2.5   ? per::ms2p5
                                               : p == 5     ? per::ms5
                                                            : per::ms10;
    tdd.pattern1.nrof_dl_slots   = (uint16_t)c.tdd_dl_slots;
    tdd.pattern1.nrof_dl_symbols = (uint8_t)c.tdd_dl_symbols;
    tdd.pattern1.nrof_ul_slots   = (uint16_t)c.tdd_ul_slots;
    tdd.pattern1.nrof_ul_symbols = (uint8_t)c.tdd_ul_symbols;
  }
  sc.ss_pbch_block_pwr = -16;

  sib1.ue_timers_and_consts_present    = true;
  sib1.ue_timers_and_consts.t300.value = ue_timers_and_consts_s::t300_opts::ms1000;
  sib1.ue_timers_and_consts.t301.value = ue_timers_and_consts_s::t301_opts::ms1000;
  sib1.ue_timers_and_consts.t310.value = ue_timers_and_consts_s::t310_opts::ms1000;
  sib1.ue_timers_and_consts.n310.value = ue_timers_and_consts_s::n310_opts::n1;
  sib1.ue_timers_and_consts.t311.value = ue_timers_and_consts_s::t311_opts::ms30000;
  sib1.ue_timers_and_consts.n311.value = ue_timers_and_consts_s::n311_opts::n1;
  sib1.ue_timers_and_consts.t319.value = ue_timers_and_consts_s::t319_opts::ms1000;

  std::vector<uint8_t> buf(1024);
  asn1::bit_ref        bref(buf.data(), (uint32_t)buf.size());
  if (msg.pack(bref) != asn1::SRSASN_SUCCESS) {
    die("cannot pack SIB1");
  }
  buf.resize(bref.distance_bytes());
  return buf;
}

/* ------------------------------------------------------------------------------------------ */
/* traffic                                                                                     */
/* ------------------------------------------------------------------------------------------ */

struct Ue {
  uint32_t idx;
  uint16_t rnti;
  bool     present     = true;
  double   absent_until = -1; // re-attach: the new RNTI appears then
  double   idle_until  = -1;
  uint32_t dl_size, ul_size;
  double   p_dl, p_ul;
  uint32_t base_al; // log2
  uint64_t sent = 0, blocked = 0;
};

struct Truth {
  FILE* f;
  FILE* ev;
};

struct Sim {
  Cell                     c = {};
  std::mt19937_64          rng;
  std::vector<Ue>          ues;
  std::set<uint16_t>       used_rntis;
  uint32_t                 fb_uss, fb_css; // 1_0/0_0 sizes: USS on CORESET 1 (BWP), CSS on CORESET#0
  std::vector<double>      al_w;           // weight per log2 AL
  uint32_t                 ss_cand[5], css_cand[5];
  bool                     free_place;
  // PHY objects
  srsran_carrier_nr_t      car_main = {}, car_cs0 = {};
  srsran_pdcch_nr_t        pdcch = {};
  srsran_pdsch_nr_t        pdsch = {};
  srsran_dmrs_sch_t        dmrs_sch = {};
  srsran_softbuffer_tx_t   sb = {};
  srsran_dci_nr_t          dci_si = {};
  srsran_ssb_t             ssb = {};
  srsran_dft_plan_t        ifft = {};
  std::vector<cpx>         grid, grid0, fbuf, tbuf;
  std::vector<uint8_t>     occ; // [14][nof_prb]
  std::vector<uint8_t>     sib1;
  std::vector<uint8_t>     sib1_tb;
  // stats
  uint64_t                 nof_dci = 0, nof_blocked = 0, nof_sib1 = 0;
};

static double urand(Sim& s)
{
  return std::uniform_real_distribution<double>(0, 1)(s.rng);
}
static uint32_t irand(Sim& s, uint32_t lo, uint32_t hi)
{
  return std::uniform_int_distribution<uint32_t>(lo, hi)(s.rng);
}

static uint16_t new_rnti(Sim& s)
{
  uint32_t lo = (uint32_t)iget("rnti_min"), hi = (uint32_t)iget("rnti_max");
  for (int t = 0; t < 100000; t++) {
    uint16_t r = (uint16_t)irand(s, lo, hi);
    if (!s.used_rntis.count(r)) {
      s.used_rntis.insert(r);
      return r;
    }
  }
  die("no free RNTI in [%#x, %#x]", lo, hi);
}

/* DCI 1_1 and 0_1 sizes of one UE: distinct from each other and from the fallback sizes
  (the size alignment of TS 38.212 7.3.1.0 keeps them apart in a real cell) */
static void draw_sizes(Sim& s, Ue& u)
{
  auto bad = [&](uint32_t v) { return v == s.fb_uss || v == s.fb_css; };
  uint32_t odl = u.dl_size, oul = u.ul_size;
  for (int t = 0; t < 1000; t++) {
    u.dl_size = irand(s, (uint32_t)iget("dl_size_min"), (uint32_t)iget("dl_size_max"));
    u.ul_size = irand(s, (uint32_t)iget("ul_size_min"), (uint32_t)iget("ul_size_max"));
    if (!bad(u.dl_size) && !bad(u.ul_size) && u.dl_size != u.ul_size && (u.dl_size != odl || u.ul_size != oul)) {
      return;
    }
  }
  die("cannot draw DCI sizes: widen the size ranges");
}

static uint32_t draw_al(Sim& s)
{
  double tot = 0;
  for (double w : s.al_w) {
    tot += w;
  }
  double r = urand(s) * tot;
  for (uint32_t l = 0; l < 5; l++) {
    if (r < s.al_w[l]) {
      return l;
    }
    r -= s.al_w[l];
  }
  return 2;
}

static void new_ue(Sim& s, Ue& u, uint16_t rnti)
{
  u.rnti = rnti;
  u.p_dl = dget("p_dl_min") + urand(s) * (dget("p_dl_max") - dget("p_dl_min"));
  u.p_ul = dget("p_ul_min") + urand(s) * (dget("p_ul_max") - dget("p_ul_min"));
  u.base_al = draw_al(s);
}

static void log_event(Sim& s, Truth& t, double ts, const char* what, const Ue& u, uint16_t old_rnti)
{
  fprintf(t.ev, "%.6f,%s,%u,0x%04x,0x%04x,%u,%u\n", ts, what, u.idx, old_rnti, u.rnti, u.dl_size, u.ul_size);
  printf("t=%.3f s %-8s ue %u: rnti 0x%04x -> 0x%04x, sizes DL %u UL %u\n", ts, what, u.idx, old_rnti, u.rnti, u.dl_size,
         u.ul_size);
}

/* ------------------------------------------------------------------------------------------ */
/* PDCCH placement                                                                             */
/* ------------------------------------------------------------------------------------------ */

/* CRBs of a candidate (CORESET-relative RB mask mapped to the carrier) */
static void cand_crbs(const Sim& s, uint32_t coreset_id, const srsran_dci_location_t& loc, std::vector<uint32_t>& out)
{
  const srsran_coreset_t* cs = coreset_id == 0 ? &s.c.cs0 : &s.c.cs1;
  bool                    mask[SRSRAN_MAX_PRB_NR] = {};
  srsran_pdcch_nr_cce_to_reg_mapping(cs, &loc, mask);
  out.clear();
  uint32_t bw = srsran_coreset_get_bw(cs);
  for (uint32_t rb = 0; rb < bw; rb++) {
    if (mask[rb]) {
      out.push_back(coreset_id == 0 ? s.c.cs0_crb + rb : s.c.cs1_crb[rb]);
    }
  }
}

static bool crbs_free(const Sim& s, const std::vector<uint32_t>& crbs, uint32_t dur)
{
  for (uint32_t l = 0; l < dur; l++) {
    for (uint32_t r : crbs) {
      if (s.occ[l * s.c.nof_prb + r]) {
        return false;
      }
    }
  }
  return true;
}

static void crbs_take(Sim& s, const std::vector<uint32_t>& crbs, uint32_t dur)
{
  for (uint32_t l = 0; l < dur; l++) {
    for (uint32_t r : crbs) {
      s.occ[l * s.c.nof_prb + r] = 1;
    }
  }
}

/* Encodes one DCI message at msg.ctx.location in CORESET msg.ctx.coreset_id */
static void put_pdcch(Sim& s, uint32_t slot, srsran_dci_msg_nr_t& msg)
{
  const bool              cs0  = msg.ctx.coreset_id == 0;
  srsran_carrier_nr_t*    car  = cs0 ? &s.car_cs0 : &s.car_main;
  const srsran_coreset_t* cs   = cs0 ? &s.c.cs0 : &s.c.cs1;
  cf_t*                   grid = (cf_t*)(cs0 ? s.grid0.data() : s.grid.data());
  srsran_slot_cfg_t       sl   = {slot};
  if (srsran_pdcch_nr_set_carrier(&s.pdcch, car, cs) < SRSRAN_SUCCESS ||
      srsran_dmrs_pdcch_put(car, cs, &sl, &msg.ctx.location, grid) < SRSRAN_SUCCESS ||
      srsran_pdcch_nr_encode(&s.pdcch, &msg, grid) < SRSRAN_SUCCESS) {
    die("PDCCH encode failed (coreset %u L=%u ncce=%u bits=%u)", msg.ctx.coreset_id, 1u << msg.ctx.location.L,
        msg.ctx.location.ncce, msg.nof_bits);
  }
}

/* Finds a free candidate for (coreset, AL), hashing (USS) or free placement */
static bool find_location(Sim& s, uint32_t coreset_id, srsran_search_space_type_t ss_type, uint16_t rnti, uint32_t slot,
                          uint32_t al, srsran_dci_location_t* loc, std::vector<uint32_t>& crbs)
{
  const srsran_coreset_t* cs    = coreset_id == 0 ? &s.c.cs0 : &s.c.cs1;
  uint32_t                n_cce = srsran_coreset_get_bw(cs) * cs->duration / 6;
  uint32_t                L     = 1u << al;
  if (n_cce < L) {
    return false;
  }
  std::vector<uint32_t> locs;
  const uint32_t*       ncand = coreset_id == 0 ? s.css_cand : s.ss_cand;
  if (s.free_place && coreset_id != 0) {
    for (uint32_t n = 0; n + L <= n_cce; n += L) {
      locs.push_back(n);
    }
    std::shuffle(locs.begin(), locs.end(), s.rng);
  } else {
    if (ncand[al] == 0) {
      return false;
    }
    srsran_search_space_t ss = {};
    ss.id                    = coreset_id == 0 ? 1 : 2;
    ss.coreset_id            = coreset_id;
    ss.type                  = ss_type;
    ss.nof_candidates[al]    = ncand[al];
    uint32_t l[SRSRAN_SEARCH_SPACE_MAX_NOF_CANDIDATES_NR];
    int      n = srsran_pdcch_nr_locations_coreset(cs, &ss, rnti, al, slot, l);
    for (int i = 0; i < n; i++) {
      locs.push_back(l[i]);
    }
  }
  for (uint32_t n : locs) {
    loc->L    = al;
    loc->ncce = n;
    cand_crbs(s, coreset_id, *loc, crbs);
    if (crbs_free(s, crbs, cs->duration)) {
      return true;
    }
  }
  return false;
}

/* ------------------------------------------------------------------------------------------ */
/* SIB1 DCI + PDSCH on CORESET#0                                                               */
/* ------------------------------------------------------------------------------------------ */

struct Sib1Plan {
  bool     ok = false;
  uint32_t rb_start, nof_rb, mcs;
};

static bool sib1_grant(Sim& s, uint32_t slot, uint32_t rb_start, uint32_t nof_rb, uint32_t mcs, srsran_dci_msg_nr_t* msg,
                       srsran_sch_cfg_nr_t* cfg)
{
  srsran_dci_dl_nr_t dci        = {};
  dci.ctx.location.L            = (uint32_t)std::log2((double)iget("sib1_al"));
  dci.ctx.location.ncce         = 0;
  dci.ctx.ss_type               = srsran_search_space_type_common_0;
  dci.ctx.coreset_id            = 0;
  dci.ctx.coreset_start_rb      = 0;
  dci.ctx.rnti_type             = srsran_rnti_type_si;
  dci.ctx.format                = srsran_dci_format_nr_1_0;
  dci.ctx.rnti                  = 0xffff;
  dci.coreset0_bw               = s.c.cs0_bw;
  dci.freq_domain_assigment     = riv_type1(s.c.cs0_bw, rb_start, nof_rb);
  dci.time_domain_assigment     = 0;
  dci.vrb_to_prb_mapping        = 0;
  dci.mcs                       = mcs;
  dci.rv                        = 0;
  dci.sii                       = 0;
  if (srsran_dci_nr_dl_pack(&s.dci_si, &dci, msg) < SRSRAN_SUCCESS) {
    return false;
  }
  // the grant exactly as NR-Scope's SIB decoder derives it: from the unpacked DCI
  srsran_dci_msg_nr_t m2 = *msg;
  srsran_dci_dl_nr_t  rx = {};
  if (srsran_dci_nr_dl_unpack(&s.dci_si, &m2, &rx) < SRSRAN_SUCCESS) {
    return false;
  }
  srsran_sch_hl_cfg_nr_t hl = {};
  hl.typeA_pos              = s.c.typeA_pos;
  *cfg                      = {};
  cfg->dmrs.typeA_pos       = s.c.typeA_pos;
  srsran_slot_cfg_t sl      = {slot};
  return srsran_ra_dl_dci_to_grant_nr(&s.car_cs0, &sl, &hl, &rx, cfg, &cfg->grant) == SRSRAN_SUCCESS;
}

/* PRBs (CORESET#0-relative) for the SIB1, outside the SSB when it is in the same slot */
static Sib1Plan plan_sib1(Sim& s, uint32_t slot, bool with_ssb)
{
  Sib1Plan p;
  std::vector<bool> freeb(s.c.cs0_bw, true);
  if (with_ssb) {
    for (uint32_t r = s.c.ssb_rb_first; r <= s.c.ssb_rb_last; r++) {
      if (r >= s.c.cs0_crb && r < s.c.cs0_crb + s.c.cs0_bw) {
        freeb[r - s.c.cs0_crb] = false;
      }
    }
  }
  // largest free run
  uint32_t best_s = 0, best_n = 0;
  for (uint32_t a = 0; a < s.c.cs0_bw;) {
    if (!freeb[a]) {
      a++;
      continue;
    }
    uint32_t b = a;
    while (b < s.c.cs0_bw && freeb[b]) {
      b++;
    }
    if (b - a > best_n) {
      best_s = a, best_n = b - a;
    }
    a = b;
  }
  const uint32_t need = (uint32_t)s.sib1.size() * 8;
  for (uint32_t mcs = (uint32_t)iget("sib1_mcs"); mcs <= 9; mcs++) {
    for (uint32_t n = 1; n <= best_n; n++) {
      srsran_dci_msg_nr_t m;
      srsran_sch_cfg_nr_t cfg;
      if (sib1_grant(s, slot, best_s, n, mcs, &m, &cfg) && cfg.grant.tb[0].tbs >= (int)need) {
        p.ok = true, p.rb_start = best_s, p.nof_rb = n, p.mcs = mcs;
        return p;
      }
    }
  }
  return p;
}

static void put_sib1(Sim& s, Truth& t, uint64_t abs_slot, double ts, uint32_t sfn, uint32_t slot, const Sib1Plan& p)
{
  srsran_dci_msg_nr_t msg;
  srsran_sch_cfg_nr_t cfg;
  if (!sib1_grant(s, slot, p.rb_start, p.nof_rb, p.mcs, &msg, &cfg)) {
    die("SIB1 grant");
  }
  std::vector<uint32_t> crbs;
  cand_crbs(s, 0, msg.ctx.location, crbs);
  if (!crbs_free(s, crbs, s.c.cs0_dur)) {
    return; // cannot happen: the SIB1 is placed first
  }
  crbs_take(s, crbs, s.c.cs0_dur);
  put_pdcch(s, slot, msg);

  // PDSCH, on the CORESET#0 grid
  s.sib1_tb.assign(cfg.grant.tb[0].tbs / 8 + 1, 0);
  memcpy(s.sib1_tb.data(), s.sib1.data(), s.sib1.size());
  cfg.grant.tb[0].softbuffer.tx = &s.sb;
  srsran_softbuffer_tx_reset(&s.sb);
  uint8_t*          data[SRSRAN_MAX_TB]  = {s.sib1_tb.data()};
  cf_t*             syms[SRSRAN_MAX_PORTS] = {(cf_t*)s.grid0.data()};
  srsran_slot_cfg_t sl                   = {slot};
  if (srsran_dmrs_sch_put_sf(&s.dmrs_sch, &sl, &cfg, &cfg.grant, (cf_t*)s.grid0.data()) < SRSRAN_SUCCESS ||
      srsran_pdsch_nr_encode(&s.pdsch, &cfg, &cfg.grant, data, syms) < SRSRAN_SUCCESS) {
    die("SIB1 PDSCH encode");
  }
  for (uint32_t l = cfg.grant.S; l < cfg.grant.S + cfg.grant.L; l++) {
    for (uint32_t r = 0; r < s.c.cs0_bw; r++) {
      if (cfg.grant.prb_idx[r]) {
        s.occ[l * s.c.nof_prb + s.c.cs0_crb + r] = 1;
      }
    }
  }
  s.nof_sib1++;
  fprintf(t.f, "%lu,%.6f,%u,%u,0xffff,DL,1_0,%u,%u,%u,%u,%zu,0,css0,%u,0,-1,\n", (unsigned long)abs_slot, ts, sfn, slot,
          msg.nof_bits, 1u << msg.ctx.location.L, msg.ctx.location.ncce, crbs.empty() ? 0 : crbs[0], crbs.size(),
          s.c.pci);
}

/* ------------------------------------------------------------------------------------------ */
/* OFDM                                                                                        */
/* ------------------------------------------------------------------------------------------ */

/* CP of symbol l of the subframe, as NR-Scope's receiver computes it (ofdm.c) */
static uint32_t cp_len(const Cell& c, uint32_t l_sf)
{
  uint32_t cpc = (l_sf % (7u << c.mu)) == 0 ? 144 + (16u << c.mu) : 144;
  return (uint32_t)std::ceil((double)cpc * c.nfft / 2048.0);
}

static uint32_t slot_len(const Cell& c, uint32_t slot)
{
  uint32_t n = 0, l0 = (slot % (1u << c.mu)) * 14;
  for (uint32_t l = 0; l < 14; l++) {
    n += cp_len(c, l0 + l) + c.nfft;
  }
  return n;
}

static void modulate(Sim& s, uint32_t slot, cpx* out)
{
  const Cell& c      = s.c;
  const double norm  = 1.0 / std::sqrt((double)c.nfft);
  const int    K     = (int)(c.nof_prb * 12);
  const double half  = std::min(c.srate / 2.0 * dget("passband"), c.srate / 2.0);
  uint32_t     l0    = (slot % (1u << c.mu)) * 14;
  // samples from the subframe start to this slot's start
  uint64_t     count = 0;
  for (uint32_t l = 0; l < l0; l++) {
    count += cp_len(c, l) + c.nfft;
  }
  const bool partial = (K / 2 + std::abs(c.rx_shift_sc)) * (double)c.scs_hz > c.srate / 2.0;
  for (uint32_t l = 0; l < 14; l++) {
    const cpx* g = &s.grid[l * K];
    std::fill(s.fbuf.begin(), s.fbuf.end(), cpx(0, 0));
    for (int k = 0; k < K; k++) {
      int b = k - K / 2 + c.rx_shift_sc;
      if (partial && std::fabs((double)b * c.scs_hz) > half) {
        continue;
      }
      s.fbuf[(b + (int)c.nfft * 4) % (int)c.nfft] = g[k];
    }
    srsran_dft_run_c(&s.ifft, (cf_t*)s.fbuf.data(), (cf_t*)s.tbuf.data());
    uint32_t cp = cp_len(c, l0 + l);
    count += cp;
    double phase = -2.0 * M_PI * std::fmod(c.fc_rx * (double)count / c.srate, 1.0);
    cpx    ph    = cpx((float)(norm * std::cos(phase)), (float)(norm * std::sin(phase)));
    for (uint32_t n = 0; n < c.nfft; n++) {
      out[cp + n] = s.tbuf[n] * ph;
    }
    for (uint32_t n = 0; n < cp; n++) {
      out[n] = out[c.nfft + n];
    }
    out += cp + c.nfft;
    count += c.nfft;
  }
}

/* ------------------------------------------------------------------------------------------ */
/* main                                                                                        */
/* ------------------------------------------------------------------------------------------ */

static void usage()
{
  printf("gs_cellsim [-c file.ini] [--key value | --key=value]...\n\nkeys (default):\n");
  for (const KeyDoc* d = keys_doc; d->key; d++) {
    printf("  %-22s %-24s %s\n", d->key, d->def[0] ? d->def : "\"\"", d->doc);
  }
}

struct Ev {
  double   t;
  int      kind; // 0 resize, 1 churn
  uint32_t ue;
};

int main(int argc, char** argv)
{
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      usage();
      return 0;
    }
    if (a == "-c" && i + 1 < argc) {
      load_ini(argv[++i]);
    } else if (a.rfind("--", 0) == 0) {
      size_t eq = a.find('=');
      if (eq != std::string::npos) {
        set_key(a.substr(2, eq - 2), a.substr(eq + 1));
      } else if (i + 1 < argc) {
        set_key(a.substr(2), argv[++i]);
      } else {
        die("missing value for %s", a.c_str());
      }
    } else {
      die("bad argument %s", a.c_str());
    }
  }
  set_srsran_verbose_level(SRSRAN_VERBOSE_NONE);

  Sim s;
  Cell& c = s.c;
  make_cell(c);
  s.rng.seed((uint64_t)iget("seed"));
  s.sib1 = make_sib1(c);

  // sizes: DCI 1_0 = 28 + FDRA (TS 38.212 7.3.1.2.1), 0_0 aligned to it
  auto fdra = [](uint32_t n) { return (uint32_t)std::ceil(std::log2(n * (n + 1) / 2.0)); };
  s.fb_uss  = 28 + fdra(c.nof_prb);
  s.fb_css  = 28 + fdra(c.cs0_bw);
  {
    std::vector<std::string> a = split(sget("ss_candidates"), ','), b = split(sget("css_candidates"), ',');
    for (int i = 0; i < 5; i++) {
      s.ss_cand[i]  = i < (int)a.size() ? (uint32_t)atoi(a[i].c_str()) : 0;
      s.css_cand[i] = i < (int)b.size() ? (uint32_t)atoi(b[i].c_str()) : 0;
    }
    s.al_w.assign(5, 0);
    for (auto& e : split(sget("al_weights"), ',')) {
      std::vector<std::string> p = split(e, ':');
      int                      L = atoi(p[0].c_str());
      int                      l = L == 1 ? 0 : L == 2 ? 1 : L == 4 ? 2 : L == 8 ? 3 : L == 16 ? 4 : -1;
      if (l < 0 || p.size() != 2) {
        die("al_weights: AL:weight with AL in 1,2,4,8,16");
      }
      s.al_w[l] = atof(p[1].c_str());
    }
    s.free_place = sget("placement") == "free";
  }

  // PHY objects
  s.car_main                 = {};
  s.car_main.pci             = c.pci;
  s.car_main.nof_prb         = c.nof_prb;
  s.car_main.scs             = c.scs;
  s.car_main.max_mimo_layers = 1;
  s.car_main.dl_center_frequency_hz = c.fc;
  s.car_cs0                  = s.car_main;
  s.car_cs0.nof_prb          = c.cs0_bw;
  srsran_pdcch_nr_args_t pa  = {};
  if (srsran_pdcch_nr_init_tx(&s.pdcch, &pa) < SRSRAN_SUCCESS) {
    die("pdcch init");
  }
  srsran_pdsch_nr_args_t pda = {};
  pda.max_prb                = c.cs0_bw;
  pda.max_layers             = 1;
  if (srsran_pdsch_nr_init_enb(&s.pdsch, &pda) < SRSRAN_SUCCESS || srsran_pdsch_nr_set_carrier(&s.pdsch, &s.car_cs0) ||
      srsran_dmrs_sch_init(&s.dmrs_sch, false) || srsran_dmrs_sch_set_carrier(&s.dmrs_sch, &s.car_cs0) ||
      srsran_softbuffer_tx_init_guru(&s.sb, SRSRAN_SCH_NR_MAX_NOF_CB_LDPC, SRSRAN_LDPC_MAX_LEN_ENCODED_CB)) {
    die("pdsch init");
  }
  srsran_dci_cfg_nr_t dcfg = {};
  dcfg.bwp_dl_initial_bw   = 275;
  dcfg.bwp_ul_initial_bw   = 275;
  dcfg.bwp_dl_active_bw    = 275;
  dcfg.bwp_ul_active_bw    = 275;
  dcfg.monitor_common_0_0  = true;
  dcfg.monitor_0_0_and_1_0 = true;
  dcfg.monitor_0_1_and_1_1 = true;
  dcfg.coreset0_bw         = c.cs0_bw;
  if (srsran_dci_nr_set_cfg(&s.dci_si, &dcfg) < SRSRAN_SUCCESS) {
    die("dci cfg");
  }
  srsran_ssb_args_t sa = {};
  sa.enable_encode     = true;
  sa.max_srate_hz      = c.srate;
  sa.min_scs           = c.scs;
  srsran_ssb_cfg_t sc  = {};
  sc.srate_hz          = c.srate;
  sc.center_freq_hz    = c.fc_rx;
  sc.ssb_freq_hz       = c.ssb_freq;
  sc.scs               = c.scs;
  sc.pattern           = c.ssb_pattern;
  sc.duplex_mode       = c.duplex;
  sc.periodicity_ms    = (uint32_t)iget("ssb_period_ms");
  if (srsran_ssb_init(&s.ssb, &sa) < SRSRAN_SUCCESS || srsran_ssb_set_cfg(&s.ssb, &sc) < SRSRAN_SUCCESS) {
    die("ssb init");
  }
  if (srsran_dft_plan_c(&s.ifft, (int)c.nfft, SRSRAN_DFT_BACKWARD) < SRSRAN_SUCCESS) {
    die("dft");
  }
  srsran_dft_plan_set_norm(&s.ifft, false);
  srsran_dft_plan_set_mirror(&s.ifft, false);
  srsran_dft_plan_set_dc(&s.ifft, false);
  s.grid.assign(14 * c.nof_prb * 12, 0);
  s.grid0.assign(14 * c.cs0_bw * 12, 0);
  s.fbuf.assign(c.nfft, 0);
  s.tbuf.assign(c.nfft, 0);
  s.occ.assign(14 * c.nof_prb, 0);

  // SIB1 plans: in a slot with and without the SSB
  Sib1Plan plan_plain = plan_sib1(s, c.n0, false), plan_ssb = plan_sib1(s, c.n0, true);
  if (!plan_plain.ok) {
    die("the SIB1 (%zu bytes) does not fit in CORESET#0's %u RB", s.sib1.size(), c.cs0_bw);
  }

  // UEs
  std::vector<std::string> rl = split(sget("rnti_list"), ',');
  uint32_t                 nue = rl.empty() ? (uint32_t)iget("nof_ues") : (uint32_t)rl.size();
  s.ues.resize(nue);
  for (uint32_t i = 0; i < nue; i++) {
    Ue& u = s.ues[i];
    u.idx = i;
    uint16_t r;
    if (!rl.empty()) {
      r = (uint16_t)strtol(rl[i].c_str(), nullptr, 0);
      s.used_rntis.insert(r);
    } else {
      r = new_rnti(s);
    }
    new_ue(s, u, r);
    u.dl_size = u.ul_size = 0;
    draw_sizes(s, u);
  }

  // events
  const double         dur = dget("duration_s");
  std::vector<Ev>      evs;
  auto                 parse_ev = [&](const char* key, int kind) {
    for (auto& e : split(sget(key), ',')) {
      std::vector<std::string> p = split(e, ':');
      if (p.size() != 2 || (uint32_t)atoi(p[1].c_str()) >= nue) {
        die("%s: t_s:ue_index with ue_index < nof_ues", key);
      }
      evs.push_back({atof(p[0].c_str()), kind, (uint32_t)atoi(p[1].c_str())});
    }
  };
  parse_ev("resize_at", 0);
  parse_ev("churn_at", 1);
  for (int kind = 0; kind < 2; kind++) {
    double rate = dget(kind == 0 ? "resize_rate" : "churn_rate");
    if (rate > 0 && nue > 0) {
      double t = 0;
      while (true) {
        t += std::exponential_distribution<double>(rate)(s.rng);
        if (t >= dur) {
          break;
        }
        evs.push_back({t, kind, irand(s, 0, nue - 1)});
      }
    }
  }
  std::sort(evs.begin(), evs.end(), [](const Ev& a, const Ev& b) { return a.t < b.t; });

  // outputs
  std::string out_path = sget("out"), truth_path = sget("truth"), yaml_path = sget("yaml");
  std::string out_abs  = out_path;
  if (!out_path.empty() && out_path[0] != '/') {
    char cwd[4096];
    if (getcwd(cwd, sizeof(cwd))) {
      out_abs = std::string(cwd) + "/" + out_path;
    }
  }
  if (yaml_path.empty()) {
    yaml_path = out_path + ".yaml";
  }
  Truth t;
  t.f  = fopen(truth_path.c_str(), "w");
  t.ev = fopen((truth_path + ".events.csv").c_str(), "w");
  if (!t.f || !t.ev) {
    die("cannot write %s", truth_path.c_str());
  }
  fprintf(t.f, "abs_slot,t_s,sfn,slot,rnti,dir,format,size,al,ncce,first_crb,nof_rb,coreset_id,ss,n_id,n_rnti,ue,"
               "payload_hex\n");
  fprintf(t.ev, "t_s,event,ue,old_rnti,rnti,dl_size,ul_size\n");
  for (auto& u : s.ues) {
    log_event(s, t, 0, "start", u, 0);
  }

  // the NR-Scope config for this cell
  {
    FILE* y = fopen(yaml_path.c_str(), "w");
    if (!y) {
      die("cannot write %s", yaml_path.c_str());
    }
    std::string mode = sget("gs_mode");
    fprintf(y,
            "# NR-Scope on gs_cellsim output (generated): %s\n"
            "nof_usrp_dev: 1\n"
            "usrp_setting_0:\n"
            "  device_name: file\n"
            "  rf_args: \"rx_file=%s,base_srate=%.0f\"\n"
            "  ssb_freq: %.0f\n"
            "  rx_center_freq: %.0f\n"
            "  band_list: [%u]\n"
            "  rx_gain: 0\n  min_rx_gain: 0\n  max_rx_gain: 0\n"
            "  srate_hz: %.0f\n"
            "  srsran_srate_hz: %.0f\n"
            "  nof_carriers: 1\n  nof_antennas: 1\n"
            "  scs_index: %u\n"
            "  rf_log_level: \"info\"\n"
            "  nof_rnti_worker_groups: 1\n  nof_bwps: 1\n"
            "  nof_workers: %ld\n"
            "  cpu_affinity: false\n  disable_cfo: false\n  rnti_expiry_s: 30\n  agc: false\n"
            "  log_name: \"%s.nrscope.csv\"\n\n"
            "log_config:\n  recording_mode: false\n  quiet: true\n  local_log: false\n  push_to_google: false\n"
            "  enable_gui: false\n\n"
            "golden_sniffer:\n  enable: true\n  max_duration: %u\n  css_sizes: auto\n  mode: %s\n"
            "  feed_rntis: false\n  report_s: 5\n  csv: \"%s\"\n",
            out_path.c_str(), out_abs.c_str(), c.srate, c.ssb_freq, c.fc_rx, c.band, c.srate, c.srate, c.mu,
            iget("nof_workers"), out_path.c_str(), std::max(2u, std::max(c.cs0_dur, c.cs1.duration)), mode.c_str(),
            sget("gs_csv").c_str());
    fclose(y);
  }

  FILE* fo = fopen(out_path.c_str(), "wb");
  if (!fo) {
    die("cannot write %s", out_path.c_str());
  }

  printf("cell: band n%u %s, %u kHz, %u PRB, PCI %u; carrier centre %.6f MHz, point A %.6f MHz, SSB %.6f MHz "
         "(offsetToPointA %u, k_SSB %u, CRB %u..%u, pattern %s)\n",
         c.band, c.tdd ? "TDD" : "FDD", c.scs_hz / 1000, c.nof_prb, c.pci, c.fc / 1e6, c.pointA / 1e6, c.ssb_freq / 1e6,
         c.offset_to_point_a, c.k_ssb, c.ssb_rb_first, c.ssb_rb_last, srsran_ssb_pattern_to_str(c.ssb_pattern));
  printf("capture: centre %.6f MHz, %.3f Msps (FFT %u)%s\n", c.fc_rx / 1e6, c.srate / 1e6, c.nfft,
         (c.nof_prb * 6 + std::abs(c.rx_shift_sc)) * (double)c.scs_hz > c.srate / 2 ? ", PARTIAL: the carrier is wider"
                                                                                     : "");
  printf("CORESET#0: index %u, CRB %u..%u, %u symbol(s); search space #0 index %u: slot n0 %u, SFN %% 2 == %u; "
         "dmrs-TypeA-Position %u\n",
         c.cs0_idx, c.cs0_crb, c.cs0_crb + c.cs0_bw - 1, c.cs0_dur, c.ss0_idx, c.n0, c.sfn_c,
         c.typeA_pos == srsran_dmrs_sch_typeA_pos_3 ? 3 : 2);
  printf("SIB1: %zu bytes, %u RB at CORESET#0 RB %u, MCS %u (with the SSB in the slot: %s)\n", s.sib1.size(),
         plan_plain.nof_rb, plan_plain.rb_start, plan_plain.mcs,
         plan_ssb.ok ? "fits beside it" : "moved to slot n0+1");
  printf("CORESET 1: CRB %u..%u, %u symbol(s), %s, DM-RS N_ID %u%s; DCI 1_0/0_0 sizes: USS %u, CSS %u\n",
         c.cs1_crb[0], c.cs1_crb[c.cs1_bw - 1], c.cs1.duration,
         c.cs1.mapping_type == srsran_coreset_mapping_type_interleaved ? "interleaved" : "non-interleaved",
         c.cs1.dmrs_scrambling_id_present ? c.cs1.dmrs_scrambling_id : c.pci,
         c.cs1.dmrs_scrambling_id_present ? " (pdcch-DMRS-ScramblingID, n_RNTI = C-RNTI)" : " (PCI)", s.fb_uss, s.fb_css);
  if (c.tdd) {
    printf("TDD: %.3g ms, %u DL slots, %u DL + %u UL symbols, %u UL slots\n", c.tdd_period_ms, c.tdd_dl_slots,
           c.tdd_dl_symbols, c.tdd_ul_symbols, c.tdd_ul_slots);
  }
  printf("UEs: %u; %zu events; %.2f s -> %s (%.0f MB)\n", nue, evs.size(), dur, out_path.c_str(),
         dur * c.srate * 8 / 1e6);
  fflush(stdout);

  // channel state
  std::vector<std::pair<uint32_t, cpx>> taps;
  {
    double tot = 0;
    for (auto& e : split(sget("taps"), ',')) {
      std::vector<std::string> p = split(e, ':');
      double                   g = std::pow(10.0, (p.size() > 1 ? atof(p[1].c_str()) : 0) / 20.0);
      double                   a = taps.empty() ? 0 : 2 * M_PI * urand(s);
      taps.push_back({(uint32_t)atoi(p[0].c_str()), cpx((float)(g * std::cos(a)), (float)(g * std::sin(a)))});
      tot += g * g;
    }
    for (auto& tp : taps) {
      tp.second /= (float)std::sqrt(tot);
    }
  }
  uint32_t max_d = 0;
  for (auto& tp : taps) {
    max_d = std::max(max_d, tp.first);
  }
  
  const double      snr   = dget("snr_db");
  const double      sigma = std::isinf(snr) ? 0 : std::sqrt(std::pow(10.0, -snr / 10.0) / 2.0);
  const double      cfo   = dget("cfo_hz");
  const float       gain  = (float)std::pow(10.0, dget("gain_db") / 20.0);
  std::normal_distribution<float> gauss(0.0f, 1.0f);
  std::mt19937_64                 nrng((uint64_t)iget("seed") * 7919 + 13);
  uint64_t                        n_abs    = 0;
  uint64_t                        to_drop  = (uint64_t)iget("drop_samples");
  const uint32_t                  max_slot = slot_len(c, 0);
  std::vector<cpx>                buf(max_slot + max_d), out(max_slot);

  const uint64_t nof_slots = (uint64_t)std::ceil(dur * 1000.0 * (1u << c.mu));
  const uint32_t sib1_per  = (uint32_t)iget("sib1_period_frames");
  const uint32_t ssb_per   = (uint32_t)iget("ssb_period_ms") / 10;
  const uint32_t max_dl = (uint32_t)iget("max_dl_per_slot"), max_ul = (uint32_t)iget("max_ul_per_slot");
  const double   fb_p = dget("fallback_prob"), css_p = dget("css_fallback_prob"), sticky = dget("al_sticky");
  const double   idle_rate = dget("idle_rate"), idle_min = dget("idle_min_s"), idle_max = dget("idle_max_s");
  const double   gap = dget("reattach_gap_s");
  const double   slot_s = 1e-3 / (1u << c.mu);
  size_t         ev_i   = 0;
  uint64_t       samples_written = 0;

  for (uint64_t abs = 0; abs < nof_slots; abs++) {
    const uint32_t sfn  = (uint32_t)(((uint64_t)iget("sfn0") + abs / c.spf) % 1024);
    const uint32_t slot = (uint32_t)(abs % c.spf);
    const double   ts   = abs * slot_s - (double)to_drop / c.srate;
    const double   tnow = abs * slot_s;

    // events
    while (ev_i < evs.size() && evs[ev_i].t <= tnow) {
      Ev& e = evs[ev_i++];
      Ue& u = s.ues[e.ue];
      if (e.kind == 0) {
        if (!u.present) {
          continue;
        }
        draw_sizes(s, u);
        log_event(s, t, tnow, "resize", u, u.rnti);
      } else if (u.present) {
        uint16_t old = u.rnti;
        u.present    = false;
        u.absent_until = tnow + gap;
        log_event(s, t, tnow, "leave", u, old);
      }
    }
    for (auto& u : s.ues) {
      if (!u.present && u.absent_until >= 0 && tnow >= u.absent_until) {
        uint16_t old = u.rnti;
        new_ue(s, u, new_rnti(s));
        if (iget("churn_new_sizes")) {
          draw_sizes(s, u);
        }
        u.present      = true;
        u.absent_until = -1;
        u.idle_until   = -1;
        log_event(s, t, tnow, "attach", u, old);
      }
      if (u.present && idle_rate > 0 && tnow >= u.idle_until && urand(s) < idle_rate * slot_s) {
        u.idle_until = tnow + idle_min + urand(s) * (idle_max - idle_min);
        fprintf(t.ev, "%.6f,idle,%u,0x%04x,0x%04x,%u,%u\n", tnow, u.idx, u.rnti, u.rnti, u.dl_size, u.ul_size);
      }
    }

    std::fill(s.grid.begin(), s.grid.end(), cpx(0, 0));
    std::fill(s.grid0.begin(), s.grid0.end(), cpx(0, 0));
    std::fill(s.occ.begin(), s.occ.end(), 0);
    const uint32_t ndl     = dl_symbols(c, slot);
    const bool     ssb_now = slot == 0 && sfn % ssb_per == 0;
    if (ssb_now) {
      for (uint32_t l = c.ssb_l0; l < c.ssb_l0 + 4; l++) {
        for (uint32_t r = c.ssb_rb_first; r <= c.ssb_rb_last; r++) {
          s.occ[l * c.nof_prb + r] = 1;
        }
      }
    }
    // SIB1 in the first Type0 slot, or the second when the SSB leaves no room
    if (sfn % sib1_per == (c.sfn_c % sib1_per) || sib1_per == 1) {
      bool ssb_in_n0 = c.n0 == 0 && sfn % ssb_per == 0;
      if (slot == c.n0 && (!ssb_in_n0 || plan_ssb.ok)) {
        put_sib1(s, t, abs, ts, sfn, slot, ssb_in_n0 ? plan_ssb : plan_plain);
      } else if (slot == (c.n0 + 1) % c.spf && ssb_in_n0 && !plan_ssb.ok) {
        put_sib1(s, t, abs, ts, sfn, slot, plan_plain);
      }
    }

    // UE DCIs
    if (ndl >= c.cs1.duration) {
      struct Req {
        Ue*  u;
        bool dl;
      };
      std::vector<Req> req;
      for (auto& u : s.ues) {
        if (!u.present || tnow < u.idle_until) {
          continue;
        }
        if (urand(s) < u.p_dl) {
          req.push_back({&u, true});
        }
        if (urand(s) < u.p_ul) {
          req.push_back({&u, false});
        }
      }
      std::shuffle(req.begin(), req.end(), s.rng);
      uint32_t ndl_dci = 0, nul_dci = 0;
      for (auto& r : req) {
        if ((r.dl && ndl_dci >= max_dl) || (!r.dl && nul_dci >= max_ul)) {
          continue;
        }
        Ue&  u        = *r.u;
        bool fallback = urand(s) < fb_p;
        bool css      = fallback && urand(s) < css_p && ndl >= c.cs0_dur;
        uint32_t cs_id = css ? 0 : 1;
        uint32_t al    = urand(s) < sticky ? u.base_al : draw_al(s);
        srsran_dci_location_t loc = {};
        std::vector<uint32_t> crbs;
        srsran_search_space_type_t sst = css ? srsran_search_space_type_common_3 : srsran_search_space_type_ue;
        bool                       ok  = false;
        // the AL first, then the others nearest to it
        for (int d : {0, -1, 1, -2, 2, -3, 3, -4, 4}) {
          int a = (int)al + d;
          if (a < 0 || a > 4) {
            continue;
          }
          if (find_location(s, cs_id, sst, u.rnti, slot, (uint32_t)a, &loc, crbs)) {
            ok = true;
            break;
          }
          if (d == 0 && !s.free_place) {
            // the hashed candidates of this AL are taken: a gNB would try another AL
            continue;
          }
        }
        if (!ok) {
          u.blocked++;
          s.nof_blocked++;
          continue;
        }
        srsran_dci_msg_nr_t msg = {};
        msg.ctx.location        = loc;
        msg.ctx.ss_type         = sst;
        msg.ctx.coreset_id      = cs_id;
        msg.ctx.rnti_type       = srsran_rnti_type_c;
        msg.ctx.rnti            = u.rnti;
        msg.ctx.format          = r.dl ? (fallback ? srsran_dci_format_nr_1_0 : srsran_dci_format_nr_1_1)
                                       : (fallback ? srsran_dci_format_nr_0_0 : srsran_dci_format_nr_0_1);
        msg.nof_bits            = fallback ? (css ? s.fb_css : s.fb_uss) : (r.dl ? u.dl_size : u.ul_size);
        if (msg.nof_bits > sizeof(msg.payload)) {
          die("DCI size %u > %zu", msg.nof_bits, sizeof(msg.payload));
        }
        msg.payload[0] = r.dl ? 1 : 0;
        for (uint32_t b = 1; b < msg.nof_bits; b++) {
          msg.payload[b] = (uint8_t)(s.rng() & 1);
        }
        const srsran_coreset_t* cs = css ? &c.cs0 : &c.cs1;
        crbs_take(s, crbs, cs->duration);
        put_pdcch(s, slot, msg);
        (r.dl ? ndl_dci : nul_dci)++;
        u.sent++;
        s.nof_dci++;
        char hex[32] = {};
        for (uint32_t b = 0; b < (msg.nof_bits + 3) / 4; b++) {
          uint32_t v = 0;
          for (uint32_t k = 0; k < 4; k++) {
            v = v * 2 + (4 * b + k < msg.nof_bits ? msg.payload[4 * b + k] : 0);
          }
          hex[b] = "0123456789abcdef"[v];
        }
        bool     nr   = !css && c.cs1.dmrs_scrambling_id_present;
        uint32_t n_id = (!css && c.cs1.dmrs_scrambling_id_present) ? c.cs1.dmrs_scrambling_id : c.pci;
        fprintf(t.f, "%lu,%.6f,%u,%u,0x%04x,%s,%s,%u,%u,%u,%u,%zu,%u,%s,%u,%d,%u,%s\n", (unsigned long)abs, ts, sfn, slot,
                u.rnti, r.dl ? "DL" : "UL", srsran_dci_format_nr_string(msg.ctx.format), msg.nof_bits, 1u << loc.L,
                loc.ncce, crbs[0], crbs.size(), cs_id, css ? "css" : "uss", n_id, nr ? 1 : 0, u.idx, hex);
      }
    }

    // CORESET#0 grid into the carrier
    for (uint32_t l = 0; l < 14; l++) {
      cpx*       dst = &s.grid[(l * c.nof_prb + c.cs0_crb) * 12];
      const cpx* src = &s.grid0[l * c.cs0_bw * 12];
      for (uint32_t k = 0; k < c.cs0_bw * 12; k++) {
        dst[k] += src[k];
      }
    }
    const uint32_t len = slot_len(c, slot);
    cpx*           x   = buf.data() + max_d;
    if (ndl > 0) {
      modulate(s, slot, x);
    } else {
      std::fill(x, x + len, cpx(0, 0));
    }
    if (ssb_now) {
      srsran_mib_nr_t mib      = {};
      mib.sfn                  = sfn;
      mib.ssb_idx              = 0;
      mib.hrf                  = false;
      mib.scs_common           = c.scs;
      mib.ssb_offset           = c.k_ssb;
      mib.dmrs_typeA_pos       = c.typeA_pos;
      mib.coreset0_idx         = c.cs0_idx;
      mib.ss0_idx              = c.ss0_idx;
      mib.cell_barred          = false;
      mib.intra_freq_reselection = true;
      srsran_pbch_msg_nr_t pbch = {};
      if (srsran_pbch_msg_nr_mib_pack(&mib, &pbch) < SRSRAN_SUCCESS ||
          srsran_ssb_add(&s.ssb, c.pci, &pbch, (cf_t*)x, (cf_t*)x) < SRSRAN_SUCCESS) {
        die("SSB");
      }
    }
    // channel: multipath, CFO, AWGN
    for (uint32_t n = 0; n < len; n++) {
      cpx y = 0;
      if (taps.empty()) {
        y = x[n];
      } else {
        for (auto& tp : taps) {
          y += tp.second * x[(int)n - (int)tp.first];
        }
      }
      if (cfo != 0) {
        double ph = 2 * M_PI * std::fmod(cfo * (double)(n_abs + n) / c.srate, 1.0);
        y *= cpx((float)std::cos(ph), (float)std::sin(ph));
      }
      if (sigma > 0) {
        y += cpx((float)sigma * gauss(nrng), (float)sigma * gauss(nrng));
      }
      out[n] = y * gain;
    }
    if (max_d) {
      std::copy(x + len - max_d, x + len, buf.begin());
    }
    n_abs += len;
    uint32_t skip = (uint32_t)std::min<uint64_t>(to_drop, len);
    to_drop -= skip;
    if (len > skip) {
      if (fwrite(out.data() + skip, sizeof(cpx), len - skip, fo) != len - skip) {
        fprintf(stderr, "gs_cellsim: write stopped at %.2f s (reader gone?)\n", tnow);
        break;
      }
      samples_written += len - skip;
    }
  }
  fclose(fo);
  fclose(t.f);
  fclose(t.ev);
  printf("done: %.2f s, %lu DCIs placed (+%lu SIB1), %lu blocked (no free candidate)\n",
         samples_written / c.srate, (unsigned long)s.nof_dci, (unsigned long)s.nof_sib1, (unsigned long)s.nof_blocked);
  srsran_pdcch_nr_free(&s.pdcch);
  srsran_pdsch_nr_free(&s.pdsch);
  srsran_dmrs_sch_free(&s.dmrs_sch);
  srsran_softbuffer_tx_free(&s.sb);
  srsran_ssb_free(&s.ssb);
  srsran_dft_plan_free(&s.ifft);
  return 0;
}
