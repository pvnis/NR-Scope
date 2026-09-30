/* Regression test for dmrs_check_pdsch: NR-Scope's own TS 38.211 DM-RS
 * generator must find the pilots srsRAN's mapper (srsran_dmrs_sch_put_sf) puts
 * in a grid, and must not find them where they are not.
 *
 * Each case builds one slot of a 273-PRB carrier: random QPSK everywhere, the
 * PDSCH DM-RS written by srsRAN over it, then a two-path channel with a timing
 * offset and white noise. The claimed configuration must score high, the
 * controls (a data symbol, the wrong scrambling ID) low. A two-layer case mixes
 * ports 1000 and 1001, which differ only by the frequency OCC.
 *
 * Exit status 0 when every case passes. */
#include "nrscope/hdr/dmrs_check.h"

#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <vector>

extern "C" {
#include "srsran/phy/ch_estimation/dmrs_sch.h"
}

namespace {

const uint32_t NOF_PRB = 273;
const uint32_t NOF_RE  = NOF_PRB * SRSRAN_NRE;
const uint32_t NSYMB   = 14;

/* srsRAN's cf_t is GCC's C _Complex float even in C++; std::complex<float>
  converts from it implicitly but not back. */
cf_t to_cf(std::complex<float> z)
{
  cf_t v;
  __real__ v = z.real();
  __imag__ v = z.imag();
  return v;
}

struct Case {
  const char* name;
  uint32_t    prb_start, nof_prb; // BWP-relative allocation
  uint32_t    S, L;
  uint32_t    slot;
  bool        scrambling_id0_present;
  uint32_t    scrambling_id0;
  uint32_t    nof_layers; // 1 or 2 (ports 1000, 1001)
  float       snr_db;
  float       min_coherence; // claimed must reach this
  float       max_control;   // both controls must stay under this
};

bool run(const Case& c, std::mt19937& rng)
{
  const uint32_t pci = 1;

  srsran_carrier_nr_t carrier = SRSRAN_DEFAULT_CARRIER_NR;
  carrier.pci                 = pci;
  carrier.nof_prb             = NOF_PRB;
  carrier.scs                 = srsran_subcarrier_spacing_30kHz;

  srsran_sch_cfg_nr_t cfg          = {};
  cfg.dmrs.type                    = srsran_dmrs_sch_type_1;
  cfg.dmrs.typeA_pos               = srsran_dmrs_sch_typeA_pos_2;
  cfg.dmrs.additional_pos          = srsran_dmrs_sch_add_pos_2;
  cfg.dmrs.length                  = srsran_dmrs_sch_len_1;
  cfg.dmrs.scrambling_id0_present  = c.scrambling_id0_present;
  cfg.dmrs.scrambling_id0          = c.scrambling_id0;
  srsran_sch_grant_nr_t& g         = cfg.grant;
  g.mapping                        = srsran_sch_mapping_type_A;
  g.S                              = c.S;
  g.L                              = c.L;
  g.nof_prb                        = c.nof_prb;
  g.nof_layers                     = 1;
  g.nof_dmrs_cdm_groups_without_data = 1;
  g.n_scid                         = false;
  g.beta_dmrs                      = 1.0f;
  for (uint32_t p = c.prb_start; p < c.prb_start + c.nof_prb; p++) {
    g.prb_idx[p] = true;
  }

  // Random QPSK everywhere, then srsRAN's DM-RS over it (port 1000)
  std::vector<cf_t>                     grid(NOF_RE * NSYMB);
  std::uniform_int_distribution<int>    bit(0, 1);
  const float                           a = (float)M_SQRT1_2;
  for (cf_t& v : grid) {
    v = to_cf({bit(rng) ? a : -a, bit(rng) ? a : -a});
  }
  std::vector<cf_t> tx = grid;
  srsran_dmrs_sch_t dmrs = {};
  srsran_slot_cfg_t slot = {};
  slot.idx               = c.slot;
  if (srsran_dmrs_sch_init(&dmrs, false) < SRSRAN_SUCCESS || srsran_dmrs_sch_set_carrier(&dmrs, &carrier) < SRSRAN_SUCCESS ||
      srsran_dmrs_sch_put_sf(&dmrs, &slot, &cfg, &g, tx.data()) < SRSRAN_SUCCESS) {
    printf("%-34s FAIL: srsRAN could not place the DM-RS\n", c.name);
    return false;
  }
  srsran_dmrs_sch_free(&dmrs);

  uint32_t symbols[SRSRAN_DMRS_SCH_MAX_SYMBOLS] = {};
  const int nof_symbols = srsran_dmrs_sch_get_symbols_idx(&cfg.dmrs, &g, symbols);

  // Two-path channel with a timing offset, per layer; port 1001 = port 1000 x OCC [+1 -1] on alternate pilots
  std::normal_distribution<float>       gauss(0.0f, 1.0f);
  std::uniform_real_distribution<float> phase(0.0f, 2.0f * (float)M_PI);
  const float                           scs = 30e3f, delay0 = 0.3e-6f, delay1 = 1.1e-6f;
  std::complex<float>                   h0[2], h1[2];
  for (int l = 0; l < 2; l++) {
    h0[l] = std::polar(1.0f, phase(rng));
    h1[l] = std::polar(0.4f, phase(rng));
  }
  const float noise_std = sqrtf(powf(10.0f, -c.snr_db / 10.0f) / 2.0f) * sqrtf((float)c.nof_layers);
  std::vector<cf_t> rx(grid.size());
  for (uint32_t l = 0; l < NSYMB; l++) {
    for (uint32_t k = 0; k < NOF_RE; k++) {
      std::complex<float> x = tx[l * NOF_RE + k];
      std::complex<float> y = 0;
      for (uint32_t layer = 0; layer < c.nof_layers; layer++) {
        const float         f = scs * (float)k;
        std::complex<float> h = h0[layer] * std::polar(1.0f, -2.0f * (float)M_PI * f * delay0) +
                                h1[layer] * std::polar(1.0f, -2.0f * (float)M_PI * f * delay1);
        std::complex<float> xl = x;
        if (layer == 1) {
          // Port 1001: same pilots with w_f(k') = -1 for k' = 1; data REs get independent symbols
          bool is_dmrs = false;
          for (int i = 0; i < nof_symbols; i++) {
            is_dmrs |= symbols[i] == l;
          }
          if (is_dmrs && (k % 2 == 0)) {
            xl = ((k / 2) % 2 == 1) ? -x : x;
          } else {
            xl = std::complex<float>(bit(rng) ? a : -a, bit(rng) ? a : -a);
          }
        }
        y += h * xl;
      }
      y += std::complex<float>(noise_std * gauss(rng), noise_std * gauss(rng));
      rx[l * NOF_RE + k] = to_cf(y);
    }
  }

  const uint32_t   n_id = c.scrambling_id0_present ? c.scrambling_id0 : pci;
  DmrsCheckResult  r    = dmrs_check_pdsch(rx.data(), NOF_PRB, 0, 0, g, symbols, (uint32_t)nof_symbols, n_id, c.slot);
  const bool ok = r.valid && r.coherence >= c.min_coherence && r.coherence_data_symbol <= c.max_control &&
                  r.coherence_wrong_nid <= c.max_control;
  printf("%-34s %s  coherence %.3f (>= %.2f)  data-symbol %.3f  wrong-N_ID %.3f (<= %.2f)  SNR %+.1f dB, %u pilots\n",
         c.name, ok ? "pass" : "FAIL", r.coherence, c.min_coherence, r.coherence_data_symbol, r.coherence_wrong_nid,
         c.max_control, r.snr_db, r.nof_pilots);
  return ok;
}

} // namespace

int main()
{
  std::mt19937 rng(7);
  const Case cases[] = {
      {"1 layer, 2 PRB, 13 symbols", 0, 2, 1, 13, 6, false, 0, 1, 20.0f, 0.80f, 0.60f},
      {"1 layer, 50 PRB, 13 symbols", 10, 50, 1, 13, 15, false, 0, 1, 20.0f, 0.95f, 0.20f},
      {"1 layer, 50 PRB, 5 symbols", 10, 50, 1, 5, 3, false, 0, 1, 20.0f, 0.95f, 0.25f},
      {"1 layer, scramblingID0 = 517", 100, 40, 1, 13, 9, true, 517, 1, 20.0f, 0.95f, 0.20f},
      {"1 layer, top of carrier", 240, 33, 1, 13, 19, false, 0, 1, 20.0f, 0.95f, 0.20f},
      {"2 layers, 50 PRB, 13 symbols", 0, 50, 1, 13, 6, false, 0, 2, 20.0f, 0.90f, 0.20f},
      {"1 layer, 50 PRB, 0 dB SNR", 10, 50, 1, 13, 4, false, 0, 1, 0.0f, 0.40f, 0.20f},
  };
  bool all = true;
  for (const Case& c : cases) {
    all &= run(c, rng);
  }
  printf("%s\n", all ? "all DM-RS check cases passed" : "DM-RS check FAILED");
  return all ? 0 : 1;
}
