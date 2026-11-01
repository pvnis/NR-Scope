#include "nrscope/hdr/ul_dmrs.h"
#include "nrscope/hdr/sensing/nr_ue_dmrs_despread.h"

#include <atomic>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

bool nrscope_ul_dmrs = false;

namespace {

struct Pending {
  uint16_t            rnti;
  srsran_sch_cfg_nr_t cfg;
  uint16_t            ports; // DM-RS port mask
};

std::mutex                         store_mx;
std::multimap<uint32_t, Pending>   store;      // by absolute slot, modulo the SFN period
FILE*                              out = nullptr;
std::mutex                         out_mx;
std::atomic<uint64_t>              n_found{0}, n_done{0}, n_missed{0}, n_skipped{0};
uint64_t                           wraps = 0;  // SFN periods seen, for an unwrapped slot number
uint32_t                           last_key = 0;

/* DCI 0_1 antenna ports, one layer, DM-RS type 1, one symbol, transform precoding off
  (TS 38.212 Table 7.3.1.1.2-8): the port and the CDM groups without data. */
bool ul_port_type1_len1(uint32_t value, uint16_t* ports, uint32_t* cdm)
{
  static const uint8_t port[] = {0, 1, 0, 1, 2, 3}, groups[] = {1, 1, 2, 2, 2, 2};
  if (value >= sizeof(port)) {
    return false;
  }
  *ports = (uint16_t)(1u << port[value]);
  *cdm   = groups[value];
  return true;
}

void report()
{
  printf("ul_dmrs: %lu uplink grants found, %lu estimated, %lu missed (slot gone first), %lu skipped "
         "(unsupported DM-RS)\n",
         (unsigned long)n_found.load(), (unsigned long)n_done.load(), (unsigned long)n_missed.load(),
         (unsigned long)n_skipped.load());
  if (out) {
    fclose(out);
  }
}

} // namespace

void ul_dmrs_open(const char* map_dump)
{
  std::string p = map_dump ? map_dump : "map2d.csv";
  const size_t dot = p.rfind(".csv");
  p = (dot == std::string::npos ? p : p.substr(0, dot)) + "_ul.csv";
  out = fopen(p.c_str(), "w");
  if (!out) {
    printf("ul_dmrs: cannot open %s\n", p.c_str());
    return;
  }
  // R: the chains' covariance, upper triangle by rows (r00, r01, .. r33), real and imaginary parts
  fprintf(out, "sfn,slot,slot_abs,rnti,prb_first,nof_prb,nof_symb,nof_re,power_db");
  for (int m = 0; m < 4; m++)
    for (int n = m; n < 4; n++)
      fprintf(out, ",r%d%d_re,r%d%d_im", m, n, m, n);
  fprintf(out, "\n");
  printf("ul_dmrs: uplink DM-RS covariances to %s\n", p.c_str());
  atexit(report);
}

void ul_dmrs_add(uint32_t sfn, uint32_t slot_idx, uint32_t slots_per_frame, uint16_t rnti,
                 const srsran_sch_cfg_nr_t& cfg, uint32_t dci_ports)
{
  n_found++;
  uint16_t ports = 0;
  uint32_t cdm   = 0;
  if (cfg.enable_transform_precoder || cfg.grant.nof_layers != 1 || cfg.dmrs.type != srsran_dmrs_sch_type_1 ||
      cfg.dmrs.length != srsran_dmrs_sch_len_1 || !ul_port_type1_len1(dci_ports, &ports, &cdm) ||
      cdm != cfg.grant.nof_dmrs_cdm_groups_without_data) {
    n_skipped++;
    return;
  }
  const uint32_t period = 1024 * slots_per_frame;
  const uint32_t key    = (sfn * slots_per_frame + slot_idx + cfg.grant.k) % period;
  std::lock_guard<std::mutex> l(store_mx);
  store.emplace(key, Pending{rnti, cfg, ports});
}

void ul_dmrs_process(DCIDecoder& d, uint32_t sfn, srsran_slot_cfg_t* slot, uint32_t slots_per_frame)
{
  const uint32_t period = 1024 * slots_per_frame;
  const uint32_t key    = (sfn * slots_per_frame + slot->idx) % period;
  std::vector<Pending> due;
  uint64_t             slot_abs;
  {
    std::lock_guard<std::mutex> l(store_mx);
    auto r = store.equal_range(key);
    for (auto it = r.first; it != r.second; ++it)
      due.push_back(it->second);
    store.erase(r.first, r.second);
    // grants whose slot went by unprocessed (a worker took it before the DCI came): drop them
    for (auto it = store.begin(); it != store.end();) {
      const uint32_t ahead = (it->first + period - key) % period;
      if (ahead > period / 2) {
        n_missed++;
        it = store.erase(it);
      } else {
        ++it;
      }
    }
    if (key + period / 2 < last_key) {
      wraps++;
    }
    last_key = key;
    slot_abs = wraps * period + key;
  }
  if (due.empty() || d.dci_cfg.enable_transform_precoding) {
    n_skipped += d.dci_cfg.enable_transform_precoding ? due.size() : 0;
    return;
  }

  const uint32_t n_ant = std::min<uint32_t>(d.ue_dl_grid.nof_rx_antennas, 4);
  const int      n_sc  = (int)(d.grid_carrier.nof_prb * SRSRAN_NRE);
  // the same placement as the PDSCH's: the uplink BWP taken where the downlink one is (TDD)
  const nr_dmrs_placement_t place = {d.crb_offset, d.bwp_start_crb, 0};
  std::vector<cf_t>         H((size_t)n_ant * n_sc);
  uint32_t                  done[SRSRAN_MAX_PORTS] = {};

  for (const Pending& p : due) {
    uint32_t  sym[SRSRAN_DMRS_SCH_MAX_SYMBOLS];
    const int n_sym = srsran_dmrs_sch_get_symbols_idx(&p.cfg.dmrs, &p.cfg.grant, sym);
    nr_dmrs_layout_t lay;
    if (n_sym <= 0 || !nr_ue_dmrs_layout(p.ports, p.cfg.dmrs.type, &lay)) {
      n_skipped++;
      continue;
    }
    uint32_t need[SRSRAN_MAX_PORTS] = {};
    for (int i = 0; i < n_sym; i++)
      for (uint32_t a = 0; a < n_ant; a++)
        need[a] |= 1u << sym[i];
    srsran_ue_dl_nr_estimate_fft_nrscope_symbols(&d.ue_dl_grid, slot, d.arg_scs_grid, need, done);

    double   R_re[4][4] = {}, R_im[4][4] = {};
    uint64_t n_re       = 0;
    for (int i = 0; i < n_sym; i++) {
      const uint32_t cinit = nr_ue_dmrs_seed(&d.carrier_ul, &p.cfg.dmrs, &p.cfg.grant, slot->idx, sym[i]);
      std::unique_ptr<bool[]> valid(new bool[n_sc]());
      for (uint32_t a = 0; a < n_ant; a++) {
        bool* v = a == 0 ? valid.get() : nullptr;   // the same REs on every chain
        nr_ue_dmrs_estimate_symbol(&lay, 0, cinit, &p.cfg.grant, &place,
                                   d.ue_dl_grid.sf_symbols[a] + (size_t)sym[i] * n_sc, n_sc, &H[(size_t)a * n_sc], v);
      }
      for (int k = 0; k < n_sc; k++) {
        if (!valid[k])
          continue;
        n_re++;
        for (uint32_t m = 0; m < n_ant; m++) {
          const cf_t hm = H[(size_t)m * n_sc + k];
          for (uint32_t n = m; n < n_ant; n++) {
            const cf_t  hn = H[(size_t)n * n_sc + k];   // hm conj(hn)
            const float ar = __real__ hm, ai = __imag__ hm, br = __real__ hn, bi = __imag__ hn;
            R_re[m][n] += ar * br + ai * bi;
            R_im[m][n] += ai * br - ar * bi;
          }
        }
      }
    }
    if (n_re == 0) {
      n_skipped++;
      continue;
    }
    double tr = 0;
    for (uint32_t m = 0; m < n_ant; m++)
      tr += R_re[m][m];
    uint32_t prb_first = 0;
    while (prb_first < SRSRAN_MAX_PRB_NR && !p.cfg.grant.prb_idx[prb_first])
      prb_first++;
    std::lock_guard<std::mutex> l(out_mx);
    if (out) {
      fprintf(out, "%u,%u,%lu,%u,%u,%u,%u,%lu,%.2f", sfn, slot->idx, (unsigned long)slot_abs, p.rnti, prb_first,
              p.cfg.grant.nof_prb, p.cfg.grant.L, (unsigned long)n_re, 10 * log10(tr / n_re / n_ant + 1e-30));
      for (int m = 0; m < 4; m++)
        for (int n = m; n < 4; n++)
          fprintf(out, ",%.4g,%.4g", R_re[m][n] / n_re, R_im[m][n] / n_re);
      fprintf(out, "\n");
      fflush(out);
    }
    n_done++;
  }
}
