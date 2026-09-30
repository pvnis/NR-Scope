#include "nrscope/hdr/dmrs_check.h"

#include <cmath>
#include <complex>
#include <vector>

namespace {

const uint32_t NSYMB_PER_SLOT = 14;
const uint32_t NC             = 1600; // TS 38.211 5.2.1

/* TS 38.211 5.2.1 pseudo-random sequence, len bits from c_init. */
std::vector<uint8_t> gold_sequence(uint32_t c_init, uint32_t len)
{
  const uint32_t       total = NC + len;
  std::vector<uint8_t> x1(total + 31, 0), x2(total + 31, 0);
  x1[0] = 1;
  for (uint32_t i = 0; i < 31; i++) {
    x2[i] = (c_init >> i) & 1;
  }
  for (uint32_t n = 0; n < total; n++) {
    x1[n + 31] = (x1[n + 3] + x1[n]) & 1;
    x2[n + 31] = (x2[n + 3] + x2[n + 2] + x2[n + 1] + x2[n]) & 1;
  }
  std::vector<uint8_t> c(len);
  for (uint32_t n = 0; n < len; n++) {
    c[n] = (x1[n + NC] + x2[n + NC]) & 1;
  }
  return c;
}

/* TS 38.211 7.4.1.1.1, n_SCID = 0 (CDM group 0, lambda-bar = 0). */
uint32_t dmrs_c_init(uint32_t slot, uint32_t symbol, uint32_t n_id, uint32_t n_scid)
{
  uint64_t v = ((uint64_t)(NSYMB_PER_SLOT * slot + symbol + 1) * (2ULL * n_id + 1ULL)) << 17;
  v += 2ULL * n_id + n_scid;
  return (uint32_t)(v & 0x7fffffffULL);
}

struct Pilot {
  uint32_t            k_abs; // subcarrier index from subcarrier 0 of CRB 0 (point A)
  std::complex<float> y;     // received RE
};

/* The grant's pilots of DM-RS configuration type 1, CDM group 0: the even
  subcarriers of each allocated PRB (k = 4n + 2k', delta = 0). */
std::vector<Pilot> grant_pilots(const cf_t*                  grid,
                                uint32_t                     nof_carrier_prb,
                                uint32_t                     crb_offset,
                                uint32_t                     bwp_start_crb,
                                const srsran_sch_grant_nr_t& grant,
                                uint32_t                     symbol)
{
  std::vector<Pilot> p;
  const uint32_t     nof_re = nof_carrier_prb * SRSRAN_NRE;
  for (uint32_t prb = 0; prb < SRSRAN_MAX_PRB_NR; prb++) {
    if (!grant.prb_idx[prb]) {
      continue;
    }
    const uint32_t crb = bwp_start_crb + prb;
    if (crb < crb_offset || crb - crb_offset >= nof_carrier_prb) {
      continue; // outside the carrier grid
    }
    for (uint32_t sc = 0; sc < SRSRAN_NRE; sc += 2) {
      const uint32_t k_grid = (crb - crb_offset) * SRSRAN_NRE + sc;
      p.push_back({crb * SRSRAN_NRE + sc, grid[symbol * nof_re + k_grid]}); // cf_t is std::complex<float> in C++
    }
  }
  return p;
}

struct Coherence {
  double   num   = 0; // sum over symbols of |sum_pairs z_m conj(z_{m+2})|
  double   den   = 0; // sum over paired pilots of (|z_m|^2 + |z_{m+2}|^2) / 2
  uint32_t pairs = 0;
};

/* Adds one symbol's pilots, divided by the pilots regenerated for c_init. */
void accumulate(Coherence& acc, const std::vector<Pilot>& pilots, uint32_t c_init)
{
  if (pilots.size() < 3) {
    return;
  }
  // r(m) uses bits 2m and 2m+1, with m = k_abs / 2 for CDM group 0
  const uint32_t             m_max = pilots.back().k_abs / 2;
  const std::vector<uint8_t> c     = gold_sequence(c_init, 2 * m_max + 2);
  const float                a     = (float)M_SQRT1_2;

  std::vector<std::complex<float> > z(pilots.size());
  for (size_t i = 0; i < pilots.size(); i++) {
    const uint32_t            m = pilots[i].k_abs / 2;
    const std::complex<float> r(a * (1 - 2 * (int)c[2 * m]), a * (1 - 2 * (int)c[2 * m + 1]));
    z[i] = pilots[i].y * std::conj(r); // |r| = 1
  }

  std::complex<double> s = 0;
  for (size_t i = 0; i + 2 < pilots.size(); i++) {
    if (pilots[i + 2].k_abs - pilots[i].k_abs != 4) {
      continue; // not neighbours on the same OCC weight (gap in the allocation)
    }
    s += std::complex<double>(z[i] * std::conj(z[i + 2]));
    acc.den += 0.5 * ((double)std::norm(z[i]) + (double)std::norm(z[i + 2]));
    acc.pairs++;
  }
  acc.num += std::abs(s);
}

float ratio(const Coherence& c)
{
  return c.den > 0 ? (float)(c.num / c.den) : 0.0f;
}

} // namespace

DmrsCheckResult dmrs_check_pdsch(const cf_t*                  grid,
                                 uint32_t                     nof_carrier_prb,
                                 uint32_t                     crb_offset,
                                 uint32_t                     bwp_start_crb,
                                 const srsran_sch_grant_nr_t& grant,
                                 const uint32_t*              dmrs_symbols,
                                 uint32_t                     nof_dmrs_symbols,
                                 uint32_t                     n_id,
                                 uint32_t                     slot_idx_in_frame)
{
  DmrsCheckResult res;
  if (grid == nullptr || nof_dmrs_symbols == 0) {
    return res;
  }
  const uint32_t n_scid = grant.n_scid ? 1 : 0;

  Coherence claimed, wrong_nid, data_sym;
  double    energy = 0;
  uint32_t  nof_re = 0;
  for (uint32_t i = 0; i < nof_dmrs_symbols; i++) {
    const uint32_t     l      = dmrs_symbols[i];
    std::vector<Pilot> pilots = grant_pilots(grid, nof_carrier_prb, crb_offset, bwp_start_crb, grant, l);
    accumulate(claimed, pilots, dmrs_c_init(slot_idx_in_frame, l, n_id, n_scid));
    accumulate(wrong_nid, pilots, dmrs_c_init(slot_idx_in_frame, l, n_id + 1, n_scid));
    for (const Pilot& p : pilots) {
      energy += std::norm(p.y);
      nof_re++;
    }
  }

  /* Control on the data symbols: every PDSCH symbol that is not DM-RS, tested
    with the sequence DM-RS would have there. They carry data, not pilots. All of
    them, so a small grant still has enough pairs for the control to mean
    something. */
  for (uint32_t l = grant.S; l < grant.S + grant.L && l < NSYMB_PER_SLOT; l++) {
    bool is_dmrs = false;
    for (uint32_t i = 0; i < nof_dmrs_symbols; i++) {
      is_dmrs |= (dmrs_symbols[i] == l);
    }
    if (!is_dmrs) {
      std::vector<Pilot> pilots = grant_pilots(grid, nof_carrier_prb, crb_offset, bwp_start_crb, grant, l);
      accumulate(data_sym, pilots, dmrs_c_init(slot_idx_in_frame, l, n_id, n_scid));
    }
  }

  if (claimed.pairs == 0) {
    return res;
  }
  res.valid                 = true;
  res.nof_pilots            = nof_re;
  res.coherence             = ratio(claimed);
  res.coherence_wrong_nid   = ratio(wrong_nid);
  res.coherence_data_symbol = ratio(data_sym);
  res.pilot_epre_db         = 10.0f * log10f((float)(energy / nof_re) + 1e-30f);
  /* For a channel flat over 4 subcarriers, coherence = S / (S + N) */
  const float rho = std::min(res.coherence, 0.9999f);
  res.snr_db      = 10.0f * log10f(rho / (1.0f - rho) + 1e-30f);
  return res;
}
