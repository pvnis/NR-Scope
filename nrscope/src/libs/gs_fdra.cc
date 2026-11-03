#include "nrscope/hdr/gs_fdra.h"
#include "nrscope/hdr/dmrs_pilots.h"
#include <cmath>
#include <complex>
#include <cstring>

/* RBG size P, TS 38.214 Table 5.1.2.2.1-1 */
static uint32_t gs_rbg_p(uint32_t nbwp, uint32_t cfg)
{
  if (nbwp <= 36) {
    return cfg == 1 ? 2 : 4;
  }
  if (nbwp <= 72) {
    return cfg == 1 ? 4 : 8;
  }
  if (nbwp <= 144) {
    return cfg == 1 ? 8 : 16;
  }
  return 16;
}

static uint32_t gs_type0_width(uint32_t nbwp, uint32_t bwp_start, uint32_t cfg)
{
  const uint32_t P = gs_rbg_p(nbwp, cfg);
  return (nbwp + bwp_start % P + P - 1) / P;
}

static uint32_t gs_type1_width(uint32_t nbwp)
{
  const double n = (double)nbwp * (nbwp + 1) / 2.0;
  return (uint32_t)std::ceil(std::log2(n));
}

uint32_t gs_fdra_hypotheses(uint32_t carrier_start_crb, uint32_t carrier_nprb, uint32_t init_bwp_start_crb,
                            uint32_t init_bwp_nprb, GsFdraHyp* out, uint32_t max_out)
{
  uint32_t n = 0;
  struct {
    uint32_t start, size;
  } bwps[2]     = {{carrier_start_crb, carrier_nprb}, {init_bwp_start_crb, init_bwp_nprb}};
  const uint32_t nb = (init_bwp_nprb == carrier_nprb && init_bwp_start_crb == carrier_start_crb) || init_bwp_nprb == 0
                          ? 1
                          : 2;
  /* the simplest first: on a tie of scores, the first one wins */
  for (uint32_t b = 0; b < nb; b++) {
    for (uint32_t off = 1; off <= 6; off++) {
      for (uint32_t t = 0; t < 5; t++) {
        GsFdraHyp h = {};
        h.off       = (uint8_t)off;
        h.nbwp      = (uint16_t)bwps[b].size;
        h.bwp_start = (uint16_t)bwps[b].start;
        if (t == 0) {
          h.type  = 1;
          h.width = (uint8_t)gs_type1_width(h.nbwp);
        } else if (t <= 2) {
          h.type    = 0;
          h.rbg_cfg = (uint8_t)t;
          h.width   = (uint8_t)gs_type0_width(h.nbwp, h.bwp_start, t);
        } else {
          h.type    = 2;
          h.rbg_cfg = (uint8_t)(t - 2);
          uint32_t w0 = gs_type0_width(h.nbwp, h.bwp_start, h.rbg_cfg), w1 = gs_type1_width(h.nbwp);
          h.width   = (uint8_t)((w0 > w1 ? w0 : w1) + 1);
        }
        if (n < max_out) {
          out[n++] = h;
        }
      }
    }
  }
  return n;
}

static uint32_t gs_bits(const uint8_t* p, uint32_t from, uint32_t n)
{
  uint32_t v = 0;
  for (uint32_t i = 0; i < n; i++) {
    v = (v << 1) | (p[from + i] & 1);
  }
  return v;
}

static bool gs_parse_type1(uint32_t riv, uint32_t N, bool* prb, uint32_t* nof)
{
  /* TS 38.214 5.1.2.2.2: RIV = N (L - 1) + start if L - 1 <= floor(N/2), else
     N (N - L + 1) + (N - 1 - start) */
  uint32_t L = riv / N + 1, start = riv % N;
  if (start + L > N) {
    L     = N - riv / N + 1;
    start = N - 1 - riv % N;
  }
  if (L < 1 || start + L > N || riv >= N * (N + 1) / 2) {
    return false;
  }
  for (uint32_t i = 0; i < L; i++) {
    prb[start + i] = true;
  }
  *nof = L;
  return true;
}

static bool gs_parse_type0(const uint8_t* bits, uint32_t width, uint32_t N, uint32_t bwp_start, uint32_t cfg, bool* prb,
                           uint32_t* nof)
{
  const uint32_t P = gs_rbg_p(N, cfg), first = P - bwp_start % P;
  uint32_t       n = 0;
  for (uint32_t g = 0; g < width; g++) {
    if (!bits[g]) {
      continue;
    }
    const uint32_t lo = g == 0 ? 0 : first + (g - 1) * P;
    const uint32_t hi = g == 0 ? (first < N ? first : N) : (lo + P < N ? lo + P : N);
    for (uint32_t k = lo; k < hi; k++) {
      prb[k] = true;
      n++;
    }
  }
  *nof = n;
  return n > 0;
}

bool gs_fdra_parse(const GsFdraHyp& h, const uint8_t* payload, uint32_t nof_bits, bool* prb, uint32_t* nof_prb)
{
  if ((uint32_t)h.off + h.width > nof_bits || h.nbwp == 0) {
    return false;
  }
  memset(prb, 0, sizeof(bool) * h.nbwp);
  *nof_prb = 0;
  const uint8_t* f = payload + h.off;
  if (h.type == 1) {
    return gs_parse_type1(gs_bits(f, 0, h.width), h.nbwp, prb, nof_prb);
  }
  if (h.type == 0) {
    return gs_parse_type0(f, h.width, h.nbwp, h.bwp_start, h.rbg_cfg, prb, nof_prb);
  }
  /* dynamic switch: MSB 0 = type 0 on the LSBs, 1 = type 1 on the LSBs */
  const uint32_t w0 = gs_type0_width(h.nbwp, h.bwp_start, h.rbg_cfg), w1 = gs_type1_width(h.nbwp);
  if (f[0] == 0) {
    return gs_parse_type0(f + h.width - w0, w0, h.nbwp, h.bwp_start, h.rbg_cfg, prb, nof_prb);
  }
  return gs_parse_type1(gs_bits(f, h.width - w1, w1), h.nbwp, prb, nof_prb);
}

void gs_fdra_rb_coherence(const cf_t* sym, uint32_t first_crb, uint32_t nof_prb, uint32_t slot_idx, uint32_t l,
                          uint32_t n_id, cf_t* s, float* d)
{
  /* the pilots of the whole carrier from point A: m = k / 2 (type 1, CDM group 0) */
  static thread_local cf_t r[NRSCOPE_DMRS_MAX_PILOTS];
  const uint32_t           m_max = (first_crb + nof_prb) * 6;
  if (m_max > NRSCOPE_DMRS_MAX_PILOTS) {
    return;
  }
  nrscope_dmrs_sequence(nrscope_dmrs_c_init(slot_idx, l, n_id, false), m_max, (float)M_SQRT1_2, r);
  for (uint32_t c = 0; c < nof_prb; c++) {
    std::complex<float> z[6];
    float               pw = 0;
    for (uint32_t j = 0; j < 6; j++) {
      const uint32_t           k  = (first_crb + c) * 12 + 2 * j; // from point A
      const std::complex<float> y  = std::complex<float>(sym[c * 12 + 2 * j]);
      const std::complex<float> rr = std::complex<float>(r[k / 2]);
      z[j]                         = y * std::conj(rr) * 2.0f; // |r|^2 = 1/2
    }
    std::complex<float> acc = 0;
    for (uint32_t j = 0; j + 2 < 6; j++) { // pairs 4 subcarriers apart: the same OCC weight
      acc += z[j] * std::conj(z[j + 2]);
      pw += 0.5f * (std::norm(z[j]) + std::norm(z[j + 2]));
    }
    __real__ s[c] = acc.real(); // cf_t is the GNU complex type here
    __imag__ s[c] = acc.imag();
    d[c] = pw;
  }
}

void gs_dmrs_shape(const cf_t* sym, uint32_t first_crb, const bool* prb, uint32_t nof_grid_prb, uint32_t slot_idx,
                   uint32_t l, uint32_t n_id, GsDmrsShape* out)
{
  *out = {0, 1, 0, 0};
  static thread_local cf_t r[NRSCOPE_DMRS_MAX_PILOTS];
  const uint32_t           m_max = (first_crb + nof_grid_prb) * 6;
  if (m_max > NRSCOPE_DMRS_MAX_PILOTS) {
    return;
  }
  nrscope_dmrs_sequence(nrscope_dmrs_c_init(slot_idx, l, n_id, false), m_max, (float)M_SQRT1_2, r);
  double              pp = 0, pm = 0, e0 = 0, e1 = 0, d0 = 0, d1 = 0;
  std::complex<float> s0 = 0, s1 = 0;
  for (uint32_t c = 0; c < nof_grid_prb; c++) {
    if (!prb[c]) {
      continue;
    }
    std::complex<float> z0[6], z1[6];
    for (uint32_t j = 0; j < 6; j++) {
      const uint32_t            m  = (first_crb + c) * 6 + j; // the same index on both combs
      const std::complex<float> rr = std::conj(std::complex<float>(r[m])) * 2.0f;
      z0[j]                        = std::complex<float>(sym[c * 12 + 2 * j]) * rr;
      z1[j]                        = std::complex<float>(sym[c * 12 + 2 * j + 1]) * rr;
      e0 += std::norm(z0[j]);
      e1 += std::norm(z1[j]);
    }
    for (uint32_t n = 0; n < 3; n++) { // CDM pairs: k = 4n (k' = 0) and 4n + 2 (k' = 1)
      pp += std::norm(z0[2 * n] + z0[2 * n + 1]);
      pm += std::norm(z0[2 * n] - z0[2 * n + 1]);
    }
    for (uint32_t j = 0; j + 2 < 6; j++) { // 4 subcarriers apart: the same cover weight
      s0 += z0[j] * std::conj(z0[j + 2]);
      s1 += z1[j] * std::conj(z1[j + 2]);
      d0 += 0.5 * (std::norm(z0[j]) + std::norm(z0[j + 2]));
      d1 += 0.5 * (std::norm(z1[j]) + std::norm(z1[j + 2]));
    }
  }
  out->occ  = pp > 0 ? (float)(pm / pp) : 0.0f;
  out->e1   = e0 > 0 ? (float)(e1 / e0) : 1.0f;
  out->coh0 = d0 > 0 ? (float)(std::abs(s0) / d0) : 0.0f;
  out->coh1 = d1 > 0 ? (float)(std::abs(s1) / d1) : 0.0f;
}

void gs_dmrs_decide(const GsDmrsShape& s, uint16_t* ports, bool* cdm2, int* dci_ports)
{
  // a second layer 11 dB under the first, or more, is left out (it is noise to a map)
  *ports     = s.occ < 0.08f ? 0x1 : s.occ > 12.0f ? 0x2 : 0x3;
  *cdm2      = s.e1 < 0.3f || s.coh1 > 0.5f;
  *dci_ports = *cdm2 ? 7 : 2;
}

float gs_fdra_score(const GsFdraHyp& h, const uint8_t* payload, uint32_t nof_bits, uint32_t first_crb,
                    uint32_t nof_grid_prb, const cf_t* s, const float* d, float* coherence)
{
  bool     prb[SRSRAN_MAX_PRB_NR];
  uint32_t n = 0;
  if (!gs_fdra_parse(h, payload, nof_bits, prb, &n) || n == 0) {
    return -1.0f;
  }
  uint32_t            in = 0, occ = 0;
  std::complex<float> sum = 0;
  float               den = 0;
  for (uint32_t i = 0; i < h.nbwp; i++) {
    if (!prb[i]) {
      continue;
    }
    const int32_t c = (int32_t)(h.bwp_start + i) - (int32_t)first_crb;
    if (c < 0 || c >= (int32_t)nof_grid_prb) {
      continue;
    }
    in++;
    const std::complex<float> sc = std::complex<float>(s[c]);
    if (d[c] > 0 && std::abs(sc) / d[c] > 0.5f) {
      occ++;
    }
    sum += sc;
    den += d[c];
  }
  if (coherence != nullptr) {
    *coherence = den > 0 ? std::abs(sum) / den : 0.0f;
  }
  return in > 0 ? (float)occ / (float)in : -1.0f;
}
