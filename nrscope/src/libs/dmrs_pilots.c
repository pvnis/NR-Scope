#include "nrscope/hdr/dmrs_pilots.h"

#include <math.h>
#include <string.h>

#define NC 1600 // TS 38.211 5.2.1
#define NSYMB_PER_SLOT 14

uint32_t nrscope_dmrs_n_id(uint32_t pci, const srsran_dmrs_sch_cfg_t* dmrs_cfg, bool n_scid)
{
  if (!n_scid && dmrs_cfg != NULL && dmrs_cfg->scrambling_id0_present) {
    return dmrs_cfg->scrambling_id0;
  }
  if (n_scid && dmrs_cfg != NULL && dmrs_cfg->scrambling_id1_present) {
    return dmrs_cfg->scrambling_id1;
  }
  return pci;
}

uint32_t nrscope_dmrs_c_init(uint32_t slot_idx, uint32_t symbol_idx, uint32_t n_id, bool n_scid)
{
  uint64_t v = ((uint64_t)(NSYMB_PER_SLOT * slot_idx + symbol_idx + 1) * (2ULL * n_id + 1ULL)) << 17;
  v += 2ULL * n_id + (n_scid ? 1 : 0);
  return (uint32_t)(v & 0x7fffffffULL);
}

void nrscope_dmrs_sequence(uint32_t c_init, uint32_t n, float amplitude, cf_t* r)
{
  if (n > NRSCOPE_DMRS_MAX_PILOTS) {
    n = NRSCOPE_DMRS_MAX_PILOTS;
  }
  /* Gold sequence, TS 38.211 5.2.1: two m-sequences of length 31 run NC steps
    before the first output bit. Bits are kept one per byte; 2n + NC + 31 of them
    is at most ~5000, cheap next to the FFT of the same symbol. */
  const uint32_t total = NC + 2 * n;
  uint8_t        x1[NC + 2 * NRSCOPE_DMRS_MAX_PILOTS + 31];
  uint8_t        x2[NC + 2 * NRSCOPE_DMRS_MAX_PILOTS + 31];
  memset(x1, 0, 31);
  x1[0] = 1;
  for (uint32_t i = 0; i < 31; i++) {
    x2[i] = (c_init >> i) & 1;
  }
  for (uint32_t k = 0; k < total; k++) {
    x1[k + 31] = (x1[k + 3] + x1[k]) & 1;
    x2[k + 31] = (x2[k + 3] + x2[k + 2] + x2[k + 1] + x2[k]) & 1;
  }
  for (uint32_t m = 0; m < n; m++) {
    const uint8_t c0 = (x1[2 * m + NC] + x2[2 * m + NC]) & 1;
    const uint8_t c1 = (x1[2 * m + 1 + NC] + x2[2 * m + 1 + NC]) & 1;
    __real__ r[m]    = amplitude * (1 - 2 * (int)c0);
    __imag__ r[m]    = amplitude * (1 - 2 * (int)c1);
  }
}

uint32_t nrscope_dmrs_k_of_index(uint32_t m, srsran_dmrs_sch_type_t type, uint32_t delta)
{
  const uint32_t n = m / 2, kp = m % 2;
  return (type == srsran_dmrs_sch_type_1) ? 4 * n + 2 * kp + delta : 6 * n + kp + delta;
}

int nrscope_dmrs_index_of_k(uint32_t k, srsran_dmrs_sch_type_t type, uint32_t delta)
{
  if (k < delta) {
    return -1;
  }
  const uint32_t off = k - delta;
  if (type == srsran_dmrs_sch_type_1) {
    return (off % 2 == 0) ? (int)(off / 2) : -1; // k' = (off / 2) % 2, n = off / 4
  }
  const uint32_t n = off / 6, kp = off % 6;
  return (kp < 2) ? (int)(2 * n + kp) : -1;
}
