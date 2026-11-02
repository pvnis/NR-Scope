/*
 * SPDX-License-Identifier: MIT
 *
 * GoldenSniffer port: Gold sequence inversion of the PDCCH DM-RS.
 * Ports scrambling.m, scrambling_x1.m, DMRS_process.m, c_init_PDCCH_DMRS_decode.m of
 * https://github.com/5G-GoldenSniffer/5G-GoldenSniffer (Copyright 2024-2026 the
 * authors, MIT License), and replaces the encrypted PDCCH_DMRS_solve_underdet.p.
 */
#include "gs_internal.h"
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#define GS_NC 1600U

/* Longest sequence offset needed: point A may be offsetToCarrier (<= 2199) RBs below
 * the grid, plus a full 275 RB carrier; 6 sequence bits per RB, then the run itself. */
#define GS_SEQ_LEN (GS_NC + 6U * (2199U + 275U) + 2U * (GS_MAX_SC / 4U) + 128U)
#define GS_SEQ_WORDS ((GS_SEQ_LEN + 63U) / 64U)

static uint8_t  gs_x1[GS_SEQ_LEN];                /* x1(n), one bit per byte */
static uint64_t gs_x1_packed[GS_SEQ_WORDS];
static uint64_t gs_x2_unit[31][GS_SEQ_WORDS];     /* x2(n) for c_init = 2^b */
static pthread_once_t gs_gold_once = PTHREAD_ONCE_INIT;

static void gs_pack_set(uint64_t* w, uint32_t i, uint8_t b)
{
  if (b) {
    w[i >> 6] |= 1ULL << (i & 63U);
  }
}

static void gs_gold_build(void)
{
  memset(gs_x1, 0, sizeof(gs_x1));
  gs_x1[0] = 1;
  for (uint32_t n = 0; n + 31 < GS_SEQ_LEN; n++) {
    gs_x1[n + 31] = gs_x1[n + 3] ^ gs_x1[n];
  }
  memset(gs_x1_packed, 0, sizeof(gs_x1_packed));
  for (uint32_t n = 0; n < GS_SEQ_LEN; n++) {
    gs_pack_set(gs_x1_packed, n, gs_x1[n]);
  }

  static uint8_t x2[GS_SEQ_LEN];
  memset(gs_x2_unit, 0, sizeof(gs_x2_unit));
  for (uint32_t b = 0; b < 31; b++) {
    memset(x2, 0, sizeof(x2));
    x2[b] = 1;
    for (uint32_t n = 0; n + 31 < GS_SEQ_LEN; n++) {
      x2[n + 31] = x2[n + 3] ^ x2[n + 2] ^ x2[n + 1] ^ x2[n];
    }
    for (uint32_t n = 0; n < GS_SEQ_LEN; n++) {
      gs_pack_set(gs_x2_unit[b], n, x2[n]);
    }
  }
}

void gs_gold_init(void)
{
  pthread_once(&gs_gold_once, gs_gold_build);
}

/* bits [start, start + len) of a packed sequence, len <= 64 */
static inline uint64_t gs_window(const uint64_t* w, uint32_t start, uint32_t len)
{
  uint32_t i  = start >> 6;
  uint32_t sh = start & 63U;
  uint64_t v  = w[i] >> sh;
  if (sh != 0 && i + 1 < GS_SEQ_WORDS) {
    v |= w[i + 1] << (64U - sh);
  }
  return len >= 64 ? v : (v & ((1ULL << len) - 1ULL));
}

void gs_gold_bits(uint32_t c_init, uint32_t offset, uint32_t len, uint8_t* out)
{
  gs_gold_init();
  for (uint32_t m = 0; m < len; m++) {
    uint32_t n = GS_NC + offset + m;
    uint8_t  v = gs_x1[n];
    for (uint32_t b = 0; b < 31; b++) {
      if ((c_init >> b) & 1U) {
        v ^= (uint8_t)((gs_x2_unit[b][n >> 6] >> (n & 63U)) & 1U);
      }
    }
    out[m] = v;
  }
}

uint32_t gs_pdcch_dmrs_cinit(uint32_t slot_idx, uint32_t l, uint32_t n_id)
{
  uint64_t s = 14ULL * slot_idx + l + 1ULL;
  return (uint32_t)((((s * (2ULL * n_id + 1ULL)) << 17U) + 2ULL * n_id) & 0x7fffffffULL);
}

bool gs_pdcch_cinit_check(uint32_t c_init, uint32_t slot_idx, uint32_t l, uint32_t* n_id)
{
  if (c_init & 1U) {
    return false;
  }
  uint32_t n = (c_init >> 1) & 0xffffU;
  if (c_init != gs_pdcch_dmrs_cinit(slot_idx, l, n)) {
    return false;
  }
  *n_id = n;
  return true;
}

/* hard bits of c * exp(j(pi/4 + rot pi/2)): c(2m) = Re < 0, c(2m+1) = Im < 0 */
static void gs_qpsk_bits(const cf_t* c, uint32_t n, uint32_t rot, uint8_t* bits)
{
  const cf_t r = cexpf(I * (float)(M_PI / 4.0 + rot * M_PI / 2.0));
  for (uint32_t m = 0; m < n; m++) {
    cf_t v          = c[m] * r;
    bits[2 * m]     = crealf(v) < 0.0f;
    bits[2 * m + 1] = cimagf(v) < 0.0f;
  }
}

int gs_dmrs_process(const cf_t* c_est, uint32_t n, uint32_t D, uint32_t c_init[4], uint32_t rot[4])
{
  gs_gold_init();
  const uint32_t M = 2 * n;
  if (M <= 31 || GS_NC + D + M > GS_SEQ_LEN) {
    return 0;
  }
  uint8_t* x2 = malloc(GS_NC + D + M);
  uint8_t  c[2 * GS_MAX_SC / 4];
  if (x2 == NULL) {
    return 0;
  }
  int found = 0;
  for (uint32_t r = 0; r < 4; r++) {
    gs_qpsk_bits(c_est, n, r, c);
    uint32_t base = GS_NC + D;
    for (uint32_t m = 0; m < M; m++) {
      x2[base + m] = c[m] ^ gs_x1[base + m]; /* xor with the known x1 */
    }
    /* does x2 satisfy its recursion? */
    bool ok = true;
    for (uint32_t k = base + M - 32 + 1; k-- > base;) {
      if (x2[k] != (x2[k + 31] ^ x2[k + 3] ^ x2[k + 2] ^ x2[k + 1])) {
        ok = false;
        break;
      }
    }
    if (!ok) {
      continue;
    }
    /* run it backwards to the initial state */
    for (uint32_t k = base; k-- > 0;) {
      x2[k] = x2[k + 31] ^ x2[k + 3] ^ x2[k + 2] ^ x2[k + 1];
    }
    uint32_t ci = 0;
    for (uint32_t b = 0; b < 31; b++) {
      ci |= (uint32_t)x2[b] << b;
    }
    c_init[found] = ci;
    rot[found]    = r;
    found++;
  }
  free(x2);
  return found;
}

/*
 * Replacement for PDCCH_DMRS_solve_underdet.p.
 *
 * A short run gives fewer than 31 bits per symbol, too few to invert x2 alone. But
 * c_init = 2^17 * ((s_l (2N+1)) mod 2^14) + 2N, s_l = 14 slot + l + 1, so all
 * symbols share N, and the 14 upper bits only depend on N mod 2^13. For each of the
 * 8192 values of N mod 2^13 the upper bits and c_init bits 1..13 are known; the last
 * three bits of N (c_init bits 14..16) enter x2 linearly, so they are matched against
 * the 8 combinations of their x2 columns. All symbols must agree on one N and one of
 * the 4 rotations of the common QPSK ambiguity.
 */
bool gs_dmrs_solve_underdet(const cf_t* c_est, uint32_t n, uint32_t dur, uint32_t D, uint32_t slot_idx,
                            uint32_t* n_id, uint32_t* rot)
{
  gs_gold_init();
  const uint32_t M = 2 * n;
  if (M > 64 || M < 8 || dur < 1 || dur > GS_MAX_DURATION || GS_NC + D + M > GS_SEQ_LEN) {
    return false;
  }
  const uint32_t start = GS_NC + D;
  const uint64_t x1w   = gs_window(gs_x1_packed, start, M);

  /* x2 window as a linear function of c_init: 4 byte-indexed tables */
  uint64_t col[31];
  for (uint32_t b = 0; b < 31; b++) {
    col[b] = gs_window(gs_x2_unit[b], start, M);
  }
  uint64_t lut[4][256];
  for (uint32_t t = 0; t < 4; t++) {
    for (uint32_t v = 0; v < 256; v++) {
      uint64_t acc = 0;
      for (uint32_t b = 0; b < 8; b++) {
        uint32_t bit = 8 * t + b;
        if (bit < 31 && ((v >> b) & 1U)) {
          acc ^= col[bit];
        }
      }
      lut[t][v] = acc;
    }
  }
  uint64_t comb[8];
  for (uint32_t h = 0; h < 8; h++) {
    comb[h] = ((h & 1U) ? col[14] : 0) ^ ((h & 2U) ? col[15] : 0) ^ ((h & 4U) ? col[16] : 0);
  }

  /* observations per rotation and symbol, x1 removed */
  uint64_t obs[4][GS_MAX_DURATION];
  uint8_t  bits[64];
  for (uint32_t r = 0; r < 4; r++) {
    for (uint32_t l = 0; l < dur; l++) {
      gs_qpsk_bits(c_est + l * n, n, r, bits);
      uint64_t w = 0;
      for (uint32_t m = 0; m < M; m++) {
        w |= (uint64_t)bits[m] << m;
      }
      obs[r][l] = w ^ x1w;
    }
  }

  uint32_t nof_sol = 0, sol_n = 0, sol_r = 0;
  for (uint32_t nlo = 0; nlo < 8192; nlo++) {
    uint64_t part[GS_MAX_DURATION];
    for (uint32_t l = 0; l < dur; l++) {
      uint32_t s  = 14 * slot_idx + l + 1;
      uint32_t hi = (s * (2 * nlo + 1)) & 0x3fffU;
      uint32_t c  = (hi << 17) | (nlo << 1);
      part[l]     = lut[0][c & 0xff] ^ lut[1][(c >> 8) & 0xff] ^ lut[2][(c >> 16) & 0xff] ^ lut[3][(c >> 24) & 0x7f];
    }
    for (uint32_t r = 0; r < 4; r++) {
      uint64_t res0 = obs[r][0] ^ part[0];
      for (uint32_t h = 0; h < 8; h++) {
        if (res0 != comb[h]) {
          continue;
        }
        bool ok = true;
        for (uint32_t l = 1; l < dur && ok; l++) {
          ok = (obs[r][l] ^ part[l]) == comb[h];
        }
        if (ok) {
          uint32_t cand = nlo | (h << 13);
          if (nof_sol == 0 || cand != sol_n || r != sol_r) {
            nof_sol++;
            sol_n = cand;
            sol_r = r;
          }
        }
      }
    }
  }
  if (nof_sol != 1) {
    return false;
  }
  *n_id = sol_n;
  *rot  = sol_r;
  return true;
}
