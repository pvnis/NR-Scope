/*
 * SPDX-License-Identifier: MIT
 *
 * GoldenSniffer port: DCI decoding by linear algebra over GF(2).
 * Ports make_G_dci.m, make_E_dci.m, DCIdecode_try.m and the dci_encode.m model of
 * https://github.com/5G-GoldenSniffer/5G-GoldenSniffer (Copyright 2024-2026 the
 * authors, MIT License).
 *
 * With the scrambling known, the PDCCH bits are an affine function of the unknowns
 * x = [payload a_0..a_{N-1} | RNTI b_15..b_0]:
 *
 *    b = b0 + x G + scr(N_ID)          (common search space, c_init = N_ID)
 *    b = b0 + x G_B + scr(N_ID)        (n_RNTI = RNTI: the RNTI also enters c_init)
 *
 * b0 is the encoding of all zeros (it carries the 24 ones before the CRC and the x1
 * part of the scrambling). G is built row by row from srsRAN's own encoder chain
 * (the encode path of lib/src/phy/phch/pdcch_nr.c), so it is exact. K = N + 16
 * unknowns, E = 108 L observations: x is read from K pivot bits through the inverse
 * of G restricted to them, and the remaining E - K bits must agree (CRC and polar
 * frozen bits checked at once). The MATLAB keeps an E x E elimination matrix; the
 * pivot form here stores K x E bits instead.
 */
#include "gs_internal.h"
#include "srsran/phy/common/phy_common.h"
#include "srsran/phy/common/sequence.h"
#include "srsran/phy/fec/crc.h"
#include "srsran/phy/fec/polar/polar_chanalloc.h"
#include "srsran/phy/fec/polar/polar_code.h"
#include "srsran/phy/fec/polar/polar_encoder.h"
#include "srsran/phy/fec/polar/polar_interleaver.h"
#include "srsran/phy/fec/polar/polar_rm.h"
#include "srsran/phy/utils/bit.h"
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#define GS_K_MAX (GS_MAX_DCI_BITS + 16)
#define GS_KW ((GS_K_MAX + 63) / 64)
#define GS_E_MAX (108 << GS_MAX_L_LOG2)
#define GS_EW ((GS_E_MAX + 63) / 64)
#define GS_PARITY_MIN 8 /* bits beyond the unknowns left to check (MATLAB: N <= E - 24) */

typedef struct {
  bool      valid;
  uint32_t  K;
  uint64_t* G;                       /* K rows x W words */
  uint16_t  piv[GS_K_MAX];           /* pivot bit positions */
  uint64_t  Pinv[GS_K_MAX][GS_KW];   /* x = ybar * Pinv */
} gs_sys_t;

typedef struct {
  int       ready; /* 0 not built, 1 built (sys[i].valid tells usability) */
  uint64_t* b0;
  gs_sys_t  sys[2];
} gs_entry_t;

static gs_entry_t      gs_tab[GS_MAX_DCI_BITS + 1][GS_MAX_L_LOG2 + 1];
static uint64_t        gs_scr[GS_MAX_L_LOG2 + 1][31][GS_EW]; /* x2 part of the scrambling for c_init = 2^b */
static pthread_mutex_t gs_gf2_mutex = PTHREAD_MUTEX_INITIALIZER;
static int             gs_gf2_inited = 0;

static srsran_polar_code_t    gs_code;
static srsran_polar_encoder_t gs_enc;
static srsran_polar_rm_t      gs_rm;
static srsran_crc_t           gs_crc;

static inline uint32_t gs_words(uint32_t E)
{
  return (E + 63U) / 64U;
}

static inline int gs_bit(const uint64_t* w, uint32_t i)
{
  return (int)((w[i >> 6] >> (i & 63U)) & 1U);
}

static void gs_pack(const uint8_t* bits, uint32_t n, uint64_t* w)
{
  memset(w, 0, gs_words(n) * sizeof(uint64_t));
  for (uint32_t i = 0; i < n; i++) {
    if (bits[i] & 1U) {
      w[i >> 6] |= 1ULL << (i & 63U);
    }
  }
}

/* srsran_pdcch_nr_encode() up to and including the scrambling, without modulation */
static int gs_encode(const uint8_t* payload, uint32_t N, uint16_t rnti, uint32_t cinit, uint32_t E, uint8_t* f)
{
  uint8_t  c_full[24 + GS_K_MAX + 24 + 24];
  uint8_t  c_prime[SRSRAN_POLAR_INTERLEAVER_K_MAX_IL];
  uint8_t  allocated[NMAX];
  uint8_t  d[NMAX];
  uint32_t K = N + 24;

  if (srsran_polar_code_get(&gs_code, (uint16_t)K, (uint16_t)E, 9U) < SRSRAN_SUCCESS) {
    return SRSRAN_ERROR;
  }
  uint8_t* c = c_full;
  srsran_bit_unpack(UINT32_MAX, &c, 24U);      /* 24 ones, c now at offset 24 */
  memcpy(c, payload, N);
  srsran_crc_attach(&gs_crc, c_full, (int)K); /* CRC over ones + payload */

  uint8_t  unpacked_rnti[16];
  uint8_t* ptr = unpacked_rnti;
  srsran_bit_unpack(rnti, &ptr, 16);
  for (uint32_t i = 0; i < 16; i++) {
    c[K - 16 + i] ^= unpacked_rnti[i];
  }
  srsran_polar_interleaver_run_u8(c, c_prime, K, true);
  srsran_polar_chanalloc_tx(c_prime, allocated, gs_code.N, gs_code.K, gs_code.nPC, gs_code.K_set, gs_code.PC_set);
  if (srsran_polar_encoder_encode(&gs_enc, allocated, d, gs_code.n) < SRSRAN_SUCCESS) {
    return SRSRAN_ERROR;
  }
  srsran_polar_rm_tx(&gs_rm, d, f, gs_code.n, E, K, 0);
  srsran_sequence_apply_bit(f, f, E, cinit);
  return SRSRAN_SUCCESS;
}

int gs_gf2_init(void)
{
  pthread_mutex_lock(&gs_gf2_mutex);
  if (gs_gf2_inited) {
    pthread_mutex_unlock(&gs_gf2_mutex);
    return SRSRAN_SUCCESS;
  }
  int ret = SRSRAN_ERROR;
  if (srsran_crc_init(&gs_crc, SRSRAN_LTE_CRC24C, 24) < SRSRAN_SUCCESS ||
      srsran_polar_code_init(&gs_code) < SRSRAN_SUCCESS ||
      srsran_polar_encoder_init(&gs_enc, SRSRAN_POLAR_ENCODER_PIPELINED, NMAX_LOG) < SRSRAN_SUCCESS ||
      srsran_polar_rm_tx_init(&gs_rm) < SRSRAN_SUCCESS) {
    goto out;
  }
  /* scrambling columns: sequence(2^b) + sequence(0) = x2 part, per AL */
  static uint8_t zeros[GS_E_MAX], s0[GS_E_MAX], sb[GS_E_MAX];
  for (uint32_t Ll = 0; Ll <= GS_MAX_L_LOG2; Ll++) {
    uint32_t E = 108U << Ll;
    srsran_sequence_apply_bit(zeros, s0, E, 0);
    for (uint32_t b = 0; b < 31; b++) {
      srsran_sequence_apply_bit(zeros, sb, E, 1U << b);
      for (uint32_t i = 0; i < E; i++) {
        sb[i] ^= s0[i];
      }
      gs_pack(sb, E, gs_scr[Ll][b]);
    }
  }
  gs_gf2_inited = 1;
  ret           = SRSRAN_SUCCESS;
out:
  pthread_mutex_unlock(&gs_gf2_mutex);
  return ret;
}

bool gs_gf2_size_valid(uint32_t nof_bits, uint32_t L_log2)
{
  if (L_log2 > GS_MAX_L_LOG2 || nof_bits < 1 || nof_bits > GS_MAX_DCI_BITS) {
    return false;
  }
  uint32_t E = 108U << L_log2;
  return nof_bits + 16 + GS_PARITY_MIN <= E && nof_bits + 24 < E;
}

/* Pivots of G and the inverse of G restricted to them */
static bool gs_sys_solve_setup(gs_sys_t* s, uint32_t E)
{
  const uint32_t W = gs_words(E);
  const uint32_t K = s->K;
  uint64_t*      work = malloc((size_t)K * W * sizeof(uint64_t));
  if (work == NULL) {
    return false;
  }
  memcpy(work, s->G, (size_t)K * W * sizeof(uint64_t));
  uint32_t rank = 0;
  for (uint32_t col = 0; col < E && rank < K; col++) {
    uint32_t r = rank;
    while (r < K && !gs_bit(&work[(size_t)r * W], col)) {
      r++;
    }
    if (r == K) {
      continue;
    }
    if (r != rank) {
      for (uint32_t w = 0; w < W; w++) {
        uint64_t t                  = work[(size_t)r * W + w];
        work[(size_t)r * W + w]     = work[(size_t)rank * W + w];
        work[(size_t)rank * W + w]  = t;
      }
    }
    for (uint32_t r2 = rank + 1; r2 < K; r2++) {
      if (gs_bit(&work[(size_t)r2 * W], col)) {
        for (uint32_t w = 0; w < W; w++) {
          work[(size_t)r2 * W + w] ^= work[(size_t)rank * W + w];
        }
      }
    }
    s->piv[rank++] = (uint16_t)col;
  }
  free(work);
  if (rank < K) {
    return false;
  }

  /* P[u][k] = G[u][piv[k]]; Gauss-Jordan [P | I] -> [I | P^-1] */
  uint64_t A[GS_K_MAX][GS_KW];
  memset(A, 0, sizeof(A));
  memset(s->Pinv, 0, sizeof(s->Pinv));
  for (uint32_t u = 0; u < K; u++) {
    for (uint32_t k = 0; k < K; k++) {
      if (gs_bit(&s->G[(size_t)u * W], s->piv[k])) {
        A[u][k >> 6] |= 1ULL << (k & 63U);
      }
    }
    s->Pinv[u][u >> 6] |= 1ULL << (u & 63U);
  }
  for (uint32_t c = 0; c < K; c++) {
    uint32_t r = c;
    while (r < K && !gs_bit(A[r], c)) {
      r++;
    }
    if (r == K) {
      return false;
    }
    if (r != c) {
      for (uint32_t w = 0; w < GS_KW; w++) {
        uint64_t t = A[r][w];
        A[r][w]    = A[c][w];
        A[c][w]    = t;
        t          = s->Pinv[r][w];
        s->Pinv[r][w] = s->Pinv[c][w];
        s->Pinv[c][w] = t;
      }
    }
    for (uint32_t r2 = 0; r2 < K; r2++) {
      if (r2 != c && gs_bit(A[r2], c)) {
        for (uint32_t w = 0; w < GS_KW; w++) {
          A[r2][w] ^= A[c][w];
          s->Pinv[r2][w] ^= s->Pinv[c][w];
        }
      }
    }
  }
  return true;
}

/* Builds the tables of one (size, AL). Caller holds the mutex. */
static void gs_entry_build(uint32_t N, uint32_t Ll)
{
  gs_entry_t*    e = &gs_tab[N][Ll];
  const uint32_t E = 108U << Ll;
  const uint32_t W = gs_words(E);
  const uint32_t K = N + 16;
  static uint8_t payload[GS_MAX_DCI_BITS], f[GS_E_MAX];

  memset(e, 0, sizeof(*e));
  e->b0 = calloc(W, sizeof(uint64_t));
  for (uint32_t i = 0; i < 2; i++) {
    e->sys[i].K = K;
    e->sys[i].G = calloc((size_t)K * W, sizeof(uint64_t));
  }
  if (e->b0 == NULL || e->sys[0].G == NULL || e->sys[1].G == NULL) {
    e->ready = 1;
    return;
  }
  memset(payload, 0, sizeof(payload));
  if (gs_encode(payload, N, 0, 0, E, f) < SRSRAN_SUCCESS) {
    e->ready = 1;
    return;
  }
  gs_pack(f, E, e->b0);

  uint64_t row[GS_EW];
  for (uint32_t u = 0; u < K; u++) {
    uint16_t rnti = 0;
    if (u < N) {
      payload[u] = 1;
    } else {
      rnti = (uint16_t)(1U << (15U - (u - N))); /* RNTI MSB first */
    }
    gs_encode(payload, N, rnti, 0, E, f);
    if (u < N) {
      payload[u] = 0;
    }
    gs_pack(f, E, row);
    for (uint32_t w = 0; w < W; w++) {
      row[w] ^= e->b0[w];
      e->sys[0].G[(size_t)u * W + w] = row[w];
      e->sys[1].G[(size_t)u * W + w] = row[w];
    }
    /* n_RNTI = RNTI: RNTI bit b (b <= 14) is also c_init bit 16 + b */
    if (u > N) {
      uint32_t b = 15U - (u - N);
      for (uint32_t w = 0; w < W; w++) {
        e->sys[1].G[(size_t)u * W + w] ^= gs_scr[Ll][16 + b][w];
      }
    }
  }
  for (uint32_t i = 0; i < 2; i++) {
    e->sys[i].valid = gs_sys_solve_setup(&e->sys[i], E);
  }
  __atomic_store_n(&e->ready, 1, __ATOMIC_RELEASE);
}

static const gs_entry_t* gs_entry_get(uint32_t N, uint32_t Ll)
{
  gs_entry_t* e = &gs_tab[N][Ll];
  if (__atomic_load_n(&e->ready, __ATOMIC_ACQUIRE)) {
    return e;
  }
  pthread_mutex_lock(&gs_gf2_mutex);
  if (!e->ready) {
    gs_entry_build(N, Ll);
  }
  pthread_mutex_unlock(&gs_gf2_mutex);
  return e;
}

void gs_gf2_prepare(uint32_t N)
{
  for (uint32_t Ll = 0; Ll <= GS_MAX_L_LOG2; Ll++) {
    if (gs_gf2_size_valid(N, Ll)) {
      gs_entry_get(N, Ll);
    }
  }
}

bool gs_dci_try(const uint8_t* bits, uint32_t E, uint32_t nof_bits, uint32_t n_id, bool n_rnti, uint8_t* payload,
                uint16_t* rnti)
{
  uint32_t Ll = 0;
  while (Ll <= GS_MAX_L_LOG2 && (108U << Ll) != E) {
    Ll++;
  }
  if (Ll > GS_MAX_L_LOG2 || !gs_gf2_size_valid(nof_bits, Ll) || !gs_gf2_inited) {
    return false;
  }
  const gs_entry_t* e = gs_entry_get(nof_bits, Ll);
  const gs_sys_t*   s = &e->sys[n_rnti ? 1 : 0];
  if (!s->valid) {
    return false;
  }
  const uint32_t W = gs_words(E);
  const uint32_t K = s->K;

  /* y = b + b0 + scr(N_ID) */
  uint64_t y[GS_EW];
  gs_pack(bits, E, y);
  for (uint32_t w = 0; w < W; w++) {
    y[w] ^= e->b0[w];
  }
  for (uint32_t b = 0; b < 16; b++) {
    if ((n_id >> b) & 1U) {
      for (uint32_t w = 0; w < W; w++) {
        y[w] ^= gs_scr[Ll][b][w];
      }
    }
  }

  /* x = y[piv] * Pinv */
  uint64_t x[GS_KW] = {0};
  for (uint32_t k = 0; k < K; k++) {
    if (gs_bit(y, s->piv[k])) {
      for (uint32_t w = 0; w < GS_KW; w++) {
        x[w] ^= s->Pinv[k][w];
      }
    }
  }
  /* consistency: x G == y on all E bits */
  uint64_t chk[GS_EW] = {0};
  for (uint32_t u = 0; u < K; u++) {
    if (gs_bit(x, u)) {
      const uint64_t* g = &s->G[(size_t)u * W];
      for (uint32_t w = 0; w < W; w++) {
        chk[w] ^= g[w];
      }
    }
  }
  for (uint32_t w = 0; w < W; w++) {
    if (chk[w] != y[w]) {
      return false;
    }
  }
  for (uint32_t i = 0; i < nof_bits; i++) {
    payload[i] = (uint8_t)gs_bit(x, i);
  }
  uint16_t r = 0;
  for (uint32_t j = 0; j < 16; j++) {
    r = (uint16_t)((r << 1) | (uint16_t)gs_bit(x, nof_bits + j));
  }
  *rnti = r;
  return true;
}
