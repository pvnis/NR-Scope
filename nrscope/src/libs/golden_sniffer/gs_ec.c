/*
 * SPDX-License-Identifier: MIT
 *
 * GoldenSniffer port: soft-information error correction of the GF(2) DCI decoder.
 * Port of the LGS error correction of 5GDescrambler (https://github.com/kit-ps/5GDescrambler,
 * pdcch/src/lgs_error_correction.rs, pdcch_macros/src/decode_mat.rs, the AL1 rule of
 * phys/src/pipeline/candidate/candidate_finding.rs). See gs_ec.h for the algorithm and
 * how it maps onto GoldenSniffer's formulation.
 *
 * Tables, per (DCI size N, AL, scrambling case), built once from srsRAN's own encoder
 * like gs_gf2.c builds its generator (5GDescrambler's precomputed "decode matrices"):
 *   G   (N + 16) x R affine generator restricted to the first R bits (R = E, or the
 *       polar length N_p when E > N_p: one copy of every coded bit);
 *   its reduced row echelon form M = T G with pivots P: M[:, P] = I, so for any word y
 *   of the code x = y_P T, and the syndrome of a word is
 *        s = y_Q + y_P M_Q        (Q = the R - K non-pivot bits, one check each).
 * Flipping bit j changes s by row j of the parity-check matrix: M_Q[k] if j is pivot
 * k, the unit vector of j if j is in Q.
 *
 * Per call, when the hard decision fails (lgs()):
 *   S = the t least reliable of the R bits (smallest_abs_indices), ascending |llr|;
 *   rows h_j, j in S, are reduced in that order (gauss_elimination over the error
 *   unknowns: a row that reduces to zero is a free unknown, kept at its hard value);
 *   the syndrome is reduced by the pivot rows; a nonzero rest = inconsistent (errors
 *   outside S, or no DCI): reject (5GDescrambler's final syndrome_match);
 *   the flips are the combination of S rows that gives s.
 * 5GDescrambler sets the error bits to 0 and solves for their values (with a 1 column
 * for the negation bit); solving for the flips of the hard decision is the same
 * system written relative to the hard decision.
 */
#include "gs_ec.h"
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
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#define EC_K_MAX (GS_MAX_DCI_BITS + 16)
#define EC_KW ((EC_K_MAX + 63) / 64)
#define EC_E_MAX (108 << GS_MAX_L_LOG2)
#define EC_R_MAX 512 /* PDCCH polar length N_p <= 2^9 */
#define EC_RW (EC_R_MAX / 64)
#define EC_TW ((GS_EC_T_MAX + 63) / 64)

typedef struct {
  bool     valid;
  uint32_t K, R, Q;
  uint16_t piv[EC_K_MAX];
  int16_t  pos_k[EC_R_MAX]; /* pivot index of a bit, -1 for a check bit */
  int16_t  pos_q[EC_R_MAX]; /* check index of a bit, -1 for a pivot bit */
  uint64_t Mq[EC_K_MAX][EC_RW]; /* syndrome row of pivot bit k */
  uint64_t T[EC_K_MAX][EC_KW];  /* x = y_P T */
  /* generator rows of the RNTI bits that do not enter the scrambling (all 16 with
   * n_RNTI = 0, the MSB with n_RNTI = RNTI): low-weight codewords, the neighbours a
   * corrected word is checked against */
  uint32_t nof_nb;
  uint64_t nb[16][EC_RW];
} ec_sys_t;

typedef struct {
  int      ready;
  uint32_t R;
  uint64_t b0[EC_RW]; /* encoding of all zeros, first R bits */
  ec_sys_t sys[2];    /* [n_rnti] */
} ec_entry_t;

static ec_entry_t*     ec_tab[GS_MAX_DCI_BITS + 1][GS_MAX_L_LOG2 + 1];
static uint64_t        ec_scr[GS_MAX_L_LOG2 + 1][31][EC_RW]; /* x2 part of the scrambling of c_init = 2^b */
static pthread_mutex_t ec_mutex  = PTHREAD_MUTEX_INITIALIZER;
static int             ec_inited = 0;

static srsran_polar_code_t    ec_code;
static srsran_polar_encoder_t ec_enc;
static srsran_polar_rm_t      ec_rm;
static srsran_crc_t           ec_crc;

static inline int ec_bit(const uint64_t* w, uint32_t i)
{
  return (int)((w[i >> 6] >> (i & 63U)) & 1U);
}

static void ec_pack(const uint8_t* bits, uint32_t n, uint64_t* w, uint32_t nw)
{
  memset(w, 0, nw * sizeof(uint64_t));
  for (uint32_t i = 0; i < n; i++) {
    if (bits[i] & 1U) {
      w[i >> 6] |= 1ULL << (i & 63U);
    }
  }
}

/* srsran_pdcch_nr_encode() up to the scrambling, as gs_encode() in gs_gf2.c */
static int ec_encode(const uint8_t* payload, uint32_t N, uint16_t rnti, uint32_t cinit, uint32_t E, uint8_t* f)
{
  uint8_t  c_full[24 + EC_K_MAX + 24 + 24];
  uint8_t  c_prime[SRSRAN_POLAR_INTERLEAVER_K_MAX_IL];
  uint8_t  allocated[NMAX];
  uint8_t  d[NMAX];
  uint32_t K = N + 24;

  if (srsran_polar_code_get(&ec_code, (uint16_t)K, (uint16_t)E, 9U) < SRSRAN_SUCCESS) {
    return SRSRAN_ERROR;
  }
  uint8_t* c = c_full;
  srsran_bit_unpack(UINT32_MAX, &c, 24U);
  memcpy(c, payload, N);
  srsran_crc_attach(&ec_crc, c_full, (int)K);
  uint8_t  unpacked_rnti[16];
  uint8_t* ptr = unpacked_rnti;
  srsran_bit_unpack(rnti, &ptr, 16);
  for (uint32_t i = 0; i < 16; i++) {
    c[K - 16 + i] ^= unpacked_rnti[i];
  }
  srsran_polar_interleaver_run_u8(c, c_prime, K, true);
  srsran_polar_chanalloc_tx(c_prime, allocated, ec_code.N, ec_code.K, ec_code.nPC, ec_code.K_set, ec_code.PC_set);
  if (srsran_polar_encoder_encode(&ec_enc, allocated, d, ec_code.n) < SRSRAN_SUCCESS) {
    return SRSRAN_ERROR;
  }
  srsran_polar_rm_tx(&ec_rm, d, f, ec_code.n, E, K, 0);
  srsran_sequence_apply_bit(f, f, E, cinit);
  return SRSRAN_SUCCESS;
}

static int ec_init(void)
{
  if (__atomic_load_n(&ec_inited, __ATOMIC_ACQUIRE)) {
    return SRSRAN_SUCCESS;
  }
  if (gs_gf2_init() < SRSRAN_SUCCESS) { /* the hard path is gs_dci_try */
    return SRSRAN_ERROR;
  }
  pthread_mutex_lock(&ec_mutex);
  int ret = SRSRAN_SUCCESS;
  if (!ec_inited) {
    if (srsran_crc_init(&ec_crc, SRSRAN_LTE_CRC24C, 24) < SRSRAN_SUCCESS ||
        srsran_polar_code_init(&ec_code) < SRSRAN_SUCCESS ||
        srsran_polar_encoder_init(&ec_enc, SRSRAN_POLAR_ENCODER_PIPELINED, NMAX_LOG) < SRSRAN_SUCCESS ||
        srsran_polar_rm_tx_init(&ec_rm) < SRSRAN_SUCCESS) {
      ret = SRSRAN_ERROR;
    } else {
      static uint8_t zeros[EC_E_MAX], s0[EC_E_MAX], sb[EC_E_MAX];
      for (uint32_t Ll = 0; Ll <= GS_MAX_L_LOG2; Ll++) {
        uint32_t R = (108U << Ll) < EC_R_MAX ? (108U << Ll) : EC_R_MAX;
        srsran_sequence_apply_bit(zeros, s0, R, 0);
        for (uint32_t b = 0; b < 31; b++) {
          srsran_sequence_apply_bit(zeros, sb, R, 1U << b);
          for (uint32_t i = 0; i < R; i++) {
            sb[i] ^= s0[i];
          }
          ec_pack(sb, R, ec_scr[Ll][b], EC_RW);
        }
      }
      __atomic_store_n(&ec_inited, 1, __ATOMIC_RELEASE);
    }
  }
  pthread_mutex_unlock(&ec_mutex);
  return ret;
}

/* Gauss-Jordan of the K x R generator: pivots, M_Q rows and T. */
static bool ec_sys_setup(ec_sys_t* s, uint64_t (*G)[EC_RW])
{
  const uint32_t K = s->K, R = s->R;
  memset(s->T, 0, sizeof(s->T));
  for (uint32_t u = 0; u < K; u++) {
    s->T[u][u >> 6] |= 1ULL << (u & 63U);
  }
  uint32_t rank = 0;
  for (uint32_t col = 0; col < R && rank < K; col++) {
    uint32_t r = rank;
    while (r < K && !ec_bit(G[r], col)) {
      r++;
    }
    if (r == K) {
      continue;
    }
    if (r != rank) {
      for (uint32_t w = 0; w < EC_RW; w++) {
        uint64_t t = G[r][w];
        G[r][w]    = G[rank][w];
        G[rank][w] = t;
      }
      for (uint32_t w = 0; w < EC_KW; w++) {
        uint64_t t    = s->T[r][w];
        s->T[r][w]    = s->T[rank][w];
        s->T[rank][w] = t;
      }
    }
    for (uint32_t r2 = 0; r2 < K; r2++) {
      if (r2 != rank && ec_bit(G[r2], col)) {
        for (uint32_t w = 0; w < EC_RW; w++) {
          G[r2][w] ^= G[rank][w];
        }
        for (uint32_t w = 0; w < EC_KW; w++) {
          s->T[r2][w] ^= s->T[rank][w];
        }
      }
    }
    s->piv[rank++] = (uint16_t)col;
  }
  if (rank < K) {
    return false;
  }
  /* bit roles, and the reduced rows restricted to the check bits */
  for (uint32_t j = 0; j < R; j++) {
    s->pos_k[j] = -1;
    s->pos_q[j] = -1;
  }
  for (uint32_t k = 0; k < K; k++) {
    s->pos_k[s->piv[k]] = (int16_t)k;
  }
  uint32_t q = 0;
  for (uint32_t j = 0; j < R; j++) {
    if (s->pos_k[j] < 0) {
      s->pos_q[j] = (int16_t)q++;
    }
  }
  s->Q = q;
  memset(s->Mq, 0, sizeof(s->Mq));
  for (uint32_t k = 0; k < K; k++) {
    for (uint32_t j = 0; j < R; j++) {
      if (s->pos_q[j] >= 0 && ec_bit(G[k], j)) {
        s->Mq[k][s->pos_q[j] >> 6] |= 1ULL << (s->pos_q[j] & 63);
      }
    }
  }
  return true;
}

/* Builds the tables of one (size, AL). Caller holds ec_mutex. */
static ec_entry_t* ec_entry_build(uint32_t N, uint32_t Ll)
{
  const uint32_t E = 108U << Ll;
  const uint32_t K = N + 16;
  ec_entry_t*    e = calloc(1, sizeof(ec_entry_t));
  uint64_t(*G)[EC_RW] = calloc(EC_K_MAX * 2, sizeof(*G));
  static uint8_t payload[GS_MAX_DCI_BITS], f[EC_E_MAX];
  if (e == NULL || G == NULL) {
    free(e);
    free(G);
    return NULL;
  }
  e->ready = 1;
  memset(payload, 0, sizeof(payload));
  if (ec_encode(payload, N, 0, 0, E, f) < SRSRAN_SUCCESS) {
    free(G);
    return e; /* sys[].valid false */
  }
  const uint32_t R = E < ec_code.N ? E : ec_code.N; /* get_optimal_rate_matching_limit() */
  e->R             = R;
  ec_pack(f, R, e->b0, EC_RW);

  uint64_t(*G1)[EC_RW] = G + EC_K_MAX;
  for (uint32_t u = 0; u < K; u++) {
    uint16_t rnti = 0;
    if (u < N) {
      payload[u] = 1;
    } else {
      rnti = (uint16_t)(1U << (15U - (u - N))); /* RNTI MSB first, as gs_gf2.c */
    }
    ec_encode(payload, N, rnti, 0, E, f);
    if (u < N) {
      payload[u] = 0;
    }
    ec_pack(f, R, G[u], EC_RW);
    for (uint32_t w = 0; w < EC_RW; w++) {
      G[u][w] ^= e->b0[w];
      G1[u][w] = G[u][w];
    }
    if (u > N) { /* n_RNTI = RNTI: RNTI bit b <= 14 is also c_init bit 16 + b */
      uint32_t b = 15U - (u - N);
      for (uint32_t w = 0; w < EC_RW; w++) {
        G1[u][w] ^= ec_scr[Ll][16 + b][w];
      }
    }
  }
  for (uint32_t i = 0; i < 2; i++) {
    ec_sys_t* s = &e->sys[i];
    s->K        = K;
    s->R        = R;
    s->nof_nb   = i ? 1 : 16;
    for (uint32_t b = 0; b < s->nof_nb; b++) {
      memcpy(s->nb[b], (i ? G1 : G)[N + b], sizeof(s->nb[b]));
    }
    s->valid = ec_sys_setup(s, i ? G1 : G);
  }
  free(G);
  return e;
}

/* Once a background build has started (gs_gf2_prebuild_start), tables are built there
   only: a caller never builds nor waits, a table not ready means no correction. */
int gs_ec_no_lazy = 0;

static const ec_entry_t* ec_entry_build_locked(uint32_t N, uint32_t Ll);

static const ec_entry_t* ec_entry_get(uint32_t N, uint32_t Ll)
{
  ec_entry_t* e = __atomic_load_n(&ec_tab[N][Ll], __ATOMIC_ACQUIRE);
  if (e != NULL) {
    return e;
  }
  if (__atomic_load_n(&gs_ec_no_lazy, __ATOMIC_ACQUIRE)) {
    return NULL;
  }
  return ec_entry_build_locked(N, Ll);
}

static const ec_entry_t* ec_entry_build_locked(uint32_t N, uint32_t Ll)
{
  ec_entry_t* e;
  pthread_mutex_lock(&ec_mutex);
  e = ec_tab[N][Ll];
  if (e == NULL) {
    e = ec_entry_build(N, Ll);
    __atomic_store_n(&ec_tab[N][Ll], e, __ATOMIC_RELEASE);
  }
  pthread_mutex_unlock(&ec_mutex);
  return e;
}

void gs_ec_prepare(uint32_t N)
{
  if (ec_init() < SRSRAN_SUCCESS) {
    return;
  }
  for (uint32_t Ll = 0; Ll <= GS_MAX_L_LOG2; Ll++) {
    if (gs_gf2_size_valid(N, Ll)) {
      ec_entry_build_locked(N, Ll);
    }
  }
}

void gs_ec_args_default(gs_ec_args_t* a)
{
  memset(a, 0, sizeof(*a));
  a->frac_ue    = 0.10f;
  a->frac_css   = 0.05f;
  a->min_parity = 16;
  a->max_flips  = 0;
  a->min_metric = -INFINITY;
  a->min_nb_margin = 0.0f;
  a->hard_first = true;
}

/* indices of the t smallest |llr[0..R)|, ascending (smallest_abs_indices) */
static void ec_least_reliable(const float* llr, uint32_t R, uint32_t t, uint16_t* S)
{
  /* max-heap of the t smallest so far */
  float    hv[GS_EC_T_MAX];
  uint16_t hi[GS_EC_T_MAX];
  uint32_t n = 0;
  for (uint32_t j = 0; j < R; j++) {
    float v = fabsf(llr[j]);
    if (n < t) {
      uint32_t c = n++;
      while (c > 0 && hv[(c - 1) / 2] < v) {
        hv[c] = hv[(c - 1) / 2];
        hi[c] = hi[(c - 1) / 2];
        c     = (c - 1) / 2;
      }
      hv[c] = v;
      hi[c] = (uint16_t)j;
    } else if (v < hv[0]) {
      uint32_t c = 0;
      for (;;) {
        uint32_t l = 2 * c + 1, r = l + 1, m = c;
        float    mv = v;
        if (l < n && hv[l] > mv) {
          m  = l;
          mv = hv[l];
        }
        if (r < n && hv[r] > mv) {
          m = r;
        }
        if (m == c) {
          break;
        }
        hv[c] = hv[m];
        hi[c] = hi[m];
        c     = m;
      }
      hv[c] = v;
      hi[c] = (uint16_t)j;
    }
  }
  /* pop: largest first, written from the back */
  while (n > 0) {
    S[n - 1] = hi[0];
    float    v = hv[n - 1];
    uint16_t x = hi[n - 1];
    n--;
    uint32_t c = 0;
    for (;;) {
      uint32_t l = 2 * c + 1, r = l + 1, m = c;
      float    mv = v;
      if (l < n && hv[l] > mv) {
        m  = l;
        mv = hv[l];
      }
      if (r < n && hv[r] > mv) {
        m = r;
      }
      if (m == c) {
        break;
      }
      hv[c] = hv[m];
      hi[c] = hi[m];
      c     = m;
    }
    if (n > 0) {
      hv[c] = v;
      hi[c] = x;
    }
  }
}

/* log2 of the number of patterns of weight <= w on r unknowns */
static double ec_log2_ball(uint32_t r, uint32_t w)
{
  double c = 1.0, v = 1.0;
  for (uint32_t i = 1; i <= w && i <= r; i++) {
    c = c * (double)(r - i + 1) / (double)i;
    v += c;
  }
  return log2(v);
}

/* 5GDescrambler's confidence sum over the t bits, divided by t times the mean |llr|
 * of the R bits: independent of the LLR scale; ~0 for a random solution (half the t
 * bits flipped), positive when few, unreliable bits were corrected */
static float ec_norm_metric(const float* llr, uint32_t R, uint32_t t, float sum)
{
  float a = 0;
  for (uint32_t j = 0; j < R; j++) {
    a += fabsf(llr[j]);
  }
  return t > 0 && a > 0 ? sum * (float)R / ((float)t * a) : 0.0f;
}

static uint32_t ec_lsb(const uint64_t* w, uint32_t nw)
{
  for (uint32_t i = 0; i < nw; i++) {
    if (w[i]) {
      return i * 64U + (uint32_t)__builtin_ctzll(w[i]);
    }
  }
  return UINT32_MAX;
}

bool gs_dci_try_soft_ex(const gs_ec_args_t* a, const float* llr, uint32_t E, uint32_t nof_bits, uint32_t n_id,
                        bool n_rnti, uint32_t effort, uint8_t* payload, uint16_t* rnti, gs_ec_info_t* info)
{
  gs_ec_info_t dummy;
  if (info == NULL) {
    info = &dummy;
  }
  memset(info, 0, sizeof(*info));
  uint32_t Ll = 0;
  while (Ll <= GS_MAX_L_LOG2 && (108U << Ll) != E) {
    Ll++;
  }
  if (Ll > GS_MAX_L_LOG2 || !gs_gf2_size_valid(nof_bits, Ll) || nof_bits + 16 + a->min_parity > E ||
      ec_init() < SRSRAN_SUCCESS) {
    return false;
  }

  /* 1. hard decision: exactly gs_dci_try */
  uint8_t hb[EC_E_MAX];
  for (uint32_t i = 0; i < E; i++) {
    hb[i] = llr[i] < 0.0f;
  }
  bool hard_ok = a->hard_first && gs_dci_try(hb, E, nof_bits, n_id, n_rnti, payload, rnti);
  if (!hard_ok && effort == 0) {
    return false;
  }

  const ec_entry_t* e = ec_entry_get(nof_bits, Ll);
  if (e == NULL || !e->sys[n_rnti ? 1 : 0].valid) {
    return hard_ok;
  }
  const ec_sys_t* s = &e->sys[n_rnti ? 1 : 0];
  const uint32_t  R = s->R, Q = s->Q, K = s->K;
  const uint32_t  QW = (Q + 63) / 64;

  /* error budget (decode_mat.rs), leaving min_parity checks */
  uint32_t t = effort;
  if (effort == GS_EC_EFFORT_AUTO || (effort == GS_EC_EFFORT_MAX && !n_rnti)) {
    t = (uint32_t)((float)E * (n_rnti ? a->frac_ue : a->frac_css));
  } else if (effort == GS_EC_EFFORT_MAX) {
    t = GS_EC_T_MAX; /* capped to Q - min_parity below */
  }
  if (t > GS_EC_T_MAX) {
    t = GS_EC_T_MAX;
  }
  t = Q > a->min_parity ? (t < Q - a->min_parity ? t : Q - a->min_parity) : 0;
  if (t > R) {
    t = R;
  }
  info->R = R;
  info->Q = Q;
  info->t = t;

  uint16_t S[GS_EC_T_MAX];
  if (t > 0) {
    ec_least_reliable(llr, R, t, S);
  }

  if (hard_ok) {
    float m = 0;
    for (uint32_t i = 0; i < t; i++) {
      m += fabsf(llr[S[i]]);
    }
    info->hard       = true;
    info->rank       = 0;
    info->eff_parity = (float)(E - K);
    info->metric_raw = m / (2.0f * (float)R);
    info->metric     = ec_norm_metric(llr, R, t, m);
    return true;
  }
  if (t == 0) {
    return false;
  }

  /* 2. y = b + b0 + scr(N_ID) over the R bits, and its syndrome */
  uint64_t y[EC_RW];
  ec_pack(hb, R, y, EC_RW);
  for (uint32_t w = 0; w < EC_RW; w++) {
    y[w] ^= e->b0[w];
  }
  for (uint32_t b = 0; b < 16; b++) {
    if ((n_id >> b) & 1U) {
      for (uint32_t w = 0; w < EC_RW; w++) {
        y[w] ^= ec_scr[Ll][b][w];
      }
    }
  }
  uint64_t syn[EC_RW] = {0};
  for (uint32_t j = 0; j < R; j++) {
    if (!ec_bit(y, j)) {
      continue;
    }
    if (s->pos_q[j] >= 0) {
      syn[s->pos_q[j] >> 6] ^= 1ULL << (s->pos_q[j] & 63);
    } else {
      const uint64_t* m = s->Mq[s->pos_k[j]];
      for (uint32_t w = 0; w < QW; w++) {
        syn[w] ^= m[w];
      }
    }
  }

  /* 3. elimination over the error unknowns, least reliable first */
  uint64_t rows[GS_EC_T_MAX][EC_RW];
  uint64_t comb[GS_EC_T_MAX][EC_TW];
  uint16_t pcol[GS_EC_T_MAX];
  bool     is_free[GS_EC_T_MAX];
  uint32_t nr = 0;
  for (uint32_t i = 0; i < t; i++) {
    const uint32_t j = S[i];
    uint64_t*      h = rows[nr];
    uint64_t*      c = comb[nr];
    if (s->pos_q[j] >= 0) {
      memset(h, 0, QW * sizeof(uint64_t));
      h[s->pos_q[j] >> 6] = 1ULL << (s->pos_q[j] & 63);
    } else {
      memcpy(h, s->Mq[s->pos_k[j]], QW * sizeof(uint64_t));
    }
    memset(c, 0, EC_TW * sizeof(uint64_t));
    c[i >> 6] = 1ULL << (i & 63U);
    for (uint32_t p = 0; p < nr; p++) {
      if (ec_bit(h, pcol[p])) {
        for (uint32_t w = 0; w < QW; w++) {
          h[w] ^= rows[p][w];
        }
        for (uint32_t w = 0; w < EC_TW; w++) {
          c[w] ^= comb[p][w];
        }
      }
    }
    uint32_t lsb = ec_lsb(h, QW);
    is_free[i]   = lsb == UINT32_MAX;
    if (!is_free[i]) {
      pcol[nr++] = (uint16_t)lsb;
    }
  }
  uint64_t flip[EC_TW] = {0};
  for (uint32_t p = 0; p < nr; p++) {
    if (ec_bit(syn, pcol[p])) {
      for (uint32_t w = 0; w < QW; w++) {
        syn[w] ^= rows[p][w];
      }
      for (uint32_t w = 0; w < EC_TW; w++) {
        flip[w] ^= comb[p][w];
      }
    }
  }
  info->rank = nr;
  if (ec_lsb(syn, QW) != UINT32_MAX) {
    return false; /* inconsistent: errors outside S, or not a DCI of this size */
  }

  /* 4. corrected word, confidence and acceptance */
  uint32_t nflip = 0;
  float    m = 0, fl = 0;
  for (uint32_t i = 0; i < t; i++) {
    float v = fabsf(llr[S[i]]);
    if (ec_bit(flip, i)) {
      y[S[i] >> 6] ^= 1ULL << (S[i] & 63U);
      nflip++;
      m -= v;
      fl += v;
    } else if (is_free[i]) {
      m -= v;
    } else {
      m += v;
    }
  }
  info->flips      = nflip;
  info->flip_llr   = fl;
  info->metric_raw = m / (2.0f * (float)R);
  info->metric     = ec_norm_metric(llr, R, t, m);
  info->eff_parity = (float)((double)Q - ec_log2_ball(nr, nflip));
  if (info->eff_parity < (float)a->min_parity || (a->max_flips > 0 && nflip > a->max_flips) ||
      info->metric < a->min_metric) {
    return false;
  }

  /* neighbours (not in 5GDescrambler): the corrected word plus one RNTI-only row is a
   * codeword too, at another RNTI. Its extra discrepancy with the hard decision,
   * sum over the row of |llr| (agreeing bits) - |llr| (corrected bits), must be positive
   * (the decode is the likelier of the two) by min_nb_margin mean |llr| */
  {
    uint64_t fm[EC_RW] = {0};
    for (uint32_t i = 0; i < t; i++) {
      if (ec_bit(flip, i)) {
        fm[S[i] >> 6] |= 1ULL << (S[i] & 63U);
      }
    }
    float mean = 0;
    for (uint32_t j = 0; j < R; j++) {
      mean += fabsf(llr[j]);
    }
    mean = mean > 0 ? mean / (float)R : 1.0f;
    float worst = INFINITY;
    for (uint32_t b = 0; b < s->nof_nb; b++) {
      float d = 0;
      for (uint32_t w = 0; w < EC_RW; w++) {
        for (uint64_t v = s->nb[b][w]; v != 0; v &= v - 1) {
          uint32_t j = w * 64U + (uint32_t)__builtin_ctzll(v);
          d += (fm[w] >> (j & 63U)) & 1U ? -fabsf(llr[j]) : fabsf(llr[j]);
        }
      }
      worst = d < worst ? d : worst;
    }
    info->nb_margin = worst / mean;
    if (info->nb_margin < a->min_nb_margin) {
      return false;
    }
  }

  /* x = y_P T */
  uint64_t x[EC_KW] = {0};
  for (uint32_t k = 0; k < K; k++) {
    if (ec_bit(y, s->piv[k])) {
      for (uint32_t w = 0; w < EC_KW; w++) {
        x[w] ^= s->T[k][w];
      }
    }
  }
  uint16_t r = 0;
  for (uint32_t j = 0; j < 16; j++) {
    r = (uint16_t)((r << 1) | (uint16_t)ec_bit(x, nof_bits + j));
  }
  if (n_rnti && (r == 0 || r > 0xFFF2)) {
    return false; /* not a UE-specific RNTI (candidate_finding.rs plausibility check) */
  }
  for (uint32_t i = 0; i < nof_bits; i++) {
    payload[i] = (uint8_t)ec_bit(x, i);
  }
  *rnti = r;
  return true;
}

bool gs_dci_try_soft(gs_t* q, const float* llr, uint32_t E, uint32_t nof_bits, uint32_t n_id, bool n_rnti,
                     uint32_t effort, uint8_t* payload, uint16_t* rnti, float* metric)
{
  gs_ec_args_t a;
  gs_ec_args_default(&a);
  if (q != NULL) {
    a.min_parity = q->args.min_parity;
  }
  gs_ec_info_t info;
  bool         ok = gs_dci_try_soft_ex(&a, llr, E, nof_bits, n_id, n_rnti, effort, payload, rnti, &info);
  if (ok && metric != NULL) {
    *metric = info.metric;
  }
  return ok;
}
