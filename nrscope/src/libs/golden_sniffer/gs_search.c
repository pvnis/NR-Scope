/*
 * SPDX-License-Identifier: MIT
 *
 * GoldenSniffer port: blind PDCCH search of one slot.
 * Ports blind_search_for_DCIs.m, power_detect_REGs.m, enum_candidates.m,
 * update_candidates.m, unwrap_QPSK.m and PDCCH_equalize.m of
 * https://github.com/5G-GoldenSniffer/5G-GoldenSniffer (Copyright 2024-2026 the
 * authors, MIT License).
 *
 * Differences from the MATLAB, all deliberate:
 *  - interleaved pairs of runs holding n > 1 candidates take Nsc_cand/2 from each run
 *    per candidate (the MATLAB takes Nsc_cand/(2n), a candidate of the wrong size);
 *  - after a DCI is found the search restarts from the largest AL the remaining runs
 *    allow, not above the current one (the MATLAB mixes stale candidates and ALs);
 *  - DCI sizes: polar codes of different sizes are nested, so with n_RNTI = 0 a
 *    codeword also solves a larger size with probability 2^-8 (the unmasked CRC bits)
 *    and a smaller one with 2^-(8 + size difference). The MATLAB stops at the most
 *    frequent size that solves; here every size is tried and, among those that solve,
 *    the most frequent so far wins. With n_RNTI = RNTI the CRC must match the RNTI of
 *    the scrambling, and the first size that solves is taken;
 *  - interleaved candidates are also cut from merged runs at a spacing learned from
 *    earlier interleaved DCIs (the MATLAB needs runs of equal length, which adjacent
 *    DCIs of similar power never give);
 *  - both scrambling hypotheses (n_RNTI = 0 and n_RNTI = RNTI) are tried for every
 *    N_ID: a common search space in a CORESET with pdcch-DMRS-ScramblingID uses that ID
 *    with n_RNTI = 0, and a UE search space may be configured with an ID equal to the PCI;
 *  - the DM-RS offset is relative to point A from SIB1 (offsetToCarrier), not to the
 *    first subcarrier of the capture;
 *  - runs must start on an RB, not on a CCE boundary from point A: a CORESET can start
 *    off the 6-RB grid (rb-Offset-r16; srsRAN's coreset1_rb_start);
 *  - a decode needs min_parity (16) observed bits beyond the unknowns, not 8.
 */
#include "gs_internal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GS_CORR_THRES 0.975f

/* GS_DEBUG=1 in the environment: one line per stage, for replaying dumped grids */
static int gs_debug(void)
{
  static int v = -1;
  if (v < 0) {
    v = getenv("GS_DEBUG") != NULL;
  }
  return v;
}
#define GS_AMPLITUDE_TOL 0.4f

void gs_args_default(gs_args_t* a)
{
  memset(a, 0, sizeof(*a));
  a->max_duration = 2;
  a->size_min     = 20;
  a->size_max     = 84;
  a->max_L_log2   = GS_MAX_L_LOG2;
  a->min_parity   = 16;
}

gs_t* gs_create(const gs_args_t* args)
{
  gs_gold_init();
  if (gs_gf2_init() != 0) {
    return NULL;
  }
  gs_t* q = calloc(1, sizeof(gs_t));
  if (q == NULL) {
    return NULL;
  }
  gs_args_default(&q->args);
  if (args != NULL) {
    q->args = *args;
  }
  if (q->args.max_duration < 1 || q->args.max_duration > GS_MAX_DURATION) {
    q->args.max_duration = 2;
  }
  if (q->args.max_L_log2 > GS_MAX_L_LOG2) {
    q->args.max_L_log2 = GS_MAX_L_LOG2;
  }
  if (q->args.nof_sizes > 0 && q->args.sizes != NULL) {
    for (uint32_t i = 0; i < q->args.nof_sizes && q->nof_sizes <= GS_MAX_DCI_BITS; i++) {
      uint32_t n = q->args.sizes[i];
      if (n >= 1 && n <= GS_MAX_DCI_BITS) {
        q->sizes[q->nof_sizes++] = n;
      }
    }
  } else {
    uint32_t lo = q->args.size_min < 1 ? 1 : q->args.size_min;
    uint32_t hi = q->args.size_max > GS_MAX_DCI_BITS ? GS_MAX_DCI_BITS : q->args.size_max;
    for (uint32_t n = lo; n <= hi; n++) {
      q->sizes[q->nof_sizes++] = n;
    }
  }
  /* ascending */
  for (uint32_t i = 1; i < q->nof_sizes; i++) {
    for (uint32_t j = i; j > 0 && q->sizes[j] < q->sizes[j - 1]; j--) {
      uint32_t t      = q->sizes[j];
      q->sizes[j]     = q->sizes[j - 1];
      q->sizes[j - 1] = t;
    }
  }
  for (uint32_t n = 0; n <= GS_MAX_DCI_BITS; n++) {
    q->css_ok[n] = q->args.nof_css_sizes == 0;
  }
  for (uint32_t i = 0; q->args.css_sizes != NULL && i < q->args.nof_css_sizes; i++) {
    if (q->args.css_sizes[i] <= GS_MAX_DCI_BITS) {
      q->css_ok[q->args.css_sizes[i]] = true;
    }
  }
  q->args.sizes         = NULL; /* not kept */
  q->args.nof_sizes     = q->nof_sizes;
  q->args.css_sizes     = NULL;
  for (uint32_t i = 0; i < q->nof_sizes; i++) {
    gs_gf2_prepare(q->sizes[i]);
  }
  return q;
}

void gs_destroy(gs_t* q)
{
  free(q);
}

/* ---------------------------------------------------------------- unwrap_QPSK.m */

void gs_unwrap_qpsk(const cf_t* y, uint32_t n, cf_t* yc)
{
  if (n == 0 || n > GS_MAX_SC) {
    return;
  }
  float phi[n];
  yc[0]  = y[0];
  phi[0] = cargf(y[0]);
  for (uint32_t i = 1; i < n; i++) {
    float t = cargf(y[i] * conjf(yc[i - 1]));
    if (fabsf(t) > 3.0f * (float)M_PI / 4.0f) {
      yc[i]  = -y[i];
      phi[i] = phi[i - 1] + (t < 0 ? t + (float)M_PI : t - (float)M_PI);
    } else if (t < -(float)M_PI / 4.0f) {
      yc[i]  = I * y[i];
      phi[i] = phi[i - 1] + t + (float)M_PI / 2.0f;
    } else if (t > (float)M_PI / 4.0f) {
      yc[i]  = -I * y[i];
      phi[i] = phi[i - 1] + t - (float)M_PI / 2.0f;
    } else {
      yc[i]  = y[i];
      phi[i] = phi[i - 1] + t;
    }
  }
  if (n < 2) {
    return;
  }
  /* linear phase fit: the total rotation is only known modulo pi/2 */
  float J_min = INFINITY;
  int   i_min = 0;
  for (int i = -2; i <= 2; i++) {
    float inc = phi[n - 1] - phi[0] + (float)i * (float)M_PI / 2.0f;
    float J   = 0;
    for (uint32_t m = 0; m < n; m++) {
      float a = cargf(yc[m]) - (phi[0] + inc * (float)m / (float)(n - 1));
      float d = a - (float)M_PI / 2.0f * roundf(a / ((float)M_PI / 2.0f));
      J += fabsf(d);
    }
    if (J < J_min) {
      J_min = J;
      i_min = i;
    }
  }
  for (uint32_t m = 0; m < n; m++) {
    float t = phi[0] + (phi[n - 1] - phi[0] + (float)i_min * (float)M_PI / 2.0f) * (float)m / (float)(n - 1) - phi[m];
    if (fabsf(t) > 3.0f * (float)M_PI / 4.0f) {
      yc[m] = -yc[m];
    } else if (t > (float)M_PI / 4.0f) {
      yc[m] = I * yc[m];
    } else if (t < -(float)M_PI / 4.0f) {
      yc[m] = -I * yc[m];
    }
  }
}

/* ---------------------------------------------------------- power_detect_REGs.m */

static bool gs_aligned(int32_t k, int32_t ref, int32_t period)
{
  return ref != GS_NO_REF && k - ref >= 0 && (k - ref) % period == 0;
}

static int gs_run_cmp(const void* a, const void* b)
{
  const gs_run_t* x = a;
  const gs_run_t* y = b;
  return (x->first > y->first) - (x->first < y->first);
}

static uint32_t gs_power_detect(gs_t* q, const gs_slot_t* s)
{
  const uint32_t N   = s->nof_sc;
  const cf_t*    Y0  = s->grid;
  uint32_t       dur = 1;

  if (q->args.max_duration > 1) {
    const cf_t* Y1 = s->grid + N;
    double      c = 0, n0 = 0, n1 = 0;
    for (uint32_t k = 0; k < N; k++) {
      float a0 = cabsf(Y0[k]), a1 = cabsf(Y1[k]);
      c += a0 * a1;
      n0 += a0 * a0;
      n1 += a1 * a1;
    }
    if (n0 > 0 && n1 > 0 && c / sqrt(n0 * n1) > GS_CORR_THRES) {
      dur = 2;
      if (q->args.max_duration > 2) {
        const cf_t* Y2 = s->grid + 2 * N;
        double      c3 = 0, n01 = 0, n2 = 0;
        for (uint32_t k = 0; k < N; k++) {
          float a01 = cabsf(Y0[k]) + cabsf(Y1[k]), a2 = cabsf(Y2[k]);
          c3 += a01 * a2;
          n01 += a01 * a01;
          n2 += a2 * a2;
        }
        if (n01 > 0 && n2 > 0 && c3 / sqrt(n01 * n2) > GS_CORR_THRES) {
          dur = 3;
        }
      }
    }
  }

  float* A = q->amp;
  for (uint32_t k = 0; k < N; k++) {
    A[k] = 0;
    for (uint32_t l = 0; l < dur; l++) {
      A[k] += cabsf(s->grid[l * N + k]);
    }
  }

  /* region growing around the strongest subcarrier, until 3 runs that are not a
   * whole number of REGs */
  gs_run_t runs[GS_MAX_RUNS];
  uint32_t nof_runs = 0;
  int      loops    = 3;
  while (loops > 0 && nof_runs < GS_MAX_RUNS) {
    uint32_t pos = 0;
    for (uint32_t k = 1; k < N; k++) {
      if (A[k] > A[pos]) {
        pos = k;
      }
    }
    uint32_t r = pos, l = pos;
    while (r + 1 < N && fabsf(A[r + 1] - A[r]) < GS_AMPLITUDE_TOL * (A[r + 1] + A[r])) {
      r++;
    }
    while (l > 0 && fabsf(A[l - 1] - A[l]) < GS_AMPLITUDE_TOL * (A[l - 1] + A[l])) {
      l--;
    }
    double p = 0;
    for (uint32_t k = l; k <= r; k++) {
      p += (double)A[k] * A[k];
    }
    runs[nof_runs].first = (int32_t)l;
    runs[nof_runs].last  = (int32_t)r;
    runs[nof_runs].pow   = (float)(p / (r - l + 1));
    nof_runs++;
    if ((r - l + 1) % GS_NSC_REG != 0) {
      loops--;
    }
    for (uint32_t k = l; k <= r; k++) {
      A[k] = 0;
    }
  }

  qsort(runs, nof_runs, sizeof(gs_run_t), gs_run_cmp);
  uint32_t i = 0;
  while (nof_runs > 0 && i + 1 < nof_runs) {
    float p0 = runs[i].pow, p1 = runs[i + 1].pow;
    if (runs[i + 1].first - runs[i].last <= 2 && p0 > 0 && p1 > 0 && p0 / p1 + p1 / p0 < 20.0f) {
      float len0 = (float)(runs[i].last - runs[i].first + 1);
      float len1 = (float)(runs[i + 1].last - runs[i + 1].first + 1);
      runs[i].pow  = (p0 * len0 + p1 * len1) / (float)(runs[i + 1].last - runs[i].first + 1);
      runs[i].last = runs[i + 1].last;
      memmove(&runs[i + 1], &runs[i + 2], (nof_runs - i - 2) * sizeof(gs_run_t));
      nof_runs--;
    } else {
      i++;
    }
  }

  const int32_t nsc_cce = GS_NSC_CCE / (int32_t)dur;
  q->nof_runs           = 0;
  for (i = 0; i < nof_runs; i++) {
    int32_t len = runs[i].last - runs[i].first + 1;
    /* RB-aligned only, not CCE-aligned: a CORESET may start off the 6-RB grid
      (rb-Offset-r16), and the CCEs are counted from its start */
    if (len % nsc_cce == 0 &&
        (gs_aligned(runs[i].first, s->ref_crb0, GS_NSC_REG) || gs_aligned(runs[i].first, s->ref_cs0, GS_NSC_REG))) {
      q->runs[q->nof_runs++] = runs[i];
    }
  }
  return dur;
}

/* ------------------------------------------------------------ enum_candidates.m */

static void gs_add_cand(gs_t* q, uint32_t nof_iv, const gs_iv_t* iv)
{
  if (q->nof_cand >= GS_MAX_CAND) {
    return;
  }
  gs_cand_t* c = &q->cand[q->nof_cand++];
  c->nof_iv    = nof_iv;
  for (uint32_t i = 0; i < nof_iv; i++) {
    c->iv[i] = iv[i];
  }
  /* resource elements are mapped in frequency order */
  for (uint32_t i = 1; i < nof_iv; i++) {
    for (uint32_t j = i; j > 0 && c->iv[j].first < c->iv[j - 1].first; j--) {
      gs_iv_t t    = c->iv[j];
      c->iv[j]     = c->iv[j - 1];
      c->iv[j - 1] = t;
    }
  }
}

static void gs_enum_candidates(gs_t* q, uint32_t dur, uint32_t Ll)
{
  const int32_t nsc_cce  = GS_NSC_CCE / (int32_t)dur;
  const int32_t nsc_cand = (1 << Ll) * nsc_cce;
  const gs_run_t* r      = q->runs;
  const uint32_t  n      = q->nof_runs;
  gs_iv_t         iv[GS_MAX_IV];

  q->nof_cand = 0;
  for (uint32_t i = 0; i < n; i++) {
    int32_t len = r[i].last - r[i].first + 1;
    if (len <= 0) {
      continue;
    }
    /* disjoint non-interleaved candidates */
    if (len % nsc_cand == 0) {
      for (int32_t c = 0; c < len / nsc_cand; c++) {
        iv[0].first = r[i].first + c * nsc_cand;
        iv[0].last  = iv[0].first + nsc_cand - 1;
        gs_add_cand(q, 1, iv);
      }
    }
    /* candidates overlapping by half */
    if (len > nsc_cand && len % nsc_cand == nsc_cand / 2) {
      int32_t nc = 1 + (len - nsc_cand) / (nsc_cand / 2);
      for (int32_t c = 0; c < nc; c++) {
        iv[0].first = r[i].first + c * (nsc_cand / 2);
        iv[0].last  = iv[0].first + nsc_cand - 1;
        gs_add_cand(q, 1, iv);
      }
    }
  }
  /* interleaved: pairs of runs of equal length */
  for (uint32_t i = 0; i + 1 < n; i++) {
    for (uint32_t j = i + 1; j < n; j++) {
      int32_t li = r[i].last - r[i].first + 1, lj = r[j].last - r[j].first + 1;
      if (li <= 0 || li != lj || (li + lj) % nsc_cand != 0 || (r[i].first - r[j].first) % nsc_cce != 0) {
        continue;
      }
      int32_t half = nsc_cand / 2;
      for (int32_t c = 0; c < (li + lj) / nsc_cand; c++) {
        iv[0].first = r[i].first + c * half;
        iv[0].last  = iv[0].first + half - 1;
        iv[1].first = r[j].first + c * half;
        iv[1].last  = iv[1].first + half - 1;
        gs_add_cand(q, 2, iv);
      }
    }
  }
  /* interleaved, learned spacing: halves cut out of any run, partner inside a run */
  const int32_t half = nsc_cand / 2;
  for (uint32_t si = 0; si < q->nof_il && half > 0 && half % nsc_cce == 0; si++) {
    int32_t sp = q->il_spacing[si];
    for (uint32_t i = 0; i < n; i++) {
      for (int32_t a = r[i].first; a + half - 1 <= r[i].last; a += nsc_cce) {
        int32_t b = a + sp;
        bool    in_run = false;
        for (uint32_t j = 0; j < n && !in_run; j++) {
          in_run = b >= r[j].first && b + half - 1 <= r[j].last;
        }
        if (!in_run) {
          continue;
        }
        bool dup = false;
        for (uint32_t c = 0; c < q->nof_cand && !dup; c++) {
          dup = q->cand[c].nof_iv == 2 && q->cand[c].iv[0].first == a && q->cand[c].iv[1].first == b &&
                q->cand[c].iv[0].last == a + half - 1;
        }
        if (!dup) {
          iv[0].first = a, iv[0].last = a + half - 1;
          iv[1].first = b, iv[1].last = b + half - 1;
          gs_add_cand(q, 2, iv);
        }
      }
    }
  }
  /* interleaved: triples of equal length, or the longest is the sum of the others */
  for (uint32_t i = 0; i + 2 < n; i++) {
    for (uint32_t j = i + 1; j + 1 < n; j++) {
      for (uint32_t k = j + 1; k < n; k++) {
        int32_t li = r[i].last - r[i].first + 1, lj = r[j].last - r[j].first + 1, lk = r[k].last - r[k].first + 1;
        if (li <= 0 || lj <= 0 || lk <= 0 || li + lj + lk != nsc_cand) {
          continue;
        }
        bool shape = (li == lj && li == lk) || li == lj + lk || lj == li + lk || lk == li + lj;
        bool align = (r[i].first - r[j].first) % nsc_cce == 0 && (r[i].first - r[k].first) % nsc_cce == 0;
        if (shape && align) {
          iv[0].first = r[i].first, iv[0].last = r[i].last;
          iv[1].first = r[j].first, iv[1].last = r[j].last;
          iv[2].first = r[k].first, iv[2].last = r[k].last;
          gs_add_cand(q, 3, iv);
        }
      }
    }
  }
}

/* ---------------------------------------------------------- update_candidates.m */

static void gs_remove_cand(gs_t* q, const gs_cand_t* c)
{
  uint32_t n = q->nof_runs;
  for (uint32_t i = 0; i < n; i++) {
    for (uint32_t k = 0; k < c->nof_iv; k++) {
      gs_run_t*      r  = &q->runs[i];
      const gs_iv_t* iv = &c->iv[k];
      if (iv->first < r->first || iv->last > r->last || r->first > r->last) {
        continue;
      }
      if (iv->first == r->first && iv->last == r->last) {
        r->first = r->last + 1;
      } else if (iv->first == r->first) {
        r->first = iv->last + 1;
      } else if (iv->last == r->last) {
        r->last = iv->first - 1;
      } else if (q->nof_runs < GS_MAX_RUNS) {
        q->runs[q->nof_runs].first = iv->last + 1;
        q->runs[q->nof_runs].last  = r->last;
        q->runs[q->nof_runs].pow   = r->pow;
        q->nof_runs++;
        r->last = iv->first - 1;
      }
    }
  }
  uint32_t m = 0;
  for (uint32_t i = 0; i < q->nof_runs; i++) {
    if (q->runs[i].last >= q->runs[i].first) {
      q->runs[m++] = q->runs[i];
    }
  }
  q->nof_runs = m;
}

/* ------------------------------------------------------------- PDCCH_equalize.m */

/* Channel from the DM-RS of symbol l_est (dm: the transmitted DM-RS, one per 4
 * subcarriers from first + 1), smoothed by keeping the first N/2 taps after removing
 * the mean phase and the delay slope; zero forcing on the data REs of every symbol.
 * yeq holds 3 len / 4 REs per symbol. */
static void gs_equalize(const gs_slot_t* s, int32_t first, uint32_t len, uint32_t dur, const cf_t* dm, uint32_t l_est,
                        cf_t* yeq)
{
  const uint32_t Nd = len / 4, Nf = len, half = Nd / 2;
  const uint32_t N  = s->nof_sc;
  cf_t  hd[Nd], w[Nf], h[Nf], Hs[Nf];
  float ph[Nd];

  const cf_t* Yl = s->grid + l_est * N + first;
  for (uint32_t i = 0; i < Nd; i++) {
    hd[i] = Yl[1 + 4 * i] / dm[i];
  }
  ph[0] = cargf(hd[0]);
  for (uint32_t i = 1; i < Nd; i++) {
    float d = cargf(hd[i]) - cargf(hd[i - 1]);
    d -= 2.0f * (float)M_PI * roundf(d / (2.0f * (float)M_PI));
    ph[i] = ph[i - 1] + d;
  }
  float slope = Nd > 1 ? (ph[Nd - 1] - ph[0]) / (float)(Nd - 1) : 0.0f;
  float mean  = 0;
  for (uint32_t i = 0; i < Nd; i++) {
    mean += ph[i];
  }
  mean /= (float)Nd;

  for (uint32_t m = 0; m < Nf; m++) {
    w[m] = cexpf(-2.0f * I * (float)M_PI * (float)m / (float)Nf);
  }
  /* kept taps: [0, half - 1] and [Nf - half - 1, Nf - 1], the two edge ones halved */
  for (uint32_t t = 0; t < Nf; t++) {
    h[t] = 0;
    if (!(t < half || t >= Nf - half - 1)) {
      continue;
    }
    cf_t acc = 0;
    for (uint32_t i = 0; i < Nd; i++) {
      uint32_t k = 1 + 4 * i;
      cf_t     H = hd[i] * cexpf(-I * (mean + (float)k * slope / 4.0f));
      acc += H * conjf(w[(uint64_t)k * t % Nf]);
    }
    h[t] = acc / (float)Nf;
    if (t == half - 1 || t == Nf - half - 1) {
      h[t] *= 0.5f;
    }
  }
  for (uint32_t k = 0; k < Nf; k++) {
    if (k % 4 == 1) {
      continue;
    }
    cf_t acc = 0;
    for (uint32_t t = 0; t < Nf; t++) {
      if (t < half || t >= Nf - half - 1) {
        acc += h[t] * w[(uint64_t)k * t % Nf];
      }
    }
    Hs[k] = 4.0f * acc * cexpf(I * (mean + (float)k * slope / 4.0f));
  }
  for (uint32_t l = 0; l < dur; l++) {
    const cf_t* Y = s->grid + l * N + first;
    uint32_t    m = 0;
    for (uint32_t k = 0; k < Nf; k++) {
      if (k % 4 != 1) {
        yeq[l * (3 * len / 4) + m++] = Y[k] / Hs[k];
      }
    }
  }
}

/* ------------------------------------------------------- blind_search_for_DCIs.m */

typedef struct {
  bool     ok;
  uint32_t n_id;
  bool     cs0;
  uint32_t rot;
  uint32_t l_est;
} gs_dmrs_res_t;

/* N_ID of one interval of a candidate from its DM-RS, against one reference point */
static bool gs_interval_dmrs(const gs_slot_t* s, const gs_iv_t* iv, uint32_t dur, int32_t ref, cf_t cdm[][GS_MAX_SC / 4],
                             gs_dmrs_res_t* res)
{
  const int32_t len = iv->last - iv->first + 1;
  const int32_t k0  = iv->first - ref;
  if (ref == GS_NO_REF || k0 < 0 || k0 % 4 != 0) {
    return false;
  }
  const uint32_t D  = (uint32_t)k0 / 2; /* 2 sequence bits per DM-RS, one DM-RS per 4 subcarriers */
  const uint32_t Nd = (uint32_t)len / 4;

  if (len >= GS_NSC_CCE) {
    /* direct method: one symbol is enough */
    for (uint32_t l = 0; l < dur; l++) {
      uint32_t c_init[4], rot[4];
      int      nf = gs_dmrs_process(cdm[l], Nd, D, c_init, rot);
      for (int i = 0; i < nf; i++) {
        uint32_t n_id;
        if (gs_pdcch_cinit_check(c_init[i], s->slot_idx, l, &n_id)) {
          res->ok    = true;
          res->n_id  = n_id;
          res->rot   = rot[i];
          res->l_est = l;
          return true;
        }
      }
    }
    return false;
  }
  if ((uint32_t)len * dur >= GS_NSC_CCE && dur > 1 && 2 * Nd <= 64) {
    /* underdetermined: all symbols together */
    cf_t all[GS_MAX_DURATION * 32];
    for (uint32_t l = 0; l < dur; l++) {
      memcpy(&all[l * Nd], cdm[l], Nd * sizeof(cf_t));
    }
    uint32_t n_id, rot;
    if (gs_dmrs_solve_underdet(all, Nd, dur, D, s->slot_idx, &n_id, &rot)) {
      res->ok    = true;
      res->n_id  = n_id;
      res->rot   = rot;
      res->l_est = 0;
      return true;
    }
  }
  return false;
}

static void gs_learn_spacing(gs_t* q, const gs_cand_t* c)
{
  if (c->nof_iv != 2 || c->iv[0].last - c->iv[0].first != c->iv[1].last - c->iv[1].first) {
    return;
  }
  int32_t sp = c->iv[1].first - c->iv[0].first;
  for (uint32_t i = 0; i < q->nof_il; i++) {
    if (q->il_spacing[i] == sp) {
      q->il_hits[i]++;
      return;
    }
  }
  uint32_t slot = q->nof_il < 4 ? q->nof_il++ : 3; /* replace the last (least seen after sorting) */
  q->il_spacing[slot] = sp;
  q->il_hits[slot]    = 1;
  for (uint32_t i = slot; i > 0 && q->il_hits[i] > q->il_hits[i - 1]; i--) {
    int32_t  t = q->il_spacing[i];
    uint32_t h = q->il_hits[i];
    q->il_spacing[i] = q->il_spacing[i - 1], q->il_hits[i] = q->il_hits[i - 1];
    q->il_spacing[i - 1] = t, q->il_hits[i - 1] = h;
  }
}

static bool gs_try_candidate(gs_t* q, const gs_slot_t* s, const gs_cand_t* c, uint32_t dur, uint32_t Ll, gs_dci_t* out)
{
  cf_t (*cdm)[GS_MAX_SC / 4]                          = q->cdm;
  cf_t (*yeq_iv)[GS_MAX_DURATION * GS_MAX_SC * 3 / 4] = q->yeq_iv;
  uint8_t* bits                                       = q->bits;
  const uint32_t          N = s->nof_sc;

  uint32_t total = 0;
  for (uint32_t i = 0; i < c->nof_iv; i++) {
    total += (uint32_t)(c->iv[i].last - c->iv[i].first + 1);
  }
  if (total * dur != ((uint32_t)GS_NSC_CCE << Ll)) {
    return false;
  }

  int32_t refs[2];
  uint32_t nof_refs = 0;
  if (s->ref_crb0 != GS_NO_REF) {
    refs[nof_refs++] = s->ref_crb0;
  }
  if (s->ref_cs0 != GS_NO_REF && s->ref_cs0 != s->ref_crb0) {
    refs[nof_refs++] = s->ref_cs0;
  }

  int64_t  n_id   = -1;
  bool     cs0    = false;
  double   epre   = 0;
  for (uint32_t i = 0; i < c->nof_iv; i++) {
    const gs_iv_t* iv  = &c->iv[i];
    const uint32_t len = (uint32_t)(iv->last - iv->first + 1);
    if (len % GS_NSC_REG != 0 || iv->first < 0 || (uint32_t)iv->last >= N) {
      return false;
    }
    /* DM-RS REs relative to the ambiguous channel: on the QPSK axes */
    for (uint32_t l = 0; l < dur; l++) {
      for (uint32_t m = 0; m < len / 4; m++) {
        uint32_t k = (uint32_t)iv->first + 1 + 4 * m;
        cdm[l][m]  = s->grid[l * N + k] / q->h_amb[l][k];
      }
      for (uint32_t k = (uint32_t)iv->first; k <= (uint32_t)iv->last; k++) {
        float a = cabsf(s->grid[l * N + k]);
        epre += a * a;
      }
    }
    gs_dmrs_res_t res = {0};
    uint32_t      r;
    for (r = 0; r < nof_refs; r++) {
      if (gs_interval_dmrs(s, iv, dur, refs[r], cdm, &res)) {
        break;
      }
    }
    if (!res.ok) {
      if (gs_debug()) {
        printf("GS   cand L%u [%d..%d]: no DM-RS c_init\n", 1U << Ll, iv->first, iv->last);
      }
      return false;
    }
    if (gs_debug()) {
      printf("GS   cand L%u [%d..%d]: N_ID %u (%s ref, symbol %u)\n", 1U << Ll, iv->first, iv->last, res.n_id,
             refs[r] == s->ref_cs0 ? "CORESET#0" : "point A", res.l_est);
    }
    bool iv_cs0 = refs[r] == s->ref_cs0;
    if (i == 0) {
      n_id = res.n_id;
      cs0  = iv_cs0;
    } else if (n_id != (int64_t)res.n_id || cs0 != iv_cs0) {
      return false;
    }
    /* transmitted DM-RS of the symbol the channel is estimated on */
    cf_t       dm[GS_MAX_SC / 4];
    const cf_t rot = cexpf(I * (float)(M_PI / 4.0 + res.rot * M_PI / 2.0));
    for (uint32_t m = 0; m < len / 4; m++) {
      dm[m] = cdm[res.l_est][m] * rot;
    }
    gs_equalize(s, iv->first, len, dur, dm, res.l_est, yeq_iv[i]);
  }

  /* hard bits, REs in frequency order within a symbol, then the next symbol */
  uint32_t E = 0;
  double   err = 0, ref_pow = 0;
  for (uint32_t l = 0; l < dur; l++) {
    for (uint32_t i = 0; i < c->nof_iv; i++) {
      uint32_t len = (uint32_t)(c->iv[i].last - c->iv[i].first + 1);
      uint32_t nre = 3 * len / 4;
      for (uint32_t m = 0; m < nre; m++) {
        cf_t v      = yeq_iv[i][l * nre + m];
        bits[E++]   = crealf(v) < 0.0f;
        bits[E++]   = cimagf(v) < 0.0f;
        cf_t d      = ((crealf(v) < 0 ? -1.0f : 1.0f) + I * (cimagf(v) < 0 ? -1.0f : 1.0f)) * (float)M_SQRT1_2;
        float e     = cabsf(v - d);
        err += e * e;
        ref_pow += 1.0;
      }
    }
  }
  if (E != (108U << Ll)) {
    return false;
  }

  /* Could this candidate be in CORESET#0? Its DM-RS says so (or cannot tell, the two
    references coinciding), and it lies inside CORESET#0's band when that is known. A
    CORESET 1 DCI next to it must not be held to CORESET#0's DCI sizes. */
  bool in_cs0 = cs0 || (s->ref_cs0 != GS_NO_REF && s->ref_cs0 == s->ref_crb0);
  if (in_cs0 && s->cs0_nof_sc > 0) {
    for (uint32_t i = 0; i < c->nof_iv; i++) {
      in_cs0 = in_cs0 && c->iv[i].first >= s->ref_cs0 && c->iv[i].last < s->ref_cs0 + (int32_t)s->cs0_nof_sc;
    }
  }

  /* n_RNTI = RNTI decodes are strong (the CRC must match the RNTI the scrambling
    used): the first one wins. n_RNTI = 0 decodes are weak across sizes (nested polar
    codes): every size that solves is kept, and the one seen most often wins, the
    smallest on a tie. */
  bool     have_weak = false;
  uint32_t weak_nb = 0, weak_count = 0;
  uint16_t weak_rnti = 0;
  uint8_t  weak_payload[GS_MAX_DCI_BITS];
  uint8_t  payload[GS_MAX_DCI_BITS];
  bool     strong   = false;
  uint32_t final_nb = 0;
  uint16_t final_rnti = 0;
  bool     final_n_rnti = false;
  for (uint32_t si = 0; si < q->nof_sizes && !strong; si++) {
    uint32_t nb = q->sizes[si];
    if (!gs_gf2_size_valid(nb, Ll) || nb + 16 + q->args.min_parity > E) {
      continue;
    }
    /* n_RNTI = RNTI is the likely case when the ID is not the PCI */
    bool first_hyp = n_id != (int64_t)s->pci;
    for (uint32_t h = 0; h < 2; h++) {
      bool     n_rnti = h == 0 ? first_hyp : !first_hyp;
      uint16_t rnti;
      /* the size list of the common search spaces applies where their DM-RS is: CORESET#0 */
      if (!n_rnti && in_cs0 && !q->css_ok[nb]) {
        continue;
      }
      if (!gs_dci_try(bits, E, nb, (uint32_t)n_id, n_rnti, payload, &rnti)) {
        continue;
      }
      if (n_rnti) {
        strong       = true;
        final_nb     = nb;
        final_rnti   = rnti;
        final_n_rnti = true;
        memcpy(out->payload, payload, nb);
        break;
      }
      if (!have_weak || q->size_count[nb] > weak_count) {
        have_weak  = true;
        weak_nb    = nb;
        weak_count = q->size_count[nb];
        weak_rnti  = rnti;
        memcpy(weak_payload, payload, nb);
      }
      break; /* both hypotheses cannot hold */
    }
  }
  if (!strong && !have_weak) {
    if (gs_debug()) {
      printf("GS   cand L%u: no DCI size solves (evm %.2f)\n", 1U << Ll, sqrt(err / (ref_pow > 0 ? ref_pow : 1)));
    }
    return false;
  }
  if (!strong) {
    final_nb     = weak_nb;
    final_rnti   = weak_rnti;
    final_n_rnti = false;
    memcpy(out->payload, weak_payload, weak_nb);
  }
  out->slot_idx      = s->slot_idx;
  out->duration      = dur;
  out->L             = 1U << Ll;
  out->cand          = *c;
  out->coreset0      = in_cs0;
  out->ref_ambiguous = s->ref_cs0 != GS_NO_REF && s->ref_cs0 == s->ref_crb0;
  out->n_id          = (uint32_t)n_id;
  out->rnti          = final_rnti;
  out->n_rnti        = final_n_rnti;
  out->nof_bits      = final_nb;
  out->evm           = (float)sqrt(err / (ref_pow > 0 ? ref_pow : 1));
  out->epre          = (float)(epre / ((double)total * dur));
  q->size_count[final_nb]++;
  gs_learn_spacing(q, c);
  return true;
}

static uint32_t gs_max_L_log2(const gs_t* q, uint32_t dur)
{
  uint32_t sc = 0;
  for (uint32_t i = 0; i < q->nof_runs; i++) {
    sc += (uint32_t)(q->runs[i].last - q->runs[i].first + 1);
  }
  uint32_t ncce = sc * dur / GS_NSC_CCE;
  if (ncce == 0) {
    return UINT32_MAX;
  }
  uint32_t Ll = 0;
  while ((2U << Ll) <= ncce) {
    Ll++;
  }
  return Ll > q->args.max_L_log2 ? q->args.max_L_log2 : Ll;
}

int gs_search_slot(gs_t* q, const gs_slot_t* s, gs_dci_t* out, uint32_t max_out)
{
  if (q == NULL || s == NULL || s->grid == NULL || s->nof_sc == 0 || s->nof_sc > GS_MAX_SC) {
    return -1;
  }
  const uint32_t dur = gs_power_detect(q, s);
  if (q->nof_runs == 0) {
    return 0;
  }
  const uint32_t N = s->nof_sc;
  if (gs_debug()) {
    printf("GS slot %u: duration %u, %u runs:", s->slot_idx, dur, q->nof_runs);
    for (uint32_t i = 0; i < q->nof_runs; i++) {
      printf(" [%d..%d]", q->runs[i].first, q->runs[i].last);
    }
    printf("\n");
  }

  /* channel estimates up to j^k per run, aligned across symbols */
  cf_t* tmp = q->tmp;
  for (uint32_t i = 0; i < q->nof_runs; i++) {
    uint32_t f = (uint32_t)q->runs[i].first, n = (uint32_t)(q->runs[i].last - q->runs[i].first + 1);
    gs_unwrap_qpsk(s->grid + f, n, &q->h_amb[0][f]);
    for (uint32_t l = 1; l < dur; l++) {
      gs_unwrap_qpsk(s->grid + l * N + f, n, tmp);
      cf_t  num = 0;
      float den = 0;
      for (uint32_t k = 0; k < n; k++) {
        num += conjf(tmp[k]) * q->h_amb[0][f + k];
        den += crealf(tmp[k] * conjf(tmp[k]));
      }
      cf_t rho = den > 0 ? num / den : 1.0f;
      cf_t rot = cexpf(I * (float)M_PI / 2.0f * roundf(cargf(rho) * 2.0f / (float)M_PI));
      for (uint32_t k = 0; k < n; k++) {
        q->h_amb[l][f + k] = tmp[k] * rot;
      }
    }
  }

  uint32_t nof_dci = 0;
  uint32_t Ll      = gs_max_L_log2(q, dur);
  while (Ll != UINT32_MAX && nof_dci < max_out) {
    gs_enum_candidates(q, dur, Ll);
    bool hit = false;
    for (uint32_t c = 0; c < q->nof_cand; c++) {
      if (gs_try_candidate(q, s, &q->cand[c], dur, Ll, &out[nof_dci])) {
        gs_cand_t used = q->cand[c];
        nof_dci++;
        gs_remove_cand(q, &used);
        hit = true;
        break;
      }
    }
    if (hit) {
      uint32_t m = gs_max_L_log2(q, dur);
      if (m == UINT32_MAX) {
        break;
      }
      Ll = m < Ll ? m : Ll;
      continue;
    }
    Ll = Ll == 0 ? UINT32_MAX : Ll - 1;
  }
  return (int)nof_dci;
}
