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
#include "gs_ec.h"
#include "srsran/phy/common/sequence.h"
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
  a->soft         = 1;
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
  /* DM-RS channel with the mean phase and the slope taken out, once */
  cf_t Hn[Nd];
  for (uint32_t i = 0; i < Nd; i++) {
    Hn[i] = hd[i] * cexpf(-I * (mean + (float)(1 + 4 * i) * slope / 4.0f));
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
      acc += Hn[i] * conjf(w[(uint64_t)k * t % Nf]);
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

/* Channel of one symbol from its own DM-RS, in O(len): the delay slope from the pilot
 * autocorrelation (noise-robust, no unwrapping), the mean phase, then each pilot's
 * derotated channel averaged over +-GS_EQ_HALF pilots (about 2 RBs), linearly
 * interpolated to the data REs; zero forcing. Unlike gs_equalize it does not ring at the
 * edges of wide candidates (AL8/16) and costs a few operations per RE. yeq: 3 len / 4. */
#define GS_EQ_HALF 3
static void gs_equalize_fast(const gs_slot_t* s, int32_t first, uint32_t len, uint32_t l, const cf_t* dm, cf_t* yeq,
                             float* w)
{
  const uint32_t Nd = len / 4;
  const cf_t*    Y  = s->grid + l * s->nof_sc + first;
  cf_t           hd[GS_MAX_SC / 4], hs[GS_MAX_SC / 4];
  cf_t           ac = 0;
  for (uint32_t i = 0; i < Nd; i++) {
    hd[i] = Y[1 + 4 * i] * conjf(dm[i]) * 2.0f; /* |dm|^2 = 1/2 */
    if (i > 0) {
      ac += hd[i] * conjf(hd[i - 1]);
    }
  }
  const float slope = Nd > 1 ? cargf(ac) / 4.0f : 0.0f; /* per subcarrier */
  for (uint32_t i = 0; i < Nd; i++) {
    hd[i] *= cexpf(-I * slope * (float)(1 + 4 * i));
  }
  /* sliding window over the derotated pilots */
  cf_t     run = 0;
  uint32_t lo = 0, hi = 0; /* window [lo, hi) */
  for (uint32_t i = 0; i < Nd; i++) {
    const uint32_t a = i > GS_EQ_HALF ? i - GS_EQ_HALF : 0, b = i + GS_EQ_HALF + 1 < Nd ? i + GS_EQ_HALF + 1 : Nd;
    while (hi < b) {
      run += hd[hi++];
    }
    while (lo < a) {
      run -= hd[lo++];
    }
    hs[i] = run / (float)(hi - lo);
  }
  uint32_t m = 0;
  for (uint32_t k = 0; k < len; k++) {
    if (k % 4 == 1) {
      continue;
    }
    /* pilot positions 1 + 4 i: interpolate between the two around k */
    float    x  = ((float)k - 1.0f) / 4.0f;
    int32_t  i0 = (int32_t)floorf(x);
    if (i0 < 0) {
      i0 = 0;
    }
    if (i0 > (int32_t)Nd - 2) {
      i0 = Nd >= 2 ? (int32_t)Nd - 2 : 0;
    }
    const float f = Nd >= 2 ? x - (float)i0 : 0.0f;
    cf_t        h = Nd >= 2 ? hs[i0] * (1.0f - f) + hs[i0 + 1] * f : hs[0];
    h *= cexpf(I * slope * (float)k);
    if (w != NULL) {
      w[m] = crealf(h * conjf(h)); /* the equalized RE's noise is 1/|h|^2 of the received one */
    }
    yeq[m++] = cabsf(h) > 0 ? Y[k] / h : 0;
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
    codes): a DCI also solves each larger size with probability 2^-8, a smaller one with
    2^-(8 + difference). The smallest size that solves wins, unless a larger one gives an
    RNTI already known (gs_set_known_rnti) and the smallest does not. (The most frequent
    size used to win: with many UEs of different sizes the histogram is flat and a wrong
    larger size won, 1.9% of the DCIs of a 100-UE synthetic cell had a false RNTI.) */
  bool     have_weak = false;
  uint32_t weak_nb = 0, weak_count = 0;
  uint16_t weak_rnti = 0;
  bool     weak_known = false;
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
      {
        const bool known = gs_known_rnti(rnti);
        /* sizes are tried in ascending order (gs_create): the first is the smallest */
        if (!have_weak || (known && !weak_known)) {
          have_weak  = true;
          weak_nb    = nb;
          weak_count = q->size_count[nb];
          weak_rnti  = rnti;
          weak_known = known;
          memcpy(weak_payload, payload, nb);
        }
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

/* RNTIs GoldenSniffer knows (learned, converged or from the RACH): a bitset written by
 * the scheduler, read by the searchers without a lock (a stale bit only changes which of
 * two weak decodes is kept) */
static uint64_t gs_known_bits[65536 / 64];

void gs_set_known_rnti(uint16_t rnti, bool known)
{
  const uint64_t m = 1ULL << (rnti & 63);
  if (known) {
    __atomic_fetch_or(&gs_known_bits[rnti >> 6], m, __ATOMIC_RELAXED);
  } else {
    __atomic_fetch_and(&gs_known_bits[rnti >> 6], ~m, __ATOMIC_RELAXED);
  }
}

bool gs_known_rnti(uint16_t rnti)
{
  return (__atomic_load_n(&gs_known_bits[rnti >> 6], __ATOMIC_RELAXED) >> (rnti & 63)) & 1;
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

/* --------------------------------------------------------- tracked UE, no search */

#define GS_TRACK_GATE 0.6f
#define GS_SOFT_MAX 8 /* soft decisions per group and slot, at most */ /* DM-RS correlation a candidate needs before it is decoded */

int gs_decode_tracked(gs_t* q, const gs_slot_t* s, const gs_track_t* t, gs_dci_t* out, uint32_t max_out)
{
  if (q == NULL || s == NULL || t == NULL || s->grid == NULL || t->duration < 1 || t->duration > GS_MAX_DURATION ||
      s->nof_sc == 0 || s->nof_sc > GS_MAX_SC) {
    return -1;
  }
  const int32_t ref = t->coreset0 ? s->ref_cs0 : s->ref_crb0;
  if (ref == GS_NO_REF) {
    return 0;
  }
  const uint32_t N   = s->nof_sc;
  const uint32_t dur = t->duration;
  const int32_t  cce = GS_NSC_CCE / (int32_t)dur;
  cf_t(*dm)[GS_MAX_SC / 4] = q->cdm; /* expected DM-RS per symbol */

  int32_t  used_lo[GS_MAX_DCI_PER_SLOT], used_hi[GS_MAX_DCI_PER_SLOT];
  uint32_t nof = 0;

  for (int Ll = GS_MAX_L_LOG2; Ll >= 0 && nof < max_out; Ll--) {
    if (!((t->al_mask >> Ll) & 1U)) {
      continue;
    }
    const int32_t span = (1 << Ll) * cce;
    for (int32_t first = t->sc_lo; first + span - 1 <= t->sc_hi && nof < max_out; first += cce) {
      if (first < 0 || first + span > (int32_t)N || span % GS_NSC_REG != 0) {
        continue;
      }
      bool overlap = false;
      for (uint32_t i = 0; i < nof && !overlap; i++) {
        overlap = first <= used_hi[i] && first + span - 1 >= used_lo[i];
      }
      const int32_t k0 = first - ref;
      if (overlap || k0 < 0 || k0 % 4 != 0) {
        continue;
      }
      const uint32_t nd = (uint32_t)span / 4;

      /* expected DM-RS, and how well the REs match it, coherently over 2 RBs at a time */
      double num = 0, den = 0;
      for (uint32_t l = 0; l < dur; l++) {
        /* as srsRAN generates it (dmrs_pdcch.c): sequence from c_init, k0/4 pilots in */
        srsran_sequence_state_t seq;
        srsran_sequence_state_init(&seq, gs_pdcch_dmrs_cinit(s->slot_idx, l, t->n_id));
        srsran_sequence_state_advance(&seq, (uint32_t)k0 / 2);
        srsran_sequence_state_gen_f(&seq, (float)M_SQRT1_2, (float*)dm[l], 2 * nd);
        const cf_t* Y = s->grid + l * N + first;
        for (uint32_t g = 0; g < nd; g += 6) {
          cf_t   acc = 0;
          double pw  = 0;
          uint32_t n = 0;
          for (uint32_t m = g; m < g + 6 && m < nd; m++, n++) {
            cf_t y = Y[1 + 4 * m];
            acc += y * conjf(dm[l][m]);
            pw += crealf(y * conjf(y));
          }
          num += cabsf(acc);
          den += sqrt(pw * n);
        }
      }
      if (gs_debug()) {
        printf("GS track 0x%04x L%u [%d..%d] gate %.2f\n", t->rnti, 1U << Ll, first, first + span - 1,
               den > 0 ? num / den : 0.0);
      }
      if (den <= 0 || num / den < GS_TRACK_GATE) {
        continue;
      }

      /* each symbol with its own pilots: a CFO turns the channel from one symbol to the next */
      cf_t*          yeq = q->yeq_iv[0];
      const uint32_t nre = 3 * (uint32_t)span / 4;
      for (uint32_t l = 0; l < dur; l++) {
        gs_equalize_fast(s, first, (uint32_t)span, l, dm[l], yeq + l * nre, NULL);
      }
      uint32_t       E   = 0;
      double         err = 0, epre = 0;
      for (uint32_t l = 0; l < dur; l++) {
        for (uint32_t m = 0; m < nre; m++) {
          cf_t v         = yeq[l * nre + m];
          q->bits[E++]   = crealf(v) < 0.0f;
          q->bits[E++]   = cimagf(v) < 0.0f;
          cf_t  d        = ((crealf(v) < 0 ? -1.0f : 1.0f) + I * (cimagf(v) < 0 ? -1.0f : 1.0f)) * (float)M_SQRT1_2;
          float e        = cabsf(v - d);
          err += e * e;
        }
        for (int32_t k = first; k < first + span; k++) {
          float a = cabsf(s->grid[l * N + k]);
          epre += a * a;
        }
      }
      if (E != (108U << Ll)) {
        continue;
      }
      for (uint32_t z = 0; z < t->nof_sizes; z++) {
        const uint32_t nb = t->sizes[z];
        uint16_t       rnti;
        if (!gs_gf2_size_valid(nb, (uint32_t)Ll) || nb + 16 + q->args.min_parity > E) {
          continue;
        }
        if (!gs_dci_try(q->bits, E, nb, t->n_id, t->n_rnti, out[nof].payload, &rnti) || rnti != t->rnti) {
          continue;
        }
        gs_dci_t* o      = &out[nof];
        o->slot_idx      = s->slot_idx;
        o->duration      = dur;
        o->L             = 1U << Ll;
        o->cand.nof_iv   = 1;
        o->cand.iv[0].first = first;
        o->cand.iv[0].last  = first + span - 1;
        o->coreset0      = t->coreset0;
        o->ref_ambiguous = s->ref_cs0 != GS_NO_REF && s->ref_cs0 == s->ref_crb0;
        o->n_id          = t->n_id;
        o->rnti          = rnti;
        o->n_rnti        = t->n_rnti;
        o->nof_bits      = nb;
        o->evm           = (float)sqrt(err / (double)(nre * dur));
        o->epre          = (float)(epre / ((double)span * dur));
        used_lo[nof]     = first;
        used_hi[nof]     = first + span - 1;
        nof++;
        break;
      }
    }
  }
  return (int)nof;
}


/* ------------------------------------------------------------------ many tracked UEs */

static bool gs_track_same_group(const gs_track_t* a, const gs_track_t* b, int32_t ref, int32_t cce)
{
  return a->n_id == b->n_id && a->n_rnti == b->n_rnti && a->duration == b->duration &&
         a->coreset0 == b->coreset0 && ((a->sc_lo - ref) % cce + cce) % cce == ((b->sc_lo - ref) % cce + cce) % cce;
}

static bool gs_track_has_size(const gs_track_t* t, uint32_t nb)
{
  for (uint32_t z = 0; z < t->nof_sizes; z++) {
    if (t->sizes[z] == nb) {
      return true;
    }
  }
  return false;
}

/* Soft decision for one candidate of a tracked group that passed the DM-RS gate but whose
 * hard decision solved no size: LLRs from the equalized REs weighted by |h|^2, the group's
 * sizes smallest first, and only an RNTI of a tracked UE at one of its sizes is accepted
 * (the soft decoder bounds false decodes at 2^-min_parity, the RNTI check adds 2^-16). */
static int gs_tracked_soft(gs_t*             q,
                           const gs_slot_t*  s,
                           const gs_track_t* tracks,
                           uint32_t          nt,
                           const gs_track_t* t0,
                           int32_t           ref,
                           int32_t           cce,
                           int32_t           first,
                           int               Ll,
                           const uint32_t*   sizes,
                           uint32_t          nsz,
                           gs_dci_t*         o)
{
  const uint32_t dur = t0->duration, span = (uint32_t)((1 << Ll) * cce), nre = 3 * span / 4;
  const int32_t  k0  = first - ref;
  static __thread float llr[108 << GS_MAX_L_LOG2];
  static __thread float w[GS_MAX_DURATION * GS_MAX_SC * 3 / 4];
  cf_t*                 yeq = q->yeq_iv[0];
  for (uint32_t l = 0; l < dur; l++) {
    gs_equalize_fast(s, first, span, l, q->cdm[l] + k0 / 4, yeq + l * nre, w + l * nre);
  }
  uint32_t E = 0;
  for (uint32_t m = 0; m < dur * nre; m++) {
    llr[E++] = crealf(yeq[m]) * w[m];
    llr[E++] = cimagf(yeq[m]) * w[m];
  }
  if (E != (108U << Ll)) {
    return 0;
  }
  for (uint32_t z = 0; z < nsz; z++) {
    const uint32_t nb = sizes[z];
    uint16_t       rnti;
    float          metric;
    if (!gs_gf2_size_valid(nb, (uint32_t)Ll) || nb + 16 + q->args.min_parity > E) {
      continue;
    }
    if (!gs_dci_try_soft(q, llr, E, nb, t0->n_id, t0->n_rnti, t0->n_rnti ? GS_EC_EFFORT_MAX : GS_EC_EFFORT_AUTO,
                         o->payload, &rnti, &metric)) {
      continue;
    }
    for (uint32_t i = 0; i < nt; i++) {
      const gs_track_t* t = &tracks[i];
      if (t->rnti == rnti && gs_track_same_group(t0, t, ref, cce) && gs_track_has_size(t, nb)) {
        o->slot_idx         = s->slot_idx;
        o->duration         = dur;
        o->L                = 1U << Ll;
        o->cand.nof_iv      = 1;
        o->cand.iv[0].first = first;
        o->cand.iv[0].last  = first + (int32_t)span - 1;
        o->coreset0         = t0->coreset0;
        o->ref_ambiguous    = s->ref_cs0 != GS_NO_REF && s->ref_cs0 == s->ref_crb0;
        o->n_id             = t0->n_id;
        o->rnti             = rnti;
        o->n_rnti           = t0->n_rnti;
        o->nof_bits         = nb;
        o->evm              = -1.0f; /* found by the soft decision */
        o->epre             = 0;
        return 1;
      }
    }
    return 0; /* the smallest size that solves is the DCI: another RNTI, not ours */
  }
  return 0;
}

int gs_decode_tracked_multi(gs_t*             q,
                            const gs_slot_t*  s,
                            const gs_track_t* tracks,
                            uint32_t          nof_tracks,
                            bool              full_band,
                            bool              want_untracked,
                            gs_dci_t*         out,
                            uint8_t*          tracked,
                            uint32_t          max_out)
{
  if (q == NULL || s == NULL || tracks == NULL || out == NULL || tracked == NULL || s->grid == NULL ||
      s->nof_sc == 0 || s->nof_sc > GS_MAX_SC) {
    return -1;
  }
  const uint32_t N   = s->nof_sc;
  uint32_t       nof = 0;
  static __thread bool done[GS_MAX_TRACKS];
  const uint32_t nt = nof_tracks < GS_MAX_TRACKS ? nof_tracks : GS_MAX_TRACKS;
  memset(done, 0, sizeof(bool) * nt);

  for (uint32_t g = 0; g < nt && nof < max_out; g++) {
    const gs_track_t* t0 = &tracks[g];
    if (done[g] || t0->duration < 1 || t0->duration > GS_MAX_DURATION) {
      continue;
    }
    const int32_t ref = t0->coreset0 ? s->ref_cs0 : s->ref_crb0;
    const uint32_t dur = t0->duration;
    const int32_t  cce = GS_NSC_CCE / (int32_t)dur;
    if (ref == GS_NO_REF || ref < 0 || t0->sc_lo < ref) {
      done[g] = true;
      continue;
    }

    /* the group: its tracks, the union of their sizes, ALs and extents */
    /* every size of the group (with many UEs, each with its own sizes, 16 were too few:
      the UEs whose sizes did not fit were never decoded) */
    uint32_t sizes[128], nsz = 0, al_mask = 0;
    int32_t  lo = t0->sc_lo, hi = t0->sc_hi;
    for (uint32_t i = g; i < nt; i++) {
      const gs_track_t* t = &tracks[i];
      if (done[i] || !gs_track_same_group(t0, t, ref, cce)) {
        continue;
      }
      done[i] = true;
      al_mask |= t->al_mask;
      lo = t->sc_lo < lo ? t->sc_lo : lo;
      hi = t->sc_hi > hi ? t->sc_hi : hi;
      for (uint32_t z = 0; z < t->nof_sizes; z++) {
        bool have = false;
        for (uint32_t k = 0; k < nsz; k++) {
          have = have || sizes[k] == t->sizes[z];
        }
        if (!have && nsz < 128) {
          sizes[nsz++] = t->sizes[z];
        }
      }
    }
    /* ascending: polar codes are nested, a DCI of K bits is also a valid codeword at every
      larger size (with another RNTI), so the smallest size that solves is the real one */
    for (uint32_t x = 1; x < nsz; x++) {
      for (uint32_t y = x; y > 0 && sizes[y - 1] > sizes[y]; y--) {
        const uint32_t t = sizes[y];
        sizes[y]         = sizes[y - 1];
        sizes[y - 1]     = t;
      }
    }
    if (full_band) {
      lo      = ref + ((t0->sc_lo - ref) % cce);
      hi      = (int32_t)N - 1;
      al_mask = (1U << (GS_MAX_L_LOG2 + 1)) - 1;
    }

    /* DM-RS of the whole band from the reference, and per RB: the coherent sum of its 3
      pilots against it and their power. A candidate's gate is then a sum over its RBs. */
    const uint32_t nrb_tot = ((uint32_t)N - (uint32_t)ref) / 12;
    const uint32_t npil    = nrb_tot * 3;
    if (nrb_tot == 0 || npil > GS_MAX_SC / 4) {
      continue;
    }
    cf_t  acc[GS_MAX_DURATION][GS_MAX_SC / 12];
    float pw[GS_MAX_DURATION][GS_MAX_SC / 12];
    for (uint32_t l = 0; l < dur; l++) {
      srsran_sequence_state_t seq;
      srsran_sequence_state_init(&seq, gs_pdcch_dmrs_cinit(s->slot_idx, l, t0->n_id));
      srsran_sequence_state_gen_f(&seq, (float)M_SQRT1_2, (float*)q->cdm[l], 2 * npil);
      const cf_t* Y = s->grid + l * N + ref;
      for (uint32_t r = 0; r < nrb_tot; r++) {
        cf_t  a = 0;
        float p = 0;
        for (uint32_t j = 0; j < 3; j++) {
          const cf_t y = Y[12 * r + 1 + 4 * j];
          a += y * conjf(q->cdm[l][3 * r + j]);
          p += crealf(y * conjf(y));
        }
        acc[l][r] = a;
        pw[l][r]  = p;
      }
    }

    int32_t used_lo[GS_MAX_DCI_PER_SLOT], used_hi[GS_MAX_DCI_PER_SLOT];
    uint32_t nused = 0;
    struct {
      int32_t first;
      int     Ll;
      float   gate;
    } soft_cand[64];
    uint32_t nsoft = 0;
    for (int Ll = GS_MAX_L_LOG2; Ll >= 0 && nof < max_out; Ll--) {
      if (!((al_mask >> Ll) & 1U)) {
        continue;
      }
      const int32_t span = (1 << Ll) * cce;
      for (int32_t first = lo; first + span - 1 <= hi && first + span <= (int32_t)N && nof < max_out; first += cce) {
        const int32_t k0 = first - ref;
        if (k0 < 0 || k0 % 12 != 0 || span % GS_NSC_REG != 0) {
          continue;
        }
        bool overlap = false;
        for (uint32_t i = 0; i < nused && !overlap; i++) {
          overlap = first <= used_hi[i] && first + span - 1 >= used_lo[i];
        }
        if (overlap) {
          continue;
        }
        const uint32_t r0 = (uint32_t)k0 / 12, nrb = (uint32_t)span / 12;
        if (r0 + nrb > nrb_tot) {
          continue;
        }
        /* coherent over 2 RBs (6 pilots) at a time, as gs_decode_tracked */
        double num = 0, den = 0;
        for (uint32_t l = 0; l < dur; l++) {
          for (uint32_t r = 0; r < nrb; r += 2) {
            cf_t     a = acc[l][r0 + r];
            double   p = pw[l][r0 + r];
            uint32_t n = 3;
            if (r + 1 < nrb) {
              a += acc[l][r0 + r + 1];
              p += pw[l][r0 + r + 1];
              n = 6;
            }
            num += cabsf(a);
            den += sqrt(p * n);
          }
        }
        if (den <= 0 || num / den < GS_TRACK_GATE) {
          continue;
        }

        cf_t*          yeq = q->yeq_iv[0];
        const uint32_t nre = 3 * (uint32_t)span / 4;
        for (uint32_t l = 0; l < dur; l++) {
          gs_equalize_fast(s, first, (uint32_t)span, l, q->cdm[l] + k0 / 4, yeq + l * nre, NULL);
        }
        uint32_t E   = 0;
        double   err = 0, epre = 0;
        for (uint32_t l = 0; l < dur; l++) {
          for (uint32_t m = 0; m < nre; m++) {
            cf_t v       = yeq[l * nre + m];
            q->bits[E++] = crealf(v) < 0.0f;
            q->bits[E++] = cimagf(v) < 0.0f;
            cf_t  d      = ((crealf(v) < 0 ? -1.0f : 1.0f) + I * (cimagf(v) < 0 ? -1.0f : 1.0f)) * (float)M_SQRT1_2;
            float e      = cabsf(v - d);
            err += e * e;
          }
          for (int32_t k = first; k < first + span; k++) {
            float a = cabsf(s->grid[l * N + k]);
            epre += a * a;
          }
        }
        if (E != (108U << Ll)) {
          continue;
        }
        /* the smallest size that solves (nested codes, above); then tracked or not */
        int      best = -1; /* 1: tracked, 0: untracked */
        uint32_t best_nb = 0;
        uint16_t best_rnti = 0;
        uint8_t  pay[GS_MAX_DCI_BITS];
        for (uint32_t z = 0; z < nsz && best < 0; z++) {
          const uint32_t nb = sizes[z];
          uint16_t       rnti;
          if (!gs_gf2_size_valid(nb, (uint32_t)Ll) || nb + 16 + q->args.min_parity > E) {
            continue;
          }
          if (!gs_dci_try(q->bits, E, nb, t0->n_id, t0->n_rnti, pay, &rnti)) {
            continue;
          }
          int is_tracked = 0;
          for (uint32_t i = 0; i < nt && !is_tracked; i++) {
            const gs_track_t* t = &tracks[i];
            is_tracked = t->rnti == rnti && gs_track_same_group(t0, t, ref, cce) && gs_track_has_size(t, nb);
          }
          if (best < 0) {
            best      = is_tracked;
            best_nb   = nb;
            best_rnti = rnti;
            memcpy(out[nof].payload, pay, nb);
          }
        }
        if (best < 0) {
          if (q->args.soft && nsoft < 64) {
            soft_cand[nsoft].first = first;
            soft_cand[nsoft].Ll    = Ll;
            soft_cand[nsoft].gate  = (float)(num / den);
            nsoft++;
          }
          continue;
        }
        if (best == 0 && !want_untracked) {
          continue;
        }
        gs_dci_t* o         = &out[nof];
        o->slot_idx         = s->slot_idx;
        o->duration         = dur;
        o->L                = 1U << Ll;
        o->cand.nof_iv      = 1;
        o->cand.iv[0].first = first;
        o->cand.iv[0].last  = first + span - 1;
        o->coreset0         = t0->coreset0;
        o->ref_ambiguous    = s->ref_cs0 != GS_NO_REF && s->ref_cs0 == s->ref_crb0;
        o->n_id             = t0->n_id;
        o->rnti             = best_rnti;
        o->n_rnti           = t0->n_rnti;
        o->nof_bits         = best_nb;
        o->evm              = (float)sqrt(err / (double)(nre * dur));
        o->epre             = (float)(epre / ((double)span * dur));
        tracked[nof]        = (uint8_t)best;
        if (nused < GS_MAX_DCI_PER_SLOT) {
          used_lo[nused]   = first;
          used_hi[nused++] = first + span - 1;
        }
        nof++;
      }
    }
    /* soft pass: the candidates with the best DM-RS gate first, at most GS_SOFT_MAX, none
      overlapping a DCI already found */
    for (uint32_t done_soft = 0; done_soft < GS_SOFT_MAX && nof < max_out;) {
      int bi = -1;
      for (uint32_t i = 0; i < nsoft; i++) {
        if (soft_cand[i].gate > 0 && (bi < 0 || soft_cand[i].gate > soft_cand[bi].gate)) {
          bi = (int)i;
        }
      }
      if (bi < 0) {
        break;
      }
      const int32_t first = soft_cand[bi].first, span = (1 << soft_cand[bi].Ll) * cce;
      soft_cand[bi].gate  = 0;
      bool overlap        = false;
      for (uint32_t i = 0; i < nused && !overlap; i++) {
        overlap = first <= used_hi[i] && first + span - 1 >= used_lo[i];
      }
      if (overlap) {
        continue;
      }
      done_soft++;
      if (gs_tracked_soft(q, s, tracks, nt, t0, ref, cce, first, soft_cand[bi].Ll, sizes, nsz, &out[nof])) {
        tracked[nof] = 1;
        if (nused < GS_MAX_DCI_PER_SLOT) {
          used_lo[nused]   = first;
          used_hi[nused++] = first + span - 1;
        }
        nof++;
      }
    }
  }
  return (int)nof;
}
