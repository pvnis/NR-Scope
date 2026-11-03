/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
/*
 * \brief Tests of sensing.join_grants (nrscope_sensing_join.h).
 *
 * 1. The join on its own, rows in and rows out: a full-band grant comes out bit for
 *    bit; the same band split over 2..8 UEs (any order, from 4 threads) too; overlaps
 *    keep the first grant; the half-carrier rule is applied after joining; different
 *    ports or DM-RS power never share a row; a precoder step at a junction cuts it; a
 *    gap is cut, or filled within join_max_gap_prb.
 * 2. The continuity check against noise and multipath: how often it cuts two grants
 *    that share one channel, and how often it lets through two that do not (a random
 *    precoder per UE, as a massive-MIMO cell would give).
 * 3. End to end, through the whole sensing pipeline (alignment, history, maps), each in
 *    its own process because the sensing context is a singleton:
 *    a. one full-band UE per slot (the old path) and the same band split over several
 *       UEs per slot (joined): the maps must be byte for byte the same;
 *    b. UEs covering 70% of the carrier, no single one half of it: the old path builds
 *       no map at all, the joined one finds the target.
 */
#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "nrscope/hdr/sensing/nr_ue_map.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"
#include "nrscope/hdr/sensing/nrscope_sensing.h"
#include "nrscope/hdr/sensing/nrscope_sensing_join.h"
#include "srsran/phy/ch_estimation/dmrs_sch.h"
#include "srsran/phy/phch/phch_cfg_nr.h"

#define N_RB 273
#define N_SC (N_RB * 12)
#define NS 14

static int failures = 0;
static void check(bool ok, const char* what, double got, double want)
{
  printf("  [%s] %s: %g (want %g)\n", ok ? " ok " : "FAIL", what, got, want);
  failures += !ok;
}

static uint64_t rng = 88172645463325252ULL;
static double   urand(void)
{
  rng ^= rng << 13;
  rng ^= rng >> 7;
  rng ^= rng << 17;
  return (double)(rng >> 11) / 9007199254740992.0;
}
static cf_t gauss(void)
{
  const double u = urand() + 1e-300, v = urand();
  const double r = sqrt(-log(u));
  return (float)(r * cos(2 * M_PI * v)) + I * (float)(r * sin(2 * M_PI * v)); // E|x|^2 = 1
}

/* A channel of n paths with delays up to tau_max seconds, unit mean power. */
typedef struct {
  int    n;
  double tau[8];
  cf_t   g[8];
} chan_t;
static chan_t chan_make(int n, double tau_max)
{
  chan_t c = {.n = n};
  double p = 0;
  for (int i = 0; i < n; i++) {
    c.tau[i] = i == 0 ? 0.2e-6 : 0.2e-6 + urand() * tau_max;
    c.g[i]   = gauss() * (float)(i == 0 ? 1.0 : exp(-2.0 * urand()));
    p += crealf(c.g[i] * conjf(c.g[i]));
  }
  for (int i = 0; i < n; i++) {
    c.g[i] /= (float)sqrt(p);
  }
  return c;
}
static cf_t chan_at(const chan_t* c, int k)
{
  cf_t h = 0;
  for (int i = 0; i < c->n; i++) {
    h += c->g[i] * cexpf(-I * (float)(2 * M_PI * 30e3 * k * c->tau[i]));
  }
  return h;
}

/* One grant's rows [NS][N_SC]: comb 0, step 2, PRBs [p0, p1), on symbols of mask. */
static cf_t  gH[NS * N_SC];
static bool  gV[NS * N_SC];
static void grant_rows(const chan_t* c, int p0, int p1, cf_t u, double snr_db, uint16_t mask, int comb)
{
  memset(gV, 0, sizeof(gV));
  const float sn = (float)pow(10.0, -snr_db / 20.0);
  for (int l = 0; l < NS; l++) {
    if (!(mask & (1u << l))) {
      continue;
    }
    for (int k = p0 * 12 + comb; k < p1 * 12; k += 2) {
      gH[l * N_SC + k] = u * chan_at(c, k) + (snr_db < 200 ? sn * gauss() : 0);
      gV[l * N_SC + k] = true;
    }
  }
}

/* What a flush handed over */
typedef struct {
  int      n;
  uint16_t ports[8];
  int      layer[8];
  uint16_t mask[8];
  cf_t     H[8][NS * N_SC];
  bool     V[8][NS * N_SC];
} sink_t;
static sink_t out;
static void sink(void* ctx, uint32_t aarx, uint16_t ports, int layer, uint16_t mask, const cf_t* H, const bool* V,
                 uint32_t n_sc)
{
  (void)ctx, (void)aarx;
  sink_t* s = &out;
  if (s->n >= 8) {
    return;
  }
  s->ports[s->n] = ports, s->layer[s->n] = layer, s->mask[s->n] = mask;
  memcpy(s->H[s->n], H, sizeof(cf_t) * NS * n_sc);
  memcpy(s->V[s->n], V, sizeof(bool) * NS * n_sc);
  s->n++;
}

static int add(nrscope_sensing_join_t* j, uint16_t ports, float beta, uint16_t rnti, uint16_t mask)
{
  return nrscope_sensing_join_add(j, 0, ports, 0, 1, beta, rnti, mask, gH, gV, N_SC, 0, 100, 3, 12345);
}

/* Row l of the flush against the reference rows: same mask, same values, bit for bit. */
static bool same_row(int o, int l, const cf_t* rH, const bool* rV)
{
  for (int k = 0; k < N_SC; k++) {
    const bool v = out.V[o][l * N_SC + k];
    if (v != rV[l * N_SC + k] || (v && out.H[o][l * N_SC + k] != rH[l * N_SC + k])) {
      return false;
    }
  }
  return true;
}
static int row_pilots(int o, int l, int* lo, int* hi)
{
  int n = 0;
  *lo = -1, *hi = -1;
  for (int k = 0; k < N_SC; k++) {
    if (out.V[o][l * N_SC + k]) {
      n++;
      *lo = *lo < 0 ? k : *lo;
      *hi = k;
    }
  }
  return n;
}

static cf_t refH[NS * N_SC];
static bool refV[NS * N_SC];

typedef struct {
  nrscope_sensing_join_t* j;
  const chan_t*           c;
  int                     n, first;
  const int*              cut;
  const int*              order;
} thr_arg_t;
static void* thr_add(void* p)
{
  thr_arg_t* a = (thr_arg_t*)p;
  static __thread cf_t H[NS * N_SC];
  static __thread bool V[NS * N_SC];
  for (int i = a->first; i < a->n; i += 4) {
    const int u = a->order[i];
    memset(V, 0, sizeof(V));
    for (int k = a->cut[u] * 12; k < a->cut[u + 1] * 12; k += 2) {
      for (int l = 2; l <= 11; l += 9) {
        H[l * N_SC + k] = refH[l * N_SC + k];
        V[l * N_SC + k] = true;
      }
    }
    nrscope_sensing_join_add(a->j, 0, 1, 0, 1, 1.0f, (uint16_t)(0x4600 + u), (1u << 2) | (1u << 11), H, V, N_SC, 0,
                             100, 3, 12345);
  }
  return NULL;
}

static void test_rows(void)
{
  printf("\nJoin, rows in and out\n");
  nrscope_sensing_join_args.max_gap_prb = 0;
  nrscope_sensing_join_args.check       = true;
  nrscope_sensing_join_t* j             = nrscope_sensing_join_alloc(0, N_SC);
  const uint16_t          mask          = (1u << 2) | (1u << 11);
  chan_t                  c             = chan_make(4, 0.5e-6);

  // 1. one full-band grant
  grant_rows(&c, 0, N_RB, 1.0f, 30, mask, 0);
  memcpy(refH, gH, sizeof(gH)), memcpy(refV, gV, sizeof(gV));
  nrscope_sensing_join_begin(j, 100, 3);
  add(j, 1, 1.0f, 0x4601, mask);
  out.n = 0;
  nrscope_sensing_join_flush_to(j, sink, NULL);
  check(out.n == 1 && out.mask[0] == mask && same_row(0, 2, refH, refV) && same_row(0, 11, refH, refV),
        "full-band grant comes out unchanged (bit exact)", out.n, 1);

  // 2. the same band over 2..8 UEs, random cut points and order, 200 trials
  int exact = 0;
  for (int t = 0; t < 200; t++) {
    const int n = 2 + t % 7;
    int       cut[10];
    cut[0] = 0, cut[n] = N_RB;
    for (int i = 1; i < n; i++) {
      cut[i] = cut[i - 1] + 1 + (int)(urand() * (N_RB - cut[i - 1] - (n - i)) / 2);
    }
    int order[10];
    for (int i = 0; i < n; i++) {
      order[i] = i;
    }
    for (int i = n - 1; i > 0; i--) {
      const int r = (int)(urand() * (i + 1));
      const int x = order[i];
      order[i] = order[r], order[r] = x;
    }
    nrscope_sensing_join_begin(j, 100, 3);
    for (int i = 0; i < n; i++) {
      const int u = order[i];
      memset(gV, 0, sizeof(gV));
      for (int l = 2; l <= 11; l += 9) {
        for (int k = cut[u] * 12; k < cut[u + 1] * 12; k += 2) {
          gH[l * N_SC + k] = refH[l * N_SC + k];
          gV[l * N_SC + k] = true;
        }
      }
      add(j, 1, 1.0f, (uint16_t)(0x4600 + u), mask);
    }
    out.n = 0;
    nrscope_sensing_join_flush_to(j, sink, NULL);
    exact += out.n == 1 && same_row(0, 2, refH, refV) && same_row(0, 11, refH, refV);
  }
  check(exact == 200, "band split over 2..8 UEs, any order: rejoined bit exact (of 200)", exact, 200);

  // 3. the same from 4 threads at once
  exact = 0;
  for (int t = 0; t < 50; t++) {
    const int n = 8;
    int       cut[9]   = {0, 20, 41, 90, 133, 170, 200, 251, N_RB};
    int       order[8] = {3, 7, 0, 5, 2, 6, 1, 4};
    nrscope_sensing_join_begin(j, 100, 3);
    pthread_t th[4];
    thr_arg_t a[4];
    for (int i = 0; i < 4; i++) {
      a[i] = (thr_arg_t){j, &c, n, i, cut, order};
      pthread_create(&th[i], NULL, thr_add, &a[i]);
    }
    for (int i = 0; i < 4; i++) {
      pthread_join(th[i], NULL);
    }
    out.n = 0;
    nrscope_sensing_join_flush_to(j, sink, NULL);
    exact += out.n == 1 && same_row(0, 2, refH, refV) && same_row(0, 11, refH, refV);
  }
  check(exact == 50, "8 UEs added from 4 threads: rejoined bit exact (of 50)", exact, 50);

  // 4. overlap: the first grant keeps the pilots both carry
  nrscope_sensing_join_begin(j, 100, 3);
  grant_rows(&c, 0, 150, 1.0f, 200, mask, 0);
  add(j, 1, 1.0f, 0x4601, mask);
  grant_rows(&c, 100, N_RB, 2.0f, 200, mask, 0); // twice the value: visible if it won
  add(j, 1, 1.0f, 0x4602, mask);
  out.n = 0;
  nrscope_sensing_join_flush_to(j, sink, NULL);
  {
    const int k = 120 * 12;
    check(out.n == 1 && cabsf(out.H[0][2 * N_SC + k] - chan_at(&c, k)) < 1e-6f, "overlap: the first grant is kept",
          out.n, 1);
  }

  // 5. half the carrier after joining: 4 x 30 PRB = 120 dropped, 5 x 30 = 150 kept
  for (int nu = 4; nu <= 5; nu++) {
    nrscope_sensing_join_begin(j, 100, 3);
    for (int u = 0; u < nu; u++) {
      grant_rows(&c, 40 + 30 * u, 70 + 30 * u, 1.0f, 30, mask, 0);
      add(j, 1, 1.0f, (uint16_t)(0x4600 + u), mask);
    }
    out.n = 0;
    nrscope_sensing_join_flush_to(j, sink, NULL);
    int lo, hi;
    const int np = out.n ? row_pilots(0, 2, &lo, &hi) : 0;
    if (nu == 4) {
      check(out.n == 0, "4 UEs x 30 PRB (120 < 136.5): dropped", out.n, 0);
    } else {
      check(out.n == 1 && np == 150 * 6, "5 UEs x 30 PRB: one row of 150 PRB", np / 6.0, 150);
    }
  }

  // 6. keys: other ports, other DM-RS power never join
  nrscope_sensing_join_begin(j, 100, 3);
  grant_rows(&c, 0, 140, 1.0f, 30, mask, 0);
  add(j, 1, 1.0f, 0x4601, mask);
  grant_rows(&c, 140, N_RB, 1.0f, 30, mask, 0);
  add(j, 2, 1.0f, 0x4602, mask); // port 1001
  out.n = 0;
  nrscope_sensing_join_flush_to(j, sink, NULL);
  check(out.n == 1 && out.ports[0] == 1, "port 1000 + port 1001 grants: not joined (only the 140 PRB one kept)",
        out.n, 1);
  nrscope_sensing_join_begin(j, 100, 3);
  grant_rows(&c, 0, 100, 1.0f, 30, mask, 0);
  add(j, 1, 1.0f, 0x4601, mask);
  grant_rows(&c, 100, N_RB, 1.0f, 30, mask, 0);
  add(j, 1, 1.41421356f, 0x4602, mask); // other CDM groups without data
  out.n = 0;
  nrscope_sensing_join_flush_to(j, sink, NULL);
  check(out.n == 1, "beta 1 and beta sqrt(2) grants: not joined (the 173 PRB one alone)", out.n, 1);

  // 7. precoder step at the junction: 90 deg cuts, 10 deg and 1 dB joins
  for (int s = 0; s < 2; s++) {
    nrscope_sensing_join_begin(j, 100, 3);
    grant_rows(&c, 0, 100, 1.0f, 30, mask, 0);
    add(j, 1, 1.0f, 0x4601, mask);
    const cf_t u = s == 0 ? cexpf(I * (float)(M_PI / 2)) : 1.12f * cexpf(I * (float)(10 * M_PI / 180));
    grant_rows(&c, 100, N_RB, u, 30, mask, 0);
    add(j, 1, 1.0f, 0x4602, mask);
    out.n = 0;
    nrscope_sensing_join_flush_to(j, sink, NULL);
    int lo, hi;
    const int np = out.n ? row_pilots(0, 2, &lo, &hi) : 0;
    if (s == 0) {
      check(out.n == 1 && np == 173 * 6 && lo == 100 * 12, "90 deg precoder step: cut, the wider side (173 PRB) kept",
            np / 6.0, 173);
    } else {
      check(out.n == 1 && np == N_RB * 6, "10 deg / 1 dB step: joined (273 PRB)", np / 6.0, N_RB);
    }
  }

  // 8. a 2 PRB gap: cut at max_gap 0, filled at max_gap 2 (fill error against the channel,
  // one the lines can bridge: 0.15 us of delay spread)
  c = chan_make(3, 0.15e-6);
  for (int g = 0; g <= 2; g += 2) {
    nrscope_sensing_join_args.max_gap_prb = g;
    nrscope_sensing_join_begin(j, 100, 3);
    grant_rows(&c, 0, 120, 1.0f, 200, mask, 0);
    add(j, 1, 1.0f, 0x4601, mask);
    grant_rows(&c, 122, N_RB, 1.0f, 200, mask, 0);
    add(j, 1, 1.0f, 0x4602, mask);
    out.n = 0;
    nrscope_sensing_join_flush_to(j, sink, NULL);
    int lo, hi;
    const int np = out.n ? row_pilots(0, 2, &lo, &hi) : 0;
    if (g == 0) {
      check(out.n == 1 && np == 151 * 6, "2 PRB gap, join_max_gap_prb 0: cut, 151 PRB side kept", np / 6.0, 151);
    } else {
      double e = 0, p = 0;
      for (int k = 120 * 12; k < 122 * 12; k += 2) {
        const cf_t d = out.H[0][2 * N_SC + k] - chan_at(&c, k);
        e += crealf(d * conjf(d));
        p += crealf(chan_at(&c, k) * conjf(chan_at(&c, k)));
      }
      const double nmse = 10 * log10(e / p);
      printf("  (%d PRB in the row)\n", np / 6);
      check(out.n == 1 && np == N_RB * 6 && nmse < -10, "2 PRB gap, join_max_gap_prb 2: filled, NMSE of the fill (dB)",
            nmse, -10);
    }
  }
  nrscope_sensing_join_args.max_gap_prb = 0;
  nrscope_sensing_join_free(j);
}

/* The continuity check over many channels: false cuts (one channel, the check says
   no) and false joins (one random precoder per UE, the check says yes). The row's slope
   comes from 200 subcarriers around the junction, as a grant gives it. */
static void test_check_stats(void)
{
  printf("\nContinuity check at a junction (1 PRB each side, chi2 99%% + %.2f |H|)\n",
         nrscope_sensing_join_args.max_step);
  printf("  %-9s %-9s %-10s %-12s %-24s\n", "SNR/pilot", "spread", "false cut", "false join", "informative (H0 / H1)");
  const double snrs[]   = {0, 5, 10, 20, 30};
  const double spread[] = {0.1e-6, 0.3e-6, 1.0e-6, 3.0e-6};
  double       worst_fc = 0, fj_inf_20 = 0;
  for (int si = 0; si < 5; si++) {
    for (int di = 0; di < 4; di++) {
      int         fc = 0, fj = 0, inf0 = 0, inf1 = 0, fj_inf = 0, n = 2000;
      static cf_t H[N_SC];
      static bool V[N_SC];
      for (int t = 0; t < n; t++) {
        chan_t      c  = chan_make(6, spread[di]);
        const int   kb = (int)(urand() * (N_RB - 40) + 20) * 12;
        const float sn = (float)pow(10.0, -snrs[si] / 20.0);
        cf_t        nz[400];
        for (int i = 0; i < 400; i++) {
          nz[i] = sn * gauss() * (float)M_SQRT1_2;
        }
        for (int pass = 0; pass < 2; pass++) {
          const cf_t u = pass == 0 ? 1.0f : gauss() * (float)M_SQRT1_2;
          memset(V, 0, sizeof(V));
          cf_t acc = 0;
          for (int k = kb - 200, i = 0; k < kb + 200; k += 2, i++) {
            H[k] = (k >= kb ? u : 1.0f) * chan_at(&c, k) + nz[i];
            V[k] = true;
            if (k > kb - 200 && k != kb) {
              acc += H[k] * conjf(H[k - 2]);
            }
          }
          float      db, deg;
          bool       inf;
          const int  r  = nrscope_sensing_join_test(H, V, N_SC, kb - 2, kb, 2, 6, cargf(acc), &db, &deg, &inf);
          if (pass == 0) {
            fc += r != 1;
            inf0 += inf;
          } else {
            fj += r == 1;
            inf1 += inf;
            fj_inf += inf && r == 1;
          }
        }
      }
      printf("  %5.0f dB  %5.1f us  %6.2f%%    %6.2f%%      %5.1f%% / %5.1f%% (false join among them %.1f%%)\n",
             snrs[si], spread[di] * 1e6, 100.0 * fc / n, 100.0 * fj / n, 100.0 * inf0 / n, 100.0 * inf1 / n,
             inf1 ? 100.0 * fj_inf / inf1 : 0.0);
      if (spread[di] <= 1.0e-6 && (double)fc / n > worst_fc) {
        worst_fc = (double)fc / n;
      }
      if (si == 3 && di == 1) {
        fj_inf_20 = inf1 ? (double)fj_inf / inf1 : 1;
      }
    }
  }
  check(worst_fc < 0.03, "worst false-cut rate, delay spread up to 1 us", worst_fc, 0.03);
  check(fj_inf_20 < 0.10, "false joins among informative junctions, 20 dB, 0.3 us", fj_inf_20, 0.10);
}

/* The cell's verdict: two UEs per slot side by side at 20 dB, each with its own random
   precoder (a beamforming cell), then with none (an unprecoded one). */
static void test_learn(void)
{
  printf("\nThe cell's verdict (join_learn)\n");
  nrscope_sensing_join_t* j    = nrscope_sensing_join_alloc(0, N_SC);
  const uint16_t          mask = 1u << 2;
  for (int mode = 0; mode < 2; mode++) {
    nrscope_sensing_join_reset();
    int joined_late = 0, late = 0;
    for (int sl = 0; sl < 600; sl++) {
      chan_t c = chan_make(4, 0.3e-6);
      nrscope_sensing_join_begin(j, 100, 3);
      const int cut = 100 + (int)(urand() * 73);
      grant_rows(&c, 0, cut, 1.0f, 20, mask, 0);
      add(j, 1, 1.0f, 0x4601, mask);
      grant_rows(&c, cut, N_RB, mode == 0 ? gauss() * (float)M_SQRT1_2 : 1.0f, 20, mask, 0);
      add(j, 1, 1.0f, 0x4602, mask);
      out.n = 0;
      nrscope_sensing_join_flush_to(j, sink, NULL);
      if (sl >= 300) {
        int lo, hi;
        late++;
        joined_late += out.n == 1 && row_pilots(0, 2, &lo, &hi) == N_RB * 6;
      }
    }
    if (mode == 0) {
      check(nrscope_sensing_join_cell() == 2, "random precoder per UE: cell judged precoded", nrscope_sensing_join_cell(),
            2);
      check(joined_late == 0, "  and no slot joined after the verdict (of 300)", joined_late, 0);
    } else {
      check(nrscope_sensing_join_cell() == 1, "no precoder: cell judged coherent", nrscope_sensing_join_cell(), 1);
      check(joined_late >= 294, "  and the UEs joined (of 300)", joined_late, 294);
    }
  }
  nrscope_sensing_join_reset();
  nrscope_sensing_join_free(j);
}

/* ---- end to end, one scene per process ---- */

static void scene(const char* mode, const char* dir)
{
  const bool split = strcmp(mode, "single") != 0;
  const bool narrow = strcmp(mode, "narrow") == 0 || strcmp(mode, "narrow_old") == 0;
  const bool use_join = strcmp(mode, "split") == 0 || strcmp(mode, "narrow") == 0;
  const uint32_t n_sc = N_SC;
  nrscope_sensing_default_args(&nrscope_sensing_args);
  nrscope_sensing_args.enable         = true;
  nrscope_sensing_args.symbols        = NR_SENSING_HISTORY_DEPTH;
  nrscope_sensing_args.max_speed_ms   = 10.0;
  nrscope_sensing_args.clutter_kernel = true;
  nrscope_sensing_args.tdd_detect     = true;
  snprintf(nrscope_sensing_args.dump, sizeof(nrscope_sensing_args.dump), "%s/map2d.csv", dir);
  nr_sensing_tdd_period_slots = 10;
  mkdir(dir, 0755);
  nrscope_sensing_t*         s  = nrscope_sensing_get(3450000000ULL, 122.88e6, 4096, 30000, 1);
  nrscope_sensing_scratch_t* sc = nrscope_sensing_scratch_alloc(n_sc);
  nrscope_sensing_join_t*    j  = nrscope_sensing_join_alloc(0, n_sc);

  srsran_carrier_nr_t carrier = SRSRAN_DEFAULT_CARRIER_NR;
  carrier.pci = 1, carrier.nof_prb = N_RB, carrier.scs = srsran_subcarrier_spacing_30kHz;
  srsran_dmrs_sch_t dmrs = {};
  srsran_dmrs_sch_init(&dmrs, false);
  srsran_dmrs_sch_set_carrier(&dmrs, &carrier);
  static cf_t tx[N_SC * NS], tx1[N_SC * NS], rx[N_SC * NS];
  const nr_dmrs_placement_t place = {0, 0, 0};
  const double              f_d = 120.0;
  const int                 d_bins = 23, t_bins = 41;
  int                       pushed = 0;
  for (int sl = 0; sl < 1800; sl++) {
    const int phase = sl % 10;
    if (phase >= 8) {
      continue;
    }
    srsran_sch_cfg_nr_t cfg = {};
    cfg.dmrs.type = srsran_dmrs_sch_type_1, cfg.dmrs.typeA_pos = srsran_dmrs_sch_typeA_pos_2;
    cfg.dmrs.additional_pos = srsran_dmrs_sch_add_pos_2, cfg.dmrs.length = srsran_dmrs_sch_len_1;
    cfg.grant.mapping = srsran_sch_mapping_type_A, cfg.grant.S = 1, cfg.grant.L = (phase == 7) ? 5 : 13;
    cfg.grant.nof_layers = 1, cfg.grant.nof_dmrs_cdm_groups_without_data = 1, cfg.grant.beta_dmrs = 1.0f;
    const uint32_t    sfn = 100 + sl / 20, slot_idx = sl % 20;
    srsran_slot_cfg_t slot = {.idx = slot_idx};
    /* the UEs of this slot: cut points drawn per slot, contiguous; "narrow" leaves the
      lowest 30% of the carrier empty and no UE more than 60 PRB */
    int n_ue = 1, cut[12] = {0, N_RB};
    if (split) {
      const int lo = narrow ? 82 : 0;
      n_ue         = narrow ? 4 : 2 + (sl % 5);
      cut[0] = lo, cut[n_ue] = N_RB;
      for (int i = 1; i < n_ue; i++) {
        const int room = N_RB - cut[i - 1] - (n_ue - i) * 8;
        cut[i]         = cut[i - 1] + 8 + (narrow ? (N_RB - lo) / n_ue - 8 : (int)(urand() * room / 2));
      }
    }
    memset(tx, 0, sizeof(tx));
    srsran_sch_cfg_nr_t ue_cfg[12];
    for (int u = 0; u < n_ue; u++) {
      ue_cfg[u]               = cfg;
      ue_cfg[u].grant.rnti    = 0x4601 + u;
      ue_cfg[u].grant.nof_prb = cut[u + 1] - cut[u];
      for (int p = 0; p < N_RB; p++) {
        ue_cfg[u].grant.prb_idx[p] = p >= cut[u] && p < cut[u + 1];
      }
      // each UE's DM-RS written on its own, then summed: the sequence is the carrier's
      memset(tx1, 0, sizeof(tx1));
      srsran_dmrs_sch_put_sf(&dmrs, &slot, &ue_cfg[u], &ue_cfg[u].grant, tx1);
      for (int i = 0; i < N_SC * NS; i++) {
        tx[i] += tx1[i];
      }
    }
    for (int l = 0; l < NS; l++) {
      const uint32_t t_in = (l == 0) ? 352 : 352 + 4096 + (l - 1) * 4384 + 288;
      const double   t    = ((double)(100 * 20 + sl) * 61440 + t_in) / 122.88e6;
      const cf_t     dop  = 0.3f * cexpf(I * (float)(2.0 * M_PI * f_d * t));
      for (int k = 0; k < N_SC; k++) {
        const cf_t h = cexpf(-I * (float)(2.0 * M_PI * d_bins * k / 4096.0)) +
                       dop * cexpf(-I * (float)(2.0 * M_PI * t_bins * k / 4096.0));
        rx[l * N_SC + k] = h * tx[l * N_SC + k];
      }
    }
    if (use_join) {
      nrscope_sensing_join_begin(j, sfn, slot_idx);
    }
    for (int u = n_ue - 1; u >= 0; u--) { // the decoders' order is not the frequency order
      const int n = nrscope_sensing_process_grant_join(s, sc, use_join ? j : NULL, 0, rx, n_sc, &place, &ue_cfg[u], 0,
                                                       1, sfn, slot_idx);
      pushed += (!use_join && n > 0) ? n : 0;
    }
    if (use_join) {
      nrscope_sensing_join_flush(j, NULL);
    }
  }
  nrscope_sensing_wait_maps(s);
  nrscope_sensing_join_report();
  // maps peaking at the target
  char path[600];
  snprintf(path, sizeof(path), "%s/map2d.csv", dir);
  FILE* f      = fopen(path, "r");
  int   n_maps = 0, ok = 0;
  static char line[4 << 20];
  if (f != NULL) {
    fgets(line, sizeof(line), f);
    while (fgets(line, sizeof(line), f)) {
      double v[21];
      char*  p = line;
      for (int i = 0; i < 21; i++) {
        v[i] = strtod(p, &p);
        p++;
      }
      const int    n_bins = (int)v[12], n_freq = (int)v[13];
      const double f_max  = v[10];
      int          best_i = 0;
      double       best   = -1;
      for (int i = 0; i < n_bins * n_freq; i++) {
        const double pw = strtod(p, &p);
        p++;
        if (pw > best) {
          best = pw, best_i = i;
        }
      }
      const int    b    = best_i / n_freq;
      const double f_hz = -f_max + (best_i % n_freq) * 2.0 * f_max / n_freq;
      ok += abs(b - t_bins) <= 1 && fabs(f_hz - f_d) <= 4.0 * f_max / n_freq;
      n_maps++;
    }
    fclose(f);
  }
  printf("scene %s: %d maps, %d peaking at the target\n", mode, n_maps, ok);
  FILE* r = fopen(strcat(strcpy(path, dir), "/result"), "w");
  if (r) {
    fprintf(r, "%d %d\n", n_maps, ok);
    fclose(r);
  }
  srsran_dmrs_sch_free(&dmrs);
}

static bool files_equal(const char* a, const char* b)
{
  FILE* fa = fopen(a, "rb");
  FILE* fb = fopen(b, "rb");
  bool  eq = fa && fb;
  while (eq) {
    const int ca = fgetc(fa), cb = fgetc(fb);
    eq = ca == cb;
    if (ca == EOF || cb == EOF) {
      break;
    }
  }
  if (fa) fclose(fa);
  if (fb) fclose(fb);
  return eq;
}

static void run_scene(const char* self, const char* mode, int* maps, int* ok)
{
  char cmd[1024], path[600];
  snprintf(cmd, sizeof(cmd), "rm -rf join_test_%s && %s --scene %s join_test_%s > join_test_%s.log 2>&1", mode, self,
           mode, mode, mode);
  if (system(cmd) != 0) {
    printf("  scene %s failed (join_test_%s.log)\n", mode, mode);
  }
  *maps = *ok = -1;
  snprintf(path, sizeof(path), "join_test_%s/result", mode);
  FILE* r = fopen(path, "r");
  if (r) {
    if (fscanf(r, "%d %d", maps, ok) != 2) {
      *maps = *ok = -1;
    }
    fclose(r);
  }
}

int main(int argc, char** argv)
{
  if (argc == 4 && strcmp(argv[1], "--scene") == 0) {
    scene(argv[2], argv[3]);
    return 0;
  }
  test_rows();
  test_check_stats();
  test_learn();

  printf("\nEnd to end, through the sensing pipeline (each scene in its own process)\n");
  char self[512];
  const ssize_t n = readlink("/proc/self/exe", self, sizeof(self) - 1);
  self[n > 0 ? n : 0] = 0;
  int m_single, ok_single, m_split, ok_split, m_nold, ok_nold, m_nar, ok_nar;
  run_scene(self, "single", &m_single, &ok_single);
  run_scene(self, "split", &m_split, &ok_split);
  run_scene(self, "narrow_old", &m_nold, &ok_nold);
  run_scene(self, "narrow", &m_nar, &ok_nar);
  check(m_single > 0 && ok_single == m_single, "one full-band UE per slot: maps peaking at the target", ok_single,
        m_single);
  check(files_equal("join_test_single/map2d.csv", "join_test_split/map2d.csv") && m_split == m_single,
        "2..6 UEs per slot, joined: maps byte for byte those of one UE", m_split, m_single);
  {
    // every row of the split scene joined back to the full carrier (join.csv)
    FILE* f    = fopen("join_test_split/map2d.csv.join.csv", "r");
    int   rows = 0, full = 0;
    char  line[4096];
    if (f && fgets(line, sizeof(line), f)) {
      while (fgets(line, sizeof(line), f)) {
        double v[11];
        char*  p = line;
        for (int i = 0; i < 11; i++) {
          v[i] = strtod(p, &p);
          p++;
        }
        rows++;
        full += v[10] == 1 && fabs(v[9] - N_RB) < 0.01 && v[5] >= 2;
      }
    }
    if (f) {
      fclose(f);
    }
    check(rows > 0 && full == rows, "  join.csv: every row 273 PRB from 2..6 grants", full, rows);
  }
  check(m_nold <= 0, "4 UEs x 48 PRB per slot, not joined: no map (every grant under half)", m_nold, 0);
  check(m_nar > 0 && ok_nar == m_nar, "4 UEs x 48 PRB per slot, joined: maps peaking at the target", ok_nar, m_nar);

  printf("\n%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
