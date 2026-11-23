/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/* sensing.backend omp: two moving paths on snapshots with gaps (random runs of measured PRBs per snapshot), built
through the receiver's own taper and delay transform (nr_ue_sensing_delay_response), not through the kernel the
back end models them with, so the test checks the model against independent data. Full band and gappy, with and
without the slow-time mean the clutter stage removes. Each case: both paths found within half a bin and one
Doppler grid step, nothing else above threshold. */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_omp.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"

#define T_KSTEP 2
#define T_L 1638            // lattice positions of 273 PRBs at k_step 2
#define T_B 50              // delay bins kept (NR_SENSING_MAP_MAX_BINS_RANGE)
#define T_NANT 2
#define T_LAMBDA (299792458.0 / 3.75e9)

static int failures = 0;
static void check(int ok, const char *what, double got, double want, double tol)
{
  printf("  %-44s got %9.3f  want %9.3f  +-%.3f  %s\n", what, got, want, tol, ok ? "ok" : "FAIL");
  failures += !ok;
}

static unsigned rng_state = 12345u;
static double urand(void) { rng_state = rng_state * 1664525u + 1013904223u; return (rng_state >> 8) / 16777216.0; }
static double nrand(void) { double u = urand() + 1e-12, v = urand(); return sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v); }

typedef struct { double bin, f_hz, amp, ant_deg; } t_path_t;

static void run_case(const char *name, int gappy, int mean_removed, int gains)
{
  printf("%s\n", name);
  const int N = NR_SENSING_IDFT_SIZE(T_KSTEP);
  enum { NS = 240 };
  static nr_sensing_snapshot_t snap[NS];
  const nr_sensing_snapshot_t *sp[NS];
  static double t_s[NS], gain[NS];
  static float h_re[T_NANT][NS * T_B], h_im[T_NANT][NS * T_B];
  const t_path_t paths[2] = {{16.4, 60.0, 0.05, 15.0}, {33.7, -90.0, 0.03, -20.0}};
  static cf_t pil[NR_SENSING_MAX_IDFT];
  static cf_t h[NR_SENSING_MAX_IDFT];
  float w[NR_SENSING_MAX_IDFT];

  // snapshot times: 240 of the first 300 slots, one symbol each, in order
  int slot = 0;
  for (int i = 0; i < NS; i++) {
    slot += 1 + (urand() < 0.25);
    t_s[i] = slot * 0.5e-3 + 2.0 / 28000.0;
  }
  for (int i = 0; i < NS; i++) {
    // measured positions: full band, or 2-5 runs of 60-500 positions
    static unsigned char m[T_L];
    memset(m, gappy ? 0 : 1, sizeof(m));
    if (gappy) {
      const int nr = 2 + (int)(urand() * 4);
      for (int r = 0; r < nr; r++) {
        const int a = (int)(urand() * (T_L - 520)), len = 60 + (int)(urand() * 440);
        for (int l = a; l < a + len; l++) m[l] = 1;
      }
    }
    int lo = 0, hi = T_L - 1;
    while (!m[lo]) lo++;
    while (!m[hi]) hi--;
    const int span = hi - lo + 1;
    if (!nr_ue_sensing_hann_coeffs(span, w))
      for (int l = 0; l < span; l++) w[l] = 1.0f;
    nr_sensing_snapshot_t *s = &snap[i];
    memset(s, 0, sizeof(*s));
    s->stream.k_step = T_KSTEP;
    s->k_first = lo * T_KSTEP;
    s->n_pilots = span;
    s->idft_size = N;
    for (int l = lo; l <= hi; l++)
      if (m[l] && w[l - lo] != 0.0f) s->mask[l >> 6] |= 1ULL << (l & 63);
    sp[i] = s;
    gain[i] = gains ? 1.0 + 3.0 * urand() : 1.0; // as sensing.gain_equalise scales a snapshot, noise included
    for (int a = 0; a < T_NANT; a++) {
      memset(pil, 0, sizeof(pil));
      for (int l = lo; l <= hi; l++) {
        if (!m[l]) continue;
        double complex x = 0.0;
        for (int p = 0; p < 2; p++)
          x += paths[p].amp * cexp(-2.0 * M_PI * I * l * paths[p].bin / N) * cexp(2.0 * M_PI * I * paths[p].f_hz * (t_s[i] - t_s[0]))
               * cexp(I * M_PI * a * sin(paths[p].ant_deg * M_PI / 180.0));
        x += 1e-4 * (nrand() + I * nrand());
        pil[l] = (cf_t)(x * w[l - lo]);
      }
      nr_ue_sensing_delay_response(T_L, pil, N, h);
      for (int b = 0; b < T_B; b++) {
        h_re[a][i * T_B + b] = (float)(gain[i] * crealf(h[b]));
        h_im[a][i * T_B + b] = (float)(gain[i] * cimagf(h[b]));
      }
    }
  }
  static double q0[NS];
  int n_q = 0;
  if (mean_removed) { // as the clutter stage's degree-0 trend: each bin loses its mean over the window
    for (int a = 0; a < T_NANT; a++)
      for (int b = 0; b < T_B; b++) {
        double mr = 0.0, mi = 0.0;
        for (int i = 0; i < NS; i++) { mr += h_re[a][i * T_B + b]; mi += h_im[a][i * T_B + b]; }
        for (int i = 0; i < NS; i++) { h_re[a][i * T_B + b] -= (float)(mr / NS); h_im[a][i * T_B + b] -= (float)(mi / NS); }
      }
    for (int i = 0; i < NS; i++) q0[i] = 1.0 / sqrt((double)NS);
    n_q = 1;
  }
  nr_sensing_slowtime_t slow[T_NANT];
  for (int a = 0; a < T_NANT; a++) {
    memset(&slow[a], 0, sizeof(slow[a]));
    slow[a].n_snap = NS;
    slow[a].n_bins = T_B;
    slow[a].t_s = t_s;
    slow[a].h_re = h_re[a];
    slow[a].h_im = h_im[a];
    slow[a].idft_size = N;
    slow[a].m_per_bin = 2.44;
    slow[a].snap = sp;
    slow[a].gain = gain;
  }
  const double span_s = t_s[NS - 1] - t_s[0];
  const double f_max = 200.0;
  const int n_freq = 2 * (int)ceil(2.0 * f_max * span_s) + 1;
  const double df = 2.0 * f_max / (n_freq - 1);
  nr_tdd_target_t tg[32];
  nr_omp_cfg_t cfg = {0};
  const int nt = nr_ue_omp_detect(slow, T_NANT, n_freq, f_max, T_LAMBDA, mean_removed ? q0 : NULL, n_q, &cfg, tg, 32);
  check(nt == 2, "paths found", nt, 2, 0);
  if (nt != 2)
    for (int j = 0; j < nt; j++)
      printf("    path %2d: bin %7.3f  %+8.2f Hz  %5.1f dB\n", j, tg[j].bin, tg[j].f_hz, tg[j].snr_dB);
  for (int p = 0; p < 2; p++) {
    int k = -1;
    double dbest = 1e9;
    for (int j = 0; j < nt; j++) {
      const double d = fabs(tg[j].bin - paths[p].bin) + fabs(tg[j].f_hz - paths[p].f_hz) / df;
      if (d < dbest) { dbest = d; k = j; }
    }
    char what[96];
    snprintf(what, sizeof(what), "path %d delay bin", p);
    check(k >= 0 && fabs(tg[k].bin - paths[p].bin) <= 0.5, what, k >= 0 ? tg[k].bin : -1, paths[p].bin, 0.5);
    snprintf(what, sizeof(what), "path %d Doppler (Hz)", p);
    check(k >= 0 && fabs(tg[k].f_hz - paths[p].f_hz) <= df, what, k >= 0 ? tg[k].f_hz : -1, paths[p].f_hz, df);
  }
}

int main(void)
{
  run_case("Full band, two moving paths", 0, 0, 0);
  run_case("Gaps (2-5 runs per snapshot), two moving paths", 1, 0, 0);
  run_case("Gaps, slow-time mean removed (trend projection)", 1, 1, 0);
  run_case("Gaps, mean removed, snapshots scaled by gains 1-4 (gain_equalise)", 1, 1, 1);
  printf("%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
