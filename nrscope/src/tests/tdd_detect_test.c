/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * Check of nr_ue_tdd_detect.c on a synthetic TDD sampling pattern, ported from
 * the OAI reference (tag isac-reference-for-nrscope-port, test/tdd_detect_test.c)
 * and built as the tdd_detect_test target.
 *
 * One target is placed at a fractional delay bin and a known Doppler, on a
 * sampling pattern of 3 DMRS symbols in each of 4 downlink slots per 5 slot TDD
 * period. It prints, in order: the replica level the pattern produces, what a plain
 * CA-CFAR reports on that map, and what nr_tdd_detect() reports.
 *
 * What it demonstrates: the strongest peak on the map is not the target but its
 * +1 replica, so a detector without replica checking picks the wrong one.
 *
 * Pass criterion (added in the port; the original only printed): the first
 * target nr_tdd_detect() reports is the real one, within half a delay bin and one
 * Doppler cell of the truth. Picking a replica would be off by 1/T_TDD = 400 Hz.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrscope/hdr/sensing/nr_ue_tdd_detect.h"

#define SCS      30000.0
#define SLOT_S   (1e-3 / 2.0)          /* mu = 1 */
#define TDD_P    5                      /* slots per TDD period */
#define DL_SLOTS 4                      /* DL slots per period */
#define NPER     32
#define NSNAP    (NPER * DL_SLOTS * 3)  /* 3 DMRS symbols per DL slot */
#define NBINS    62
#define NFREQ    257                    /* odd -> DC on the grid */
#define IDFT     2048
#define WLEN     1638

int main(void)
{
  static double t[NSNAP];
  static float hr[NSNAP * NBINS], hi[NSNAP * NBINS];
  static float win[WLEN];

  /* Hann taper, same as nr_ue_sensing_apply_hann() */
  for (int j = 0; j < WLEN; j++)
    win[j] = (float)(0.5 * (1.0 - cos(2.0 * M_PI * j / (WLEN - 1))));

  /* TDD sampling pattern: DL slots 0..3 of every 5, 3 symbols each */
  int n = 0;
  for (int q = 0; q < NPER; q++)
    for (int s = 0; s < DL_SLOTS; s++)
      for (int k = 0; k < 3; k++) {
        const int sym[3] = {2, 7, 11};
        t[n++] = (q * TDD_P + s) * SLOT_S + sym[k] * (1.0 / (14.0 * 2.0 * 1000.0));
      }

  nr_tdd_obs_t obs = {
    .n_snap = n, .t_s = t, .h_re = hr, .h_im = hi,
    .n_bins = NBINS, .idft_size = IDFT, .win_len = WLEN, .win_start = 0, .win = win,
    .m_per_bin = 2.44, .n_freq = NFREQ,
    .f_max_hz = 600.0,                       /* wide enough to hold the +-400 Hz replicas */
    .lambda_m = 299792458.0 / 3.41e9,
    .t_tdd_s = TDD_P * SLOT_S,
  };

  /* one target: fractional delay bin, 50 Hz Doppler */
  const double D_TRUE = 20.3, F_TRUE = 50.0, AMP = 3000.0;
  static double wr[NBINS], wi[NBINS];
  for (int b = 0; b < NBINS; b++)
    nr_tdd_psf_range(&obs, (double)b - D_TRUE, &wr[b], &wi[b]);

  srand(1);
  for (int i = 0; i < n; i++) {
    const double ph = 2.0 * M_PI * F_TRUE * t[i];
    const double cr = cos(ph), ci = sin(ph);
    for (int b = 0; b < NBINS; b++) {
      const double vr = AMP * (wr[b] * cr - wi[b] * ci);
      const double vi = AMP * (wr[b] * ci + wi[b] * cr);
      const double nr_ = 6.0 * ((double)rand() / RAND_MAX - 0.5);
      const double ni_ = 6.0 * ((double)rand() / RAND_MAX - 0.5);
      hr[i * NBINS + b] = (float)(vr + nr_);
      hi[i * NBINS + b] = (float)(vi + ni_);
    }
  }

  printf("pattern: %d samples over %.1f ms, T_TDD %.2f ms -> replicas every %.0f Hz (%.1f m/s)\n",
         n, (t[n-1] - t[0]) * 1e3, obs.t_tdd_s * 1e3, 1.0 / obs.t_tdd_s,
         (1.0 / obs.t_tdd_s) * obs.lambda_m / 2.0);
  printf("truth  : bin %.2f (%.1f m), %.1f Hz (%.2f m/s)\n\n",
         D_TRUE, D_TRUE * obs.m_per_bin, F_TRUE, F_TRUE * obs.lambda_m / 2.0);

  /* how strong are the replicas, before any detection? */
  double w0r, w0i, w1r, w1i;
  nr_tdd_psf_doppler(&obs, 0.0, &w0r, &w0i);
  nr_tdd_psf_doppler(&obs, 1.0 / obs.t_tdd_s, &w1r, &w1i);
  printf("Doppler PSF: |W(0)| = %.4f, |W(1/T_TDD)| = %.4f  -> replica %.1f dB down\n\n",
         hypot(w0r, w0i), hypot(w1r, w1i),
         20.0 * log10(hypot(w1r, w1i) / hypot(w0r, w0i)));

  /* what a naive CFAR sees, with no replica checking */
  nr_tdd_cfg_t cfg;
  nr_tdd_cfg_default(&cfg);
  nr_tdd_cmap_t map = {0};
  nr_tdd_cmap_alloc(&map, NBINS, NFREQ);
  nr_tdd_periodogram(&obs, &map);
  int cb[64], cf[64];
  double cs[64];
  const int nc = nr_tdd_cfar(&map, &cfg, 64, cb, cf, cs);
  printf("CA-CFAR alone: %d candidate peak(s)\n", nc);
  for (int i = 0; i < nc && i < 8; i++) {
    const double f = -obs.f_max_hz + cf[i] * (2.0 * obs.f_max_hz / (NFREQ - 1));
    printf("   bin %2d  %+7.1f Hz (%+6.2f m/s)  %5.1f dB\n",
           cb[i], f, f * obs.lambda_m / 2.0, cs[i]);
  }
  nr_tdd_cmap_free(&map);

  /* with replica checking */
  nr_tdd_target_t tgt[NR_TDD_MAX_TARGETS];
  nr_tdd_target_t rej[NR_TDD_MAX_REJECTED];
  int nrej = 0;
  const int nt = nr_tdd_detect(&obs, &cfg, tgt, NR_TDD_MAX_TARGETS, rej, &nrej);
  printf("\nnr_tdd_detect(): %d target(s), %d rejected\n", nt, nrej);
  for (int i = 0; i < nt; i++)
    printf("   TARGET   bin %6.2f (%6.2f m)  %+7.2f Hz (%+6.3f m/s)  %5.1f dB  iter %d\n",
           tgt[i].bin, tgt[i].range_m, tgt[i].f_hz, tgt[i].speed_ms, tgt[i].snr_dB, tgt[i].iteration);
  for (int i = 0; i < nrej && i < NR_TDD_MAX_REJECTED; i++)
    printf("   %-8s bin %6.2f (%6.2f m)  %+7.2f Hz (%+6.3f m/s)  %5.1f dB  iter %d\n",
           rej[i].verdict == NR_TDD_VERDICT_REPLICA ? "REPLICA" : "DUP",
           rej[i].bin, rej[i].range_m, rej[i].f_hz, rej[i].speed_ms, rej[i].snr_dB, rej[i].iteration);

  const double cell_hz = 2.0 * obs.f_max_hz / (NFREQ - 1);
  if (nt < 1) {
    printf("\nFAILED: no target reported\n");
    return 1;
  }
  const double e_bin = tgt[0].bin - D_TRUE, e_hz = tgt[0].f_hz - F_TRUE;
  printf("\nerror: %+.3f bins, %+.3f Hz\n", e_bin, e_hz);
  if (fabs(e_bin) > 0.5 || fabs(e_hz) > cell_hz) {
    printf("FAILED: first target off by more than 0.5 bin or one Doppler cell (%.2f Hz)\n", cell_hz);
    return 1;
  }
  printf("PASSED: the target, not a replica (tolerance 0.5 bin, %.2f Hz)\n", cell_hz);
  return 0;
}
