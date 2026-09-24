/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_music.h"

/* A steering vector the slow-trend filter has taken more than this fraction of is not
scored: what is left of it is too short for the ratio in P(f) to mean anything, and the
filter has removed whatever the data had there anyway. It only happens right around 0 Hz. */
#define NR_MUSIC_MIN_STEER_FRAC 0.05

/* Relative convergence of the Jacobi sweeps, on the off-diagonal norm */
#define NR_MUSIC_JACOBI_TOL 1e-12
#define NR_MUSIC_JACOBI_MAX_SWEEPS 50

void nr_ue_music_eig_hermitian(int L, double *a_re, double *a_im, double *lambda, double *v_re, double *v_im)
{
  for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++) {
      v_re[i * L + j] = (i == j) ? 1.0 : 0.0;
      v_im[i * L + j] = 0.0;
    }

  double fro2 = 0.0;
  for (int i = 0; i < L * L; i++)
    fro2 += a_re[i] * a_re[i] + a_im[i] * a_im[i];
  const double tol2 = NR_MUSIC_JACOBI_TOL * NR_MUSIC_JACOBI_TOL * fro2;

  for (int sweep = 0; sweep < NR_MUSIC_JACOBI_MAX_SWEEPS; sweep++) {
    double off2 = 0.0;
    for (int p = 0; p < L; p++)
      for (int q = p + 1; q < L; q++)
        off2 += a_re[p * L + q] * a_re[p * L + q] + a_im[p * L + q] * a_im[p * L + q];
    if (off2 <= tol2)
      break;

    for (int p = 0; p < L; p++) {
      for (int q = p + 1; q < L; q++) {
        const double mag = hypot(a_re[p * L + q], a_im[p * L + q]);
        if (mag * mag <= tol2 / (L * L))
          continue;

        /* One rotation W = D G. D = diag(1, e^{-j phi}) at (p, q) makes a_pq real, with
        phi its phase; G is then the real Jacobi rotation of Numerical Recipes, which
        zeroes it. A <- W^H A W, V <- V W. */
        const double er = a_re[p * L + q] / mag; // e^{-j phi} = conj(a_pq) / |a_pq|
        const double ei = -a_im[p * L + q] / mag;
        const double theta = (a_re[q * L + q] - a_re[p * L + p]) / (2.0 * mag);
        const double t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
        const double c = 1.0 / sqrt(1.0 + t * t);
        const double s = t * c;

        // columns of A and V: x_p <- c x_p - s e^{-j phi} x_q, x_q <- s x_p + c e^{-j phi} x_q
        for (int k = 0; k < L; k++) {
          double *mr[2] = {a_re, v_re};
          double *mi[2] = {a_im, v_im};
          for (int m = 0; m < 2; m++) {
            const double pr = mr[m][k * L + p], pi = mi[m][k * L + p];
            const double qr = mr[m][k * L + q], qi = mi[m][k * L + q];
            const double eqr = er * qr - ei * qi;
            const double eqi = er * qi + ei * qr;
            mr[m][k * L + p] = c * pr - s * eqr;
            mi[m][k * L + p] = c * pi - s * eqi;
            mr[m][k * L + q] = s * pr + c * eqr;
            mi[m][k * L + q] = s * pi + c * eqi;
          }
        }
        // rows of A: y_p <- c y_p - s e^{j phi} y_q, y_q <- s y_p + c e^{j phi} y_q
        for (int k = 0; k < L; k++) {
          const double pr = a_re[p * L + k], pi = a_im[p * L + k];
          const double qr = a_re[q * L + k], qi = a_im[q * L + k];
          const double eqr = er * qr + ei * qi;
          const double eqi = er * qi - ei * qr;
          a_re[p * L + k] = c * pr - s * eqr;
          a_im[p * L + k] = c * pi - s * eqi;
          a_re[q * L + k] = s * pr + c * eqr;
          a_im[q * L + k] = s * pi + c * eqi;
        }
        // exact by construction, so rounding does not leave them to the next sweep
        a_re[p * L + q] = a_im[p * L + q] = 0.0;
        a_re[q * L + p] = a_im[q * L + p] = 0.0;
        a_im[p * L + p] = a_im[q * L + q] = 0.0;
      }
    }
  }

  for (int i = 0; i < L; i++)
    lambda[i] = a_re[i * L + i];

  // decreasing order, columns of V permuted along
  for (int i = 0; i < L - 1; i++) {
    int best = i;
    for (int j = i + 1; j < L; j++)
      if (lambda[j] > lambda[best])
        best = j;
    if (best == i)
      continue;
    const double tmp = lambda[i];
    lambda[i] = lambda[best];
    lambda[best] = tmp;
    for (int k = 0; k < L; k++) {
      const double r = v_re[k * L + i], im = v_im[k * L + i];
      v_re[k * L + i] = v_re[k * L + best];
      v_im[k * L + i] = v_im[k * L + best];
      v_re[k * L + best] = r;
      v_im[k * L + best] = im;
    }
  }
}

/// One column of X: observation o at range bin b
typedef struct {
  const float *re;
  const float *im;
  int stride;
  int bin;
} nr_music_col_t;

/* The columns of range bin b and their Gram matrix G = X^H X. Returns L. */
static int nr_music_gram(const nr_sensing_slowtime_t *obs,
                         const int *use,
                         int n_use,
                         int n,
                         int n_bins,
                         int b,
                         nr_music_col_t *col,
                         double *g_re,
                         double *g_im)
{
  int L = 0;
  for (int u = 0; u < n_use; u++) {
    const nr_sensing_slowtime_t *s = &obs[use[u]];
    for (int d = -NR_MUSIC_RANGE_HALF; d <= NR_MUSIC_RANGE_HALF; d++) {
      const int bb = b + d;
      if (bb < 0 || bb >= n_bins)
        continue;
      col[L++] = (nr_music_col_t){.re = s->h_re, .im = s->h_im, .stride = s->n_bins, .bin = bb};
    }
  }

  for (int k = 0; k < L; k++) {
    for (int l = k; l < L; l++) {
      double sr = 0.0;
      double si = 0.0;
      for (int i = 0; i < n; i++) {
        const double kr = col[k].re[(size_t)i * col[k].stride + col[k].bin];
        const double ki = col[k].im[(size_t)i * col[k].stride + col[k].bin];
        const double lr = col[l].re[(size_t)i * col[l].stride + col[l].bin];
        const double li = col[l].im[(size_t)i * col[l].stride + col[l].bin];
        // conj(x_k) x_l
        sr += kr * lr + ki * li;
        si += kr * li - ki * lr;
      }
      g_re[k * L + l] = sr;
      g_im[k * L + l] = si;
      g_re[l * L + k] = sr;
      g_im[l * L + k] = -si;
    }
  }
  return L;
}

static int nr_music_cmp_double(const void *a, const void *b)
{
  const double x = *(const double *)a;
  const double y = *(const double *)b;
  return (x > y) - (x < y);
}

int nr_ue_music_doppler(const nr_sensing_slowtime_t *obs, int n_obs, nr_sensing_map_t *map, nr_music_stats_t *stats)
{
  if (n_obs > NR_MUSIC_MAX_OBS)
    n_obs = NR_MUSIC_MAX_OBS;

  int ref = -1;
  for (int o = 0; o < n_obs && ref < 0; o++)
    if (obs[o].n_snap >= 8)
      ref = o;
  if (ref < 0)
    return 0;

  const int n = obs[ref].n_snap;
  const double *t = obs[ref].t_s;

  /* Observations share one steering vector only if they share the sample instants. The
  antennas always do, a layer the scheduler dropped in part of the window does not. */
  int use[NR_MUSIC_MAX_OBS];
  int n_use = 0;
  int n_bins = map->n_bins;
  for (int o = 0; o < n_obs; o++) {
    const nr_sensing_slowtime_t *s = &obs[o];
    if (s->n_snap != n)
      continue;
    bool same = true;
    for (int i = 0; i < n && same; i++)
      same = fabs(s->t_s[i] - t[i]) < 1e-7;
    if (!same) {
      DEBUG("sensing music: observation %d sampled at other instants, left out\n", o);
      continue;
    }
    use[n_use++] = o;
    if (s->n_bins < n_bins)
      n_bins = s->n_bins;
  }
  if (n_use == 0 || n_bins < 1)
    return 0;

  const int n_freq = map->n_freq;
  const double f_max = map->f_max_hz;
  const double df = (n_freq > 1) ? (2.0 * f_max) / (n_freq - 1) : 0.0;

  // the subspace the clutter stage removed, see nr_ue_music.h
  double q[NR_CLUTTER_SLOW_TREND_MAX + 1][NR_SENSING_HISTORY_DEPTH];
  const int n_q = nr_ue_sensing_slow_basis(t, n, obs[ref].trend_degree, q);

  /* Steering vectors a(f)[i] = exp(j*2*pi*f*t_i), the signal the DFT correlates against,
  and |a_p(f)|^2 = n - |Q^T a|^2 of their projection. The vectors themselves are kept
  unprojected: the u_k are built from projected data, so u_k^H a_p = u_k^H a. */
  float *st_re = malloc_or_fail((size_t)n_freq * n * sizeof(float));
  float *st_im = malloc_or_fail((size_t)n_freq * n * sizeof(float));
  double *norm2 = malloc_or_fail((size_t)n_freq * sizeof(double));
  for (int f = 0; f < n_freq; f++) {
    const double freq = -f_max + f * df;
    float *ar = &st_re[(size_t)f * n];
    float *ai = &st_im[(size_t)f * n];
    for (int i = 0; i < n; i++) {
      // reduced to one turn first, as in the DFT, to keep the trig argument small
      const double turns = freq * t[i];
      const double ph = 2.0 * M_PI * (turns - trunc(turns));
      ar[i] = (float)cos(ph);
      ai[i] = (float)sin(ph);
    }
    double lost = 0.0;
    for (int k = 0; k < n_q; k++) {
      double cr = 0.0;
      double ci = 0.0;
      for (int i = 0; i < n; i++) {
        cr += q[k][i] * ar[i];
        ci += q[k][i] * ai[i];
      }
      lost += cr * cr + ci * ci;
    }
    norm2[f] = (double)n - lost;
  }

  const int L_max = n_use * (2 * NR_MUSIC_RANGE_HALF + 1);
  nr_music_col_t *col = malloc_or_fail((size_t)L_max * sizeof(*col));
  double *g_re = malloc_or_fail((size_t)L_max * L_max * sizeof(double));
  double *g_im = malloc_or_fail((size_t)L_max * L_max * sizeof(double));
  double *v_re = malloc_or_fail((size_t)L_max * L_max * sizeof(double));
  double *v_im = malloc_or_fail((size_t)L_max * L_max * sizeof(double));
  double *lambda = malloc_or_fail((size_t)L_max * sizeof(double));
  double *u_re = malloc_or_fail((size_t)NR_MUSIC_MAX_ORDER * n * sizeof(double));
  double *u_im = malloc_or_fail((size_t)NR_MUSIC_MAX_ORDER * n * sizeof(double));

  /* Pass 1: the noise reference, as the median over bins of the largest eigenvalue. The
  few bins with a target sit above it and do not move it. */
  double lam_max[NR_SENSING_MAP_MAX_BINS_RANGE];
  for (int b = 0; b < n_bins; b++) {
    const int L = nr_music_gram(obs, use, n_use, n, n_bins, b, col, g_re, g_im);
    nr_ue_music_eig_hermitian(L, g_re, g_im, lambda, v_re, v_im);
    lam_max[b] = lambda[0];
  }
  double sorted[NR_SENSING_MAP_MAX_BINS_RANGE];
  memcpy(sorted, lam_max, (size_t)n_bins * sizeof(double));
  qsort(sorted, n_bins, sizeof(double), nr_music_cmp_double);
  const double noise_ref = sorted[n_bins / 2];
  const double thr = noise_ref * pow(10.0, NR_MUSIC_EIG_THRESH_DB / 10.0);

  // Pass 2: model order and pseudo-spectrum, bin by bin
  int bins_with_signal = 0;
  int order_max = 0;
  for (int b = 0; b < n_bins; b++) {
    const int L = nr_music_gram(obs, use, n_use, n, n_bins, b, col, g_re, g_im);
    nr_ue_music_eig_hermitian(L, g_re, g_im, lambda, v_re, v_im);

    /* L columns span at most L signal dimensions. The noise subspace lives in the n
    dimensional slow-time space, not among the columns, so K = L is still a valid model. */
    const int K_cap = L < NR_MUSIC_MAX_ORDER ? L : NR_MUSIC_MAX_ORDER;
    int K = 0;
    if (NR_MUSIC_FIXED_ORDER > 0)
      K = NR_MUSIC_FIXED_ORDER < K_cap ? NR_MUSIC_FIXED_ORDER : K_cap;
    else
      while (K < K_cap && lambda[K] > thr && lambda[K] > 0.0)
        K++;

    float *row = &map->power[b * n_freq];
    if (K == 0) {
      for (int f = 0; f < n_freq; f++)
        row[f] = 1.0f;
      continue;
    }
    bins_with_signal++;
    if (K > order_max)
      order_max = K;

    // u_k = X v_k / sqrt(lambda_k)
    for (int k = 0; k < K; k++) {
      const double inv = 1.0 / sqrt(lambda[k]);
      double *ur = &u_re[(size_t)k * n];
      double *ui = &u_im[(size_t)k * n];
      for (int i = 0; i < n; i++) {
        double sr = 0.0;
        double si = 0.0;
        for (int l = 0; l < L; l++) {
          const double xr = col[l].re[(size_t)i * col[l].stride + col[l].bin];
          const double xi = col[l].im[(size_t)i * col[l].stride + col[l].bin];
          const double vr = v_re[l * L + k];
          const double vi = v_im[l * L + k];
          sr += xr * vr - xi * vi;
          si += xr * vi + xi * vr;
        }
        ur[i] = sr * inv;
        ui[i] = si * inv;
      }
    }

    for (int f = 0; f < n_freq; f++) {
      if (norm2[f] < NR_MUSIC_MIN_STEER_FRAC * n) {
        row[f] = 1.0f;
        continue;
      }
      const float *ar = &st_re[(size_t)f * n];
      const float *ai = &st_im[(size_t)f * n];
      double proj = 0.0;
      for (int k = 0; k < K; k++) {
        const double *ur = &u_re[(size_t)k * n];
        const double *ui = &u_im[(size_t)k * n];
        // u^H a
        double cr = 0.0;
        double ci = 0.0;
        for (int i = 0; i < n; i++) {
          cr += ur[i] * ar[i] + ui[i] * ai[i];
          ci += ur[i] * ai[i] - ui[i] * ar[i];
        }
        proj += cr * cr + ci * ci;
      }
      double den = norm2[f] - proj;
      const double den_min = norm2[f] / NR_MUSIC_PSEUDO_MAX;
      if (den < den_min)
        den = den_min;
      row[f] = (float)(norm2[f] / den);
    }
  }

  // bins the observations do not all reach are left at the floor
  for (int b = n_bins; b < map->n_bins; b++)
    for (int f = 0; f < n_freq; f++)
      map->power[b * n_freq + f] = 1.0f;

  if (stats != NULL) {
    stats->n_obs = n_use;
    stats->n_cols = L_max;
    stats->noise_ref = noise_ref;
    stats->bins_with_signal = bins_with_signal;
    stats->order_max = order_max;
  }

  free(st_re);
  free(st_im);
  free(norm2);
  free(col);
  free(g_re);
  free(g_im);
  free(v_re);
  free(v_im);
  free(lambda);
  free(u_re);
  free(u_im);
  return n_use;
}
