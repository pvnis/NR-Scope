/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * \brief Target detection on a range-Doppler map with TDD induced impulsive
 *        sidelobes. See nr_ue_tdd_detect.h for the method and the references.
 *
 * NOT WIRED IN.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/sensing_defs.h"
#include "nrscope/hdr/sensing/nr_ue_tdd_detect.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------------- */
/* small helpers                                                             */
/* ------------------------------------------------------------------------- */

/* Doppler of grid index j. The grid spans [-f_max, +f_max] inclusive, matching
   nr_ue_sensing_range_doppler(). Note that this puts zero Doppler exactly on a grid
   point only when n_freq is odd; with an even n_freq the DC of a static return
   straddles two cells. */
static inline double nr_tdd_freq_of(const nr_tdd_obs_t *obs, int j)
{
  if (obs->n_freq <= 1)
    return 0.0;
  const double df = 2.0 * obs->f_max_hz / (obs->n_freq - 1);
  return -obs->f_max_hz + j * df;
}

static inline double nr_tdd_df(const nr_tdd_obs_t *obs)
{
  return (obs->n_freq > 1) ? 2.0 * obs->f_max_hz / (obs->n_freq - 1) : 0.0;
}

/* Grid index, fractional, of a Doppler in Hz. Inverse of nr_tdd_freq_of(). */
static inline double nr_tdd_index_of(const nr_tdd_obs_t *obs, double f_hz)
{
  const double df = nr_tdd_df(obs);
  return (df > 0.0) ? (f_hz + obs->f_max_hz) / df : 0.0;
}

/* cos/sin of -2*pi*turns, with the phase reduced to one turn first. The products
   frequency*time and lattice*offset reach hundreds of turns, and evaluating the trig
   that far out spends most of the mantissa on the integer part of the phase. */
static inline void nr_tdd_rot(double turns, double *re, double *im)
{
  const double frac = turns - trunc(turns);
  const double ph = -2.0 * M_PI * frac;
  *re = cos(ph);
  *im = sin(ph);
}

/* ------------------------------------------------------------------------- */
/* configuration                                                             */
/* ------------------------------------------------------------------------- */

void nr_tdd_cfg_default(nr_tdd_cfg_t *cfg)
{
  memset(cfg, 0, sizeof(*cfg));

  /* Per cell false alarm rate. The paper uses 1e-6; with a 62 x 384 map that is
  about 0.02 false cells expected per map before any peak checking. */
  cfg->pfa = 1e-6;

  /* Guard cells have to cover the mainlobe or the target leaks into its own noise
  estimate and censors itself. Two bins covers a Hann tapered mainlobe on a full
  band grant; the Doppler mainlobe is narrower relative to the grid, hence one. */
  cfg->cfar_guard_bin = 2;
  cfg->cfar_guard_freq = 1;
  cfg->cfar_train_bin = 4;
  cfg->cfar_train_freq = 8;

  /* Eq. (23). The paper uses 0 in simulation, where removal is near perfect, and
  0.5 on measurements to survive residual clutter. This pipeline refines delay only
  by interpolation, so removal is imperfect by construction and a demanding gamma
  would reject real targets: start permissive. */
  cfg->gamma = 0.2;
  cfg->n_sidelobes = 1;
  cfg->ell_bin = 2.0;
  cfg->ell_freq = 3.0;

  cfg->focus_zoom = 16;
  cfg->focus_span = 2;

  cfg->max_targets = 8;
  cfg->max_iter = 64;
}

/* ------------------------------------------------------------------------- */
/* complex map                                                               */
/* ------------------------------------------------------------------------- */

bool nr_tdd_cmap_alloc(nr_tdd_cmap_t *map, int n_bins, int n_freq)
{
  map->n_bins = n_bins;
  map->n_freq = n_freq;
  map->re = calloc((size_t)n_bins * n_freq, sizeof(*map->re));
  map->im = calloc((size_t)n_bins * n_freq, sizeof(*map->im));
  if (map->re == NULL || map->im == NULL) {
    free(map->re);
    free(map->im);
    map->re = NULL;
    map->im = NULL;
    return false;
  }
  return true;
}

void nr_tdd_cmap_free(nr_tdd_cmap_t *map)
{
  free(map->re);
  free(map->im);
  map->re = NULL;
  map->im = NULL;
  map->n_bins = 0;
  map->n_freq = 0;
}

void nr_tdd_cmap_copy(const nr_tdd_cmap_t *src, nr_tdd_cmap_t *dst)
{
  const size_t n = (size_t)src->n_bins * src->n_freq;
  dst->n_bins = src->n_bins;
  dst->n_freq = src->n_freq;
  memcpy(dst->re, src->re, n * sizeof(*src->re));
  memcpy(dst->im, src->im, n * sizeof(*src->im));
}

double nr_tdd_cmap_power(const nr_tdd_cmap_t *map, int bin, int f)
{
  const int i = bin * map->n_freq + f;
  return (double)map->re[i] * map->re[i] + (double)map->im[i] * map->im[i];
}

/* ------------------------------------------------------------------------- */
/* point spread function                                                     */
/* ------------------------------------------------------------------------- */

double nr_tdd_dirichlet(double x, int A)
{
  const double den = sin(M_PI * x);
  /* At integer x both terms vanish. L'Hopital gives cos(A*pi*x)/cos(pi*x), which is
  +-1 depending on the parity of A and of the integer. */
  if (fabs(den) < 1e-12) {
    const long k = lround(x);
    return (((A - 1) * k) % 2 == 0) ? 1.0 : -1.0;
  }
  return sin(A * M_PI * x) / (A * den);
}

double nr_tdd_psf_doppler_closed(double m_over_Mp, int M_DL, int R, int M_TDD)
{
  /* Eq. (10), normalised to 1 at the origin: the paper's M_DL * R prefactor is the
  value at m = 0, so dividing it out makes this directly comparable with
  nr_tdd_psf_doppler(). */
  return nr_tdd_dirichlet(m_over_Mp, M_DL) * nr_tdd_dirichlet(M_TDD * m_over_Mp, R);
}

void nr_tdd_psf_doppler(const nr_tdd_obs_t *obs, double x_hz, double *re, double *im)
{
  double ar = 0.0;
  double ai = 0.0;
  for (int i = 0; i < obs->n_snap; i++) {
    double c;
    double s;
    nr_tdd_rot(x_hz * obs->t_s[i], &c, &s);
    ar += c;
    ai += s;
  }
  const double g = (obs->n_snap > 0) ? 1.0 / obs->n_snap : 0.0;
  *re = ar * g;
  *im = ai * g;
}

/* The slow-time samples a unit target at f0_hz leaves after the slow-trend filter:
   s_i = exp(+j*2*pi*f0*t_i) minus its projection on the removed subspace,
       d = s - Q (Q^T s).
   Without a filter d = s. d_re, d_im receive n_snap entries. */
static void nr_tdd_filtered_tone(const nr_tdd_obs_t *obs, double f0_hz, double *d_re, double *d_im)
{
  const int n = obs->n_snap;
  for (int i = 0; i < n; i++) {
    double c;
    double s;
    nr_tdd_rot(f0_hz * obs->t_s[i], &c, &s);
    d_re[i] = c;
    d_im[i] = -s; // nr_tdd_rot() gives exp(-j...), the tone is exp(+j...)
  }
  if (obs->trend_q == NULL)
    return;
  for (int k = 0; k < obs->trend_n_q; k++) {
    const double *q = &obs->trend_q[(size_t)k * n];
    double cr = 0.0;
    double ci = 0.0;
    for (int i = 0; i < n; i++) {
      cr += q[i] * d_re[i];
      ci += q[i] * d_im[i];
    }
    for (int i = 0; i < n; i++) {
      d_re[i] -= cr * q[i];
      d_im[i] -= ci * q[i];
    }
  }
}

/* Doppler response at f_hz of the filtered tone d of nr_tdd_filtered_tone(), with the
   normalisation of nr_tdd_psf_doppler(): (1/n) sum_i d_i exp(-j*2*pi*f*t_i). Without a
   filter this is exactly W_D(f - f0); with one it is W_D(f - f0) minus the response of
   what the filter removed, and no longer 1 at f0. */
static void nr_tdd_psf_doppler_filtered(const nr_tdd_obs_t *obs,
                                        const double *d_re,
                                        const double *d_im,
                                        double f_hz,
                                        double *re,
                                        double *im)
{
  double ar = 0.0;
  double ai = 0.0;
  for (int i = 0; i < obs->n_snap; i++) {
    double c;
    double s;
    nr_tdd_rot(f_hz * obs->t_s[i], &c, &s);
    ar += d_re[i] * c - d_im[i] * s;
    ai += d_re[i] * s + d_im[i] * c;
  }
  const double g = (obs->n_snap > 0) ? 1.0 / obs->n_snap : 0.0;
  *re = ar * g;
  *im = ai * g;
}

void nr_tdd_psf_range(const nr_tdd_obs_t *obs, double u_bins, double *re, double *im)
{
  const double inv_n = 1.0 / obs->idft_size;
  double ar = 0.0;
  double ai = 0.0;
  double wsum = 0.0;

  /* sum_j win[j] * exp(+j*2*pi*j*u/N). nr_tdd_rot() carries a minus sign, so the
  imaginary part is negated back here to keep the IDFT convention of
  nr_ue_sensing_delay_response(). */
  for (int j = 0; j < obs->win_len; j++) {
    const double w = (obs->win != NULL) ? obs->win[j] : 1.0;
    double c;
    double s;
    nr_tdd_rot(j * u_bins * inv_n, &c, &s);
    ar += w * c;
    ai -= w * s;
    wsum += w;
  }

  /* the exp(j*2*pi*a*u/N) of the grant start; it cancels at u = 0, which is why the
  peak of a resolved path does not move when the scheduler moves the allocation */
  double ca;
  double sa;
  nr_tdd_rot(obs->win_start * u_bins * inv_n, &ca, &sa);
  sa = -sa;

  const double g = (wsum > 0.0) ? 1.0 / wsum : 0.0;
  *re = (ar * ca - ai * sa) * g;
  *im = (ar * sa + ai * ca) * g;
}

/* ------------------------------------------------------------------------- */
/* periodogram                                                               */
/* ------------------------------------------------------------------------- */

void nr_tdd_periodogram(const nr_tdd_obs_t *obs, nr_tdd_cmap_t *out)
{
  double *wr = malloc((size_t)obs->n_snap * sizeof(*wr));
  double *wi = malloc((size_t)obs->n_snap * sizeof(*wi));
  if (wr == NULL || wi == NULL) {
    free(wr);
    free(wi);
    ERROR("tdd detect: cannot allocate the transform kernel\n");
    return;
  }

  for (int f = 0; f < obs->n_freq; f++) {
    const double freq = nr_tdd_freq_of(obs, f);
    for (int i = 0; i < obs->n_snap; i++)
      nr_tdd_rot(freq * obs->t_s[i], &wr[i], &wi[i]);

    for (int b = 0; b < obs->n_bins; b++) {
      double ar = 0.0;
      double ai = 0.0;
      for (int i = 0; i < obs->n_snap; i++) {
        const double hr = obs->h_re[(size_t)i * obs->n_bins + b];
        const double hi = obs->h_im[(size_t)i * obs->n_bins + b];
        ar += hr * wr[i] - hi * wi[i];
        ai += hr * wi[i] + hi * wr[i];
      }
      out->re[b * out->n_freq + f] = (float)ar;
      out->im[b * out->n_freq + f] = (float)ai;
    }
  }
  free(wr);
  free(wi);
}

/* One cell of the transform at an arbitrary Doppler, used by the focused analysis
   where the grid is not fine enough. */
static void nr_tdd_cell_at(const nr_tdd_obs_t *obs, int bin, double f_hz, double *re, double *im)
{
  double ar = 0.0;
  double ai = 0.0;
  for (int i = 0; i < obs->n_snap; i++) {
    double c;
    double s;
    nr_tdd_rot(f_hz * obs->t_s[i], &c, &s);
    const double hr = obs->h_re[(size_t)i * obs->n_bins + bin];
    const double hi = obs->h_im[(size_t)i * obs->n_bins + bin];
    ar += hr * c - hi * s;
    ai += hr * s + hi * c;
  }
  *re = ar;
  *im = ai;
}

/* ------------------------------------------------------------------------- */
/* CA-CFAR                                                                   */
/* ------------------------------------------------------------------------- */

int nr_tdd_cfar(const nr_tdd_cmap_t *map,
                const nr_tdd_cfg_t *cfg,
                int max_out,
                int bin_out[max_out],
                int freq_out[max_out],
                double snr_dB_out[max_out])
{
  /* The CFAR only needs the power |C|^2. We compute it once here, then run the power
  version. The training ring reads every cell many times, so this is also faster than
  squaring on every read. */
  const size_t n_cells = (size_t)map->n_bins * map->n_freq;
  float *power = malloc(n_cells * sizeof(*power));
  if (power == NULL) {
    ERROR("tdd detect: cannot allocate the CFAR power map\n");
    return 0;
  }
  for (size_t i = 0; i < n_cells; i++)
    power[i] = (float)((double)map->re[i] * map->re[i] + (double)map->im[i] * map->im[i]);

  const int n_found = nr_tdd_cfar_power(power, map->n_bins, map->n_freq, cfg, max_out, bin_out, freq_out, snr_dB_out);
  free(power);
  return n_found;
}

int nr_tdd_cfar_power(const float *power,
                      int n_bins,
                      int n_freq,
                      const nr_tdd_cfg_t *cfg,
                      int max_out,
                      int bin_out[max_out],
                      int freq_out[max_out],
                      double snr_dB_out[max_out])
{
  const int gb = cfg->cfar_guard_bin;
  const int gf = cfg->cfar_guard_freq;
  const int hb = gb + cfg->cfar_train_bin;
  const int hf = gf + cfg->cfar_train_freq;
  int n_found = 0;

  for (int b = 0; b < n_bins; b++) {
    for (int f = 0; f < n_freq; f++) {
      const double p = power[b * n_freq + f];
      if (p <= 0.0)
        continue;

      /* Local maximum first: without it a single strong return is reported once per
      cell of its mainlobe and the CLEAN loop spends its whole budget on one target. */
      bool is_max = true;
      for (int db = -1; db <= 1 && is_max; db++) {
        for (int df = -1; df <= 1; df++) {
          if (db == 0 && df == 0)
            continue;
          const int bb = b + db;
          const int ff = f + df;
          if (bb < 0 || bb >= n_bins || ff < 0 || ff >= n_freq)
            continue;
          if (power[bb * n_freq + ff] > p) {
            is_max = false;
            break;
          }
        }
      }
      if (!is_max)
        continue;

      /* Training ring, guard band excluded. The ring is clipped at the edges rather
      than wrapped: the Doppler axis is a finite span and the delay axis is
      truncated, so neither is periodic. */
      double sum = 0.0;
      int nt = 0;
      for (int db = -hb; db <= hb; db++) {
        for (int df = -hf; df <= hf; df++) {
          if (abs(db) <= gb && abs(df) <= gf)
            continue;
          const int bb = b + db;
          const int ff = f + df;
          if (bb < 0 || bb >= n_bins || ff < 0 || ff >= n_freq)
            continue;
          sum += power[bb * n_freq + ff];
          nt++;
        }
      }
      if (nt < 8)
        continue; // too few training cells to estimate anything

      const double noise = sum / nt;
      if (noise <= 0.0)
        continue;

      // cell averaging threshold factor for the requested per cell false alarm rate
      const double alpha = nt * (pow(cfg->pfa, -1.0 / nt) - 1.0);
      if (p <= alpha * noise)
        continue;

      const double snr_dB = 10.0 * log10(p / noise);

      // insert strongest first
      int pos = 0;
      while (pos < n_found && snr_dB_out[pos] >= snr_dB)
        pos++;
      if (pos >= max_out)
        continue;
      for (int i = (n_found < max_out ? n_found : max_out - 1); i > pos; i--) {
        bin_out[i] = bin_out[i - 1];
        freq_out[i] = freq_out[i - 1];
        snr_dB_out[i] = snr_dB_out[i - 1];
      }
      bin_out[pos] = b;
      freq_out[pos] = f;
      snr_dB_out[pos] = snr_dB;
      if (n_found < max_out)
        n_found++;
    }
  }
  return n_found;
}

/* ------------------------------------------------------------------------- */
/* focused Fourier analysis                                                  */
/* ------------------------------------------------------------------------- */

void nr_tdd_focus(const nr_tdd_obs_t *obs, const nr_tdd_cfg_t *cfg, int bin0, int f0, nr_tdd_target_t *target)
{
  memset(target, 0, sizeof(*target));

  /* ---- Doppler, exact ----
  The transform can be evaluated at any frequency, so the peak is found on a grid
  focus_zoom times finer instead of being interpolated. */
  const double df = nr_tdd_df(obs);
  const double f_centre = nr_tdd_freq_of(obs, f0);
  const int zoom = (cfg->focus_zoom > 0) ? cfg->focus_zoom : 1;
  const int span = (cfg->focus_span > 0) ? cfg->focus_span : 1;
  const double step = df / zoom;

  double best_f = f_centre;
  double best_p = -1.0;
  for (int j = -span * zoom; j <= span * zoom; j++) {
    const double f = f_centre + j * step;
    if (f < -obs->f_max_hz || f > obs->f_max_hz)
      continue;
    double re;
    double im;
    nr_tdd_cell_at(obs, bin0, f, &re, &im);
    const double p = re * re + im * im;
    if (p > best_p) {
      best_p = p;
      best_f = f;
    }
  }

  /* ---- delay, parabolic ----
  Exact refinement would need the frequency domain lattice, which is not kept once
  the delay response has been formed. A parabola through the log power of the three
  bins around the peak is the standard substitute and lands within a few hundredths
  of a bin on a smooth mainlobe. The residual is what gamma in Eq. (23) absorbs. */
  double d_hat = bin0;
  if (bin0 >= 1 && bin0 + 1 < obs->n_bins) {
    double re;
    double im;
    double y[3];
    for (int k = 0; k < 3; k++) {
      nr_tdd_cell_at(obs, bin0 - 1 + k, best_f, &re, &im);
      const double p = re * re + im * im;
      y[k] = log(p > 0.0 ? p : 1e-30);
    }
    const double den = y[0] - 2.0 * y[1] + y[2];
    if (fabs(den) > 1e-12) {
      double delta = 0.5 * (y[0] - y[2]) / den;
      if (delta > 0.5)
        delta = 0.5;
      if (delta < -0.5)
        delta = -0.5;
      d_hat = bin0 + delta;
    }
  }

  /* ---- amplitude ----
  The model is C(q,f) = amp * W_R(q - d) * W_D'(f; f0), with W_D' the response of the
  tone after the slow-trend filter (nr_tdd_psf_doppler_filtered()). W_R is 1 at the
  origin; W_D' is 1 at f0 only without a filter, so it is evaluated there rather than
  assumed. amp then follows from the observed cell alone and the subtraction lands on
  zero there. */
  double cre;
  double cim;
  nr_tdd_cell_at(obs, bin0, best_f, &cre, &cim);

  double wre;
  double wim;
  nr_tdd_psf_range(obs, (double)bin0 - d_hat, &wre, &wim);

  double *d_re = malloc((size_t)obs->n_snap * sizeof(*d_re));
  double *d_im = malloc((size_t)obs->n_snap * sizeof(*d_im));
  if (d_re != NULL && d_im != NULL) {
    double dre;
    double dim;
    nr_tdd_filtered_tone(obs, best_f, d_re, d_im);
    nr_tdd_psf_doppler_filtered(obs, d_re, d_im, best_f, &dre, &dim);
    // W_R * W_D'(f0)
    const double pr = wre * dre - wim * dim;
    const double pi = wre * dim + wim * dre;
    wre = pr;
    wim = pi;
  } else {
    ERROR("tdd detect: cannot allocate the filtered tone, amplitude left unfiltered\n");
  }
  free(d_re);
  free(d_im);
  const double wmag2 = wre * wre + wim * wim;
  if (wmag2 > 1e-12) {
    target->amp_re = (cre * wre + cim * wim) / wmag2;
    target->amp_im = (cim * wre - cre * wim) / wmag2;
  } else {
    target->amp_re = cre;
    target->amp_im = cim;
  }

  target->bin = d_hat;
  target->range_m = d_hat * obs->m_per_bin;
  target->f_hz = best_f;
  target->speed_ms = best_f * obs->lambda_m / 2.0;
}

/* ------------------------------------------------------------------------- */
/* coherent removal                                                          */
/* ------------------------------------------------------------------------- */

void nr_tdd_psf_subtract(const nr_tdd_obs_t *obs, const nr_tdd_target_t *target, nr_tdd_cmap_t *map)
{
  double *rr = malloc((size_t)map->n_bins * sizeof(*rr));
  double *ri = malloc((size_t)map->n_bins * sizeof(*ri));
  double *dr = malloc((size_t)map->n_freq * sizeof(*dr));
  double *di = malloc((size_t)map->n_freq * sizeof(*di));
  double *tone_re = malloc((size_t)obs->n_snap * sizeof(*tone_re));
  double *tone_im = malloc((size_t)obs->n_snap * sizeof(*tone_im));
  if (rr == NULL || ri == NULL || dr == NULL || di == NULL || tone_re == NULL || tone_im == NULL) {
    free(rr);
    free(ri);
    free(dr);
    free(di);
    free(tone_re);
    free(tone_im);
    ERROR("tdd detect: cannot allocate the PSF axes\n");
    return;
  }

  /* Separable: the response is the outer product of a range PSF and a Doppler PSF,
  so it costs n_bins + n_freq evaluations instead of n_bins * n_freq. The replicas
  are not synthesised separately, they are already inside the Doppler PSF because
  the sampling pattern that produces them is what it is built from.

  The Doppler PSF is the one of the target's tone after the slow-trend filter, the
  same filter the samples went through, so only what the data still holds is removed.
  Without a filter it is W_D(f - f0) exactly. */
  nr_tdd_filtered_tone(obs, target->f_hz, tone_re, tone_im);
  for (int b = 0; b < map->n_bins; b++)
    nr_tdd_psf_range(obs, (double)b - target->bin, &rr[b], &ri[b]);
  for (int f = 0; f < map->n_freq; f++)
    nr_tdd_psf_doppler_filtered(obs, tone_re, tone_im, nr_tdd_freq_of(obs, f), &dr[f], &di[f]);
  free(tone_re);
  free(tone_im);

  for (int b = 0; b < map->n_bins; b++) {
    for (int f = 0; f < map->n_freq; f++) {
      // amp * W_R * W_D
      const double xr = rr[b] * dr[f] - ri[b] * di[f];
      const double xi = rr[b] * di[f] + ri[b] * dr[f];
      const double vr = target->amp_re * xr - target->amp_im * xi;
      const double vi = target->amp_re * xi + target->amp_im * xr;
      const int i = b * map->n_freq + f;
      map->re[i] -= (float)vr;
      map->im[i] -= (float)vi;
    }
  }

  free(rr);
  free(ri);
  free(dr);
  free(di);
}

/* ------------------------------------------------------------------------- */
/* power features check                                                      */
/* ------------------------------------------------------------------------- */

/* Mean power over an ellipse centred on a fractional cell. */
static bool nr_tdd_ellipse_mean(const nr_tdd_cmap_t *map,
                                double bin_c,
                                double freq_c,
                                double ell_bin,
                                double ell_freq,
                                double *mean)
{
  const int b0 = (int)floor(bin_c - ell_bin);
  const int b1 = (int)ceil(bin_c + ell_bin);
  const int f0 = (int)floor(freq_c - ell_freq);
  const int f1 = (int)ceil(freq_c + ell_freq);
  double sum = 0.0;
  int n = 0;

  for (int b = b0; b <= b1; b++) {
    if (b < 0 || b >= map->n_bins)
      continue;
    for (int f = f0; f <= f1; f++) {
      if (f < 0 || f >= map->n_freq)
        continue;
      const double u = (b - bin_c) / ell_bin;
      const double v = (f - freq_c) / ell_freq;
      if (u * u + v * v > 1.0)
        continue;
      sum += nr_tdd_cmap_power(map, b, f);
      n++;
    }
  }
  if (n == 0)
    return false;
  *mean = sum / n;
  return true;
}

bool nr_tdd_sidelobe_check(const nr_tdd_cmap_t *before,
                           const nr_tdd_cmap_t *after,
                           const nr_tdd_obs_t *obs,
                           const nr_tdd_cfg_t *cfg,
                           const nr_tdd_target_t *target)
{
  if (obs->t_tdd_s <= 0.0)
    return true; // no TDD structure declared, nothing to test against

  const double df_replica = 1.0 / obs->t_tdd_s;
  int n_tested = 0;

  for (int k = 1; k <= cfg->n_sidelobes; k++) {
    for (int sign = -1; sign <= 1; sign += 2) {
      const double f_k = target->f_hz + sign * k * df_replica;
      if (f_k < -obs->f_max_hz || f_k > obs->f_max_hz)
        continue; // replica outside the grid, cannot be looked at

      const double j_k = nr_tdd_index_of(obs, f_k);
      double p_before;
      double p_after;
      if (!nr_tdd_ellipse_mean(before, target->bin, j_k, cfg->ell_bin, cfg->ell_freq, &p_before))
        continue;
      if (!nr_tdd_ellipse_mean(after, target->bin, j_k, cfg->ell_bin, cfg->ell_freq, &p_after))
        continue;

      n_tested++;
      /* Eq. (23). Removing a true target takes its own replicas down with it.
      Removing a replica leaves the target that made it, and its other replicas,
      exactly where they were. */
      if (p_after > (1.0 - cfg->gamma) * p_before)
        return false;
    }
  }

  /* Nothing could be tested: the replicas all fell off the grid. Rejecting is the
  conservative answer, and it is also a configuration error worth seeing, since a
  grid narrower than 1/T_TDD makes this whole method inoperative. */
  if (n_tested == 0) {
    WARNING(
          "tdd detect: no replica of the peak at %.1f Hz lies inside the +-%.1f Hz grid, "
          "widen f_max to at least %.1f Hz\n",
          target->f_hz,
          obs->f_max_hz,
          df_replica);
    return false;
  }
  return true;
}

/* ------------------------------------------------------------------------- */
/* Algorithm 1                                                               */
/* ------------------------------------------------------------------------- */

/* Record a candidate the loop threw out, at its grid position. Counted even when the
   array is full, so a caller can tell that markers are missing rather than absent. */
static void nr_tdd_note_rejected(const nr_tdd_obs_t *obs,
                                 nr_tdd_target_t *rejected_out,
                                 int *n_rejected,
                                 int bin,
                                 int f,
                                 double snr_dB,
                                 int iteration,
                                 nr_tdd_verdict_t verdict)
{
  if (n_rejected == NULL)
    return;
  const int slot = (*n_rejected)++;
  if (rejected_out == NULL || slot >= NR_TDD_MAX_REJECTED)
    return;
  nr_tdd_target_t *r = &rejected_out[slot];
  memset(r, 0, sizeof(*r));
  r->bin = bin;
  r->range_m = bin * obs->m_per_bin;
  r->f_hz = nr_tdd_freq_of(obs, f);
  r->speed_ms = r->f_hz * obs->lambda_m / 2.0;
  r->snr_dB = snr_dB;
  r->iteration = iteration;
  r->verdict = verdict;
}

int nr_tdd_detect(const nr_tdd_obs_t *obs,
                  const nr_tdd_cfg_t *cfg,
                  nr_tdd_target_t *out,
                  int max_out,
                  nr_tdd_target_t *rejected_out,
                  int *n_rejected)
{
  if (n_rejected != NULL) *n_rejected = 0;
  if (obs->n_snap < 8 || obs->n_bins < 3 || obs->n_freq < 3) return 0;

  nr_tdd_cmap_t work = {0};
  nr_tdd_cmap_t saved = {0};
  
  if (!nr_tdd_cmap_alloc(&work, obs->n_bins, obs->n_freq))
    return 0;
  if (!nr_tdd_cmap_alloc(&saved, obs->n_bins, obs->n_freq)) {
    nr_tdd_cmap_free(&work);
    return 0;
  }

  nr_tdd_periodogram(obs, &work);

  int cand_bin[NR_TDD_MAX_CANDIDATES];
  int cand_freq[NR_TDD_MAX_CANDIDATES];
  double cand_snr[NR_TDD_MAX_CANDIDATES];
  int n_cand = nr_tdd_cfar(&work, cfg, NR_TDD_MAX_CANDIDATES, cand_bin, cand_freq, cand_snr);

  int n_out = 0;
  int p = 0; // index of the candidate under test, as in Algorithm 1
  int iter = 0;

  while (p < n_cand && n_out < max_out && n_out < cfg->max_targets && iter < cfg->max_iter) {
    iter++;

    /* Skip a candidate sitting on a target already accepted: it is a residue of an
    imperfect removal, not a new object. */
    bool duplicate = false;
    for (int i = 0; i < n_out; i++) {
      const double d_bin = cand_bin[p] - out[i].bin;
      const double d_f = nr_tdd_freq_of(obs, cand_freq[p]) - out[i].f_hz;
      if (fabs(d_bin) <= cfg->ell_bin && fabs(d_f) <= cfg->ell_freq * nr_tdd_df(obs)) {
        duplicate = true;
        break;
      }
    }
    if (duplicate) {
      nr_tdd_note_rejected(obs, rejected_out, n_rejected, cand_bin[p], cand_freq[p], cand_snr[p], iter,
                           NR_TDD_VERDICT_DUPLICATE);
      p++;
      continue;
    }

    nr_tdd_target_t cand;
    nr_tdd_focus(obs, cfg, cand_bin[p], cand_freq[p], &cand);
    cand.snr_dB = cand_snr[p];
    cand.iteration = iter;

    nr_tdd_cmap_copy(&work, &saved);
    nr_tdd_psf_subtract(obs, &cand, &work);

    cand.verdict = NR_TDD_VERDICT_TARGET;
    if (nr_tdd_sidelobe_check(&saved, &work, obs, cfg, &cand)) {
      /* Accepted. The removal stays, so weaker targets that were buried under this
      one, or under its replicas, can now be found. The candidate list is rebuilt
      from the updated map and the scan restarts from the strongest. */
      out[n_out++] = cand;
      n_cand = nr_tdd_cfar(&work, cfg, NR_TDD_MAX_CANDIDATES, cand_bin, cand_freq, cand_snr);
      p = 0;
    } else {
      // rejected as a replica or an artefact: put the map back and move on
      nr_tdd_note_rejected(obs, rejected_out, n_rejected, cand_bin[p], cand_freq[p], cand_snr[p], iter,
                           NR_TDD_VERDICT_REPLICA);
      nr_tdd_cmap_copy(&saved, &work);
      p++;
    }
  }

  nr_tdd_cmap_free(&work);
  nr_tdd_cmap_free(&saved);
  return n_out;
}
