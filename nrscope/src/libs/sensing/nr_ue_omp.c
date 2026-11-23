/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/* Mask-aware 2D-OMP back end, see nr_ue_omp.h */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_omp.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"
#include "nrscope/hdr/sensing/sensing_defs.h"
#include "nrscope/hdr/sensing/nr_ue_aoa.h"

typedef double complex dc_t;

static int cmp_double(const void *a, const void *b)
{
  const double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

/// lattice index of a snapshot's first pilot, as nr_ue_map.c computes it
static inline int omp_lattice_start(const nr_sensing_snapshot_t *s)
{
  return (s->k_first - s->k_first % s->stream.k_step) / s->stream.k_step;
}

/* k[b] = gain * K_s(b - u0), b = 0 .. nb-1: the snapshot's masked kernel, the exact delay response of a unit path
at fractional bin u0 in that snapshot (taper over its measured pilots only). */
static void omp_kernel(const nr_sensing_snapshot_t *s, double gain, double u0, int nb, float *kr, float *ki, dc_t *k)
{
  nr_ue_sensing_clutter_kernel_masked(s->n_pilots, omp_lattice_start(s), s->mask, u0, s->idft_size, nb, kr, ki);
  for (int b = 0; b < nb; b++)
    k[b] = gain * ((double)kr[b] + I * (double)ki[b]);
}

/* Solve G z = r for n_rhs right-hand sides (G: m x m Hermitian positive definite, row major; r: m x n_rhs, row
major, overwritten by z). Gaussian elimination with partial pivoting. False if G is singular. */
static bool omp_solve(dc_t *G, int m, dc_t *r, int n_rhs)
{
  for (int c = 0; c < m; c++) {
    int p = c;
    for (int i = c + 1; i < m; i++)
      if (cabs(G[i * m + c]) > cabs(G[p * m + c]))
        p = i;
    if (cabs(G[p * m + c]) < 1e-300)
      return false;
    if (p != c) {
      for (int j = 0; j < m; j++) {
        const dc_t tmp = G[c * m + j];
        G[c * m + j] = G[p * m + j];
        G[p * m + j] = tmp;
      }
      for (int j = 0; j < n_rhs; j++) {
        const dc_t tmp = r[c * n_rhs + j];
        r[c * n_rhs + j] = r[p * n_rhs + j];
        r[p * n_rhs + j] = tmp;
      }
    }
    for (int i = c + 1; i < m; i++) {
      const dc_t f = G[i * m + c] / G[c * m + c];
      for (int j = c; j < m; j++)
        G[i * m + j] -= f * G[c * m + j];
      for (int j = 0; j < n_rhs; j++)
        r[i * n_rhs + j] -= f * r[c * n_rhs + j];
    }
  }
  for (int c = m - 1; c >= 0; c--)
    for (int j = 0; j < n_rhs; j++) {
      dc_t s = r[c * n_rhs + j];
      for (int k = c + 1; k < m; k++)
        s -= G[c * m + k] * r[k * n_rhs + j];
      r[c * n_rhs + j] = s / G[c * m + c];
    }
  return true;
}

/// vertex offset of a parabola through three log powers around the middle one, in [-0.5, 0.5]; 0 if undefined
static double omp_vertex(double pm, double p0, double pp)
{
  if (pm <= 0.0 || p0 <= 0.0 || pp <= 0.0)
    return 0.0;
  const double ym = log(pm), y0 = log(p0), yp = log(pp);
  const double den = ym - 2.0 * y0 + yp;
  if (den >= 0.0)
    return 0.0;
  const double d = 0.5 * (ym - yp) / den;
  return (d > -0.5 && d < 0.5) ? d : 0.0;
}

/* The kernel part of an atom at fractional bin u into a ([n][B], unit gain), and z_ai = <kernel_i, r_ai> into
z ([na][n]): the correlation of the residual with a path at u, before its Doppler phase. */
static void omp_zeta(const nr_sensing_slowtime_t *s0, const dc_t *R, int na, int n, int B, double u, float *kr,
                     float *ki, dc_t *a, dc_t *z)
{
  for (int i = 0; i < n; i++)
    omp_kernel(s0->snap[i], 1.0, u, B, kr, ki, &a[(size_t)i * B]);
  for (int q = 0; q < na; q++)
    for (int i = 0; i < n; i++) {
      dc_t acc = 0.0;
      for (int b = 0; b < B; b++)
        acc += conj(a[(size_t)i * B + b]) * R[((size_t)q * n + i) * B + b];
      z[(size_t)q * n + i] = acc;
    }
}

/// sum over the antennas of |sum_i exp(-j2pi f t_i) z_ai|^2, z laid out [na][n]
static double omp_power_at(const dc_t *z, int na, int n, const double *t, double f)
{
  double p = 0.0;
  for (int q = 0; q < na; q++) {
    dc_t c = 0.0;
    for (int i = 0; i < n; i++) {
      const double turns = f * (t[i] - t[0]);
      c += cexp(-2.0 * M_PI * I * (turns - trunc(turns))) * z[(size_t)q * n + i];
    }
    p += creal(c * conj(c));
  }
  return p;
}

int nr_ue_omp_detect(const nr_sensing_slowtime_t *slow,
                     int n_ant,
                     int n_freq,
                     double f_max_hz,
                     double lambda_m,
                     const double *trend_q,
                     int trend_n_q,
                     const nr_omp_cfg_t *cfg,
                     nr_tdd_target_t *out,
                     int max_out)
{
  const nr_sensing_slowtime_t *s0 = &slow[0];
  const int n = s0->n_snap, B = s0->n_bins;
  if (n < 8 || B < 3 || n_freq < 3 || s0->snap == NULL || max_out <= 0)
    return 0;
  const double span = s0->t_s[n - 1] - s0->t_s[0];
  if (span <= 0.0)
    return 0;
  const int max_paths = cfg && cfg->max_paths > 0 ? cfg->max_paths : NR_OMP_MAX_PATHS_DEFAULT;
  const double thr = pow(10.0, (cfg && cfg->threshold_db != 0.0 ? cfg->threshold_db : NR_OMP_THRESHOLD_DB_DEFAULT) / 10.0);
  const double res = 1.0 / span; // Doppler resolution cell
  double *bin_floor = NULL;   // local_noise: per range bin, its median over the Doppler cells searched

  // antennas fitted together: same snapshots at the same times as chain 0
  int ant[NR_AOA_MAX_ANT];
  int na = 0;
  for (int a = 0; a < n_ant && a < NR_AOA_MAX_ANT; a++) {
    const nr_sensing_slowtime_t *s = &slow[a];
    if (s->n_snap != n || s->n_bins != B || s->h_re == NULL || s->h_im == NULL)
      continue;
    bool same = true;
    for (int i = 0; i < n && same; i++)
      same = fabs(s->t_s[i] - s0->t_s[i]) < 1e-9;
    if (same)
      ant[na++] = a;
  }
  if (na == 0)
    return 0;

  const int ND = 2 * B - 1; // kernel at integer shifts d = -(B-1) .. B-1, entry d + B - 1
  float *kr = malloc((size_t)ND * sizeof(float));
  float *ki = malloc((size_t)ND * sizeof(float));
  dc_t *K = malloc((size_t)n * ND * sizeof(dc_t));
  dc_t *R0 = malloc((size_t)na * n * B * sizeof(dc_t));
  dc_t *R = malloc((size_t)na * n * B * sizeof(dc_t));
  dc_t *Y = malloc((size_t)na * n * B * sizeof(dc_t));
  dc_t *E = malloc((size_t)n_freq * n * sizeof(dc_t));
  double *f = malloc((size_t)n_freq * sizeof(double));
  double *P = malloc((size_t)n_freq * B * sizeof(double));
  dc_t *A = malloc((size_t)max_paths * n * B * sizeof(dc_t));
  dc_t *G = malloc((size_t)max_paths * max_paths * sizeof(dc_t));
  dc_t *Z = malloc((size_t)max_paths * na * sizeof(dc_t));
  dc_t *tmp = malloc((size_t)B * sizeof(dc_t));
  double *sel_bin = malloc((size_t)max_paths * sizeof(double));
  double *sel_f = malloc((size_t)max_paths * sizeof(double));
  int nt = 0;
  if (!kr || !ki || !K || !R0 || !R || !Y || !E || !f || !P || !A || !G || !Z || !tmp || !sel_bin || !sel_f) {
    LOG_E(NR_PHY, "sensing omp: out of memory, map left without paths\n");
    goto done;
  }

  /* Noise weighting. The map scales snapshot i by gain[i] (sensing.gain_equalise), which scales its noise too, so
  a plain fit would be ruled by the amplified sparse snapshots, and so would the noise estimate. Every snapshot is
  divided back by its gain here, where the noise is the same in all of them: the fit below is then the least
  squares weighted by 1/gain^2, whatever the toggle says. The kernels are built at unit gain to match. */
  for (int i = 0; i < n; i++)
    omp_kernel(s0->snap[i], 1.0, (double)(B - 1), ND, kr, ki, &K[(size_t)i * ND]);
  for (int q = 0; q < na; q++) {
    const nr_sensing_slowtime_t *s = &slow[ant[q]];
    for (int i = 0; i < n; i++) {
      const double g = s0->gain && s0->gain[i] > 0.0 ? s0->gain[i] : 1.0;
      for (int b = 0; b < B; b++) {
        const size_t k = (size_t)i * B + b;
        R0[(size_t)q * n * B + k] = ((double)s->h_re[k] + I * (double)s->h_im[k]) / g;
      }
    }
  }
  memcpy(R, R0, (size_t)na * n * B * sizeof(dc_t));
  for (int j = 0; j < n_freq; j++) {
    f[j] = -f_max_hz + 2.0 * f_max_hz * j / (n_freq - 1);
    for (int i = 0; i < n; i++) {
      const double turns = f[j] * (s0->t_s[i] - s0->t_s[0]);
      E[(size_t)j * n + i] = cexp(-2.0 * M_PI * I * (turns - trunc(turns)));
    }
  }
  const double df = 2.0 * f_max_hz / (n_freq - 1);

  double noise = -1.0, e_prev = 0.0, first_best = 0.0;
  for (size_t k = 0; k < (size_t)na * n * B; k++)
    e_prev += creal(R[k] * conj(R[k]));
  int ns = 0;

  for (int it = 0; it < max_paths && nt < max_out; it++) {
    // y_ai(u) = sum_b r_ai[b] conj(K_i(b - u)): the residual against a unit path at bin u in each snapshot
    for (int q = 0; q < na; q++)
      for (int i = 0; i < n; i++) {
        const dc_t *r = &R[((size_t)q * n + i) * B];
        const dc_t *k = &K[(size_t)i * ND];
        dc_t *y = &Y[((size_t)q * n + i) * B];
        for (int u = 0; u < B; u++) {
          dc_t acc = 0.0;
          for (int b = 0; b < B; b++)
            acc += r[b] * conj(k[b - u + B - 1]);
          y[u] = acc;
        }
      }
    // correlation map over the antennas: P(f, u) = sum_a |sum_i exp(-j2pi f t_i) y_ai(u)|^2
    for (int j = 0; j < n_freq; j++)
      for (int u = 0; u < B; u++) {
        double p = 0.0;
        for (int q = 0; q < na; q++) {
          dc_t c = 0.0;
          for (int i = 0; i < n; i++)
            c += E[(size_t)j * n + i] * Y[((size_t)q * n + i) * B + u];
          p += creal(c * conj(c));
        }
        P[(size_t)j * B + u] = p;
      }
    if (noise < 0.0) { // once, from the first map: a residual's median falls as paths are taken out
      double *srt = malloc((size_t)n_freq * B * sizeof(double));
      if (srt == NULL)
        break;
      memcpy(srt, P, (size_t)n_freq * B * sizeof(double));
      qsort(srt, (size_t)n_freq * B, sizeof(double), cmp_double);
      noise = srt[(size_t)n_freq * B / 2] / M_LN2;
      if (cfg && cfg->local_noise && (bin_floor = malloc((size_t)B * sizeof(double))) != NULL)
        for (int u = 0; u < B; u++) {   // the same first correlation map: its Doppler profile per bin
          int m = 0;
          for (int j = 0; j < n_freq; j++)
            if (fabs(f[j]) > res) srt[m++] = P[(size_t)j * B + u];
          qsort(srt, (size_t)m, sizeof(double), cmp_double);
          bin_floor[u] = m > 0 ? srt[m / 2] / M_LN2 : 0.0;
        }
      free(srt);
      if (noise <= 0.0)
        break;
    }
    int jb = -1, ub = -1;
    double best = 0.0;
    for (int j = 0; j < n_freq; j++) {
      if (fabs(f[j]) <= res)
        continue; // static residue lives here
      for (int u = 0; u < B; u++) {
        /* not within one cell of a path already taken: a periodogram cannot tell two paths that close apart
        either, and what is left there is the taken path's own fitting error */
        bool near = false;
        for (int p = 0; p < ns && !near; p++)
          near = fabs(u - sel_bin[p]) <= 1.0 && fabs(f[j] - sel_f[p]) <= res;
        if (near)
          continue;
        if (bin_floor != NULL && P[(size_t)j * B + u] <= thr * bin_floor[u])
          continue;   // not above its own bin's floor: smeared over Doppler, not a point in it
        if (P[(size_t)j * B + u] > best) {
          best = P[(size_t)j * B + u];
          jb = j;
          ub = u;
        }
      }
    }
    if (jb < 0 || best <= thr * noise)
      break;
    /* below the model's accuracy floor under the strongest path: what is left there is fitting error (kernel in
    float, finite refinement), not a path */
    if (nt > 0 && best < first_best * pow(10.0, -NR_OMP_DYNAMIC_RANGE_DB / 10.0))
      break;
    if (nt == 0)
      first_best = best;

    /* Sub-bin refinement on the exact correlation, not on grid values: a refinement error leaves the derivative
    of the path's shape in the residual, which a strong path turns into extra paths at its sidelobes. Delay: a
    parabola through the exact correlation at u +- delta (0.1, then 0.01 bin). Doppler: golden section within one
    grid step. Alternated twice; each delay evaluation recomputes the masked kernels, a Doppler one only a sum. */
    double ur = ub + ((ub > 0 && ub < B - 1) ? omp_vertex(P[(size_t)jb * B + ub - 1], best, P[(size_t)jb * B + ub + 1]) : 0.0);
    double fr = f[jb] + ((jb > 0 && jb < n_freq - 1) ? omp_vertex(P[(size_t)(jb - 1) * B + ub], best, P[(size_t)(jb + 1) * B + ub]) * df : 0.0);
    dc_t *a = &A[(size_t)ns * n * B];
    for (int pass = 0; pass < 2; pass++) {
      /* delay: Newton-like parabola steps on the exact correlation, the probe spacing following the step, until
      the step is under 1e-4 bin (at most 8 steps, each step at most 0.5 bin from the grid's estimate) */
      double dl = 0.05;
      const double u_start = ur;
      for (int k = 0; k < 8; k++) {
        double pw[3];
        for (int m = 0; m < 3; m++) {
          omp_zeta(s0, R, na, n, B, ur + (m - 1) * dl, kr, ki, a, Y);
          pw[m] = omp_power_at(Y, na, n, s0->t_s, fr);
        }
        const double den = pw[0] - 2.0 * pw[1] + pw[2];
        if (den >= 0.0)
          break; // not a maximum here: keep what we have
        double step = 0.5 * (pw[0] - pw[2]) / den * dl;
        if (step > 2.0 * dl) step = 2.0 * dl;
        if (step < -2.0 * dl) step = -2.0 * dl;
        ur += step;
        if (ur > u_start + 0.5) ur = u_start + 0.5;
        if (ur < u_start - 0.5) ur = u_start - 0.5;
        if (fabs(step) < 1e-4)
          break;
        dl = fmax(fabs(step), 0.002);
      }
      omp_zeta(s0, R, na, n, B, ur, kr, ki, a, Y); // kernel part at ur in a, z_ai in Y
      double lo = fr - df, hi = fr + df;
      const double gr = 0.5 * (sqrt(5.0) - 1.0);
      double x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo);
      double p1 = omp_power_at(Y, na, n, s0->t_s, x1), p2 = omp_power_at(Y, na, n, s0->t_s, x2);
      for (int k = 0; k < 40; k++) {
        if (p1 >= p2) { hi = x2; x2 = x1; p2 = p1; x1 = hi - gr * (hi - lo); p1 = omp_power_at(Y, na, n, s0->t_s, x1); }
        else { lo = x1; x1 = x2; p1 = p2; x2 = lo + gr * (hi - lo); p2 = omp_power_at(Y, na, n, s0->t_s, x2); }
      }
      fr = 0.5 * (lo + hi);
    }
    for (int i = 0; i < n; i++) { // the atom as the map holds it: Doppler phase and the snapshot's gain
      const double turns = fr * (s0->t_s[i] - s0->t_s[0]);
      const double g = s0->gain && s0->gain[i] > 0.0 ? s0->gain[i] : 1.0;
      const dc_t ph = g * cexp(2.0 * M_PI * I * (turns - trunc(turns)));
      for (int b = 0; b < B; b++)
        a[(size_t)i * B + b] *= ph;
    }
    sel_bin[ns] = ur;
    sel_f[ns] = fr;
    for (int k = 0; k < trend_n_q && trend_q; k++) // the slow trend the clutter stage removed from the data
      for (int b = 0; b < B; b++) {
        dc_t c = 0.0;
        for (int i = 0; i < n; i++)
          c += trend_q[(size_t)k * n + i] * a[(size_t)i * B + b];
        for (int i = 0; i < n; i++)
          a[(size_t)i * B + b] -= c * trend_q[(size_t)k * n + i];
      }
    for (int i = 0; i < n; i++) { // back to the weighted domain of R0
      const double g = s0->gain && s0->gain[i] > 0.0 ? s0->gain[i] : 1.0;
      for (int b = 0; b < B; b++)
        a[(size_t)i * B + b] /= g;
    }
    ns++;

    // joint least squares over every atom so far, one amplitude per atom and antenna
    for (int p = 0; p < ns; p++)
      for (int qq = 0; qq < ns; qq++) {
        dc_t s = 0.0;
        const dc_t *ap = &A[(size_t)p * n * B], *aq = &A[(size_t)qq * n * B];
        for (size_t k = 0; k < (size_t)n * B; k++)
          s += conj(ap[k]) * aq[k];
        G[p * ns + qq] = s;
      }
    for (int p = 0; p < ns; p++)
      for (int q = 0; q < na; q++) {
        dc_t s = 0.0;
        const dc_t *ap = &A[(size_t)p * n * B], *r0 = &R0[(size_t)q * n * B];
        for (size_t k = 0; k < (size_t)n * B; k++)
          s += conj(ap[k]) * r0[k];
        Z[p * na + q] = s;
      }
    if (!omp_solve(G, ns, Z, na)) {
      LOG_W(NR_PHY, "sensing omp: singular fit after %d paths, stopped\n", ns);
      ns--;
      break;
    }
    double e = 0.0;
    for (int q = 0; q < na; q++)
      for (size_t k = 0; k < (size_t)n * B; k++) {
        dc_t v = R0[(size_t)q * n * B + k];
        for (int p = 0; p < ns; p++)
          v -= Z[p * na + q] * A[(size_t)p * n * B + k];
        R[(size_t)q * n * B + k] = v;
        e += creal(v * conj(v));
      }
    if (e > e_prev * (1.0 + 1e-9)) // least squares over a growing set can only lower the residual
      LOG_W(NR_PHY, "sensing omp: residual grew %.3e -> %.3e at path %d\n", e_prev, e, ns);
    e_prev = e;

    out[nt].bin = ur;
    out[nt].range_m = ur * s0->m_per_bin;
    out[nt].f_hz = fr;
    out[nt].speed_ms = fr * lambda_m / 2.0;
    out[nt].snr_dB = 10.0 * log10(best / noise);
    out[nt].iteration = it;
    out[nt].verdict = NR_TDD_VERDICT_TARGET;
    nt++;
  }

done:
  free(bin_floor);
  free(kr);
  free(ki);
  free(K);
  free(R0);
  free(R);
  free(Y);
  free(E);
  free(f);
  free(P);
  free(A);
  free(G);
  free(Z);
  free(tmp);
  free(sel_bin);
  free(sel_f);
  return nt;
}
