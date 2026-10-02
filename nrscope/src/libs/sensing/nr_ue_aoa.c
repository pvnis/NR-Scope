/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_aoa.h"
#include "nrscope/hdr/sensing/nr_ue_tdd_detect.h"

/* Per-chain phase offsets, see NR_AOA_CAL_DEG. */
static const double nr_aoa_cal_deg[NR_AOA_MAX_ANT] = NR_AOA_CAL_DEG;

void nr_ue_aoa_cell_vector(const nr_sensing_slowtime_t *slow,
                           int n_ant,
                           int bin,
                           double f_hz,
                           double *x_re,
                           double *x_im)
{
  for (int a = 0; a < n_ant; a++) {
    const nr_sensing_slowtime_t *s = &slow[a];
    double ar = 0.0;
    double ai = 0.0;

    for (int i = 0; i < s->n_snap; i++) {
      const double hr = s->h_re[(size_t)i * s->n_bins + bin];
      const double hi = s->h_im[(size_t)i * s->n_bins + bin];
      /* Same sign convention and the same modulo-one-turn reduction as the Doppler
      transform in nr_ue_sensing_range_doppler(): C[f] = sum_i h_i exp(-j*2*pi*f*t_i),
      with f*t reaching hundreds of turns over a long window, where cos() would spend
      most of the mantissa on the integer part of the phase. The convention has to
      match or the angle would be read off a conjugated cell. */
      const double turns = f_hz * s->t_s[i];
      const double ph = -2.0 * M_PI * (turns - trunc(turns));
      const double wr = cos(ph);
      const double wi = sin(ph);
      ar += hr * wr - hi * wi;
      ai += hr * wi + hi * wr;
    }

    /* Calibration applied after the transform rather than to every sample: it is a
    property of the chain and constant over slow time, so it commutes with the sum
    and one rotation per antenna does the same work as n_snap of them. */
    const double c = nr_aoa_cal_deg[a] * M_PI / 180.0;
    const double cr = cos(c);
    const double ci = -sin(c); // exp(-j*c)
    x_re[a] = ar * cr - ai * ci;
    x_im[a] = ar * ci + ai * cr;
  }
}

double nr_ue_aoa_estimate(const double *x_re, const double *x_im, int n_ant, double *quality_dB)
{
  double p[NR_AOA_GRID];
  double p_sum = 0.0;
  double p_peak = -1.0;
  int i_peak = 0;

  for (int g = 0; g < NR_AOA_GRID; g++) {
    const double sin_th = -1.0 + 2.0 * (double)g / (NR_AOA_GRID - 1);
    double ar = 0.0;
    double ai = 0.0;
    /* a_R^H(theta) x, with a_R(theta)[m] = exp(j*2*pi*(d/lambda)*m*sin(theta)) and
    element 0 as the phase reference, matching eps_m*d = m*d of a ULA. */
    for (int a = 0; a < n_ant; a++) {
      const double ph = -2.0 * M_PI * NR_AOA_ELEMENT_SPACING * a * sin_th;
      const double wr = cos(ph);
      const double wi = sin(ph);
      ar += x_re[a] * wr - x_im[a] * wi;
      ai += x_re[a] * wi + x_im[a] * wr;
    }
    p[g] = ar * ar + ai * ai;
    p_sum += p[g];
    if (p[g] > p_peak) {
      p_peak = p[g];
      i_peak = g;
    }
  }

  if (quality_dB != NULL) {
    const double p_mean = p_sum / NR_AOA_GRID;
    /* Peak over mean. Mr for a single plane wave, so 6 dB at four antennas; a cell
    holding two arrivals, or nothing but noise, falls short of that. */
    *quality_dB = (p_mean > 0.0 && p_peak > 0.0) ? 10.0 * log10(p_peak / p_mean) : 0.0;
  }

  /* Sub-grid refinement by a parabola through the peak and its neighbours. The grid
  is already fine enough that this moves the answer very little; it is here so the
  reported angle does not sit on grid points, which reads as quantisation when the
  same target is tracked across maps. */
  double idx = (double)i_peak;
  if (i_peak > 0 && i_peak < NR_AOA_GRID - 1) {
    const double den = p[i_peak - 1] - 2.0 * p[i_peak] + p[i_peak + 1];
    if (den != 0.0) {
      const double delta = 0.5 * (p[i_peak - 1] - p[i_peak + 1]) / den;
      if (delta > -1.0 && delta < 1.0)
        idx += delta;
    }
  }

  double sin_th = -1.0 + 2.0 * idx / (NR_AOA_GRID - 1);
  /* The refinement can step outside the unit interval at the very edge of the grid,
  where asin() would return NaN and poison the dump. */
  if (sin_th > 1.0) sin_th = 1.0;
  if (sin_th < -1.0) sin_th = -1.0;

  return asin(sin_th) * 180.0 / M_PI;
}

int nr_ue_aoa_process(const nr_sensing_map_t *map,
                      const nr_sensing_slowtime_t *slow,
                      int n_ant,
                      const nr_aoa_cell_t *cells,
                      int n_cells,
                      nr_sensing_aoa_t out[NR_SENSING_AOA_MAX])
{
  if (n_ant < 2) return 0; // one antenna carries no spatial information

  /* Every antenna must have produced samples over the same window, or the spatial
  vector would mix cells from different observations. The transform is fed the same
  reference symbols for all of them, so a mismatch means one of the histories failed
  and the map it would contribute is not there either. */
  for (int a = 0; a < n_ant; a++) {
    if (slow[a].n_snap < 1 || slow[a].n_bins != map->n_bins) {
      LOG_W(NR_PHY,
            "aoa: rx%d has %d snapshots over %d bins against %d bins on the map, no angles\n",
            a,
            slow[a].n_snap,
            slow[a].n_bins,
            map->n_bins);
      return 0;
    }
    if (slow[a].n_snap != slow[0].n_snap) {
      LOG_W(NR_PHY, "aoa: rx%d gathered %d snapshots against %d on rx0, no angles\n", a, slow[a].n_snap, slow[0].n_snap);
      return 0;
    }
    /* The same instants, not only as many. A window that starts elsewhere but holds
    the same count passes the test above, and the phase between the antennas would
    then carry each target's Doppler times the offset on top of its angle. The
    chains share the sample clock, so equal means equal to the sample. */
    bool same = slow[a].t0_sample == slow[0].t0_sample;
    for (int i = 0; same && i < slow[0].n_snap; i++)
      same = slow[a].t_s[i] == slow[0].t_s[i];
    if (!same) {
      LOG_W(NR_PHY, "aoa: rx%d sampled other instants than rx0 over its %d snapshots, no angles\n", a, slow[a].n_snap);
      return 0;
    }
  }

  if (map->n_bins * map->n_freq < 1) return 0;

  /* The speed of a cell is its Doppler times lambda/2, so a map carrying no carrier
  has no speed axis to report on. Bail rather than divide by zero and fill the dump
  with infinities. */
  if (map->carrier_hz <= 0.0) {
    LOG_W(NR_PHY, "aoa: map carries no carrier frequency, no angles\n");
    return 0;
  }

  /* First delay bin that can hold a target.

  The direct path occupies map->bin_los and the bins its mainlobe covers. Nothing at
  or below it can be a target: a reflected path is never shorter than the direct one.

  The mainlobe half width follows from the taper rather than a constant: the pilots
  are Hann windowed over win_len lattice points and transformed to idft_size, so the
  response is oversampled by idft_size/win_len and the Hann half width of two lattice
  bins becomes 2*idft_size/win_len delay bins. A narrow grant widens it, which is
  exactly when a fixed guard would be too small. */
  const int guard = (slow[0].win_len > 0) ? (int)ceil(2.0 * slow[0].idft_size / slow[0].win_len) : NR_AOA_LOS_GUARD_MIN;
  const int b_min = (int)ceil(map->bin_los) + (guard > NR_AOA_LOS_GUARD_MIN ? guard : NR_AOA_LOS_GUARD_MIN);

  /* Nothing left to search. Happens when the direct path lands near the end of the
  truncated axis, which means the window is mistimed rather than that the scene is
  empty, so it is worth saying out loud. */
  if (b_min >= map->n_bins) {
    LOG_W(NR_PHY,
          "aoa: LoS at bin %.2f leaves no searchable bins below %d, no angles\n",
          map->bin_los,
          map->n_bins);
    return 0;
  }

  // ---- 1. The cells to measure -------------------------------------------------

  nr_aoa_cell_t cand[NR_AOA_MAX_CANDIDATES];
  int n_cand = 0;

  if (cells != NULL) {
    /* --sensing-tdd-detect on: take the detector's targets as they are. They already
    carry their true Doppler, even outside the map's window. */
    n_cand = n_cells < NR_AOA_MAX_CANDIDATES ? n_cells : NR_AOA_MAX_CANDIDATES;
    if (n_cand > 0)
      memcpy(cand, cells, (size_t)n_cand * sizeof(*cand));
  } else {
    /* Detector off: find the peaks on the averaged map with CFAR. Each cell is compared
    with the average of the cells around it (not with the whole map), so a target next
    to strong clutter and one in a quiet area are both judged fairly.

    The CFAR threshold formula assumes a single map. Ours is the average of several
    (antennas x layers), which makes the noise smoother, so the real false-alarm rate
    is lower than the configured one: on the safe side. */
    nr_tdd_cfg_t cfg;
    nr_tdd_cfg_default(&cfg);
    int c_bin[NR_AOA_MAX_CANDIDATES];
    int c_freq[NR_AOA_MAX_CANDIDATES];
    double c_snr[NR_AOA_MAX_CANDIDATES];
    n_cand = nr_tdd_cfar_power(map->power, map->n_bins, map->n_freq, &cfg, NR_AOA_MAX_CANDIDATES, c_bin, c_freq, c_snr);

    /* Grid index to Hz. Same grid as nr_ue_sensing_range_doppler(), and it has to stay
    that way: an angle read at a frequency the map does not use would be the angle of
    a different cell. */
    const double df = (map->n_freq > 1) ? (2.0 * map->f_max_hz) / (map->n_freq - 1) : 0.0;
    for (int i = 0; i < n_cand; i++)
      cand[i] = (nr_aoa_cell_t){.bin = c_bin[i], .f_hz = -map->f_max_hz + c_freq[i] * df, .snr_dB = c_snr[i]};
  }

  // ---- 2. The angle of each cell, strongest first -------------------------------

  const double lambda = C_M_PER_S / map->carrier_hz;
  int kept_bin[NR_SENSING_AOA_MAX];
  int n_out = 0;

  for (int i = 0; i < n_cand && n_out < NR_SENSING_AOA_MAX; i++) {
    const int bin = (int)lround(cand[i].bin);

    // On or below the direct path, or off the map: not a target
    if (bin < b_min || bin >= map->n_bins)
      continue;

    // One complex number per antenna at this cell, then the angle from those numbers
    double x_re[NR_AOA_MAX_ANT];
    double x_im[NR_AOA_MAX_ANT];
    nr_ue_aoa_cell_vector(slow, n_ant, bin, cand[i].f_hz, x_re, x_im);

    double quality_dB = 0.0;
    const double angle_deg = nr_ue_aoa_estimate(x_re, x_im, n_ant, &quality_dB);

    /* Same delay bin and same angle as a stronger cell already kept: that is the same
    place, so the same object. Skip it. See NR_AOA_SEP_BINS. */
    bool merged = false;
    for (int j = 0; j < n_out && !merged; j++)
      merged = abs(bin - kept_bin[j]) <= NR_AOA_SEP_BINS && fabs(angle_deg - out[j].angle_deg) <= NR_AOA_SEP_ANGLE_DEG;
    if (merged)
      continue;

    kept_bin[n_out] = bin;

    /* range_m is the extra path length over the direct path, which is what the
    localisation needs. It is measured from bin_los and not from bin 0, because the
    sync puts the direct path a couple of bins away from 0. */
    out[n_out++] = (nr_sensing_aoa_t){
        .range_m = (float)((cand[i].bin - map->bin_los) * map->m_per_bin),
        .speed_ms = (float)(cand[i].f_hz * lambda / 2.0),
        .angle_deg = (float)angle_deg,
        .power_dB = (float)cand[i].snr_dB,
        .quality_dB = (float)quality_dB,
    };
  }

  return n_out;
}
