/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/* Alignment of the static scene across snapshots, and grouping of snapshots by grant */

#include <math.h>
#include <pthread.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_sensing_align.h"
#include "nrscope/hdr/sensing/nr_ue_map.h"

/* One stream's reference and tracking state. Streams, not antennas: the correction is
   estimated once per symbol and shared by every antenna, so one state serves them all. */
typedef struct {
  bool used;
  nr_sensing_stream_t key;
  bool have_ref;
  /// reference of the static scene, lattice indexed, and the positions it holds
  float ref_re[NR_SENSING_MAX_IDFT];
  float ref_im[NR_SENSING_MAX_IDFT];
  bool ref_valid[NR_SENSING_MAX_IDFT];
  /// shift found for the previous symbol: the centre of the next search
  double delay;
  /// symbols in a row that did not match the reference
  int n_bad;
  /// the last symbol handled, so the other antennas reuse its correction
  bool last_valid;
  uint64_t last_t;
  nr_sensing_align_t last;
} align_state_t;

static align_state_t align_states[NR_SENSING_MAX_STREAMS];
static pthread_mutex_t align_lock = PTHREAD_MUTEX_INITIALIZER;

/// The state of one stream, claimed on first use; NULL when the table is full
static align_state_t *align_state(const nr_sensing_stream_t *key)
{
  align_state_t *free_slot = NULL;
  for (int i = 0; i < NR_SENSING_MAX_STREAMS; i++) {
    if (!align_states[i].used) {
      if (free_slot == NULL)
        free_slot = &align_states[i];
      continue;
    }
    if (nr_ue_sensing_stream_eq(&align_states[i].key, key))
      return &align_states[i];
  }
  if (free_slot != NULL) {
    memset(free_slot, 0, sizeof(*free_slot));
    free_slot->used = true;
    free_slot->key = *key;
  }
  return free_slot;
}

/// Start the reference over from these pilots
static void align_ref_reset(align_state_t *s, const cf_t *p, int a, int n)
{
  memset(s->ref_valid, 0, sizeof(s->ref_valid));
  for (int i = a; i < a + n; i++) {
    s->ref_re[i] = crealf(p[i]);
    s->ref_im[i] = cimagf(p[i]);
    s->ref_valid[i] = true;
  }
  s->have_ref = true;
  s->delay = 0.0;
  s->n_bad = 0;
}

/// Blend aligned pilots into the reference; positions never seen before are taken as they are
static void align_ref_update(align_state_t *s, const cf_t *p, int a, int n)
{
  const float w = NR_SENSING_ALIGN_REF_WEIGHT;
  for (int i = a; i < a + n; i++) {
    if (!s->ref_valid[i]) {
      s->ref_re[i] = crealf(p[i]);
      s->ref_im[i] = cimagf(p[i]);
      s->ref_valid[i] = true;
    } else {
      s->ref_re[i] += w * ((float)crealf(p[i]) - s->ref_re[i]);
      s->ref_im[i] += w * ((float)cimagf(p[i]) - s->ref_im[i]);
    }
  }
}

/* c(d) = sum_k x_k exp(+j 2 pi (lo + k) d / N), with x_k = p conj(R) at lattice position lo + k.
   The phasor is advanced by a fixed step rather than recomputed per position. 
   INSERT HERE EQUATION NUMBER AT FINAL REPORT */
static void align_corr(const double *xr, const double *xi, int lo, int len, int N, double d, double *cr, double *ci)
{
  const double turns = (double)lo * d / N;
  const double ph0 = 2.0 * M_PI * (turns - trunc(turns));
  double rr = cos(ph0);
  double ri = sin(ph0);
  const double sr = cos(2.0 * M_PI * d / N);
  const double si = sin(2.0 * M_PI * d / N);
  double ar = 0.0;
  double ai = 0.0;
  for (int k = 0; k < len; k++) {
    ar += xr[k] * rr - xi[k] * ri;
    ai += xr[k] * ri + xi[k] * rr;
    const double nr = rr * sr - ri * si;
    ri = rr * si + ri * sr;
    rr = nr;
  }
  *cr = ar;
  *ci = ai;
}

/* Saturation to the int16 range went with the Q15 representation and is gone.
In float the rotation below is exact to within rounding, and clamping the result
to +-32767 would now be an arbitrary clip applied to values whose scale is set
by the estimator rather than by the sample type. */

/// Undo a delay of d bins, a phase of phi and a gain of amp:
///     p_i <- p_i exp(+j (2 pi i d / N - phi)) / amp
/// amp is 1 when only the delay and the phase are corrected.
/// INSERT HERE EQUATION NUMBER AT FINAL REPORT
static void align_apply(cf_t *p, int a, int n, int N, double d, double phi, double amp)
{
  const double turns = (double)a * d / N;
  const double ph0 = 2.0 * M_PI * (turns - trunc(turns)) - phi;
  /* The 1/amp is folded into the starting phasor. The recurrence below only rotates it
  (sr^2 + si^2 = 1), so its length, hence the scale, is carried unchanged to every pilot. */
  double rr = cos(ph0) / amp;
  double ri = sin(ph0) / amp;
  const double sr = cos(2.0 * M_PI * d / N);
  const double si = sin(2.0 * M_PI * d / N);
  for (int i = a; i < a + n; i++) {
    const double re = crealf(p[i]) * rr - cimagf(p[i]) * ri;
    const double im = crealf(p[i]) * ri + cimagf(p[i]) * rr;
    p[i] = (float)re + I * (float)im;
    const double nr = rr * sr - ri * si;
    ri = rr * si + ri * sr;
    rr = nr;
  }
}

bool nr_ue_sensing_align_symbol(uint64_t t_sample,
                                nr_sensing_stream_t stream,
                                const pilot_lattice_t *lat,
                                cf_t *pilots,
                                nr_sensing_align_t *out)
{
  const int N = NR_SENSING_IDFT_SIZE(lat->k_step);
  const int a = (lat->k_first - lat->k_first % lat->k_step) / lat->k_step;
  const int n = lat->n;
  nr_sensing_align_t res = {.corr = 1.0, .amp = 1.0};

  if (a < 0 || a + n > NR_SENSING_MAX_IDFT || n < 2) {
    if (out)
      *out = res;
    return true; // not a lattice this file can handle, leave it as it is
  }

  pthread_mutex_lock(&align_lock);
  const nr_sensing_stream_t key = nr_ue_sensing_stream(stream.ports, stream.layer, lat);
  align_state_t *s = align_state(&key);
  if (s == NULL) {
    pthread_mutex_unlock(&align_lock);
    if (out)
      *out = res;
    return true; // more streams than the table holds: this one is not aligned
  }

  // Another antenna of a symbol already handled: same correction, reference untouched.
  if (s->last_valid && s->last_t == t_sample) {
    res = s->last;
    pthread_mutex_unlock(&align_lock);
    if (res.applied)
      align_apply(pilots, a, n, N, res.delay_bins, res.phase_rad, res.amp);
    if (out)
      *out = res;
    return !res.dropped;
  }

  if (!s->have_ref) {
    // The first symbol of a stream is the reference: nothing to correct yet.
    align_ref_reset(s, pilots, a, n);
  } else {
    // a. cross terms p conj(R) over the positions both hold
    double xr[NR_SENSING_MAX_IDFT];
    double xi[NR_SENSING_MAX_IDFT];
    double e_p = 0.0;
    double e_r = 0.0;
    int n_overlap = 0;
    for (int k = 0; k < n; k++) {
      const int i = a + k;
      if (!s->ref_valid[i]) {
        xr[k] = xi[k] = 0.0;
        continue;
      }
      const double pr = crealf(pilots[i]), pi = cimagf(pilots[i]);
      const double rr = s->ref_re[i], ri = s->ref_im[i];
      xr[k] = pr * rr + pi * ri;
      xi[k] = pi * rr - pr * ri;
      e_p += pr * pr + pi * pi;
      e_r += rr * rr + ri * ri;
      n_overlap++;
    }

    bool good = false;
    if (n_overlap >= NR_SENSING_ALIGN_MIN_OVERLAP && e_p > 0.0 && e_r > 0.0) {
      
      // coarse search around the previous shift ...
      const int n_steps = (int)lround(NR_SENSING_ALIGN_SEARCH_BINS / NR_SENSING_ALIGN_STEP_BINS);
      double pw[2 * n_steps + 1];
      int best = -n_steps; // the first point computed, so no comparison reads an unfilled entry
      for (int k = -n_steps; k <= n_steps; k++) {
        double cr, ci;
        align_corr(xr, xi, a, n, N, s->delay + k * NR_SENSING_ALIGN_STEP_BINS, &cr, &ci);
        pw[k + n_steps] = cr * cr + ci * ci;
        if (pw[k + n_steps] > pw[best + n_steps])
          best = k;
      }
      // ... then a parabola through the best point and its neighbours
      double frac = 0.0;
      if (best > -n_steps && best < n_steps) {
        const double y0 = pw[best + n_steps - 1], y1 = pw[best + n_steps], y2 = pw[best + n_steps + 1];
        const double den = y0 - 2.0 * y1 + y2;
        if (den < 0.0) {
          const double f = 0.5 * (y0 - y2) / den;
          if (f > -1.0 && f < 1.0)
            frac = f;
        }
      }
      const double d_hat = s->delay + (best + frac) * NR_SENSING_ALIGN_STEP_BINS;
      double cr, ci;
      align_corr(xr, xi, a, n, N, d_hat, &cr, &ci);

      res.delay_bins = d_hat;
      res.phase_rad = atan2(ci, cr);

      // Compute the correlation (normalized) around the estimated bin D hat. If it was that the reference symbol is not 
      // a reference at all (issues with sender or anything else), then we drop it
      res.corr = sqrt(cr * cr + ci * ci) / sqrt(e_p * e_r);
      good = res.corr >= NR_SENSING_ALIGN_MIN_CORR;

#if NR_SENSING_ALIGN_AMPLITUDE
      /* Magnitude of the least squares gain of this symbol against the reference. With
      p_i = g R_i exp(-j 2 pi i D / N), the correlation at D is c = g * sum |R_i|^2 = g * e_r,
      so |g| = |c| / e_r. The phase of g is phase_rad above; this is the rest of it. e_r is
      summed over the same overlapping positions as c, so the two are consistent even when
      the grant covers only part of the reference. See NR_SENSING_ALIGN_AMPLITUDE in the
      header for why the amplitude has to go as well.

      Clamped rather than trusted blindly: past the clamp the scene changed rather than its
      level, see NR_SENSING_ALIGN_AMP_MIN. */
      double amp = sqrt(cr * cr + ci * ci) / e_r;
      if (amp < NR_SENSING_ALIGN_AMP_MIN)
        amp = NR_SENSING_ALIGN_AMP_MIN;
      if (amp > NR_SENSING_ALIGN_AMP_MAX)
        amp = NR_SENSING_ALIGN_AMP_MAX;
      res.amp = amp;
#endif
    } else {

      // Something bad happens with the static scene, we cannot follow it
      res.corr = 0.0;
    }

    if (good) {
      /* b. undo the shift, the phase and the gain, c. let the reference follow. The
      reference is fed the normalised symbol, so it stays at the level of the symbol it was
      started from and does not drift towards the fluctuation it is there to remove. */
      align_apply(pilots, a, n, N, res.delay_bins, res.phase_rad, res.amp);
      align_ref_update(s, pilots, a, n);
      res.applied = true;
      s->delay = res.delay_bins;
      s->n_bad = 0;
    } else if (++s->n_bad >= NR_SENSING_ALIGN_RESET_AFTER) {
      LOG_W(NR_PHY,
            "sensing: static scene no longer matches its reference (corr %.2f) for %d symbols, "
            "reference restarted (ports 0x%x layer %d)\n",
            res.corr,
            NR_SENSING_ALIGN_RESET_AFTER,
            key.ports,
            key.layer);
      align_ref_reset(s, pilots, a, n);
      res = (nr_sensing_align_t){.corr = 1.0, .amp = 1.0};
    } else {
      res.dropped = true;
    }
  }

  s->last_valid = true;
  s->last_t = t_sample;
  s->last = res;
  pthread_mutex_unlock(&align_lock);

  if (out)
    *out = res;
  return !res.dropped;
}

int nr_ue_sensing_group_by_grant(const nr_sensing_history_t *hist, int idx[], int n)
{
  if (n < 2)
    return n;

  // the most common (pilot count, first pilot); scanning newest first gives a tie to the latest grant
  int best = n - 1;
  int best_count = 0;
  for (int i = n - 1; i >= 0; i--) {
    const nr_sensing_snapshot_t *s = &hist->ring[idx[i]];
    int count = 0;
    for (int j = 0; j < n; j++) {
      const nr_sensing_snapshot_t *q = &hist->ring[idx[j]];
      count += (q->n_pilots == s->n_pilots && q->k_first == s->k_first);
    }
    if (count > best_count) {
      best_count = count;
      best = i;
    }
  }

  const int n_keep = hist->ring[idx[best]].n_pilots;
  const int k_keep = hist->ring[idx[best]].k_first;
  int kept = 0;
  for (int i = 0; i < n; i++) {
    const nr_sensing_snapshot_t *s = &hist->ring[idx[i]];
    if (s->n_pilots == n_keep && s->k_first == k_keep)
      idx[kept++] = idx[i];
  }
  LOG_D(NR_PHY, "sensing: grouped by grant, %d of %d snapshots kept (%d pilots from k=%d)\n", kept, n, n_keep, k_keep);
  return kept;
}
