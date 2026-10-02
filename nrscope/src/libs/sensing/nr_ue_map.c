/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_sensing.h"
#include "nrscope/hdr/sensing/nr_ue_map.h"
#include "nrscope/hdr/sensing/nr_ue_sensing_align.h"

bool nr_ue_sensing_history_init(nr_sensing_history_t *hist,
                                int depth,
                                int scs_hz,
                                int ofdm_symbol_size,
                                uint64_t carrier_hz)
{
  pthread_mutex_init(&hist->lock, NULL);
  hist->ring = calloc(depth, sizeof(*hist->ring));
  if (hist->ring == NULL) {
    LOG_E(NR_PHY, "sensing: cannot allocate history of %d snapshots\n", depth);
    return false;
  }
  hist->depth = depth;
  hist->head = 0;
  hist->count = 0;
  hist->scs_hz = scs_hz;
  hist->ofdm_symbol_size = ofdm_symbol_size;
  hist->carrier_hz = carrier_hz;
  return true;
}

void nr_ue_sensing_history_free(nr_sensing_history_t *hist)
{
  free(hist->ring);
  hist->ring = NULL;
  hist->depth = 0;
  hist->head = 0;
  hist->count = 0;
  pthread_mutex_destroy(&hist->lock);
}

/* The counter slot of one stream, claimed on first use. NULL when the table is full,
which only happens if the scheduler has cycled through more distinct allocations than
NR_SENSING_MAX_STREAMS; such a stream is still recorded in the ring and still
transformable, it just never triggers a map on its own. */
static int *nr_ue_sensing_counter(nr_sensing_history_t *hist, const nr_sensing_stream_t *stream, bool claim)
{
  int free_slot = -1;
  for (int i = 0; i < NR_SENSING_MAX_STREAMS; i++) {
    if (!hist->since_last_map[i].used) {
      if (free_slot < 0)
        free_slot = i;
      continue;
    }
    if (nr_ue_sensing_stream_eq(&hist->since_last_map[i].key, stream))
      return &hist->since_last_map[i].count;
  }

  if (!claim || free_slot < 0) {
    if (claim)
      LOG_W(NR_PHY,
            "sensing: more than %d streams seen, ports 0x%x layer %d k_step %d will not trigger a map\n",
            NR_SENSING_MAX_STREAMS,
            stream->ports,
            stream->layer,
            stream->k_step);
    return NULL;
  }

  hist->since_last_map[free_slot].used = true;
  hist->since_last_map[free_slot].key = *stream;
  hist->since_last_map[free_slot].count = 0;
  return &hist->since_last_map[free_slot].count;
}

void nr_ue_sensing_history_push(nr_sensing_history_t *hist,
                                uint64_t t_sample,
                                nr_sensing_stream_t stream,
                                const pilot_lattice_t *lattice,
                                int idft_size,
                                const cf_t *h_time)
{
  if (hist->ring == NULL) return;

  /* The caller names the allocation and the layer, the comb comes from the data, so
  the two halves of the key cannot drift apart. */
  stream = nr_ue_sensing_stream(stream.ports, stream.layer, lattice);

  // other DL actors push into this history concurrently, and a map may be copying it
  pthread_mutex_lock(&hist->lock);
  nr_sensing_snapshot_t *s = &hist->ring[hist->head];
  s->t_sample = t_sample;
  s->stream = stream;
  s->n_pilots = lattice->n;
  s->k_first = lattice->k_first;
  s->idft_size = idft_size;

  // The IDFT has already been computed with the 0-padded idft_size,
  // idft_size is always 2048 here, so n_bins is always NR_SENSING_DELAY_BINS
  //s->n_bins = idft_size < NR_SENSING_DELAY_BINS ? idft_size : NR_SENSING_DELAY_BINS;
  
  // always take a portion of the bins, we limit the range for computing the map
  s->n_bins = NR_SENSING_MAP_MAX_BINS_RANGE;
  
  // copy the part of the range profile we want
  memcpy(s->h, h_time, s->n_bins * sizeof(cf_t));

  hist->head = (hist->head + 1) % hist->depth;
  if (hist->count < hist->depth)
    hist->count++;
  
  // count the number of symbol we have in the stream
  int *count = nr_ue_sensing_counter(hist, &stream, true);
  if (count)
    (*count)++;
  pthread_mutex_unlock(&hist->lock);
}

/* One TDD period in seconds. The slow-time sampling pattern repeats over this, which
is what makes the Doppler spectrum repeat every 1/T and sets both the replica spacing
the detector looks for and the comb the removal below takes out. */
static double nr_ue_sensing_tdd_period_s(const nr_sensing_history_t *hist)
{
  const double slot_dur_s = 1e-3 / ((double)hist->scs_hz / 15000.0); // 0.5 ms at 30 kHz
  return NR_SENSING_TDD_PERIOD_SLOTS * slot_dur_s;
}

// Compute maximum doppler frequency for a given maximum speed and carrier frequency
// Recall f_d = 2*v/lambda = 2*v*f_c/c (and we add + 1/tdd)
static double nr_ue_sensing_axis_f_max(const nr_sensing_history_t *hist, double max_speed_ms)
{
  const double t_period_s = nr_ue_sensing_tdd_period_s(hist);

  /* f = 2*v/lambda, written with the carrier instead of lambda to keep it one line.
  A history with no carrier cannot convert a speed, so it falls back to the margin
  alone, which is the minimum the detector needs to work at all. */
  const double f_speed = (hist->carrier_hz > 0 && max_speed_ms > 0.0)
                             ? 2.0 * max_speed_ms * (double)hist->carrier_hz / C_M_PER_S
                             : 0.0;

  /* Plus one period so the TDD detector can look for a candidate's copies on both sides */
  return f_speed + 1.0 / t_period_s;
}

/* Criterion 1. Computationally, we want 2 doppler points per cell. But we
do not want to exceed the map's maximum number of bins, NR_SENSING_MAP_MAX_BINS_FREQ,
which is a hard limit on the size of the power map.

n_freq = (int)(4 * f_max * t_span) | 1 must stay at
or below the last odd index the map can hold; the | 1 can only round an even count up
by one, so requiring the product itself to fit is enough and the bound does not depend
on how the count is rounded. 

Recall the formula in .h file
*/
double nr_ue_sensing_span_grid(double f_max_hz)
{
  if (!(f_max_hz > 0.0))
    return HUGE_VAL;
  return (double)(NR_SENSING_MAP_MAX_BINS_FREQ - 1) / (4.0 * f_max_hz);
}

/* Criterion 2, Range migration. The cell is the width of a peak,
c / (n_pilots * k_step * scs), not m_per_bin. m_per_bin is c / (idft_size * k_step *
scs) and the IDFT zero pads the grant out to idft_size, so the axis is oversampled and
a bin is finer than a peak. Drifting by one bin costs nothing when the peak is over a
bin wide; drifting by one cell is what breaks coherent integration. */
double nr_ue_sensing_span_migration(int n_pilots, int k_step, int scs_hz, double max_speed_ms)
{
  if (n_pilots <= 0 || k_step <= 0 || scs_hz <= 0 || !(max_speed_ms > 0.0))
    return HUGE_VAL;
  const double cell_m = C_M_PER_S / ((double)n_pilots * k_step * scs_hz);
  return cell_m / (2.0 * max_speed_ms);
}

/* Criterion 3, Doppler migration. An acceleration "a" moves the target by 2*a*t/lambda
Hz over the window, and that has to stay inside one Doppler cell, 1/t:
2*a*t/lambda < 1/t, so t^2 < lambda / 2a and t < sqrt(lambda / 2a) */
double nr_ue_sensing_span_accel(double carrier_hz)
{
  if (!(carrier_hz > 0.0) || !(NR_SENSING_TARGET_MAX_ACCEL_MS2 > 0.0))
    return HUGE_VAL;
  const double lambda = C_M_PER_S / carrier_hz;
  return sqrt(lambda / (2.0 * NR_SENSING_TARGET_MAX_ACCEL_MS2));
}

/* Criterion 4, the lower bound. dv = lambda / (2 * t_span), so reaching a target dv
needs t_span >= lambda / (2 * dv). */
double nr_ue_sensing_span_min(double carrier_hz)
{
  if (!(carrier_hz > 0.0) || !(NR_SENSING_TARGET_DV_MS > 0.0))
    return 0.0;
  return (C_M_PER_S / carrier_hz) / (2.0 * NR_SENSING_TARGET_DV_MS);
}

nr_sensing_span_t nr_ue_sensing_span_bounds(const nr_sensing_history_t *hist,
                                            double max_speed_ms,
                                            int n_pilots,
                                            int k_step)
{
  /* Criteria 2 (range migration) and 3 (acceleration) are left out of the window:
  they only cost a target actually moving that fast a few dB, smeared over two cells,
  while bounding the window by them coarsens every map for every target. They remain
  callable on their own, to report what a window means for fast targets. So the
  window has one upper bound, the grid budget, which is computational. */
  (void)n_pilots;
  (void)k_step;
  nr_sensing_span_t s;
  s.t_max_s = nr_ue_sensing_span_grid(nr_ue_sensing_axis_f_max(hist, max_speed_ms));
  s.t_max_by = NR_SENSING_SPAN_GRID;
  s.t_min_s = nr_ue_sensing_span_min((double)hist->carrier_hz);
  s.feasible = s.t_min_s <= s.t_max_s;
  return s;
}

/* The window a map is actually held to. The short end of the feasible range, because
the clutter residue grows with t_span faster than the lines it feeds concentrate; see
the OBSERVATION WINDOW block. When the requirements contradict each other the upper
bound wins, since undersampling the Doppler peaks or smearing the target across range
cells corrupts the map, while missing the resolution target only makes it coarser. */
double nr_ue_sensing_window_s(const nr_sensing_history_t *hist,
                                     double max_speed_ms,
                                     int n_pilots,
                                     int k_step)
{
  const nr_sensing_span_t s = nr_ue_sensing_span_bounds(hist, max_speed_ms, n_pilots, k_step);
  return s.feasible ? s.t_min_s : s.t_max_s;
}

/* The most recent snapshots of one reference signal, oldest first, as indices into
the ring. Shared by the transform and the snapshot dump so the two can never
disagree about which samples a map was built from. Returns how many were found.

Two bounds, and whichever is reached first stops the walk:

  n_snap_max   how many snapshots, the caller's --sensing-symbols. Caps the cost of
               the transform, which is O(n_freq * n_bins * n_snap).
  max_span_s   how long a window they may cover, from nr_ue_sensing_span_bounds().
               This is the one that decides the resolution, the Doppler grid density
               and how much clutter residue ends up on the map, so it cannot be left
               to whatever the scheduler delivered. Zero or less disables it.

Time is measured back from the newest snapshot of this stream, which is the first
match the walk finds: it runs from hist->head backwards, so pushes are visited newest
first. Push order is only nearly time order (see the sort below), so a snapshot a few
samples out of sequence can end the walk a hair early; at hundreds of snapshots per
window that costs nothing, and cutting stale entries is the whole point. */
static int nr_ue_sensing_gather(const nr_sensing_history_t *hist,
                                const nr_sensing_stream_t *stream,
                                int n_snap_max,
                                double max_span_s,
                                int idx[])
{
  const double fs = (double)hist->ofdm_symbol_size * hist->scs_hz;
  /* The bound as a sample count, so the test below is integer and no snapshot time
  is ever converted to double. fs is the OFDM sample clock t_sample is counted on. */
  const uint64_t span_samples =
      (max_span_s > 0.0 && fs > 0.0) ? (uint64_t)(max_span_s * fs) : UINT64_MAX;
  uint64_t t_newest = 0;
  bool have_newest = false;

  int n = 0;
  for (int k = 1; k <= hist->count && n < n_snap_max; k++) {
    const int i = (hist->head - k + hist->depth) % hist->depth;
    if (!nr_ue_sensing_stream_eq(&hist->ring[i].stream, stream))
      continue;

    const uint64_t t = hist->ring[i].t_sample;
    if (!have_newest) {
      t_newest = t;
      have_newest = true;
    } else if (t < t_newest && t_newest - t > span_samples) {
      /* Older than the window. Guarded on t < t_newest so a snapshot pushed out of
      sequence, which would make the unsigned difference enormous, is kept rather
      than read as ancient. */
      break;
    }
    idx[n++] = i;
  }
  // reverse into push order, oldest first
  for (int a = 0, b = n - 1; a < b; a++, b--) {
    const int tmp = idx[a];
    idx[a] = idx[b];
    idx[b] = tmp;
  }
  /* Then into time order. Push order is not time order: the DL actors finish their slots in
  any order (--num-dl-actors > 1), so their pushes interleave, while every user of idx[]
  assumes t_sample increases and computes times as t_sample - t_first without a sign. The
  input is nearly sorted, so an insertion sort is close to linear; equal times keep their
  push order. */
  for (int a = 1; a < n; a++) {
    const int cur = idx[a];
    const uint64_t t_cur = hist->ring[cur].t_sample;
    int b = a - 1;
    while (b >= 0 && hist->ring[idx[b]].t_sample > t_cur) {
      idx[b + 1] = idx[b];
      b--;
    }
    idx[b + 1] = cur;
  }
  return n;
}

void nr_ue_sensing_dump_map(const char *path,
                            int frame,
                            int slot,
                            int aarx,
                            const nr_sensing_map_t *map,
                            const nr_sensing_marker_t *markers,
                            int n_markers,
                            const nr_sensing_aoa_t *aoa,
                            int n_aoa)
{
  if (path == NULL) return;

  static pthread_mutex_t map_lock = PTHREAD_MUTEX_INITIALIZER;
  // TODO: THIS WILL BE REMOVED WHEN WE KNOW THE NEW KERNEL L PATH FITTING IS WORTH
  /* Files already started in this run, so each one is truncated and gets its header
  once. Several paths exist when --sensing-clutter-compare writes its .L1.csv. */
  static const char *started[4];
  static int n_started = 0;

  pthread_mutex_lock(&map_lock);
  // TODO: THIS WILL BE REMOVED WHEN WE KNOW THE NEW KERNEL L PATH FITTING IS WORTH
  bool map_started = false;
  for (int i = 0; i < n_started && !map_started; i++)
    map_started = strcmp(started[i], path) == 0;
  FILE *f = fopen(path, map_started ? "a" : "w");
  if (f == NULL) {
    pthread_mutex_unlock(&map_lock);
    LOG_E(NR_PHY, "sensing: cannot open map file %s\n", path);
    return;
  }

  if (!map_started) {
    fprintf(f,
            "# frame,slot,aarx,ports,layer,k_step,k_offset,n_snapshots,t_span_s,m_per_bin,f_max_hz,carrier_hz,"
            "n_bins,n_freq,"
            "n_marker,n_aoa,bin_los,n_pilots_min,n_pilots_max,n_positions,noise_ref,"
            "power[n_bins*n_freq],marker[n_marker*4: range_m,speed_ms,snr_dB,verdict],"
            "aoa[n_aoa*5: range_m,speed_ms,angle_deg,power_dB,quality_dB]\n");
    if (n_started < (int)(sizeof(started) / sizeof(started[0])))
      started[n_started++] = strdup(path);
  }

  if (markers == NULL) n_markers = 0;
  if (aoa == NULL) n_aoa = 0;

  /* noise_ref closes the metadata block, behind bin_los, which closed it before. Every
  leading offset a reader of an older dump expects is then unchanged and only the
  payload start moves, which is the same reason the marker and AoA blocks were appended
  rather than inserted.

  %.9g rather than a fixed precision: this is the one field with no natural scale. It
  is the level the map was divided by, so it carries whatever the channel and the
  receive chain happened to be doing, and on the outdoor dumps that moved over 25 dB.
  Enough digits that dividing it back out returns the transform's own numbers. */
  fprintf(f,
          "%d,%d,%d,%d,%d,%d,%d,%d,%.9f,%.6f,%.6f,%.0f,%d,%d,%d,%d,%.4f,%d,%d,%d,%.9g",
          frame,
          slot,
          aarx,
          map->stream.ports,
          map->stream.layer,
          map->stream.k_step,
          map->stream.k_offset,
          map->n_snapshots,
          map->t_span_s,
          map->m_per_bin,
          map->f_max_hz,
          map->carrier_hz,
          map->n_bins,
          map->n_freq,
          n_markers,
          n_aoa,
          map->bin_los,
          map->n_pilots_min,
          map->n_pilots_max,
          map->n_positions,
          map->noise_ref);

  for (int i = 0; i < map->n_bins * map->n_freq; i++)
    fprintf(f, ",%.6g", map->power[i]);
  /* After the payload rather than before it, so the leading fields keep the offsets
  a reader of an older dump expects and the two layouts stay distinguishable by
  length alone. The AoA block follows the markers for the same reason: a reader that
  knows nothing of it still finds the power and the markers where it expects them. */
  for (int i = 0; i < n_markers; i++)
    fprintf(f, ",%.3f,%.4f,%.2f,%d", markers[i].range_m, markers[i].speed_ms, markers[i].snr_dB, markers[i].verdict);
  for (int i = 0; i < n_aoa; i++)
    fprintf(f,
            ",%.3f,%.4f,%.2f,%.2f,%.2f",
            aoa[i].range_m,
            aoa[i].speed_ms,
            aoa[i].angle_deg,
            aoa[i].power_dB,
            aoa[i].quality_dB);
  fprintf(f, "\n");
  fclose(f);
  pthread_mutex_unlock(&map_lock);
}

int nr_ue_sensing_since_last_map(const nr_sensing_history_t *hist, const nr_sensing_stream_t *stream)
{
  /* Read only: a stream that has never been pushed has nothing to count, and
  claiming a slot for it here would spend the table on queries. The lock is taken
  on a const history because the counter table is written by concurrent pushes. */
  nr_sensing_history_t *h = (nr_sensing_history_t *)hist;
  pthread_mutex_lock(&h->lock);
  const int *count = nr_ue_sensing_counter(h, stream, false);
  const int value = count ? *count : 0;
  pthread_mutex_unlock(&h->lock);
  return value;
}

void nr_ue_sensing_clear_since_last_map(nr_sensing_history_t *hist, const nr_sensing_stream_t *stream)
{
  pthread_mutex_lock(&hist->lock);
  int *count = nr_ue_sensing_counter(hist, stream, false);
  if (count)
    *count = 0;
  pthread_mutex_unlock(&hist->lock);
}

bool nr_ue_sensing_history_take(nr_sensing_history_t *src,
                                const nr_sensing_stream_t *stream,
                                int n_layers,
                                int min_count,
                                nr_sensing_history_t *dst)
{
  // allocated before taking the lock, so the pushes are not held up by malloc
  nr_sensing_snapshot_t *ring = malloc_or_fail((size_t)src->depth * sizeof(*ring));

  pthread_mutex_lock(&src->lock);
  if (min_count > 0) {
    const int *count = nr_ue_sensing_counter(src, stream, false);
    if (count == NULL || *count < min_count) {
      // another actor took this map since the caller looked at the count
      pthread_mutex_unlock(&src->lock);
      free(ring);
      return false;
    }
  }
  /* Clear the counters of the layers this map is built from: all of them for a
  layer-averaged map, else the stream's own. OAI cleared layers 0..n_layers-1 in
  both cases, so a layer 1 map cleared layer 0's counter and never its own: layer 1
  then mapped on every slot and layer 0 never again after its first map. */
  for (int j = 0; j < n_layers; j++) {
    nr_sensing_stream_t c = *stream;
    c.layer = n_layers > 1 ? (uint8_t)j : stream->layer;
    int *count = nr_ue_sensing_counter(src, &c, false);
    if (count)
      *count = 0;
  }
  *dst = *src;
  memcpy(ring, src->ring, (size_t)src->depth * sizeof(*ring));
  pthread_mutex_unlock(&src->lock);

  dst->ring = ring;
  // the copy must not inherit the source's mutex state; it is only read by the map task
  pthread_mutex_init(&dst->lock, NULL);
  return true;
}

/// Lattice index of a snapshot's first pilot, the a_m of the kernel
static inline int nr_ue_sensing_lattice_start(const nr_sensing_snapshot_t *s)
{
  return (s->k_first - s->k_first % s->stream.k_step) / s->stream.k_step;
}

/* The distinct grants of a window, "(n=1638 a=0 k0=6)x12, ...", one entry per
   (n_m, a_m) pair with the snapshots it covers and the k_first it came from.

   Only ever called when the kernel table overflows, so the O(n^2) dedup is not worth
   avoiding. Truncated with a trailing "..." rather than grown: the first few pairs
   already say what the scheduler is doing. */
static void nr_ue_sensing_format_grants(const nr_sensing_history_t *hist, const int *idx, int n, char *buf, size_t len)
{
  size_t off = 0;
  buf[0] = '\0';

  for (int i = 0; i < n; i++) {
    const nr_sensing_snapshot_t *s = &hist->ring[idx[i]];
    const int a_m = nr_ue_sensing_lattice_start(s);

    // already listed by an earlier snapshot
    bool seen = false;
    for (int j = 0; j < i && !seen; j++) {
      const nr_sensing_snapshot_t *p = &hist->ring[idx[j]];
      seen = p->n_pilots == s->n_pilots && nr_ue_sensing_lattice_start(p) == a_m;
    }
    if (seen)
      continue;

    int cnt = 0;
    for (int j = i; j < n; j++) {
      const nr_sensing_snapshot_t *p = &hist->ring[idx[j]];
      if (p->n_pilots == s->n_pilots && nr_ue_sensing_lattice_start(p) == a_m)
        cnt++;
    }

    const int w =
        snprintf(buf + off, len - off, "%s(n=%d a=%d k0=%d)x%d", off ? ", " : "", s->n_pilots, a_m, s->k_first, cnt);
    if (w < 0 || (size_t)w >= len - off) {
      if (len >= 4)
        memcpy(buf + len - 4, "...", 4);
      return;
    }
    off += w;
  }
}

void nr_ue_sensing_clutter_kernel(int n_pilots, int a_m, double u0, int idft_size, int n_bins, float *kre, float *kim)
{
  float w[NR_SENSING_MAX_IDFT];

  /* The same taper the data went through, from the same generator, so the model
  cannot drift from it. The fallback matches nr_ue_sensing_apply_hann(): a lattice too
  short to window was left rectangular there too. */
  if (!nr_ue_sensing_hann_coeffs(n_pilots, w))
    for (int l = 0; l < n_pilots; l++)
      w[l] = 1.0f;

  for (int b = 0; b < n_bins; b++) {
    const double u = b - u0;

    /* K(u) = sum_l w[l] exp(j*2*pi*(a_m + l)*u/N), i.e. the transform of (4) applied
    to the taper placed at a_m. Computed here rather than through
    nr_ue_sensing_delay_response() because that path is fixed point: the taper alone
    is far louder than the channel it normally multiplies, so it would saturate the
    IDFT, and the saturation is the one kind of error a common scale factor cannot
    absorb. The transform is the same, in double. */

    /* Phase of the first pilot, reduced to one turn. a_m*u/N reaches a couple of
    hundred turns at the far end of the axis, where cos() would spend most of the
    mantissa on the integer part of the argument. */
    const double turns = (double)a_m * u / idft_size;
    const double ph0 = 2.0 * M_PI * (turns - trunc(turns));
    double cr = cos(ph0);
    double ci = sin(ph0);

    /* From one pilot to the next the phase advances by a fixed step, so the sum is a
    rotating phasor rather than n_pilots trig calls. |u| <= n_bins << N keeps the step
    small, and in double the recurrence drifts far below the float it is stored in. */
    const double dph = 2.0 * M_PI * u / idft_size;
    const double dr = cos(dph);
    const double di = sin(dph);

    double sr = 0.0;
    double si = 0.0;
    for (int l = 0; l < n_pilots; l++) {
      sr += w[l] * cr;
      si += w[l] * ci;
      const double nr = cr * dr - ci * di;
      ci = cr * di + ci * dr;
      cr = nr;
    }

    kre[b] = (float)sr;
    kim[b] = (float)si;
  }
}

/* Sub-bin offset of the peak at bin b, by a parabola through it and its two neighbours.

   Fitted on the profile in dB, not in power. The pilots are Hann tapered, so the
   mainlobe of the delay response is very nearly a parabola on a log scale and only
   roughly one on a linear scale; interpolating the power biases the vertex back
   towards the bin centre, which is the very error being removed.

   Returns 0 at the edges and on a plateau, where the vertex is not defined, and when
   the vertex lands more than half a bin away, which a true local maximum cannot
   produce. */
static double nr_ue_sensing_peak_frac(const double *p, int b, int n_bins)
{
  if (b <= 0 || b >= n_bins - 1 || p[b - 1] <= 0.0 || p[b] <= 0.0 || p[b + 1] <= 0.0)
    return 0.0;
  const double ym = 10.0 * log10(p[b - 1]);
  const double y0 = 10.0 * log10(p[b]);
  const double yp = 10.0 * log10(p[b + 1]);
  const double den = ym - 2.0 * y0 + yp;
  if (den >= 0.0)
    return 0.0;
  const double d = 0.5 * (ym - yp) / den;
  return (d > -0.5 && d < 0.5) ? d : 0.0;
}

/* Fit one static path at the fractional bin u and subtract it from the residual.

   The path contributes alpha * K_m(b - u) to snapshot m, with K_m the kernel of that
   snapshot's grant. alpha is the least squares fit over the whole window,

       alpha = sum_m sum_b r_m[b] conj(K_m(b - u)) / sum_m sum_b |K_m(b - u)|^2 .

   IMPORTANT FOR PAPER: 
   One alpha and not one per snapshot: that is what keeps this a clutter filter. K_m
   carries the grant variation while alpha stays constant, so a mover at the same range
   survives, where a per-snapshot amplitude would delete everything at that delay
   whatever its Doppler.

   The residual is raw, before the pilot-count (n_max) gain: K_m already carries the height
   that gain corrects, so applying it first would correct the same thing twice.

   Only bins within NR_CLUTTER_KERNEL_HALF_SPAN of u are fitted and touched.

   res_re, res_im : [n][n_bins] residual, row major, updated in place */
static void nr_ue_sensing_remove_path(const nr_sensing_history_t *hist,
                                      const int *idx,
                                      int n,
                                      int n_bins,
                                      double u,
                                      float *res_re,
                                      float *res_im)
{
  int b_lo = (int)floor(u) - NR_CLUTTER_KERNEL_HALF_SPAN;
  int b_hi = (int)ceil(u) + NR_CLUTTER_KERNEL_HALF_SPAN;
  if (b_lo < 0)
    b_lo = 0;
  if (b_hi > n_bins - 1)
    b_hi = n_bins - 1;
  const int span = b_hi - b_lo + 1;
  if (span <= 0)
    return;

  /* One kernel per distinct (n_m, a_m). The kernel depends on b - u only, so the one
  of bins b_lo..b_hi is the kernel centred at u - b_lo, over span bins. */
  float *pool[NR_CLUTTER_MAX_KERNELS][2] = {{NULL}};
  int pool_pilots[NR_CLUTTER_MAX_KERNELS];
  int pool_start[NR_CLUTTER_MAX_KERNELS];
  int pool_n = 0;
  const float *kr[NR_SENSING_HISTORY_DEPTH] = {NULL};
  const float *ki[NR_SENSING_HISTORY_DEPTH] = {NULL};

  for (int i = 0; i < n; i++) {
    const nr_sensing_snapshot_t *s = &hist->ring[idx[i]];
    const int a_m = nr_ue_sensing_lattice_start(s);
    int j = 0;
    while (j < pool_n && (pool_pilots[j] != s->n_pilots || pool_start[j] != a_m))
      j++;

    if (j == pool_n) {
      if (pool_n == NR_CLUTTER_MAX_KERNELS) {
        char grants[512];
        nr_ue_sensing_format_grants(hist, idx, n, grants, sizeof(grants));
        LOG_W(NR_PHY,
              "sensing: more than %d distinct grants over %d snapshots, clutter kernel skipped from snapshot %d on. "
              "Grants: %s\n",
              NR_CLUTTER_MAX_KERNELS,
              n,
              i,
              grants);
        break;
      }
      pool[j][0] = malloc((size_t)span * sizeof(float));
      pool[j][1] = malloc((size_t)span * sizeof(float));
      if (pool[j][0] == NULL || pool[j][1] == NULL) {
        LOG_E(NR_PHY, "sensing: cannot allocate a clutter kernel, path at bin %.2f left in\n", u);
        break;
      }
      nr_ue_sensing_clutter_kernel(s->n_pilots, a_m, u - b_lo, s->idft_size, span, pool[j][0], pool[j][1]);
      pool_pilots[j] = s->n_pilots;
      pool_start[j] = a_m;
      pool_n++;
    }
    kr[i] = pool[j][0];
    ki[i] = pool[j][1];
  }

  // alpha = sum r conj(K) / sum |K|^2
  double num_r = 0.0;
  double num_i = 0.0;
  double den = 0.0;
  for (int i = 0; i < n; i++) {
    if (kr[i] == NULL)
      continue;
    const size_t row = (size_t)i * n_bins + b_lo;
    for (int j = 0; j < span; j++) {
      const double rr = res_re[row + j];
      const double ri = res_im[row + j];
      num_r += rr * kr[i][j] + ri * ki[i][j];
      num_i += ri * kr[i][j] - rr * ki[i][j];
      den += (double)kr[i][j] * kr[i][j] + (double)ki[i][j] * ki[i][j];
    }
  }

  if (den > 0.0) {
    const double a_re = num_r / den;
    const double a_im = num_i / den;
    // r -= alpha * K
    for (int i = 0; i < n; i++) {
      if (kr[i] == NULL)
        continue;
      const size_t row = (size_t)i * n_bins + b_lo;
      for (int j = 0; j < span; j++) {
        res_re[row + j] -= (float)(a_re * kr[i][j] - a_im * ki[i][j]);
        res_im[row + j] -= (float)(a_re * ki[i][j] + a_im * kr[i][j]);
      }
    }
  }

  for (int j = 0; j < NR_CLUTTER_MAX_KERNELS; j++) {
    free(pool[j][0]);
    free(pool[j][1]);
  }
}

/*
We first create the time axis that we want to use for the basis.
q is the output of the grahm-schmidt procedure (which is useful for the dot product and no matrix inverse).
Look that we normalize the time axis [-1,1] for numerical stability and center it.
Shared with the MUSIC stage, which has to project its steering vectors on the same subspace.
*/
int nr_ue_sensing_slow_basis(const double *t, int n, int degree, double q[][NR_SENSING_HISTORY_DEPTH])
{
  if (degree < 0 || n <= 0)
    return 0;
  if (degree > NR_CLUTTER_SLOW_TREND_MAX)
    degree = NR_CLUTTER_SLOW_TREND_MAX;

  int n_q = 0;
  double t_mid = 0.0;
  for (int i = 0; i < n; i++)
    t_mid += t[i];
  t_mid /= n;
  const double t_half = (t[n - 1] > t[0]) ? 0.5 * (t[n - 1] - t[0]) : 1.0;
  for (int k = 0; k <= degree; k++) {
    double *v = q[n_q];
    for (int i = 0; i < n; i++) {
      const double tau = (t[i] - t_mid) / t_half;
      double p = 1.0;
      for (int e = 0; e < k; e++)
        p *= tau;
      v[i] = p;
    }
    for (int j = 0; j < n_q; j++) {
      double d = 0.0;
      for (int i = 0; i < n; i++)
        d += q[j][i] * v[i];
      for (int i = 0; i < n; i++)
        v[i] -= d * q[j][i];
    }
    double norm = 0.0;
    for (int i = 0; i < n; i++)
      norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm < 1e-9)
      continue; // not independent of the lower degrees on these sample times
    for (int i = 0; i < n; i++)
      v[i] /= norm;
    n_q++;
  }
  return n_q;
}

/* Orthonormal basis of the comb tones on the sample times t[0..n-1], for
NR_CLUTTER_COMB_REMOVE.

The tones are exp(j*2*pi*k*f0*t) for k = -K..K without k = 0, f0 the fundamental. They
are made orthonormal by Gram-Schmidt, and first made orthogonal to q[0..n_q-1], the
polynomial basis the slow trend already removed: the two stages then take disjoint
subspaces and the second cannot put back what the first took. On a window holding many
cycles of f0 they are nearly orthogonal to the polynomials anyway, so this mostly
matters at the short end.

A tone that collapses under the orthogonalisation is dropped rather than normalised:
that happens when the sample times cannot tell it apart from one already in the basis,
and keeping it would amplify noise along a direction the data does not resolve.

e_re, e_im receive the basis, one vector of n entries per row, and the number of rows
written is returned. Room for 2*NR_COMB_MAX_HARMONIC rows is needed. */
static int nr_ue_sensing_comb_basis(const double *t,
                                    int n,
                                    double f0_hz,
                                    const double q[][NR_SENSING_HISTORY_DEPTH],
                                    int n_q,
                                    double *e_re,
                                    double *e_im)
{
  if (n <= 0 || !(f0_hz > 0.0))
    return 0;

  int n_e = 0;
  for (int k = 1; k <= NR_COMB_MAX_HARMONIC; k++) {
    for (int sgn = 1; sgn >= -1; sgn -= 2) {
      double *vr = &e_re[(size_t)n_e * n];
      double *vi = &e_im[(size_t)n_e * n];
      const double f = sgn * k * f0_hz;
      for (int i = 0; i < n; i++) {
        /* Same reduction as the Doppler transform: f*t reaches many turns over the
        window and cos() spends its mantissa on the integer part otherwise. */
        const double turns = f * t[i];
        const double ph = 2.0 * M_PI * (turns - trunc(turns));
        vr[i] = cos(ph);
        vi[i] = sin(ph);
      }

      // orthogonalise against the real polynomial basis, then against the tones already kept
      for (int j = 0; j < n_q; j++) {
        double cr = 0.0, ci = 0.0;
        for (int i = 0; i < n; i++) {
          cr += q[j][i] * vr[i];
          ci += q[j][i] * vi[i];
        }
        for (int i = 0; i < n; i++) {
          vr[i] -= cr * q[j][i];
          vi[i] -= ci * q[j][i];
        }
      }
      for (int j = 0; j < n_e; j++) {
        const double *ur = &e_re[(size_t)j * n];
        const double *ui = &e_im[(size_t)j * n];
        double cr = 0.0, ci = 0.0;
        for (int i = 0; i < n; i++) { // <u, v> with u conjugated
          cr += ur[i] * vr[i] + ui[i] * vi[i];
          ci += ur[i] * vi[i] - ui[i] * vr[i];
        }
        for (int i = 0; i < n; i++) {
          vr[i] -= cr * ur[i] - ci * ui[i];
          vi[i] -= cr * ui[i] + ci * ur[i];
        }
      }

      double nrm2 = 0.0;
      for (int i = 0; i < n; i++)
        nrm2 += vr[i] * vr[i] + vi[i] * vi[i];
      /* Against n, the norm a tone starts with. Below a thousandth of it the direction
      is not resolved by these sample times. */
      if (nrm2 <= 1e-3 * (double)n)
        continue;
      const double inv = 1.0 / sqrt(nrm2);
      for (int i = 0; i < n; i++) {
        vr[i] *= inv;
        vi[i] *= inv;
      }
      n_e++;
    }
  }
  return n_e;
}

static int nr_ue_sensing_cmp_double(const void *a, const void *b)
{
  const double x = *(const double *)a;
  const double y = *(const double *)b;
  return (x > y) - (x < y);
}

int nr_ue_sensing_range_doppler(const nr_sensing_history_t *hist,
                                const nr_sensing_stream_t *stream,
                                int n_snap_max,
                                double max_speed_ms,
                                nr_sensing_clutter_t clutter_mode,
                                int max_paths,
                                nr_sensing_map_t *map,
                                nr_sensing_slowtime_t *slow_out)
{
  if (hist->ring == NULL) return 0;

  /* The Doppler axis does not depend on the snapshots, so it is known before they are
  gathered, and it has to be: its half width sets criterion 1 of the window, which is
  an input to the gather rather than something discovered afterwards. */
  const double f_max = nr_ue_sensing_axis_f_max(hist, max_speed_ms);

  // The window, bounded above by the grid budget
  const double max_span_s = nr_ue_sensing_window_s(hist, max_speed_ms, 0, 0);

  int idx[NR_SENSING_HISTORY_DEPTH];
  int n = nr_ue_sensing_gather(hist, stream, n_snap_max, max_span_s, idx);
#if NR_SENSING_GROUP_BY_GRANT
  // Testing switch, see nr_ue_sensing_align.h: one grant per map, so one point spread function.
  n = nr_ue_sensing_group_by_grant(hist, idx, n);
#endif

  // We require a mininum of symbols
  if (n < NR_SENSING_MIN_SNAPSHOTS) return 0;

  // widest grant gathered, for the log below
  int n_pilots_widest = hist->ring[idx[0]].n_pilots;
  for (int i = 1; i < n; i++)
    if (hist->ring[idx[i]].n_pilots > n_pilots_widest)
      n_pilots_widest = hist->ring[idx[i]].n_pilots;

  const nr_sensing_snapshot_t *first = &hist->ring[idx[0]];

  // bandwidth 
  const double fs = (double)hist->ofdm_symbol_size * hist->scs_hz;

  // meter per bin, for the range axis
  const double m_per_bin = C_M_PER_S / ((double)first->idft_size * first->stream.k_step * hist->scs_hz);

  // Number of bins for the range
  int n_bins = first->n_bins;

  // time between the first and last symbol taken into account for the map
  const double t_span = (double)(hist->ring[idx[n - 1]].t_sample - first->t_sample) / fs;
  if (t_span <= 0.0) return 0;

  /*
  DOPPLER AXIS

  With TDD, the reference symbols only arrive during the DL part of each TDD period.
  Because of that, the Doppler spectrum repeats every 1/T_TDD: a target at f also
  shows up at f + k/T_TDD. With a 5 ms period that is every 200 Hz, so every 8.8 m/s at 3.41 GHz.

  The axis goes as far as the fastest target we want to see (max_speed_ms), plus one
  TDD period of margin. Why the margin: the TDD detector tells a real peak from a copy by looking where the
  copies of a candidate would be, one period away on each side. Without the margin a
  target at the very edge of the axis could not be checked.

  So every target appears several times on a map, 200 Hz apart. With
  --sensing-tdd-detect the detector tries to say which one is real; without it, all the copies
  are there, evenly spaced.

  A target faster than the axis still shows up, folded back inside it at a wrong
  speed, so set max_speed_ms (in m/s) above the fastest thing expected.

  As discussed in our notes, the goal is also to limit t_span such that oversampling
  is always possible, and having a coarser grid than the peaks's target themselves is
  impossible (see notes). The limit of t_span is be computed wrt to the maximum
  speed we want, and set a priori by us. Also, NR_SENSING_MAP_MAX_BINS_FREQ can be
  updated by us if we think we can computationally support more points.
  The axis itself, f_speed and f_max, is built at the top of this function: the window
  the gather is allowed to cover is derived from f_max, so it has to exist before the
  snapshots do. Only the point count, which needs t_span, is left here.

  t_span does not limit how far the axis goes. It sets the resolution: two targets
  closer than 1/t_span Hz in Doppler cannot be separated. It does limit how many
  points that axis can be drawn with, which is why nr_ue_sensing_gather() bounds it;
  see the OBSERVATION WINDOW block in nr_ue_map.h.
  */

  /* The bounds the window was held to (grid budget and resolution target), reported
  with each map. */
  const nr_sensing_span_t sp = nr_ue_sensing_span_bounds(hist, max_speed_ms, n_pilots_widest, first->stream.k_step);
  static const char *const span_by[] = {"grid budget", "range migration", "target acceleration"};
  const double lambda_m = C_M_PER_S / (double)hist->carrier_hz;
  const double dv_got = lambda_m / (2.0 * t_span);
  const double dv_err_pct = 100.0 * (dv_got / NR_SENSING_TARGET_DV_MS - 1.0);

  if (!sp.feasible)
    LOG_W(NR_PHY,
          "sensing: requirements contradict. %.1f m/s of coverage allows at most %.3f s (%s), but %.2f m/s of "
          "resolution needs at least %.3f s. This map %.3f s, %.2f m/s, %+.0f%%. "
          "Lower --sensing-max-speed or raise NR_SENSING_TARGET_DV_MS.\n",
          max_speed_ms,
          sp.t_max_s,
          span_by[sp.t_max_by],
          NR_SENSING_TARGET_DV_MS,
          sp.t_min_s,
          t_span,
          dv_got,
          dv_err_pct);
  else if (dv_got > NR_SENSING_TARGET_DV_MS * (1.0 + NR_SENSING_TARGET_DV_SLACK))
    /* The bounds are satisfiable, the data ran out. Tested on the resolution rather
    than on t_span against t_min, because the window aims at t_min and lands just under
    it on every healthy map; see NR_SENSING_TARGET_DV_SLACK. Which cap ran out first is
    the actionable part: the two the code owns, or the scheduler, which no constant
    here can fix. */
    LOG_W(NR_PHY,
          "sensing: window %.3f s short of the %.3f s that %.2f m/s needs, limited by %s. "
          "Resolution %.2f m/s, %+.0f%%. %d snapshots at %.0f/s, %d pilots.\n",
          t_span,
          sp.t_min_s,
          NR_SENSING_TARGET_DV_MS,
          (n >= n_snap_max)    ? "--sensing-symbols"
          : (n >= hist->depth) ? "NR_SENSING_HISTORY_DEPTH"
                               : "the scheduler, too few snapshots offered",
          dv_got,
          dv_err_pct,
          n,
          n / t_span,
          n_pilots_widest);
  else
    LOG_I(NR_PHY,
          "sensing: window %.3f s in [%.3f, %.3f] (upper: %s), %d snapshots at %.0f/s, %d pilots. "
          "Resolution %.2f m/s (%+.0f%%) over +-%.1f m/s, %d Doppler points at %.2f per cell.\n",
          t_span,
          sp.t_min_s,
          sp.t_max_s,
          span_by[sp.t_max_by],
          n,
          n / t_span,
          n_pilots_widest,
          dv_got,
          dv_err_pct,
          f_max * lambda_m / 2.0,
          (int)(4.0 * f_max * t_span) | 1,
          ((int)(4.0 * f_max * t_span) | 1) / (2.0 * f_max * t_span));

  /* Number of Doppler points: 2 per resolution cell (1/t_span). That is the same
  density as the TDD detector's grid, so the same CFAR settings work on both. Odd, so
  that 0 Hz falls exactly on a point. */
  int n_freq = (int)(4.0 * f_max * t_span) | 1;
  if (n_freq < 3)
    n_freq = 3;
  if (n_freq > NR_SENSING_MAP_MAX_BINS_FREQ) {
    /* Should be unreachable now that the window is bounded */
    LOG_W(NR_PHY,
          "sensing: Doppler grid clamped, %d points wanted for f_max %.0f Hz over %.3f s but only %d available. "
          "The grid is %.2f points per resolution cell instead of 2 and peak heights are understated. "
          "Lower --sensing-max-speed or the window bound (%.3f s), or raise NR_SENSING_MAP_MAX_BINS_FREQ.\n",
          n_freq,
          f_max,
          t_span,
          NR_SENSING_MAP_MAX_BINS_FREQ - 1,
          (double)((NR_SENSING_MAP_MAX_BINS_FREQ - 1) | 1) / (2.0 * f_max * t_span),
          max_span_s);
    n_freq = (NR_SENSING_MAP_MAX_BINS_FREQ - 1) | 1;
  }

  map->n_bins = n_bins;
  map->n_freq = n_freq;
  map->n_snapshots = n;
  map->stream = first->stream;
  map->m_per_bin = m_per_bin;
  map->f_max_hz = f_max;
  map->t_span_s = t_span;
  map->carrier_hz = (double)hist->carrier_hz;

  double t[NR_SENSING_HISTORY_DEPTH];
  for (int i = 0; i < n; i++) t[i] = (double)(hist->ring[idx[i]].t_sample - first->t_sample) / fs;

  /* AMPLITUDE NORMALISATION ACROSS SNAPSHOTS

  The issue is the following: across grants, the number of resource allocated for
  pilots can change. You get the full band of pilots slot 0, and slot 1 you can get only
  half of the band, the starting pilot changes.

  What changes is the amplitude of the IDFT. Since its a sum, if you sum less positive 
  elements (exponentials), the sum will compute a small value, even if a target is there.

  The Doppler transform reads a varying sequence as motion and spreads the
  target off zero Doppler, and the clutter removal below, which subtracts a single
  mean, cannot cancel something whose amplitude keeps moving. A full band slot has
  1638 pilots (constant on our config) here and the narrowest the band check admits has 
  822.

  Scaling every snapshot to a common pilot count removes it. The reference is the 
  widest grant in the window, so the full band slots are left untouched and only 
  the narrow ones are brought up.

  This corrects the height of a peak, not its width: fewer pilots also make the
  mainlobe wider. Both are arguments for feeding the history only
  snapshots of one fixed allocation; this is what keeps the map usable until then.
  
  */ 
  /* Over the snapshots that survived the migration trim, unlike n_pilots_widest above,
  which was taken before it: gain[] divides by this, so it has to be a grant the window
  really holds. n_pilots_max <= n_pilots_widest always. */
  double gain[NR_SENSING_HISTORY_DEPTH];
  int n_pilots_min = hist->ring[idx[0]].n_pilots;
  int n_pilots_max = hist->ring[idx[0]].n_pilots;
  for (int i = 1; i < n; i++) {
    const int np = hist->ring[idx[i]].n_pilots;
    if (np < n_pilots_min) n_pilots_min = np;
    if (np > n_pilots_max) n_pilots_max = np;
  }

  if (n_pilots_min <= 0) {
    LOG_E(NR_PHY, "sensing: snapshot with %d pilots in the history, cannot normalise\n", n_pilots_min);
    return 0;
  }

  // The gain factor
  for (int i = 0; i < n; i++)
    gain[i] = (double)n_pilots_max / (double)hist->ring[idx[i]].n_pilots;

  // for printing some statistics in the logs
  int n_positions = 0;
  for (int i = 0; i < n; i++) {
    const int a_m = nr_ue_sensing_lattice_start(&hist->ring[idx[i]]);
    bool seen = false;
    for (int j = 0; j < i && !seen; j++)
      seen = nr_ue_sensing_lattice_start(&hist->ring[idx[j]]) == a_m;
    if (!seen)
      n_positions++;
  }
  map->n_pilots_min = n_pilots_min;
  map->n_pilots_max = n_pilots_max;
  map->n_positions = n_positions;

  if (n_pilots_max != n_pilots_min)
    LOG_W(NR_PHY,
          "sensing: grant width varies over the map window, %d to %d pilots. "
          "Heights are normalised but the delay resolution still moves: %.2f to %.2f m.\n",
          n_pilots_min,
          n_pilots_max,
          C_M_PER_S / ((double)n_pilots_max * first->stream.k_step * hist->scs_hz),
          C_M_PER_S / ((double)n_pilots_min * first->stream.k_step * hist->scs_hz));

  /* 
  CLUTTER REMOVAL : 2 modes

  1) Subtracting each bin's mean over the window. This assumes that the clutter contribution
  is constant over the symbol m. This can be assumed if the grant size of DMRS does not change.

  Clutter at range bin r is then:
    h[r,m]= c[r] + t[r,m] + n[r,m]
  So we can remove the slow time mean per range bin:
    h[r]= (1/M) m=0 sum M−1 with h[r,m]
  and then substract this to the range profiles

  Now recall that having clutter removal implies introducing copies of the
  spectrum of the mean computed into the final 2d range-doppler map. This is because we
  do not have uniform sampling for the slow-time axis.

  --------

  2) Substracting something depending on the symbol m. This assumes the clutter contribution 
  is not constant over the symbol m, especially when dealing with different grant sizes.

  Recall that the complex gain that comes with the doppler tap (the exponential) in the range profile is
    alpha_l exp(..k_first..) A_n_m(b - u_l) exp(j 2 pi c_m (b - u_l) / N)
  
  where c_m = a_m + (n_m - 1)/2 and u_l = tau_l N p delta_f
  Recall that p is the k_step, N the idft size, a_m the lattice index of the first pilot in symbol m,
  and n_m the number of pilots measured in symbol m.

  The first 2 terms are constant. So a static path is alpha * K_m(b - u_l), with one
  alpha for the window and K_m known from the grant of snapshot m. The delays u_l are
  not known: the direct path is taken at bin_los, and further static paths are found
  one at a time as the strongest peak of the slow-time mean of what is left, until
  that peak is no longer static or no longer above the noise. See NR_CLUTTER_KERNEL
  in nr_ue_map.h, then the mean of mode 1 runs on the residual.
  */
 
  /*
  Store the bin LOS. We know this bin is the LOS, and we can use it for the localization
  The bin the kernel is centred on. The direct path is static and the shortest path,
  so it is the earliest strong arrival of the window's energy profile; in a lab it is
  also tens of dB above every reflection, outdoors it need not be
  */
  double e_prof[NR_SENSING_MAP_MAX_BINS_RANGE];
  int u0 = 0;
  double best = -1.0;
  for (int b = 0; b < n_bins; b++) {
    double e = 0.0;
    for (int i = 0; i < n; i++) {
      const cf_t *h = hist->ring[idx[i]].h;
      e += (double)crealf(h[b]) * crealf(h[b]) + (double)cimagf(h[b]) * cimagf(h[b]);
    }
    e_prof[b] = e;
    if (e > best) {
      best = e;
      u0 = b;
    }
  }
  /* The earliest local peak within NR_SENSING_LOS_FIRST_DB of the strongest, see
  there: the direct path arrives first, it does not have to be the strongest. */
  const double e_min = best * pow(10.0, -NR_SENSING_LOS_FIRST_DB / 10.0);
  for (int b = 0; b < u0; b++) {
    const bool peak = (b == 0 || e_prof[b] >= e_prof[b - 1]) && e_prof[b] >= e_prof[b + 1];
    if (peak && e_prof[b] >= e_min) {
      u0 = b;
      break;
    }
  }

  /* Sub-bin refinement of the direct path. Worth doing rather than reporting the
  integer bin: at 2.44 m/bin the rounding alone is up to 1.2 m of bias, and it lands
  directly on D = L + dR in nr_ue_target_position(). */
  map->bin_los = (double)u0 + nr_ue_sensing_peak_frac(e_prof, u0, n_bins);
  LOG_I(NR_PHY, "LOS bin to check its stability: %.2f (peak bin %d)\n", map->bin_los, u0);

  /* The residual everything below works on: the raw snapshots, from which the kernel
  mode subtracts the static paths it fits. Kept before the gain, see
  nr_ue_sensing_remove_path(). */
  float *res_re = malloc((size_t)n * n_bins * sizeof(float));
  float *res_im = malloc((size_t)n * n_bins * sizeof(float));
  if (res_re == NULL || res_im == NULL) {
    LOG_E(NR_PHY, "sensing: cannot allocate the residual of %d snapshots\n", n);
    free(res_re);
    free(res_im);
    return 0;
  }
  /* DIRECT PATH NORMALISATION, see NR_CLUTTER_LOS_NORM.

  One complex number per snapshot, c[i], and every bin of that snapshot is multiplied
  by it while the residual is filled below.

  Measuring the direct path: not h[u0] on its own. bin_los is fractional, so a single
  integer bin samples the point spread function off its peak, and one bin carries one
  bin's worth of noise and of whatever else leaks into it. Instead the snapshot is
  matched filtered against the kernel of its own grant, centred on the fractional
  bin_los:

      alpha_i = sum_b h_i[b] conj(K_i(b - u0))  /  sum_b |K_i(b - u0)|^2

  which is the least squares amplitude of the direct path in that snapshot. This is
  the same expression nr_ue_sensing_remove_path() evaluates, with one difference that
  it sums over snapshots as well and gets a single alpha for the
  window, and a single alpha is exactly what cannot represent a drifting gain.
  Evaluating it per snapshot instead is the measurement of that drift.

  Correcting with it: c_i = ref * conj(alpha_i) / |alpha_i|^2, with ref the median of
  |alpha|. Applied to the direct path itself this gives ref exactly, real and positive,
  in every snapshot: the amplitude is levelled to a common value and the phase is
  zeroed. ref is a median rather than a mean so that a blocked snapshot does not drag
  the level everything else is scaled to.

  Amplitude and phase are both taken raw, without smoothing. Smoothing the amplitude
  would keep the estimate's own noise out of the other bins, but the direct path is
  20 to 30 dB above everything and the matched filter runs over a hundred bins, so that
  noise is some 40 dB down and not worth the machinery. Smoothing the phase would be
  wrong in any case: the fast jumps, a precoder change or a timing loop correction, are
  what this is here to remove. */
  double c_re[NR_SENSING_HISTORY_DEPTH];
  double c_im[NR_SENSING_HISTORY_DEPTH];
  for (int i = 0; i < n; i++) {
    c_re[i] = 1.0;
    c_im[i] = 0.0;
  }

#if NR_CLUTTER_LOS_NORM
  {
    /* Same span as a path removal, so the filter sees the same skirt the kernel mode
    would subtract, clipped to the map. */
    int b_lo = (int)floor(map->bin_los) - NR_CLUTTER_KERNEL_HALF_SPAN;
    int b_hi = (int)ceil(map->bin_los) + NR_CLUTTER_KERNEL_HALF_SPAN;
    if (b_lo < 0) b_lo = 0;
    if (b_hi > n_bins - 1) b_hi = n_bins - 1;
    const int span = b_hi - b_lo + 1;

    double a_re[NR_SENSING_HISTORY_DEPTH];
    double a_im[NR_SENSING_HISTORY_DEPTH];
    double a_abs[NR_SENSING_HISTORY_DEPTH];
    float *kre = malloc((size_t)span * sizeof(float));
    float *kim = malloc((size_t)span * sizeof(float));

    if (span > 0 && kre != NULL && kim != NULL) {
      for (int i = 0; i < n; i++) {
        const nr_sensing_snapshot_t *s = &hist->ring[idx[i]];
        /* The kernel of this snapshot's own grant. Rebuilt per snapshot rather than
        pooled by (n_pilots, a_m) as remove_path does: it costs O(span) and the
        transform that follows is orders of magnitude more. */
        nr_ue_sensing_clutter_kernel(s->n_pilots,
                                     nr_ue_sensing_lattice_start(s),
                                     map->bin_los - b_lo,
                                     s->idft_size,
                                     span,
                                     kre,
                                     kim);
        double num_r = 0.0, num_i = 0.0, den = 0.0;
        for (int j = 0; j < span; j++) {
          const double hr = crealf(s->h[b_lo + j]);
          const double hi = cimagf(s->h[b_lo + j]);
          num_r += hr * kre[j] + hi * kim[j];
          num_i += hi * kre[j] - hr * kim[j];
          den += (double)kre[j] * kre[j] + (double)kim[j] * kim[j];
        }
        a_re[i] = (den > 0.0) ? num_r / den : 0.0;
        a_im[i] = (den > 0.0) ? num_i / den : 0.0;
        a_abs[i] = sqrt(a_re[i] * a_re[i] + a_im[i] * a_im[i]);
      }

      double sorted[NR_SENSING_HISTORY_DEPTH];
      memcpy(sorted, a_abs, (size_t)n * sizeof(double));
      qsort(sorted, n, sizeof(sorted[0]), nr_ue_sensing_cmp_double);
      const double ref = sorted[n / 2];

      const double corr_max = pow(10.0, NR_CLUTTER_LOS_MAX_CORR_DB / 20.0);
      int n_clamped = 0;
      if (ref > 0.0) {
        for (int i = 0; i < n; i++) {
          if (a_abs[i] <= 0.0) {
            n_clamped++;
            continue; // no direct path at all, leave the snapshot as it is
          }
          /* c = ref * conj(alpha) / |alpha|^2, then the magnitude clamped. The phase
          is kept whatever the clamp does to the magnitude: rotating a snapshot back
          onto the common reference is always right, it is only scaling it up that
          becomes dangerous when there is nothing left to scale. */
          double scale = ref / a_abs[i];
          if (scale > corr_max) {
            scale = corr_max;
            n_clamped++;
          } else if (scale < 1.0 / corr_max) {
            scale = 1.0 / corr_max;
            n_clamped++;
          }
          const double u = scale / a_abs[i]; // folds the 1/|alpha| of the phase removal
          c_re[i] = a_re[i] * u;
          c_im[i] = -a_im[i] * u;
        }
      }
      if (n_clamped > 0)
        LOG_D(NR_PHY,
              "sensing: direct path normalisation clamped on %d of %d snapshots, the direct path was lost there\n",
              n_clamped,
              n);
    }
    free(kre);
    free(kim);
  }
#endif

  for (int i = 0; i < n; i++) {
    const cf_t *h = hist->ring[idx[i]].h;
    for (int b = 0; b < n_bins; b++) {
      // complex multiply by c[i]: levels the amplitude and zeroes the phase of the LOS
      res_re[(size_t)i * n_bins + b] = (float)((double)crealf(h[b]) * c_re[i] - (double)cimagf(h[b]) * c_im[i]);
      res_im[(size_t)i * n_bins + b] = (float)((double)crealf(h[b]) * c_im[i] + (double)cimagf(h[b]) * c_re[i]);
    }
  }

  if (clutter_mode == NR_CLUTTER_KERNEL) {
    if (max_paths < 1)
      max_paths = 1;
    if (max_paths > NR_CLUTTER_MAX_PATHS)
      max_paths = NR_CLUTTER_MAX_PATHS;
    double u_fit[NR_CLUTTER_MAX_PATHS];
    int n_fit = 0;

    // Path 0: the direct path, always, at the refined fractional bin.
    nr_ue_sensing_remove_path(hist, idx, n, n_bins, map->bin_los, res_re, res_im);
    u_fit[n_fit++] = map->bin_los;

    /* Further static paths, strongest first, while the strongest peak left in the
    slow-time mean is both static and above the noise; see NR_CLUTTER_KERNEL. */
    while (n_fit < max_paths) {
      double coh[NR_SENSING_MAP_MAX_BINS_RANGE]; // |slow-time mean|^2
      double inc[NR_SENSING_MAP_MAX_BINS_RANGE]; // slow-time mean of |r|^2
      for (int b = 0; b < n_bins; b++) {
        double sr = 0.0;
        double si = 0.0;
        double sp = 0.0;
        for (int i = 0; i < n; i++) {
          // with the gain, so a varying grant height does not read as non-static
          const double x = gain[i] * res_re[(size_t)i * n_bins + b];
          const double y = gain[i] * res_im[(size_t)i * n_bins + b];
          sr += x;
          si += y;
          sp += x * x + y * y;
        }
        coh[b] = (sr * sr + si * si) / ((double)n * n);
        inc[b] = sp / n;
      }

      // strongest peak of the mean, away from the paths already fitted
      int b_best = -1;
      for (int b = 0; b < n_bins; b++) {
        bool near = false;
        for (int l = 0; l < n_fit && !near; l++)
          near = fabs(b - u_fit[l]) < NR_CLUTTER_MIN_SEP_BINS;
        if (!near && (b_best < 0 || coh[b] > coh[b_best]))
          b_best = b;
      }
      if (b_best < 0)
        break;

      /* Noise power per snapshot from the median bin, which the few paths do not move.
      Averaging n snapshots divides it by n in the mean. */
      double sorted[NR_SENSING_MAP_MAX_BINS_RANGE];
      memcpy(sorted, inc, n_bins * sizeof(double));
      qsort(sorted, n_bins, sizeof(double), nr_ue_sensing_cmp_double);
      const double noise_mean = sorted[n_bins / 2] / n;

      const double s_static = inc[b_best] > 0.0 ? coh[b_best] / inc[b_best] : 0.0;
      if (s_static < NR_CLUTTER_STATIC_MIN || coh[b_best] < NR_CLUTTER_SNR_MIN * noise_mean)
        break;

      const double u = b_best + nr_ue_sensing_peak_frac(coh, b_best, n_bins);
      nr_ue_sensing_remove_path(hist, idx, n, n_bins, u, res_re, res_im);
      u_fit[n_fit++] = u;
      LOG_D(NR_PHY,
            "sensing: static path %d removed at bin %.2f, S %.2f, %.1f dB above noise\n",
            n_fit - 1,
            u,
            s_static,
            10.0 * log10(coh[b_best] / noise_mean));
    }
  }

  /* Stage 2, and the whole of the mean mode: remove what varies slowly in time from each
  bin, after the gain. See NR_CLUTTER_SLOW_TREND_DEGREE; degree 0 is the plain mean. In
  the kernel mode this runs on what the fits left, so it is no longer dominated by the
  strong paths and takes the weak static scatterers the loop did not reach.

  The slow functions 1, tau, tau^2, ... of the snapshot times (tau = t centred and scaled
  to about [-1, 1], for conditioning) are made orthonormal by Gram-Schmidt into q_0..q_K,
  and each bin's slow-time sequence z loses its projection on them:
      z <- z - sum_k q_k (q_k^T z)
  For degree 0, q_0 = 1/sqrt(n) and this is exactly z - mean(z).

  The residual holds the result, gain applied, so everything below reads it as is. */

  const int trend_degree = (clutter_mode == NR_CLUTTER_NONE) ? -1 : NR_CLUTTER_SLOW_TREND_DEGREE;
  double q[NR_CLUTTER_SLOW_TREND_MAX + 1][NR_SENSING_HISTORY_DEPTH];
  const int n_q = nr_ue_sensing_slow_basis(t, n, trend_degree, q);

  /*
  For each range bin, we remove now the slow part. Apply first the gain, and then we compute
  z <- z - q(q^T z)
  */
  double zr[NR_SENSING_HISTORY_DEPTH];
  double zi[NR_SENSING_HISTORY_DEPTH];
  for (int b = 0; b < n_bins; b++) {
    // normalised first: a trend fitted to unnormalised snapshots carries the very
    // amplitude modulation it is supposed to cancel
    for (int i = 0; i < n; i++) {
      zr[i] = gain[i] * res_re[(size_t)i * n_bins + b];
      zi[i] = gain[i] * res_im[(size_t)i * n_bins + b];
    }
    for (int j = 0; j < n_q; j++) {
      double cr = 0.0;
      double ci = 0.0;
      for (int i = 0; i < n; i++) {
        cr += q[j][i] * zr[i];
        ci += q[j][i] * zi[i];
      }
      for (int i = 0; i < n; i++) {
        zr[i] -= cr * q[j][i];
        zi[i] -= ci * q[j][i];
      }
    }
    for (int i = 0; i < n; i++) {
      res_re[(size_t)i * n_bins + b] = (float)zr[i];
      res_im[(size_t)i * n_bins + b] = (float)zi[i];
    }
  }

  /* COMB REMOVAL, rank one across range. See NR_CLUTTER_COMB_REMOVE.

  Three steps:

    1. P[b][j] = <e_j, z_b>, the amount of tone j in range bin b. Subtracting all of
       P would be the per bin projection, which deletes a target that happens to sit
       on a line.
    2. Keep only the rank one part of P, sigma * u v^H, found by power iteration. It
       is the one modulation shape common to every range bin, scaled by u[b], which is
       what the clutter is and what a single bin's target is not.
    3. Subtract that part alone.

  P is n_bins by 2K, a few thousand entries, so the power iteration on it costs
  nothing next to the transform that follows. */
#if NR_CLUTTER_COMB_REMOVE
  if (clutter_mode != NR_CLUTTER_NONE && n_q > 0) {
    const double f0 = 1.0 / (NR_COMB_TDD_MULTIPLE * nr_ue_sensing_tdd_period_s(hist));
    const int n_e_max = 2 * NR_COMB_MAX_HARMONIC;

    double *e_re = malloc((size_t)n_e_max * n * sizeof(double));
    double *e_im = malloc((size_t)n_e_max * n * sizeof(double));
    double *p_re = malloc((size_t)n_bins * n_e_max * sizeof(double));
    double *p_im = malloc((size_t)n_bins * n_e_max * sizeof(double));

    if (e_re == NULL || e_im == NULL || p_re == NULL || p_im == NULL) {
      LOG_E(NR_PHY, "sensing: cannot allocate the comb basis, comb left in\n");
    } else {
      const int n_e = nr_ue_sensing_comb_basis(t, n, f0, (const double(*)[NR_SENSING_HISTORY_DEPTH])q, n_q, e_re, e_im);

      if (n_e > 0) {
        // 1. tone content of every range bin
        for (int b = 0; b < n_bins; b++) {
          for (int j = 0; j < n_e; j++) {
            const double *er = &e_re[(size_t)j * n];
            const double *ei = &e_im[(size_t)j * n];
            double cr = 0.0, ci = 0.0;
            for (int i = 0; i < n; i++) {
              const double xr = res_re[(size_t)i * n_bins + b];
              const double xi = res_im[(size_t)i * n_bins + b];
              cr += er[i] * xr + ei[i] * xi; // conj(e) . z
              ci += er[i] * xi - ei[i] * xr;
            }
            p_re[(size_t)b * n_e + j] = cr;
            p_im[(size_t)b * n_e + j] = ci;
          }
        }

        /* 2. dominant singular pair of P by power iteration, alternating
        u <- P v, v <- P^H u with a normalisation each time. v starts flat so no
        harmonic is favoured. */
        double vr[2 * NR_COMB_MAX_HARMONIC], vi[2 * NR_COMB_MAX_HARMONIC];
        double ur[NR_SENSING_MAP_MAX_BINS_RANGE], ui[NR_SENSING_MAP_MAX_BINS_RANGE];
        for (int j = 0; j < n_e; j++) {
          vr[j] = 1.0 / sqrt((double)n_e);
          vi[j] = 0.0;
        }
        double sigma = 0.0;
        for (int it = 0; it < NR_COMB_POWER_ITERS; it++) {
          double nu = 0.0;
          for (int b = 0; b < n_bins; b++) {
            double sr = 0.0, si = 0.0;
            for (int j = 0; j < n_e; j++) {
              const double pr = p_re[(size_t)b * n_e + j];
              const double pi = p_im[(size_t)b * n_e + j];
              sr += pr * vr[j] - pi * vi[j];
              si += pr * vi[j] + pi * vr[j];
            }
            ur[b] = sr;
            ui[b] = si;
            nu += sr * sr + si * si;
          }
          if (nu <= 0.0) break;
          nu = 1.0 / sqrt(nu);
          for (int b = 0; b < n_bins; b++) { ur[b] *= nu; ui[b] *= nu; }

          double nv = 0.0;
          for (int j = 0; j < n_e; j++) {
            double sr = 0.0, si = 0.0;
            for (int b = 0; b < n_bins; b++) { // P^H u, so P conjugated
              const double pr = p_re[(size_t)b * n_e + j];
              const double pi = p_im[(size_t)b * n_e + j];
              sr += pr * ur[b] + pi * ui[b];
              si += pr * ui[b] - pi * ur[b];
            }
            vr[j] = sr;
            vi[j] = si;
            nv += sr * sr + si * si;
          }
          if (nv <= 0.0) break;
          sigma = sqrt(nv);
          nv = 1.0 / sigma;
          for (int j = 0; j < n_e; j++) { vr[j] *= nv; vi[j] *= nv; }
        }

        // 3. subtract sigma * u[b] * conj(v[j]) along e_j, and nothing else
        if (sigma > 0.0) {
          double removed = 0.0, total = 0.0;
          for (int b = 0; b < n_bins; b++) {
            for (int j = 0; j < n_e; j++) {
              // c = sigma * u[b] * conj(v[j])
              const double cr = sigma * (ur[b] * vr[j] + ui[b] * vi[j]);
              const double ci = sigma * (ui[b] * vr[j] - ur[b] * vi[j]);
              const double *er = &e_re[(size_t)j * n];
              const double *ei = &e_im[(size_t)j * n];
              for (int i = 0; i < n; i++) {
                res_re[(size_t)i * n_bins + b] -= (float)(cr * er[i] - ci * ei[i]);
                res_im[(size_t)i * n_bins + b] -= (float)(cr * ei[i] + ci * er[i]);
              }
              removed += cr * cr + ci * ci;
              total += p_re[(size_t)b * n_e + j] * p_re[(size_t)b * n_e + j]
                       + p_im[(size_t)b * n_e + j] * p_im[(size_t)b * n_e + j];
            }
          }
          /* How much of the tone content was common to every range bin. Near 1 says
          the comb really is one modulation and the rank one model fits; well below it
          says several independent sources and a higher rank would be needed. */
          LOG_D(NR_PHY,
                "sensing: comb removed at %.1f Hz and %d harmonics, %d basis vectors, rank-1 share %.2f\n",
                f0,
                NR_COMB_MAX_HARMONIC,
                n_e,
                (total > 0.0) ? removed / total : 0.0);
        }
      }
    }
    free(e_re);
    free(e_im);
    free(p_re);
    free(p_im);
  }
#endif

  /*
  Hand the conditioned samples back before they are consumed. This is the last
  point at which they exist as such: the loop below keeps only |C|^2. Placed after the
  slow trend is removed, which itself ran on normalised snapshots, because a trend fitted
  to unnormalised ones would carry the amplitude modulation the gain is there to remove.
  */
  if (slow_out != NULL) {
    slow_out->n_snap = n;
    slow_out->n_bins = n_bins;
    slow_out->t0_sample = first->t_sample;
    for (int i = 0; i < n; i++) {
      slow_out->t_s[i] = t[i];
      for (int b = 0; b < n_bins; b++) {
        const size_t k = (size_t)i * n_bins + b;
        slow_out->h_re[k] = res_re[k];
        slow_out->h_im[k] = res_im[k];
      }
    }
    /* The widest grant in the window, which is what gain[] normalised everything to.
    When the grants differ this is the best single answer rather than a right one:
    each snapshot was scaled to this pilot count but still carries the mainlobe shape
    of its own, so a caller modelling the response assumes a width most snapshots do
    not have. One more reason to feed the history a fixed allocation. 
    
    This is the n_max thing where we rescale, see math doc.
    */
    slow_out->win_len = n_pilots_max;
    slow_out->win_start = (first->k_first - first->k_first % first->stream.k_step) / first->stream.k_step;
    slow_out->idft_size = first->idft_size;
    slow_out->m_per_bin = m_per_bin;
    slow_out->trend_degree = trend_degree;
  }

  // Computes the actual DFT for the doppler estimation.
  // Because the arrival times of the symbols are not uniform,
  // we cannot use FFT, but we have to compute the DFT directly.

  /* Normalisation of the transform: divide the sum by the coherent gain of the slow
  time window, so power is the squared amplitude of a tone rather than of its sum.

  Without it the sum of n terms makes the map scale as n^2, and n is not a constant:
  it is whatever the scheduler delivered inside the window. Two maps of the same scene
  taken seconds apart then differ in level by (n1/n2)^2 for no physical reason, which
  is 6 dB for a factor of two in snapshot count. Nothing inside one map notices, the
  CFAR being a ratio against its own training cells, but everything that compares maps
  over time does: the dumped power, a tracker, and the eye going from one map to the
  next.

  The divisor is the coherent gain sum(w_i), not n, so that a slow time taper can be
  introduced later without touching this: for the rectangular window in use sum(w_i)
  is exactly n. A coherent target then holds its height whatever n is, and the noise
  floor falls as 1/n, which is where the integration gain belongs and where it stays
  visible. */
  double w_sum = 0.0;
  for (int i = 0; i < n; i++)
    w_sum += 1.0; // rectangular slow-time window; becomes sum(w_i) once tapered
  const double dft_norm = (w_sum > 0.0) ? 1.0 / (w_sum * w_sum) : 1.0;

  const double df = (n_freq > 1) ? (2.0 * f_max) / (n_freq - 1) : 0.0;
  double wr[NR_SENSING_HISTORY_DEPTH];
  double wi[NR_SENSING_HISTORY_DEPTH];
  for (int f = 0; f < n_freq; f++) {
    const double freq = -f_max + f * df;
    for (int i = 0; i < n; i++) {
      /* Reduced modulo one turn before the trig call. freq * t[i] reaches a few
      hundred turns at the edge of the grid for a window of a hundred milliseconds,
      and evaluating cos() that far out spends most of the mantissa on the integer
      part of the phase. Taking the fractional turn first keeps the argument inside
      [-2pi, 2pi], where the result is accurate to full double precision. */
      const double turns = freq * t[i];
      const double ph = -2.0 * M_PI * (turns - trunc(turns));
      wr[i] = cos(ph);
      wi[i] = sin(ph);
    }
    for (int b = 0; b < n_bins; b++) {
      double ar = 0.0;
      double ai = 0.0;
      for (int i = 0; i < n; i++) {
        const double hr = res_re[(size_t)i * n_bins + b];
        const double hi = res_im[(size_t)i * n_bins + b];
        ar += hr * wr[i] - hi * wi[i];
        ai += hr * wi[i] + hi * wr[i];
      }
      map->power[b * n_freq + f] = (float)((ar * ar + ai * ai) * dft_norm);
    }
  }

  /* Scale the finished map to its own noise floor, so maps can be compared with each
  other; see noise_ref in nr_ue_map.h.

  The floor is estimated as the median cell. Noise holds the large majority of the
  cells, so the median lands in it, and the things that are not noise are exactly the
  things a mean would be dragged by: the direct path occupies a couple of range bins
  out of n_bins, and the replica comb some 5% of the Doppler axis. Neither moves a
  median.

  Taken on a stride rather than on every cell. A map is up to n_bins * n_freq cells and
  sorting all of them would cost a noticeable fraction of the transform itself, while
  the median of a few thousand of them is already far more accurate than this needs to
  be: the quantity being estimated is a noise level that only has to be right to a
  fraction of a dB. The stride is odd so it does not lock onto the Doppler axis and
  sample the same few frequencies in every range bin. */
  const int n_cells = n_bins * n_freq;
  double samp[NR_SENSING_MAP_NORM_SAMPLES];
  int n_samp = 0;
  /* Rounded up, so the stride always carries the walk across the whole map. Rounding
  down instead gives a stride of 1 whenever n_cells is just over the sample budget, and
  the walk then stops at the budget having seen only the first cells, which are the low
  range bins where the direct path sits: the "noise floor" comes out of the clutter.
  Measured that way on the short window dumps it was up to 5.8 dB out. Forced odd so it
  cannot divide n_freq and sample the same few Doppler columns in every range bin. */
  int stride = (n_cells + NR_SENSING_MAP_NORM_SAMPLES - 1) / NR_SENSING_MAP_NORM_SAMPLES;
  if (stride < 1)
    stride = 1;
  stride |= 1;
  for (int c = 0; c < n_cells && n_samp < NR_SENSING_MAP_NORM_SAMPLES; c += stride)
    samp[n_samp++] = map->power[c];

  qsort(samp, n_samp, sizeof(samp[0]), nr_ue_sensing_cmp_double);
  const double ref = (n_samp > 0) ? samp[n_samp / 2] : 0.0;

  if (ref > 0.0) {
    map->noise_ref = ref;
    const float inv = (float)(1.0 / ref);
    for (int c = 0; c < n_cells; c++)
      map->power[c] *= inv;
  } else {
    /* An all zero or degenerate map, which the clutter stage can produce if every
    snapshot was identical. Left unscaled and flagged with 1.0 so a reader can tell
    the difference between "not normalised" and "normalised by 1". */
    map->noise_ref = 0.0;
    LOG_W(NR_PHY, "sensing: map noise floor estimated at %g, left unnormalised\n", ref);
  }

  free(res_re);
  free(res_im);

  return n;
}


void nr_ue_sensing_dump_snapshots(const char *path,
                                  const nr_sensing_history_t *hist,
                                  const nr_sensing_stream_t *stream,
                                  int n_snap_max,
                                  double max_speed_ms)
{
  if (path == NULL || hist->ring == NULL) return;

  int idx[NR_SENSING_HISTORY_DEPTH];
  /* Same window as the map, derived the same way from the same max_speed_ms, so the
  file still holds exactly the snapshots the map was built from. */
  /* Criteria 1 and 3 only, as the map's first pass does. The migration trim the map
  applies afterwards needs the gathered grant, so a dump can hold a few snapshots the
  map then dropped; the file carries n_pilots per snapshot, which is what a reader
  needs to reproduce the trim. */
  int n = nr_ue_sensing_gather(hist, stream, n_snap_max, nr_ue_sensing_window_s(hist, max_speed_ms, 0, 0), idx);
#if NR_SENSING_GROUP_BY_GRANT
  // same selection as the map, so a dump shows exactly the snapshots its map used
  n = nr_ue_sensing_group_by_grant(hist, idx, n);
#endif
  if (n < 1) return;

  FILE *f = fopen(path, "w");
  if (f == NULL) {
    LOG_E(NR_PHY, "sensing: cannot open snapshot file %s\n", path);
    return;
  }

  const double fs = (double)hist->ofdm_symbol_size * hist->scs_hz;
  const uint64_t t0 = hist->ring[idx[0]].t_sample;

  fprintf(f,
          "# snap,t_sample,t_rel_s,k_step,k_first,n_pilots,idft_size,n_bins,scs_hz,fs_hz,carrier_hz,"
          "re[0],im[0],...,re[n_bins-1],im[n_bins-1]\n");

  for (int i = 0; i < n; i++) {
    const nr_sensing_snapshot_t *s = &hist->ring[idx[i]];
    fprintf(f,
            "%d,%lu,%.9f,%d,%d,%d,%d,%d,%d,%.0f,%.0f",
            i,
            (unsigned long)s->t_sample,
            (double)(s->t_sample - t0) / fs,
            s->stream.k_step,
            s->k_first,
            s->n_pilots,
            s->idft_size,
            s->n_bins,
            hist->scs_hz,
            fs,
            (double)hist->carrier_hz);
    for (int b = 0; b < s->n_bins; b++)
      /* %g rather than the %d this printed in Q15: the snapshots are complex
      float now, and an integer conversion would truncate every value below one
      to zero, which is most of a normalised delay response. */
      fprintf(f, ",%g,%g", crealf(s->h[b]), cimagf(s->h[b]));
    fprintf(f, "\n");
  }

  fclose(f);
  LOG_I(NR_PHY,
        "sensing: dumped %d raw snapshots (ports 0x%x layer %d k_step %d) to %s\n",
        n,
        stream->ports,
        stream->layer,
        stream->k_step,
        path);
}
