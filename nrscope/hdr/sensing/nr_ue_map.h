/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * \brief UE sensing pipeline, cross-slot stage: accumulating the per-symbol delay
 *        responses and transforming them into a range-Doppler map.
 *
 * Split from nr_ue_sensing.h because the two stages run on different time scales.
 * The delay response of a symbol is produced inside the receive path; the Doppler
 * axis needs a hundred milliseconds of them, so it is accumulated here and
 * transformed on a worker thread.
 */

#ifndef __NR_UE_MAP__H__
#define __NR_UE_MAP__H__

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include "nrscope/hdr/sensing/sensing_defs.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"

/* Speed of light in m/s */
#define C_M_PER_S 299792458.0

// REQUIREMENTS FOR SENSING PROPERTIES
// Remember we assume TDD, 30 KHz SCS and 3.45 GHz

#define NR_SENSING_MIN_SNAPSHOTS  32 
#define NR_SENSING_TARGET_DV_MS   0.3  // sets the minimum window t_span^min

/* How far past NR_SENSING_TARGET_DV_MS a map may land before it is worth a warning.
The window is aimed at t_span^min itself, and snapshots arrive where the scheduler put
them, so the last one before the bound leaves the achieved span a little short and the
resolution a little coarse on every healthy map. Without slack the warning would fire
at a 0% miss, every map, and stop meaning anything. */
#define NR_SENSING_TARGET_DV_SLACK 0.10 // 10%
#define NR_SENSING_TARGET_MAX_ACCEL_MS2 2.0 // maximum acceleration of targets in [m/s^2]

// REQUIREMENTS FOR THE SENSING MAP

/*
At 30 kHz, 3.45 GHz its NR_SENSING_MAP_MAX_BINS_RANGE x 2.44 m/bin in maps.
Two bins in the delay domain correspond to a distance
Because with pilots every p subcarriers, the effective grid spacing feeding the IDFT is Δf_p = p·Δf. 
A length-L IDFT of a grid with spacing Δf_p produces delay samples spaced
by Δτ = 1 / (L · p · Δf) => so Δd = c·Δτ = c / (L · p · Δf)
*/
#define NR_SENSING_MAP_MAX_BINS_RANGE 50

/*
Make sure the Doppler axis is fine enough to resolve the targets.
The maximum doppler f_max is 2v_max/lambda + 1 / TDD. The grid spans [-f_max,f_max]
so 2f_max. We want 2 grid doppler points per cell, so n_freq = 2 * 2f_max * T_span points.
So the bin spacing in m/s/bin is (lambda/2) * Delta_f = (lambda/2) 2f_max / (n_freq-1) =~ (lambda/2) 1 / 2T = lambda / 4T = delta_v/2.
So with delta_v = 0.15m/s it gives 0.075 m/s per bin. This is just because we are putting 2 grid points per resolution cell !! else it would have been just 0.15m/s/bin.
So the bin spacing for doppler depends on the time window T. 
To convert this in Hz per bin we can just use the relation f = 2v/lambda, so delta_f = 2 delta_v / lambda = 2 * 0.15 / 0.087 = 3.45 Hz/bin.

NR_SENSING_MAP_MAX_BINS_FREQ only serves when 4f_max*T_span exceeds this value. If it exceeds it
then the spacing is coarser than 2 points per cell.
It defines the size of the power map for the doppler axis.
So its a computational cost
*/
#define NR_SENSING_MAP_MAX_BINS_FREQ 2048

/* Cells sampled to estimate a map's noise floor, see noise_ref. A median over this
many is accurate to well under a dB, which is all the estimate has to be, and it keeps
the sort off the critical path: the whole map can be 128 * 1024 cells. */
#define NR_SENSING_MAP_NORM_SAMPLES 4096
#define NR_SENSING_HISTORY_DEPTH 2048 // maximum number of symbols in history


/* OBSERVATION WINDOW

Until this existed the window was bounded by a snapshot count alone
(--sensing-symbols), so its duration was whatever the scheduler happened to deliver:
on outdoor dumps it wandered between 0.41 s and 0.64 s from one map to the next, and
reached 22 s once when the ring still held stale entries. Every quantity that matters
is a function of the duration, so it has to be an input.

Nothing here is a free constant. The window is derived from the three requirements at
the top of this file and from the axis the caller asked for; see
nr_ue_sensing_span_bounds().

UPPER: only criterion 1 is applied. Criteria 2 and 3 are physical but soft: they
cost a target that really moves that fast a few dB, smeared over two cells, while
bounding the window by them coarsened every map for every target, slow ones
included. Their functions stay, to report what a window means for fast targets.

  1. Grid density, nr_ue_sensing_span_grid(). The Doppler axis needs 4 * f_max *
     t_span points to hold two per resolution cell and NR_SENSING_MAP_MAX_BINS_FREQ
     caps that. Past the cap the grid is coarser than the peaks it samples and a
     target landing between two points loses up to 13 dB, with nothing in the map to
     show for it. Computational, and the only upper bound applied.

  2. Range migration, nr_ue_sensing_span_migration(). The target must stay inside one
     range cell: t < cell_m / (2 v), as the delay axis is path length and the path of
     a target at v changes at up to 2 v. cell_m is c / (n_pilots * k_step * scs), the
     width of a peak, not m_per_bin, which is only how finely that peak is sampled.
     Above about 1.5 m/s this is tighter than the grid.

  3. Doppler migration, nr_ue_sensing_span_accel(). While accelerating, the target
     must stay inside one Doppler cell: t < sqrt(lambda / 2a), 0.148 s at 3.41 GHz
     and 2 m/s^2. Needs an assumption
     about the scene, NR_SENSING_TARGET_MAX_ACCEL_MS2, and is the only bound here that
     is not derived from the radio.

LOWER, one criterion:

  4. Velocity resolution, nr_ue_sensing_span_min(). lambda / (2 * t_span) must reach
     NR_SENSING_TARGET_DV_MS, so t_span >= lambda / (2 * dv).

The two ends can cross only through criterion 1: 4 * f_max * t_min above the grid's
1023 points, e.g. 0.1 m/s (0.44 s) with more than about 17 m/s of coverage. Then
nr_ue_sensing_span_bounds() says so rather than silently picking one end.

For reference, were criteria 2 and 3 applied, at 3.41 GHz on the full carrier 0.15
m/s of resolution (0.147 s) would be impossible above about 10 m/s of coverage, and
0.1 m/s (0.22 s) at any coverage.

NOT a bound, deliberately: clutter residue. Measured on the dumps the replica comb
grows 2.45 dB per dB of t_span, against the 1.0 that concentrating a fixed amount of
energy into narrower lines would give, so the excess is the residue itself growing
with the window: 19.9 dB over the floor at 0.64 s against 3.0 dB at 0.12 s. That is a
gradient, not a constraint. It argues for sitting near the short end of whatever range
comes out below, which is why the window returned is the lower bound and not the
upper. The price of the short end is coherent gain, 10*log10 of the snapshots dropped,
and a blind zone around zero Doppler that scales as 1/t_span; averaging consecutive
maps buys most of the gain back non-coherently. */
typedef struct {
  /// shortest window that meets NR_SENSING_TARGET_DV_MS, seconds
  double t_min_s;
  /// longest window all three upper criteria allow, seconds
  double t_max_s;
  /// the one that produced t_max_s, for the log line
  enum { NR_SENSING_SPAN_GRID, NR_SENSING_SPAN_MIGRATION, NR_SENSING_SPAN_ACCEL } t_max_by;
  /// false when t_min_s > t_max_s: the requirements cannot all be met at once
  bool feasible;
} nr_sensing_span_t;

/* One OFDM symbol's delay response

The Doppler axis is built across slots, not within one: a slot spans 0.5 ms, which
gives ~86 m/s of Doppler resolution however densely it is sampled. Only the total
observation span matters, so snapshots are accumulated over many slots.

t_sample is an absolute count on the OFDM sample clock, not a slot index, because
the reference symbols are not uniformly spaced: 3 DMRS symbols land irregularly
inside a slot, PDSCH is only scheduled in some slots, and TRS arrives every 20 ms.
Keeping the true instant lets the Doppler transform work on the real sampling
pattern instead of pretending it is uniform. */
typedef struct {
  /// absolute time of this OFDM symbol, in samples of the OFDM sample clock
  uint64_t t_sample;
  /// which measurement process produced it; snapshots of different streams never mix
  nr_sensing_stream_t stream;
  /// number of pilots the delay transform used; sets the delay resolution
  int n_pilots;
  /* IDFT size used; bin b is a delay of b / (idft_size * k_step * scs_hz). Generally,
  this is bigger than n_pilots. This implies smoothing in the delay domain. */
  int idft_size;
  /* Grid subcarrier of the first pilot. Kept because the range point spread
  function carries an exp(j*2*pi*a*u/N) factor of the grant start a, which a
  detector doing coherent removal needs and the delay response alone does not
  reveal. It cancels at the peak of a resolved path, which is why the profile
  itself never needed it. */
  int k_first;
  /// valid entries in h[]
  int n_bins;
  /// complex delay response
  cf_t h[NR_SENSING_MAP_MAX_BINS_RANGE];
} nr_sensing_snapshot_t;

/// Ring buffer of delay responses for one Rx antenna
/// We can assume 3.41 GHz for the carrier, and 30 KHz for the SCS.
struct nr_sensing_history_s {
  nr_sensing_snapshot_t *ring;
  int depth;
  /// next index to write
  int head;
  /// entries held, saturates at depth
  int count;
  /* Snapshots pushed since the last map, counted per stream rather than in total.
  Streams arrive at rates the scheduler decides, so a shared counter would let a
  burst of one trigger a map of another, built from a handful of samples and a
  correspondingly useless Doppler resolution. Slots are claimed on first use and
  never released; a stream the scheduler has abandoned simply stops counting. */
  struct {
    bool used;
    nr_sensing_stream_t key;
    int count;
  } since_last_map[NR_SENSING_MAX_STREAMS];
  int scs_hz;
  int ofdm_symbol_size;
  /// DL carrier frequency in Hz, needed to turn the Doppler axis into speeds
  uint64_t carrier_hz;
  /* Several DL actors (--num-dl-actors) push into the same history at the same time, and
  a map copies it while they do. Held by push, by the counter functions and by
  nr_ue_sensing_history_take(); a copy handed to a map task is only read by that task. */
  pthread_mutex_t lock;
};

/* 
One range-Doppler map. power[b * n_freq + f] is bin b at the f-th frequency of a
grid spanning [-f_max_hz, +f_max_hz].
Bin b is a path length of b * m_per_bin metres.
Frequency f_d converts to speed as v = f_d * lambda / 2, with lambda = c / carrier_hz.
The data is used for the python plots.

power is in units of the map's own noise floor, not of anything absolute: it is
divided by noise_ref before it is handed out, so 1.0 is the floor and a cell's value
is how far it stands above it. See noise_ref.
*/
typedef struct {
  int n_bins;
  int n_freq;
  /// number of snapshots the map was built from
  int n_snapshots;
  /* Which measurement process the map describes. Carried rather than derived so a
  dumped map says on its own which beam it is of: at rank 2 the layers produce one
  map each, in the same slot, from the same antennas, on the same delay axis, and
  the layer index is the only thing that tells them apart. */
  nr_sensing_stream_t stream;
  // The bin where the power is maximum, meaning its the LOS
  double bin_los;
  /// smallest and largest grant (pilots per snapshot) in the window
  int n_pilots_min;
  int n_pilots_max;
  /// distinct grant positions a_m (lattice index of the first pilot) in the window
  int n_positions;
  double m_per_bin;
  double f_max_hz;
  /// observation span in seconds
  double t_span_s;
  /* Absolute time of the window's first snapshot, s, on the capture's sample clock.
  Lets maps of consecutive windows be placed on one time axis (map averaging). */
  double t_start_s;
  double carrier_hz;
  /* Level power[] was divided by: the median cell of this map, which is the noise
  floor since noise holds the large majority of the cells and the median ignores the
  rest. Kept so the division can be undone and, more usefully, so the level itself can
  be read: it is the map wide gain that the outdoor dumps showed swinging 10 to 25 dB
  between consecutive windows as the direct path was blocked and cleared.

  Why the division happens at all: a CFAR does not need it, being a ratio against its
  own training cells, but everything that compares one map with another does. Averaging
  consecutive maps, which is how the SNR lost to a short window is won back, is
  dominated by the loudest map when the levels differ by 10 dB rather than improving
  anything. A tracker associating detections by strength has the same problem, and so
  does the eye going from one map to the next.

  This is the cruder half of the fix. It equalises finished maps; it cannot repair what
  the gain drift did inside the window, where it also broke the clutter cancellation.
  That needs the drift removed per snapshot, before the transform. */
  double noise_ref;
  float power[NR_SENSING_MAP_MAX_BINS_RANGE * NR_SENSING_MAP_MAX_BINS_FREQ];
} nr_sensing_map_t;

/* The four criteria of the OBSERVATION WINDOW block, one function each, so a caller
   sizing a run can ask for any one of them on its own. All return seconds.

   Criterion 1: the largest t_span with 4 * f_max * t_span <= MAX_BINS_FREQ - 1, i.e.
   the point budget at two per resolution cell. f_max is the half width of the Doppler
   axis, max_speed_ms plus the replica margin, as nr_ue_sensing_range_doppler() builds
   it. */
double nr_ue_sensing_span_grid(double f_max_hz);

/* Criterion 2: cell_m / (2 max_speed_ms), the time a target takes to cross one range
   cell. n_pilots is the widest grant in the window, so cell_m is the narrowest peak
   and the bound the tightest; a narrower grant has a wider peak and a looser bound,
   which this deliberately does not exploit. */
double nr_ue_sensing_span_migration(int n_pilots, int k_step, int scs_hz, double max_speed_ms);

/// Criterion 3: sqrt(lambda / 2a), with a = NR_SENSING_TARGET_MAX_ACCEL_MS2.
double nr_ue_sensing_span_accel(double carrier_hz);

/// Criterion 4, the lower bound: lambda / (2 * NR_SENSING_TARGET_DV_MS).
double nr_ue_sensing_span_min(double carrier_hz);

/* The bounds the transform uses: criterion 1 above and criterion 4 below.
   n_pilots and k_step are unused since criterion 2 left the window; kept so callers
   do not change.

   n_pilots is the widest grant expected in the window. It is only known after the
   snapshots are gathered, and the gather needs a bound first, so pass 0 to leave
   criterion 2 out; the caller then re-evaluates with the real grant and trims. See
   nr_ue_sensing_range_doppler().

   feasible is false when the requirements contradict each other, and the caller is
   expected to report that rather than quietly honour one end. */
nr_sensing_span_t nr_ue_sensing_span_bounds(const nr_sensing_history_t *hist,
                                            double max_speed_ms,
                                            int n_pilots,
                                            int k_step);

/* The window one map is held to: the short end of the feasible range (criterion 4),
   or the upper bound when the two cross. This is what nr_ue_sensing_range_doppler()
   gathers over, and what a caller deciding when to build a map has to use, so that
   the cadence and the window are the same number. */
double nr_ue_sensing_window_s(const nr_sensing_history_t *hist, double max_speed_ms, int n_pilots, int k_step);

/// Allocate a history able to hold depth snapshots. Returns false on failure.
bool nr_ue_sensing_history_init(nr_sensing_history_t *hist,
                                int depth,
                                int scs_hz,
                                int ofdm_symbol_size,
                                uint64_t carrier_hz);

void nr_ue_sensing_history_free(nr_sensing_history_t *hist);

/* Append one time sample. Copies NR_SENSING_MAP_MAX_BINS_RANGE bins of h_time.
   Overwrites the oldest entry once full.

   stream identifies the measurement process; k_step and k_offset are filled in from
   lattice, so the caller only has to supply the port bitmap and the layer. */
void nr_ue_sensing_history_push(nr_sensing_history_t *hist,
                                uint64_t t_sample,
                                nr_sensing_stream_t stream,
                                const pilot_lattice_t *lattice,
                                int idft_size,
                                const cf_t *h_time);

/// Snapshots of this stream pushed since its last map, 0 if it has none
int nr_ue_sensing_since_last_map(const nr_sensing_history_t *hist, const nr_sensing_stream_t *stream);
void nr_ue_sensing_clear_since_last_map(nr_sensing_history_t *hist, const nr_sensing_stream_t *stream);

/* Take what one map needs from a live history, in one step under its lock: check the
   counter, clear it, copy the ring.

   The counters of stream are cleared for layers 0..n_layers-1, and the ring is copied into
   dst, which gets its own allocation (release it with free(dst->ring)). When min_count > 0
   nothing is done and false is returned unless stream has reached min_count, so two DL
   actors that both saw the count cannot start the same map twice. */
bool nr_ue_sensing_history_take(nr_sensing_history_t *src,
                                const nr_sensing_stream_t *stream,
                                int n_layers,
                                int min_count,
                                nr_sensing_history_t *dst);

/* The slow-time samples nr_ue_sensing_range_doppler() worked from, optionally
   handed back to the caller: gathered for one reference signal, in chronological
   order, normalised for pilot count and with the clutter mean removed.

   Exported because a detector needs the samples, not the power map. Two reasons it
   cannot just re-read the ring itself: it would duplicate the gather and the
   conditioning, and it would have to stay in step with them; and the conditioning
   is exactly what makes a static target constant in slow time, which any coherent
   method depends on.

   Nothing about the Doppler grid is exported. The samples are in the time domain,
   so a caller is free to build a grid of its own, and one that wants to see the
   TDD replicas has to: the grid built here deliberately stops at the point where
   they begin.

   Buffers are caller allocated and caller sized:
     t_s          at least hist->depth doubles
     h_re, h_im   at least hist->depth * n_bins floats, n_bins as the map reports it
   n_snap and n_bins are set on return and say how much was written. */
typedef struct {
  int n_snap;
  int n_bins;
  /* t_sample of the first snapshot, the instant t_s counts from. t_s alone cannot
  tell two windows apart that sample the same pattern from different starts, which
  is what the AoA has to rule out across antennas. */
  uint64_t t0_sample;
  /// [n_snap] sample instants in seconds, relative to the first
  double *t_s;
  /// [n_snap][n_bins] conditioned delay responses, row major
  float *h_re;
  float *h_im;
  /// pilots of the widest snapshot; the count everything was normalised to
  int win_len;
  /// first lattice index carrying a pilot, (k_first - k_first % k_step) / k_step
  int win_start;
  /// IDFT length the delay responses were produced with
  int idft_size;
  /// metres of path length per delay bin
  double m_per_bin;
  /* Degree of the slow trend removed from every bin, -1 when nothing was removed. The
  samples then lie outside the span of the polynomials of that degree in t_s, which a
  model based method has to know: see nr_ue_sensing_slow_basis(). */
  int trend_degree;
} nr_sensing_slowtime_t;

/* Orthonormal basis of the slow trend of degree `degree` on the sample times t[0..n-1]:
   the polynomials 1, tau, ..., tau^degree with tau = t centred and scaled to [-1, 1],
   made orthonormal by Gram-Schmidt. This is the subspace the clutter stage removes from
   every bin. q needs room for NR_CLUTTER_SLOW_TREND_MAX + 1 vectors of n entries.
   Returns the number of vectors written, 0 when degree < 0. */
int nr_ue_sensing_slow_basis(const double *t, int n, int degree, double q[][NR_SENSING_HISTORY_DEPTH]);

/* How the static scene is removed from the slow-time samples before the Doppler
   transform. */
typedef enum {
  /// leave it in
  NR_CLUTTER_NONE = 0,
  /* Subtract each bin's slow-time mean (or its slow trend, see
  NR_CLUTTER_SLOW_TREND_DEGREE). Assumes the clutter contribution is constant, or slowly
  varying, in m, which holds only while the grant does not move */
  NR_CLUTTER_MEAN = 1,
  /* Fit the strongest static paths with each snapshot's own point spread function and
  remove them, then take the mean of what is left.

  A static path at bin u contributes alpha * K_m(b - u) to bin b of snapshot m, with
  K_m(u) = A_n_m(u) * exp(j*2*pi*c_m*u/N) the kernel of that snapshot's grant and
  alpha one complex amplitude constant over the window. Fitting that one scalar by
  least squares and subtracting removes the path with the right shape and the right
  height in every snapshot, however the grant moved.

  The paths are removed one at a time:
    1. the direct path, always, at bin_los;
    2. then, up to NR_CLUTTER_MAX_PATHS in total, the strongest peak of the slow-time
       mean of what is left (the zero-Doppler row of the map). Movers average out in
       that mean, so its peaks are static scatterers.
  The loop stops as soon as the peak fails one of two tests:
    - static:  S = |slow-time mean|^2 / slow-time mean of |r|^2 must be at least
               NR_CLUTTER_STATIC_MIN. S is 1 for a constant value, near 0 for a mover
               or for noise, so a target is never fitted as clutter;
    - strong:  the peak must be NR_CLUTTER_SNR_MIN above the noise left in the mean.
  So the number of paths is not fixed: it is however many static peaks stand out.

  Each fit removes one path. The mean that follows takes the weak static scatterers
  the loop left, and by then it is no longer dominated by the strong ones. The two are
  complementary, which is why this mode runs both.

  NR_CLUTTER_MAX_PATHS 1 gives the direct path only. A target slower than one Doppler
  cell looks static and is removed, as with the plain mean. */
  NR_CLUTTER_KERNEL = 2,
} nr_sensing_clutter_t;

/// Distinct (n_m, a_m) pairs a window is expected to hold; see NR_CLUTTER_KERNEL
#define NR_CLUTTER_MAX_KERNELS 512

/// Static paths NR_CLUTTER_KERNEL fits at most, the direct path included
#define NR_CLUTTER_MAX_PATHS 4

/// Minimum S = |mean|^2 / mean(|r|^2) for a peak to be taken as static; 1 is perfectly static
#define NR_CLUTTER_STATIC_MIN 0.8

/// Minimum power of a peak of the slow-time mean over its noise (median power / M), linear; 10 dB
#define NR_CLUTTER_SNR_MIN 10.0

/// A new path must be at least this many bins from one already fitted, so a residue is not refitted
#define NR_CLUTTER_MIN_SEP_BINS 3

/* Largest correction the direct-path normalisation may apply to one snapshot, in dB
either way.

Every snapshot is divided by the complex amplitude of its own direct path, which is
the one thing in the scene known to be static, so whatever it does from one snapshot
to the next is the gain and phase drift of the chain rather than the scene. See
NR_CLUTTER_LOS_NORM.

The clamp is there for the case the normalisation is worst at: the direct path
briefly blocked, by the operator walking through it among other things. The measured
amplitude then collapses, and dividing by it would multiply that snapshot, which is
mostly noise by then, up to the level of all the others. Clamping leaves such a
snapshot under-corrected, which costs a little cancellation, rather than letting it
dominate the window. Hitting the clamp is logged: it says the direct path was lost,
which is worth knowing on its own. */
#define NR_CLUTTER_LOS_MAX_CORR_DB 10.0

/* COMB REMOVAL

What is left after the kernel fits and the slow trend. Both of those remove things
that are constant or slow in the window: the kernel fit removes the mean of a path's
amplitude sequence, the trend removes the mean and a quadratic of each bin's. The comb
is neither. It is a line spectrum at multiples of 100 Hz, which over a 146 ms window is
some fifteen cycles, so no polynomial of degree two comes near it and nothing in either
stage has any representation of it. Measured on the outdoor dumps it carried 81% of the
map energy in 5% of the cells, 20 dB over the floor, and 90 to 95% of what the TDD
detector accepted as targets sat on one of its lines.

Where it comes from: whatever repeats in the snapshot stream. The 200 Hz family is the
TDD period itself, 5 ms, and is not removable at the source, since a UE has no downlink
to measure during the uplink slots. The 100 Hz family needs a 10 ms cycle on top of
that, one radio frame, and its origin is still open. Both are covered by one tone
family at 100 Hz, 200 Hz being its second harmonic.

The removal is the same projection the slow trend already does, with the tones added to
the subspace. What it cannot be is that projection done per range bin: a target sitting
on a comb line is then indistinguishable from the comb and goes with it, measured at
-180 dB in simulation, which is not a filter but a deletion.

What separates them is range. The modulation multiplies a whole snapshot, so every
range bin sees the same time signature scaled by its own clutter strength: one
modulation, many amplitudes, a rank one outer product. A target sits in one or two
bins and has no such structure. Measured on the dumps the comb block of a map is 87 to
99.5% rank one and its shape holds to 0.96-0.996 from map to map, so the model is the
data's own.

So the tone amplitudes are fitted jointly across range under a rank one constraint
rather than per bin. Counting the freedom says why it works: 128 bins by 12 tones is
1536 free numbers, enough to describe anything including one bin's target, while one
profile and one modulation is 128 + 12 = 140 and can only describe what is present at
every range with a common shape. The target contributes about 1/n_bins of the fit.
Simulated with a target placed exactly on a line: -0.4 dB against the -180 dB of the
per bin version, at the same clutter suppression.

Only the tones are constrained this way. The polynomial stays per bin, where it belongs:
each range bin holds its own scatterers with their own levels and phases, and there is
no common shape to exploit near zero Doppler. */
#ifndef NR_CLUTTER_COMB_REMOVE
#define NR_CLUTTER_COMB_REMOVE 1
#endif

/* Fundamental of the comb, as a multiple of the TDD period. The sampling pattern
repeats every TDD period, which alone would put lines at multiples of 1/T_TDD; the
dumps also show the family that needs twice that period, so the fundamental is
1/(2*T_TDD) = 100 Hz at 5 ms and the TDD lines are its even harmonics. */
#define NR_COMB_TDD_MULTIPLE 2

/* Harmonics of the fundamental removed, each one as a pair +k and -k. Six covers
600 Hz, past which the dumps show the lines already down in the floor. The subspace
removed is 2 * this many vectors, so it also sets what the window costs: with n
snapshots a fraction 2*NR_COMB_MAX_HARMONIC/n of any target's energy goes with it,
0.4 dB at 12 vectors and 120 snapshots. */
#define NR_COMB_MAX_HARMONIC 6

/// Power iterations for the dominant singular pair. Converges in a handful on a matrix this small.
#define NR_COMB_POWER_ITERS 24

/* Whether to divide each snapshot by its own direct path before anything else.

What it fixes: the clutter stages all assume a static scatterer gives the same
value in every snapshot. The kernel fit assumes it literally, holding one complex
alpha per path across the whole window, and the slow trend assumes it up to a
quadratic. Anything that multiplies a whole snapshot breaks that assumption without
being visible to either: transmit power, the precoder at rank 2, propagation as the
scene moves, and residual carrier frequency offset, which is the awkward one because a
phase that advances with time is a Doppler shift, so it slides the clutter off zero
where the trend subspace can no longer reach it.

How: the direct path is static and tens of dB above everything else, so its complex
amplitude per snapshot is a clean measurement of exactly that common factor. Dividing
it out puts the static scene back to constant, which is what the stages below need,
and leaves every phase referenced to the direct path, which is the natural bistatic
reference anyway.

What it does not fix: anything that is not a common factor. The point spread function
changing shape with the grant width is not, which is what K_m is for, and a scatterer
that genuinely fluctuates is not either. */
#ifndef NR_CLUTTER_LOS_NORM
#define NR_CLUTTER_LOS_NORM 0
#endif

/* How far below the strongest bin the direct path may sit, in dB. The direct path is
the shortest path, so it is the earliest arrival, not necessarily the strongest: outdoors
a cluster of reflections close to the receiver can match or beat it, and taking the
strongest bin then made bin_los jump between the two ends of that cluster from map to
map (bins 2 and 6.5 on one run). bin_los is the earliest local peak of the window's
energy profile within this margin of the strongest. */
#define NR_SENSING_LOS_FIRST_DB 6.0

/* Bins either side of a path that its kernel is computed on and subtracted from.
Beyond ~48 bins the Hann skirt is some 100 dB down, below anything a cf_t response
holds, and limiting it keeps the cost of a path independent of the map size. */
#define NR_CLUTTER_KERNEL_HALF_SPAN 48

/* Slow-trend clutter filter, a testing switch: rebuild after changing it.

   After the static paths are fitted (kernel mode), or straight away (mean mode), each range
   bin's slow-time sequence z[m] loses everything that varies slowly in time, not only its
   mean. The slow part is a polynomial of this degree in the snapshot time t_m, fitted by
   least squares for each bin:
       z[m] <- z[m] - (c_0 + c_1 t_m + ... + c_K t_m^K)

     0  the plain mean, i.e. the behaviour before this switch existed
     1  mean and linear drift
     2  mean, drift and curvature

   Why: outdoors the static scene is not exactly constant over a window (something moving
   slowly next to the UE or the direct path, a slow gain or phase drift). The mean leaves that
   slow part in, just off 0 Hz, and the slow-time sampling copies it to +-200 Hz, +-400 Hz, ...
   (+-8.8, +-17.6 m/s), where it shows up as mirror pairs of fake targets next to the direct
   path. Removing the trend removes it, and so its copies.

   Price: a target slow enough to look like a trend over the window is removed with it.
   Measured on an outdoor window of 116 ms at degree 2: about 5 dB less leak above 2 m/s next
   to the direct path, while a target at 0.5 m/s loses ~2 dB, 1 m/s ~0.3 dB, 0.25 m/s ~15 dB. */
#ifndef NR_CLUTTER_SLOW_TREND_DEGREE
#define NR_CLUTTER_SLOW_TREND_DEGREE 0
#endif
/// highest degree NR_CLUTTER_SLOW_TREND_DEGREE may take
#define NR_CLUTTER_SLOW_TREND_MAX 3
#if NR_CLUTTER_SLOW_TREND_DEGREE < 0 || NR_CLUTTER_SLOW_TREND_DEGREE > NR_CLUTTER_SLOW_TREND_MAX
#error "NR_CLUTTER_SLOW_TREND_DEGREE must be between 0 and NR_CLUTTER_SLOW_TREND_MAX"
#endif

/* ---------------------------------------------------------------------------
 Runtime clutter-removal parameters.

 The NR_CLUTTER_* / NR_COMB_* / NR_SENSING_LOS_FIRST_DB macros above stay as the
 compile-time DEFAULTS and, where an array is sized by one (NR_CLUTTER_MAX_PATHS,
 NR_CLUTTER_MAX_KERNELS, NR_COMB_MAX_HARMONIC, NR_CLUTTER_SLOW_TREND_MAX), as the
 hard CAPACITY. This struct carries the ACTIVE value nr_ue_sensing_range_doppler()
 uses, so one build can run at different settings. Every field must stay within its
 capacity; the transform
 asserts the ones that size arrays.
--------------------------------------------------------------------------- */
typedef struct {
  nr_sensing_clutter_t clutter_mode; ///< NONE / MEAN / KERNEL
  int    max_paths;          ///< kernel mode: static paths fitted, 1..NR_CLUTTER_MAX_PATHS
  int    trend_degree;       ///< slow-trend polynomial degree, 0..NR_CLUTTER_SLOW_TREND_MAX
  int    kernel_half_span;   ///< bins either side of a path its kernel touches
  double static_min;         ///< S = |mean|^2/mean(|r|^2) for a peak to count as static
  double snr_min;            ///< peak-over-noise (linear) for a path to be fitted
  int    min_sep_bins;       ///< minimum separation between two fitted paths
  double los_first_db;       ///< how far below the strongest the LOS may sit
  bool   los_norm;           ///< divide each snapshot by its own direct path first
  double los_max_corr_db;    ///< clamp, dB either way, on that division
  bool   comb_remove;        ///< remove the rank-1 TDD replica comb
  int    comb_harmonic;      ///< harmonics removed, 1..NR_COMB_MAX_HARMONIC
  int    comb_tdd_multiple;  ///< comb fundamental = 1/(this * T_TDD)
  int    comb_power_iters;   ///< power iterations for the rank-1 comb
} nr_sensing_params_t;

/// Fill p with the compile-time defaults (clutter_mode left at KERNEL; the caller sets it).
static inline void nr_sensing_params_default(nr_sensing_params_t *p)
{
  p->clutter_mode      = NR_CLUTTER_KERNEL;
  p->max_paths         = NR_CLUTTER_MAX_PATHS;
  p->trend_degree      = NR_CLUTTER_SLOW_TREND_DEGREE;
  p->kernel_half_span  = NR_CLUTTER_KERNEL_HALF_SPAN;
  p->static_min        = NR_CLUTTER_STATIC_MIN;
  p->snr_min           = NR_CLUTTER_SNR_MIN;
  p->min_sep_bins      = NR_CLUTTER_MIN_SEP_BINS;
  p->los_first_db      = NR_SENSING_LOS_FIRST_DB;
  p->los_norm          = NR_CLUTTER_LOS_NORM;
  p->los_max_corr_db   = NR_CLUTTER_LOS_MAX_CORR_DB;
  p->comb_remove       = NR_CLUTTER_COMB_REMOVE;
  p->comb_harmonic     = NR_COMB_MAX_HARMONIC;
  p->comb_tdd_multiple = NR_COMB_TDD_MULTIPLE;
  p->comb_power_iters  = NR_COMB_POWER_ITERS;
}

/* Delay response (IDFT) of one grant's Hann window, sampled at u = b - u0.

   This is K_m of NR_CLUTTER_KERNEL: the transform of the taper alone, placed where
   that snapshot's pilots were.

   Its scale is arbitrary and does not need to match the data: the amplitude fitted
   against it absorbs any constant factor, so only the shape and the phase matter.

   n_pilots  : pilots the snapshot measured, the length of the taper
   a_m       : lattice index of its first pilot, (k_first - k_first % k_step) / k_step
   u0        : delay bin the kernel is centred on, fractional. An integer centre
               leaves exp(-j*2*pi*c_m*frac/N) at the peak, which moves with the grant
   idft_size : the N the delay response was produced with
   n_bins    : entries to fill, one per delay bin of the map
   kre, kim  : receive n_bins entries each */
void nr_ue_sensing_clutter_kernel(int n_pilots, int a_m, double u0, int idft_size, int n_bins, float *kre, float *kim);

/* Range-Doppler map from the accumulated history, for one reference signal.

   The slow-time transform is a direct DFT at the recorded timestamps.
   Not an DFT over a uniform grid: the reference symbols are not uniformly spaced.

   max_speed_ms is the fastest target the Doppler axis has to cover, in m/s.
  The axis reaches that speed plus one TDD period of margin,
   which the TDD detector needs to see the replicas of a target at the edge. The grid
   holds 2 points per resolution cell, the resolution being 1/t_span Hz.

   clutter_mode says how the static scene is removed first. The direct path is
   static, far stronger than any reflection, and its sidelobes cover the delay axis,
   so without this the zero-Doppler ridge hides everything that moves.

   max_paths caps the static paths NR_CLUTTER_KERNEL fits, direct path included,
   clamped to 1..NR_CLUTTER_MAX_PATHS. Pass NR_CLUTTER_MAX_PATHS normally, 1 for the
   direct path only. Ignored by the other modes.

   n_freq_fixed, when positive, is the number of Doppler points to use instead of the
   2 per resolution cell this window's t_span asks for. The axis runs over the same
   +-f_max whatever the window, so a map built with another map's n_freq lands on that
   map's grid exactly, which is what averaging maps cell by cell needs: two Rx chains
   whose windows differ by a snapshot otherwise get 379 and 383 points. Pass 0 to size
   the grid from the window.

   slow_out receives the conditioned slow-time samples when it is not NULL; see
   nr_sensing_slowtime_t. Passing NULL skips it entirely and costs nothing.

   Returns the number of snapshots used, 0 if too few. */
int nr_ue_sensing_range_doppler(const nr_sensing_history_t *hist,
                                const nr_sensing_stream_t *stream,
                                int n_snap_max,
                                double max_speed_ms,
                                const nr_sensing_params_t *p,
                                int n_freq_fixed,
                                nr_sensing_map_t *map,
                                nr_sensing_slowtime_t *slow_out);

/* 
One decision of the TDD sidelobe detector, in map coordinates, so it can be drawn
on the plot. Verdict mirrors nr_tdd_verdict_t: 0 accepted, 1 rejected as a replica of another
peak, 2 rejected as residue of a peak already accepted. 
*/
typedef struct {
  float range_m;
  float speed_ms;
  float snr_dB;
  int verdict;
} nr_sensing_marker_t;

/// Cells per map the AoA stage reports on; see nr_ue_aoa.h
#define NR_SENSING_AOA_MAX 16

/* One angle estimate, attached to the range-Doppler cell it was measured at.

Kept here next to nr_sensing_marker_t rather than in nr_ue_aoa.h because it is a
dump payload: nr_ue_sensing_dump_map() writes it, and putting it in the algorithm
header would make this one include that one for no other reason. */
typedef struct {
  float range_m;
  float speed_ms;
  /// estimated angle of arrival in degrees, 0 at the array normal
  float angle_deg;
  /// how far the cell is above the local noise (CFAR estimate), in dB
  float power_dB;
  /* Spatial peak over the mean of the spatial spectrum, in dB. A single clean plane
  wave on Mr antennas gives 10*log10(Mr), i.e. 6 dB at four; well below that means
  the cell holds several arrivals or noise and the angle should not be trusted. */
  float quality_dB;
} nr_sensing_aoa_t;

/* Append a map to a CSV file, one line per map.

   markers and aoa are independent and either may be absent: the TDD detector and
   the AoA stage are separate options and neither is a precondition for the other.
   Pass NULL / 0 for the one that did not run. */
void nr_ue_sensing_dump_map(const char *path,
                            int frame,
                            int slot,
                            int aarx,
                            const nr_sensing_map_t *map,
                            const nr_sensing_marker_t *markers,
                            int n_markers,
                            const nr_sensing_aoa_t *aoa,
                            int n_aoa);

/* Which map triggers the raw snapshot dump, counted from the first map produced.
Debug aid for checking the delay-response model against the data, so it fires once
per run rather than on every map. */
#define NR_SENSING_SNAP_DUMP_MAP 10

/* Dump the raw slow-time snapshots a map was built from, one line per snapshot.

   Raw on purpose: neither the pilot-count gain nor the clutter mean is applied, so
   what lands in the file is h_m[b] exactly as the per-symbol IDFT produced it. That
   is what the delay-response model describes, and the conditioning can be reapplied
   offline from the per-snapshot n_pilots the file carries.

   Selection is the same as nr_ue_sensing_range_doppler(): the most recent snapshots
   of this stream, written oldest first, bounded both by n_snap_max and by the window
   the map would have used. max_speed_ms is taken for that reason alone, and has to
   be the value the map is built with or the file stops describing it. Overwrites the
   file.

   Columns: snap, t_sample, t_rel_s, k_step, k_first, n_pilots, idft_size, n_bins,
   scs_hz, fs_hz, carrier_hz, then n_bins interleaved (re, im) pairs. */
void nr_ue_sensing_dump_snapshots(const char *path,
                                  const nr_sensing_history_t *hist,
                                  const nr_sensing_stream_t *stream,
                                  int n_snap_max,
                                  double max_speed_ms);

/* test_record_samples mode: one file per map holding every (chain, layer) it used, so a
 reprocessor (scripts/sensing/replay_samples.py) can rebuild the whole map offline and
 sweep any parameter. One call per (chain, layer) block appends to the same file; the
 first block passes truncate and header true. aarx tags the chain. Columns: aarx, layer,
 snap, t_sample, t_rel_s, k_step, k_first, n_pilots, idft_size, n_bins, scs_hz, fs_hz,
 carrier_hz, max_speed_ms, tdd_slots, then n_bins interleaved (re, im) pairs. Returns rows written. */
int nr_ue_sensing_dump_snapshots_tagged(const char *path,
                                        const nr_sensing_history_t *hist,
                                        const nr_sensing_stream_t *stream,
                                        int n_snap_max,
                                        double max_speed_ms,
                                        int aarx,
                                        bool truncate,
                                        bool header);

#endif // __NR_UE_MAP__H__
