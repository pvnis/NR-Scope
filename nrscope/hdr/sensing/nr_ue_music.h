/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * \brief Doppler-only MUSIC on the conditioned slow-time samples of a range-Doppler map.
 *
 * WHERE IT SITS
 *
 * nr_ue_sensing_range_doppler() gathers the snapshots, normalises them for the pilot
 * count, removes the clutter (kernel fits, then the slow trend) and finally runs a DFT
 * at the true sample instants, bin by bin. This stage replaces only that last DFT. It
 * reads the samples exactly as the DFT did, through nr_sensing_slowtime_t, so every
 * conditioning step before it is shared and unchanged.
 *
 * THE MODEL
 *
 * In range bin b, the slow-time vector of one observation (one Rx antenna, one layer) is
 *     z[i] = sum_k beta_k exp(j*2*pi*f_k*t_i) + noise,   i = 0..n-1
 * with the t_i the real, irregular, sample instants. MUSIC needs several observations of
 * that vector in which the beta_k of different targets vary independently, so that the
 * covariance has one dimension per target. Three sources are used, stacked as columns:
 *   - the Rx antennas: a target at another angle has another phase progression across them;
 *   - the layers: each layer carries another precoder, so another beta per path;
 *   - the neighbouring range bins b-R..b+R: two targets at slightly different delays
 *     leave different amplitude patterns across them.
 * Without --sensing-antenna-avg / --sensing-layer-avg only the range bins remain, 2R+1
 * columns, and two targets at the same delay are then coherent and cannot be split.
 *
 * Slow-time smoothing (subarrays shifted in time) is not used: it needs a uniform
 * sampling, which the TDD pattern and the scheduler do not give.
 *
 * THE COMPUTATION
 *
 * With X the n x L matrix of those columns, the signal subspace is taken from the L x L
 * Gram matrix G = X^H X rather than from the n x n covariance: its non-zero eigenvalues
 * are the same, an eigenvector v of G gives the unit vector u = X v / sqrt(lambda) of the
 * covariance, and L is a few tens where n is a few hundred. The eigen decomposition is a
 * cyclic complex Jacobi, exact enough at this size and free of any LAPACK dependency.
 *
 * The pseudo-spectrum is
 *     P(f) = |a(f)|^2 / (|a(f)|^2 - sum_k |u_k^H a(f)|^2)
 * so 1 where nothing is (the floor, 0 dB) and large at a target.
 *
 * STEERING VECTORS AND THE CLUTTER FILTER
 *
 * The slow-trend stage projects every bin out of span{q}, the polynomials of degree
 * trend_degree in t. The data therefore never has a component there. The steering vector
 * is projected the same way, a_p = a - Q Q^T a, and |a_p|^2 is what enters P(f). Without it
 * the removed subspace would read as "noise subspace" and dig an artificial hole around
 * 0 Hz on top of the one the filter really made. Near 0 Hz, where the filter leaves almost
 * nothing of a(f), P is set to the floor.
 *
 * MODEL ORDER
 *
 * The number of signals K is chosen per bin. The largest eigenvalue of a bin that holds
 * only noise sets the reference, measured as the median over bins of the largest
 * eigenvalue (after clutter removal most bins hold no target). An eigenvalue counts as a
 * signal when it is NR_MUSIC_EIG_THRESH_DB above that reference, up to
 * NR_MUSIC_MAX_ORDER and L. NR_MUSIC_FIXED_ORDER forces K instead.
 *
 * WHAT THE OUTPUT IS NOT
 *
 * P(f) is not a power. Peak heights say how well a steering vector fits the signal
 * subspace, not how strong the target is, so CFAR thresholds and SNRs taken on the DFT
 * map do not carry over. The eigenvalue test also integrates less than the DFT does (it
 * works on the L observations, not coherently along the n samples), so a target the DFT
 * shows only a few dB above the floor may get no signal dimension at all here.
 */

#ifndef __NR_UE_MUSIC__H__
#define __NR_UE_MUSIC__H__

#include "nrscope/hdr/sensing/nr_ue_map.h"

/// Neighbouring range bins taken on each side as extra observations (R above)
#define NR_MUSIC_RANGE_HALF 1

/// Largest number of signals MUSIC may assume in one range bin
#define NR_MUSIC_MAX_ORDER 4

/// Force this number of signals in every bin (capped at L); 0 selects it from the eigenvalues
#define NR_MUSIC_FIXED_ORDER 0

/// An eigenvalue is a signal when it is this far above the noise reference, in dB
#define NR_MUSIC_EIG_THRESH_DB 3.0

/// Observations (Rx antennas x layers) MUSIC accepts for one map
#define NR_MUSIC_MAX_OBS 16

/// Ceiling on the pseudo-spectrum, 60 dB above its floor
#define NR_MUSIC_PSEUDO_MAX 1e6

/// What one MUSIC map was built from, for the logs
typedef struct {
  /// observations whose time axis matched the first one and were used
  int n_obs;
  /// columns per range bin away from the edges, n_obs * (2R+1)
  int n_cols;
  /// the noise reference, largest eigenvalue of a noise-only bin
  double noise_ref;
  /// range bins given at least one signal dimension
  int bins_with_signal;
  /// highest model order used in any bin
  int order_max;
} nr_music_stats_t;

/* Doppler MUSIC map of one window.

   obs      : n_obs conditioned slow-time sets, as nr_ue_sensing_range_doppler() exports
              them, one per (Rx antenna, layer). An entry with n_snap 0 is skipped, and so
              is one whose sample instants differ from those of the first usable entry.
   map      : on input carries the grid of the DFT map of the same window (n_bins,
              n_freq, f_max_hz); its power[] is overwritten with the pseudo-spectrum, on
              the same grid, so the dump and the plot scripts read it unchanged.
   stats    : optional, filled when not NULL.

   Returns the number of observations used, 0 when nothing could be computed (power[]
   is then left untouched). */
int nr_ue_music_doppler(const nr_sensing_slowtime_t *obs, int n_obs, nr_sensing_map_t *map, nr_music_stats_t *stats);

/* Eigen decomposition of a Hermitian L x L matrix by cyclic Jacobi rotations.

   a_re, a_im : [L*L] row major, destroyed
   lambda     : [L] eigenvalues, sorted in decreasing order
   v_re, v_im : [L*L] row major, column k is the unit eigenvector of lambda[k]

   Exposed for the unit check only. */
void nr_ue_music_eig_hermitian(int L, double *a_re, double *a_im, double *lambda, double *v_re, double *v_im);

#endif // __NR_UE_MUSIC__H__
