/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * \brief Target detection on a range-Doppler map whose slow-time sampling is
 *        interrupted by the TDD pattern.
 *
 * Implements the method of
 *   M. Henninger, L. Giroto de Oliveira, S. Saur, A. Grudnitsky, T. Wild and
 *   S. Mandelli, "Target Detection for ISAC with TDD Transmission",
 *   arXiv:2504.19260v2, 2025.
*
 *
 * THE PROBLEM
 *
 * PDSCH only exists in downlink slots, so the slow-time samples feeding the
 * Doppler transform come in bursts separated by the uplink gaps, and that pattern
 * repeats with the TDD period T_TDD. A periodic sampling pattern cannot tell
 * Doppler f from Doppler f + k/T_TDD apart: writing every sample time as
 * t = a_r + q*T_TDD, the q dependence of exp(-j2*pi*(k/T_TDD)*t) is exp(-j2*pi*k*q),
 * which is 1 for every q. Only the within-period positions a_r separate the two,
 * and there are far too few of them to do it well.
 *
 * A single target therefore paints a comb of replicas along the Doppler axis at
 * multiples of 1/T_TDD, typically only 10 to 15 dB down. The paper calls them
 * impulsive sidelobes. They are not a windowing artefact and no taper removes
 * them: a taper shapes the response inside one period, while these come from the
 * period itself.
 *
 *
 * THE IDEA
 *
 * The replicas are deterministic. Their shape is the point spread function of the
 * sampling pattern, which is known, so a candidate peak can be tested: coherently
 * subtract the full modelled response of the candidate, including its replicas,
 * and see whether the replica positions lose power.
 *
 *   - a real target sits at the centre of its own replica comb, so removing it
 *     takes the whole comb down with it;
 *   - a replica is not the centre of anything. Removing it leaves the target that
 *     produced it untouched, and the other replicas of that target stay where they
 *     are.
 *
 *
 * WHAT DIFFERS FROM THE PAPER
 *
 * 1. Sampling model. The paper keeps a uniform grid of M OFDM symbols and zeroes
 *    the uplink ones, H = H_full * diag(d), which makes the Doppler PSF a closed
 *    form: a Dirichlet kernel times a Dirac train, their Eq. (10). This pipeline
 *    instead keeps only the samples it has, each with its true timestamp, and runs
 *    a direct transform at those instants.
 *
 *    So the closed form does not apply here, and the PSF is computed numerically
 *    from the actual sample times as
 *
 *        W_D(x) = sum_i exp(-j*2*pi*x*t_i) .
 *
 *    That is exact for any pattern, including one with dropped slots or scheduling
 *    holes, where Eq. (10) would already be wrong. nr_tdd_psf_doppler_closed()
 *    provides the paper's form anyway, for cross-checking against a synthetic
 *    strictly periodic pattern.
 *
 * 2. Range PSF. The paper assumes a rectangular band limitation, giving a
 *    Dirichlet kernel of order N. Here the frequency window is whatever taper was
 *    applied to the pilot lattice, so the range PSF is likewise computed from the
 *    window samples. It carries the exp(j*2*pi*a*u/N) factor of the grant start a,
 *    which the paper does not have because it always starts at subcarrier 0.
 *
 * 3. Removal domain. Only the PSF removal of their Sec. IV-B2b is implemented, not
 *    the CSI removal of IV-B2a. CSI removal needs the frequency domain matrix,
 *    which this pipeline no longer holds by the time the map exists: the history
 *    stores the delay response, not the lattice it came from. The paper reports
 *    the two as equivalent in performance, with PSF removal cheaper.
 *
 *
 * DOPPLER SPAN, AND WHY IT MUST BE WIDE HERE
 *
 * A map meant for a human should stop at +-1/(2*T_TDD), beyond which every feature
 * is a replica. This detector needs the opposite: it can only test a candidate by
 * looking at where its replicas would be, so it needs a grid that reaches at least
 * one replica spacing either side of the peak, i.e. f_max_hz >= 1/T_TDD, and more
 * if several replicas are to be checked.
 * 
 * => CHANGE THE F_MAX IN THE SENSING DUMP OF THE 2 MAP
 *
 * The two are not in conflict. Detect on the wide grid, where the artefacts are
 * visible and can be rejected on purpose; display the narrow one.
 */

#ifndef __NR_UE_TDD_DETECT__H__
#define __NR_UE_TDD_DETECT__H__

#include <stdbool.h>
#include <stdint.h>

/// Hard cap on reported targets, so callers can use plain arrays
#define NR_TDD_MAX_TARGETS 16

/// Candidate peaks the CFAR stage may return
#define NR_TDD_MAX_CANDIDATES 128

/// Rejected candidates reported back, so a plot can show what was thrown away
#define NR_TDD_MAX_REJECTED 32

/* What the CLEAN loop decided about a candidate. Reported rather than kept private
   because a rejection is the interesting output: a map showing only survivors looks
   the same whether the detector ran or not, and gives no way to see it working. */
typedef enum {
  /// accepted: removing it took its own replicas down with it
  NR_TDD_VERDICT_TARGET = 0,
  /// rejected: removing it left the replica positions untouched, so it is one
  NR_TDD_VERDICT_REPLICA = 1,
  /// rejected: inside the resolution of a target already accepted, so it is residue
  NR_TDD_VERDICT_DUPLICATE = 2,
  /* set aside: its mirror at -f on the same range bin is within NR_TDD_MIRROR_DB of it.
  One moving reflector has one Doppler sign; a balanced +-f pair is what a static path
  whose amplitude flickers looks like (someone shadowing it), but also what the swinging
  limbs of a person walking with no radial speed look like. Which one it is cannot be
  told from the map, so it is neither accepted nor thrown out, and not localised. */
  NR_TDD_VERDICT_UNCERTAIN = 3,
} nr_tdd_verdict_t;

/* How close in power, in dB, a target's mirror cell at -f must be for the target to be
   set aside as NR_TDD_VERDICT_UNCERTAIN. A target standing this far above its mirror
   keeps its verdict even with a modulation underneath it. */
#define NR_TDD_MIRROR_DB 6.0

/* One detected target, with off-grid range and Doppler from nr_tdd_focus(). */
typedef struct {
  /// delay in bins, fractional; range_m = bin * m_per_bin
  double bin;
  /// bistatic path length in metres
  double range_m;
  /// Doppler in Hz
  double f_hz;
  /// speed in m/s, f_hz * lambda / 2 (monostatic convention, see nr_ue_map.h)
  double speed_ms;
  /// complex amplitude at the peak of the periodogram
  double amp_re;
  double amp_im;
  /// peak power over the CFAR noise estimate, in dB
  double snr_dB;
  /// CLEAN iteration that examined it
  int iteration;
  /// what the loop decided; see nr_tdd_verdict_t
  int verdict;
} nr_tdd_target_t;

/* Complex range-Doppler map. Complex, not power: PSF removal is coherent, so the
   phase of every bin is needed. Indexed [bin * n_freq + f], matching the layout of
   nr_sensing_map_t::power. */
typedef struct {
  int n_bins;
  int n_freq;
  float *re;
  float *im;
} nr_tdd_cmap_t;

/* Everything the detector reads. It owns none of it.

   The slow-time samples are the complex delay responses already produced by
   nr_ue_sensing_delay_response(), one per reference symbol, with their timestamps.
   They must already be normalised for pilot count and, if wanted, clutter removed:
   this module does not repeat that work.

   Sign conventions follow the existing pipeline, and both matter for the PSF:
     - delay:   the response is an IDFT, X[q] = sum_p x[p] exp(+j*2*pi*p*q/N)
     - Doppler: the transform is C[f] = sum_i h_i exp(-j*2*pi*f*t_i) */
typedef struct {
  /* ---- slow time ---- */
  /// number of snapshots
  int n_snap;
  /// [n_snap] sample instants in seconds, increasing, origin arbitrary
  const double *t_s;
  /// [n_snap][n_bins] real and imaginary parts of the delay responses
  const float *h_re;
  const float *h_im;

  /* ---- delay axis ---- */
  /// bins held per snapshot
  int n_bins;
  /// IDFT length the delay responses were produced with
  int idft_size;
  /// pilots actually measured, i.e. the length of win[]
  int win_len;
  /// first lattice index carrying a pilot, the a of exp(j*2*pi*a*u/N)
  int win_start;
  /// [win_len] frequency taper, NULL for rectangular
  const float *win;
  /// metres of path length per delay bin
  double m_per_bin;

  /* ---- Doppler grid ---- */
  /// points on the Doppler grid
  int n_freq;
  /// grid spans [-f_max_hz, +f_max_hz]; see "Doppler span" above
  double f_max_hz;
  /// carrier wavelength, only used to report speeds
  double lambda_m;

  /* ---- TDD ---- */
  /// TDD period in seconds; replicas sit at multiples of 1/t_tdd_s
  double t_tdd_s;

  /* ---- slow-trend filter ----
  The clutter stage projects every bin out of the span of a few slow functions of t
  (see nr_ue_sensing_slow_basis()), so the samples here hold a target's tone minus its
  projection on them. The Doppler model has to be filtered the same way, or removing a
  target subtracts the part the filter already took out and leaves it behind as a
  ghost near 0 Hz. That part is large when a TDD replica of the target lands near 0 Hz.
  NULL / 0 when no trend was removed. */
  /// [trend_n_q][n_snap] orthonormal basis of the removed subspace, row major
  const double *trend_q;
  int trend_n_q;
} nr_tdd_obs_t;

/* Tuning. nr_tdd_cfg_default() fills in values that work for a 273 RB, 30 kHz,
   5 slot TDD carrier; anything scenario dependent is called out there. */
typedef struct {
  /// CA-CFAR false alarm probability per cell
  double pfa;
  /// guard cells either side of the cell under test, in bins and in Doppler
  int cfar_guard_bin;
  int cfar_guard_freq;
  /// training cells either side, beyond the guard band
  int cfar_train_bin;
  int cfar_train_freq;

  /// required fractional power drop at a replica for a peak to be accepted, Eq. (23)
  double gamma;
  /// replica orders +-1..n_sidelobes to test
  int n_sidelobes;
  /// half axes of the ellipse averaged at each replica, in bins and Doppler cells
  double ell_bin;
  double ell_freq;

  /// fine grid points per Doppler cell in the focused analysis
  int focus_zoom;
  /// half width of the focused Doppler search, in Doppler cells
  int focus_span;

  /// stop after this many accepted targets
  int max_targets;
  /// stop after this many CLEAN iterations, accepted or not
  int max_iter;
} nr_tdd_cfg_t;

void nr_tdd_cfg_default(nr_tdd_cfg_t *cfg);

/* ---- map helpers ---- */
bool nr_tdd_cmap_alloc(nr_tdd_cmap_t *map, int n_bins, int n_freq);
void nr_tdd_cmap_free(nr_tdd_cmap_t *map);
void nr_tdd_cmap_copy(const nr_tdd_cmap_t *src, nr_tdd_cmap_t *dst);
/// power of one cell, |C|^2
double nr_tdd_cmap_power(const nr_tdd_cmap_t *map, int bin, int f);

/* ---- point spread function ---- */

/* Dirichlet kernel of order A, Eq. (8): sin(A*pi*x) / (A*sin(pi*x)), continued to
   1 at the integers where both terms vanish. Only needed by the closed form. */
double nr_tdd_dirichlet(double x, int A);

/* The paper's closed form Doppler PSF, Eq. (10), for a strictly periodic pattern of
   R repetitions of M_TDD symbols of which the first M_DL carry data.
   Provided to validate nr_tdd_psf_doppler() against a synthetic uniform pattern;
   it is not what the detector uses. */
double nr_tdd_psf_doppler_closed(double m_over_Mp, int M_DL, int R, int M_TDD);

/* Doppler PSF of the actual sampling pattern, sum_i exp(-j*2*pi*x*t_i), evaluated
   at a frequency offset x in Hz. Normalised so that W(0) = 1. */
void nr_tdd_psf_doppler(const nr_tdd_obs_t *obs, double x_hz, double *re, double *im);

/* Range PSF of the actual frequency window, evaluated at a delay offset u in bins:
     W_R(u) = exp(j*2*pi*a*u/N) * sum_j win[j] * exp(j*2*pi*j*u/N)
   Normalised so that W_R(0) = 1. */
void nr_tdd_psf_range(const nr_tdd_obs_t *obs, double u_bins, double *re, double *im);

/* ---- pipeline stages ---- */

/* Complex range-Doppler map, the slow-time transform of obs->h onto the Doppler
   grid, i.e. Eq. (13) with the range transform already done upstream.
   O(n_bins * n_freq * n_snap). */
void nr_tdd_periodogram(const nr_tdd_obs_t *obs, nr_tdd_cmap_t *out);

/* Two dimensional cell averaging CFAR over |C|^2, Sec. IV-A.
   Local maxima whose power exceeds alpha * (mean of the training ring) are
   returned strongest first, where alpha = Nt * (pfa^(-1/Nt) - 1) for Nt training
   cells. Edges shrink the ring rather than wrapping: the Doppler axis is a finite
   span and the delay axis is truncated, so neither is periodic here.
   Returns the number written. */
int nr_tdd_cfar(const nr_tdd_cmap_t *map,
                const nr_tdd_cfg_t *cfg,
                int max_out,
                int bin_out[max_out],
                int freq_out[max_out],
                double snr_dB_out[max_out]);

/* Same CFAR, on a map that already holds power, laid out power[bin * n_freq + f].
   The AoA uses it on the antenna-averaged range-Doppler map, which only has power.
   nr_tdd_cfar() computes |C|^2 and calls this. */
int nr_tdd_cfar_power(const float *power,
                      int n_bins,
                      int n_freq,
                      const nr_tdd_cfg_t *cfg,
                      int max_out,
                      int bin_out[max_out],
                      int freq_out[max_out],
                      double snr_dB_out[max_out]);

/* Focused Fourier analysis around a grid peak, Sec. IV-B1.

   Doppler is refined exactly: the transform is re-evaluated on a grid
   focus_zoom times finer over +-focus_span cells, which costs nothing beyond a
   small direct sum and needs no interpolation.

   Delay is refined by fitting a parabola to the log power of the three bins around
   the peak. Refining it exactly would need the frequency domain lattice, which the
   history does not keep; the residual delay error leaves a small uncancelled
   remainder after removal, which the gamma of Eq. (23) has to tolerate.

   The amplitude is then the value the model must have at its own peak to reproduce
   the observed cell, so that subtracting it lands on zero there. */
void nr_tdd_focus(const nr_tdd_obs_t *obs,
                  const nr_tdd_cfg_t *cfg,
                  int bin0,
                  int f0,
                  nr_tdd_target_t *target);

/* Coherently subtract a target's full modelled response, replicas included,
   Eqs. (21) and (22).

   The response is separable, C_peak[q][f] = amp * W_R(q - d) * W_D(f - f0), so it
   costs n_bins + n_freq transform evaluations and one outer product rather than a
   full 2D synthesis. The replicas are not added by hand: they are already in
   W_D, because the sampling pattern that creates them is what W_D is built from. */
void nr_tdd_psf_subtract(const nr_tdd_obs_t *obs, const nr_tdd_target_t *target, nr_tdd_cmap_t *map);

/* Power features check, Eq. (23).

   For each replica order k in +-1..n_sidelobes, average the power over an ellipse
   centred where that replica of the target would be, before and after the removal.
   The peak is accepted when every tested replica dropped by at least gamma.
   Replicas falling outside the grid are skipped; if none is left to test the peak
   is rejected, since nothing was verified. */
bool nr_tdd_sidelobe_check(const nr_tdd_cmap_t *before,
                           const nr_tdd_cmap_t *after,
                           const nr_tdd_obs_t *obs,
                           const nr_tdd_cfg_t *cfg,
                           const nr_tdd_target_t *target);

/* Algorithm 1 end to end: periodogram, CFAR, then a CLEAN loop of refine, remove
   and check. Accepted peaks are kept removed so weaker targets can surface;
   rejected ones are restored and skipped.

   rejected_out, when not NULL, receives the candidates the loop threw out, each
   tagged with why. It must hold NR_TDD_MAX_REJECTED entries; anything beyond that is
   counted in *n_rejected but not written. Rejections carry grid coordinates only,
   not refined ones: a candidate that failed is not worth a focused analysis, and a
   marker on a map does not need sub-bin accuracy.

   Returns the number of accepted targets written to out[]. */
int nr_tdd_detect(const nr_tdd_obs_t *obs,
                  const nr_tdd_cfg_t *cfg,
                  nr_tdd_target_t *out,
                  int max_out,
                  nr_tdd_target_t *rejected_out,
                  int *n_rejected);

/*
 * WIRING, when the time comes
 *
 * nr_ue_sensing_range_doppler() already gathers the snapshots, normalises them for
 * pilot count, removes the clutter mean and runs the slow-time transform, but it
 * keeps only |C|^2 and it evaluates a deliberately narrow Doppler grid. To feed
 * this module it would have to
 *
 *   1. widen its grid to at least +-1/t_tdd_s while keeping the narrow one for the
 *      map it hands to the plotting path,
 *   2. fill an nr_tdd_obs_t with the same normalised, clutter removed snapshots it
 *      already has in hand, plus the taper from nr_ue_sensing_apply_hann() and the
 *      lattice start,
 *   3. call nr_tdd_detect() and report the survivors instead of, or alongside,
 *      nr_ue_sensing_find_peaks(), which has no notion of replicas and reports
 *      every one of them as a target.
 *
 * The taper is the one thing not currently kept anywhere: nr_ue_sensing_apply_hann()
 * applies it in place and forgets it. Either regenerate it from win_len or have that
 * function hand back the coefficients.
 */

#endif // __NR_UE_TDD_DETECT__H__
