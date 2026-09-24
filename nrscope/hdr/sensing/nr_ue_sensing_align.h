/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 *
 * What we were doing before this file
 *
 * 1. The UE decides where each OFDM symbol starts: the FFT window. Bin 0 of a delay
 *    profile means "at the window start", not "at the gNB". Every path is measured from it.
 * 2. The gNB and UE clocks never agree exactly, so the signal slowly drifts inside the
 *    window. OAI corrects that once per frame on the SSB (nr_adjust_synch_ue()): it reads
 *    one sample more or less, and the window moves by one sample.
 * 3. A map adds up hundreds of symbols over a short time. After a correction
 *    every path sits one bin (2.44 m) closer or further, and the whole symbol is rotated
 *    by one constant phase (2*pi*1638/4096 = 144 deg per sample at 273 RB, because the
 *    pilots are counted from the carrier edge, not from DC). Inside one map the static
 *    scene then exists in two versions. Subtracting the mean cannot remove the difference,
 *    and because it switches on and off in time it spreads over every speed.
 *
 *    Now: for each symbol before the taper and IDFT:
 *   a. Compare the pilots p_i with a reference R_i of the static scene, trying delays d a
 *      few bins around the last shift found:
 *          c(d) = sum_i p_i conj(R_i) exp(+j 2 pi i d / N)
 *      If the symbol is the reference delayed by D bins and rotated by phi, then c(d)
 *      peaks exactly at d = D and its phase there is phi.
 *   b. Undo it: p_i <- p_i exp(+j 2 pi i D / N) exp(-j phi). On the pilots a delay is an
 *      exact phase ramp, so nothing is interpolated and nothing is lost at the edges.
 *   c. Update the reference slowly with the aligned symbol.
 *
 *   Amplitude (NR_SENSING_ALIGN_AMPLITUDE). The correlation at the estimated delay is not
 *   only a phase, it is a complex gain. Model the aligned symbol as the reference times one
 *   complex scalar, p_i = g R_i (plus movers and noise). Then
 *          c(D) = sum_i p_i conj(R_i) = g sum_i |R_i|^2 = g e_r
 *   so the least squares gain is g = c(D) / e_r: its phase is phi above, and its magnitude
 *          |g| = |c(D)| / e_r
 *   says how much louder or quieter the whole static scene is in this symbol than in the
 *   reference. Step b divides by it too:
 *          p_i <- p_i exp(+j 2 pi i D / N) exp(-j phi) / |g|
 *
 *   Why it is needed. Measured on an outdoor capture (357 snapshots, 131 ms, map 267.5): after
 *   delay and phase were aligned, the direct path still changed its amplitude by about 14 %
 *   (std/mean) from symbol to symbol, mostly between slots and across scheduling gaps, while
 *   the delay jitter left was 0.01 bin and the phase jitter 0.013 rad. The clutter stage in
 *   nr_ue_map.c removes a constant or a slow polynomial per range bin, and a gain that jumps
 *   is neither: the static scene times (1 + eps_m) leaves eps_m * scene in every snapshot,
 *   and eps_m, being irregular in time, has a flat spectrum. It lands on every Doppler bin
 *   at a level proportional to the static power of that range bin, which is the high floor
 *   around the direct path. On that capture the floor was 33-39 dB above thermal at 0-5 m
 *   of excess range, 11-19 dB at 15-49 m; dividing each snapshot by its gain brought it to
 *   12-16 dB and 7-12 dB.
 *
 *   What it costs. |g| is fitted on the whole symbol, which the direct path dominates by tens
 *   of dB, so a mover barely moves it and keeps its own Doppler. It is only scaled along with
 *   the scene, by the same few percent. But if the mover is what changes the direct path (a
 *   person walking through the link), that modulation is now moved onto the mover and the
 *   rest of the symbol instead of staying on the direct path: a small effect on the target,
 *   a large one on the floor.
 *
 *   The gain is one scalar per symbol, estimated once and applied to every Rx antenna, like
 *   the delay and the phase, so the amplitude ratios between antennas are kept for AoA.
 *   Since the reference is built from normalised symbols, it keeps the level of the symbol it
 *   started from: every snapshot of a stream is brought to that level.
 *
 *   The shift and phase are estimated on the first Rx antenna that sees a symbol and
 *   reused unchanged for the other antennas, so the phase differences between antennas,
 *   which AoA relies on, are kept.
 * 
 *   Consequence: range and Doppler become relative to the static scene, whose strongest
 *   part is the LOS. The LOS stays at the same bin and at 0 m/s, which is the bistatic
 *   convention anyway. A symbol that does not match the reference at all is left out;
 *   if that goes on, the reference restarts from the current symbol.
 *
 * If the normalized correlation gives a bad value, it means that the pilots
 * are not "a reference", something bad happened in the scene. What could cause a low correlation are: (i) the precoding change
 * when the gNB changes beam, so channel estimation is for a new Hw' and now Hw. (ii) A symbol with a very poor estimate
 * such as with low SNR. (iii) LOS obstructated and so the dominant path is not the LOS. In the end, with enough symbols
 * over time, its better to drop some of them than take them because they will make doppler spread worse.
 * 
 * What is this reference: suppose we start a stream (port, layer, step, offset), then the reference becomes the first symbol.
 * If not, we do the correlation and we correct accordingly. Now R is a vector of the same size of the symbols. One R per stream,
 * shared with all Rx antennas. R is updated with a exponential moving average (EMA) as 
 * R <- R + w (p_corrected - R). R is a weighted average of the last approx 64 symbols, most recent weights more.
 * Why its good: Static paths have the same value symbol after symbol, so they add up and survive the averaging. 
 * Movers rotate in phase, so they average toward zero. After a few dozen symbols R is an estimate of the static scene, 
 * cleaner than any single symbol, which is why it's a good thing to correlate against.
 * 
 * --------------------------------------------------------------------------------------------------
 * 
 * Grouping by grant (NR_SENSING_GROUP_BY_GRANT). When a map is built, only the snapshots
 * with the most common grant (same size, same start) are kept, so all of them have the
 * same point spread function. The Doppler resolution is unchanged; the SNR drops by the
 * share of snapshots left out.
 */

#ifndef __NR_UE_SENSING_ALIGN__H__
#define __NR_UE_SENSING_ALIGN__H__

#include <stdbool.h>
#include <stdint.h>
#include "nrscope/hdr/sensing/sensing_defs.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"

/// 1: align every symbol's pilots to the reference of the static scene (see above)
#define NR_SENSING_ALIGN 1
/// 1: build each map from the snapshots of its most common grant only (see above)
#define NR_SENSING_GROUP_BY_GRANT 0

/// how far around the previous shift the delay is searched, in bins (a timing correction is 1 sample)
#define NR_SENSING_ALIGN_SEARCH_BINS 3.0
/// coarse search step, in bins; the peak is then refined below this
#define NR_SENSING_ALIGN_STEP_BINS 0.25
/// weight of a new aligned symbol in the reference; 1/64 follows the scene over ~20 slots
#define NR_SENSING_ALIGN_REF_WEIGHT (1.0f / 64.0f)
/// normalised correlation (1 = identical scene) below which a symbol is left out
#define NR_SENSING_ALIGN_MIN_CORR 0.5
/// symbols in a row below that before the reference restarts from the current one
#define NR_SENSING_ALIGN_RESET_AFTER 64
/// fewest pilots shared with the reference for an estimate to be trusted
#define NR_SENSING_ALIGN_MIN_OVERLAP 64
/// 1: also divide every aligned symbol by the magnitude of its gain against the reference (see above)
#define NR_SENSING_ALIGN_AMPLITUDE 1
/* Range |g| is clamped to before dividing by it, as a factor on the reference level (0.5 and
2 are -6 and +6 dB). The fluctuation it corrects is a few tens of percent. A gain outside
this range with a correlation still above NR_SENSING_ALIGN_MIN_CORR means the scene itself
changed, not its level: dividing fully would amplify the noise of a faded symbol, or push a
strong one into the int16 saturation of the pilots. */
#define NR_SENSING_ALIGN_AMP_MIN 0.5
#define NR_SENSING_ALIGN_AMP_MAX 2.0

/// What the alignment found for one symbol
typedef struct {
  /// delay of the symbol relative to the reference, in bins
  double delay_bins;
  /// phase of the symbol relative to the reference, in radians
  double phase_rad;
  /* magnitude of the gain relative to the reference, |c(D)| / e_r, after the clamp; the
  pilots were divided by it when NR_SENSING_ALIGN_AMPLITUDE is set, 1 otherwise */
  double amp;
  /// normalised correlation with the reference, 1 for an identical scene
  double corr;
  /// the correction was applied to the pilots
  bool applied;
  /// the symbol did not match the reference and should be left out
  bool dropped;
} nr_sensing_align_t;

/* Align one symbol's lattice vector to its stream's reference. 
Call it after nr_ue_sensing_extract_lattice() and before the taper.

   t_sample : the symbol's time, from nr_ue_sensing_symbol_time(). A second call with the
              same time and stream (another Rx antenna) reuses the first correction.
   stream   : only ports and layer are read; the comb comes from lat
   lat      : the symbol's lattice, as extracted
   pilots   : lattice indexed vector, modified in place
   out      : what was found, may be NULL
   Returns false when the symbol should be left out. */
bool nr_ue_sensing_align_symbol(uint64_t t_sample,
                                nr_sensing_stream_t stream,
                                const pilot_lattice_t *lat,
                                cf_t *pilots,
                                nr_sensing_align_t *out);

/* Keep, among the n snapshots idx[] points to in hist, only those whose grant (pilot count
   and first pilot) is the most common one. Order is kept. A tie goes to the most recent
   grant. Returns how many are left. */
int nr_ue_sensing_group_by_grant(const nr_sensing_history_t *hist, int idx[], int n);

#endif
