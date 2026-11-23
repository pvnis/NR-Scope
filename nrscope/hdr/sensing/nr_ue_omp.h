/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/* Mask-aware 2D-OMP back end (sensing.backend: omp), an alternative to the TDD detector's CLEAN on the periodogram.

Why: on a commercial cell every snapshot measures only the PRBs the gNB scheduled. A periodogram treats the holes
as zeros, so each moving target spreads its own gap sidelobes over range and Doppler, where they become false
detections or hide weak targets. OMP models a path exactly as the data hold it: in snapshot m its delay response
is the masked clutter kernel K_m(b - u) (the transform of that snapshot's taper over the pilots it really measured,
nr_ue_sensing_clutter_kernel_masked) times the snapshot's gain and exp(j 2 pi f t_m). Greedily: the strongest
atom of the residual's correlation map, sub-bin refinement, joint least squares over every atom found, repeat,
until the strongest cell falls under threshold x noise (noise fixed once, from the first map). The antennas are
fitted jointly (one atom, one amplitude per antenna), so their phases, which the AoA reads, are untouched.

Measured on the 5 Nov commercial-cell recording (OPEN_ISSUES 4b4, Python prototype): at an equal false-detection
rate it finds 11-22 points more injected targets without a heavy user, 3-5 points more with one. */

#ifndef NR_UE_OMP_H
#define NR_UE_OMP_H

#include "nrscope/hdr/sensing/nr_ue_map.h"
#include "nrscope/hdr/sensing/nr_ue_tdd_detect.h"

/// paths an OMP map looks for at most, when sensing.omp_max_paths is 0
#define NR_OMP_MAX_PATHS_DEFAULT 16
/* strongest cell over the map's noise (median / ln 2, exponential cells) that still counts as a path, dB, when
sensing.omp_threshold_db is 0. 11.4 dB = a factor 13.8: a false alarm probability of 1e-6 per cell. */
#define NR_OMP_THRESHOLD_DB_DEFAULT 11.4
/* a candidate this far below the strongest path of the map is the model's fitting error, not a path: the kernels
are float and the refinement finite, which leaves roughly -60 dB of a path behind */
#define NR_OMP_DYNAMIC_RANGE_DB 50.0

typedef struct {
  int max_paths;       ///< 0: NR_OMP_MAX_PATHS_DEFAULT
  double threshold_db; ///< 0: NR_OMP_THRESHOLD_DB_DEFAULT
  int local_noise;     ///< 1: each range bin's own floor too (its Doppler median): a bin smeared over Doppler (phase noise on the direct path) is not a path
} nr_omp_cfg_t;

/* Paths of one map.

   slow      : n_ant slow-time blocks from nr_ue_sensing_range_doppler(), the residual after clutter removal, with
               snap[] and gain[] filled (see nr_sensing_slowtime_t). Antennas whose times differ from slow[0]'s
               are left out of the fit.
   n_freq, f_max_hz : the map's Doppler grid, f_j = -f_max + j 2 f_max / (n_freq - 1)
   lambda_m  : carrier wavelength, for speed = f lambda / 2
   trend_q   : the slow-trend basis the clutter stage removed ([trend_n_q][n_snap], may be NULL), so each atom is
               projected out of it exactly as the data were
   out       : up to max_out paths, strongest first: bin, range_m, f_hz, speed_ms, snr_dB (over the noise),
               iteration, verdict NR_TDD_VERDICT_TARGET. Paths within one Doppler resolution cell of zero are not
               looked for (static residue).
   Returns the number of paths. */
int nr_ue_omp_detect(const nr_sensing_slowtime_t *slow,
                     int n_ant,
                     int n_freq,
                     double f_max_hz,
                     double lambda_m,
                     const double *trend_q,
                     int trend_n_q,
                     const nr_omp_cfg_t *cfg,
                     nr_tdd_target_t *out,
                     int max_out);

#endif
