/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * WHY LAST, AND NOT BEFORE THE DELAY-DOPPLER TRANSFORM
 *
 * The joint two-domain framework estimates the AoA first, beamforms to separate the
 * angles, and only then transforms delay and Doppler. That ordering rests on the
 * angles being resolvable: the beamformer has to push the arrivals from other
 * directions down into the noise, or the delay-Doppler model it hands on is not the
 * one the transform assumes.
 *
 * With four antennas at half-wavelength spacing that condition does not hold. The
 * angular resolution is 30 degrees at the array normal,
 * against about 4.6 m in range and 0.2 m/s
 * in Doppler over a 200 ms window. Angle is by far the coarsest axis here, so it is
 * the last one that should be asked to separate anything.
 *
 * Run last, the same four antennas are worth much more:
 *   - the delay and Doppler transforms have already integrated the cell coherently,
 *     so the spatial snapshot is taken at the best SNR the pipeline ever reaches;
 *   - clutter removal has already taken the direct path out of the cell, so the
 *     vector is not dominated by the LoS steering;
 *   - the angle attaches to the cell by construction, so there is no parameter
 *     grouping step and no chance of pairing an angle with the wrong target.
 */

#ifndef __NR_UE_AOA__H__
#define __NR_UE_AOA__H__

#include "nrscope/hdr/sensing/nr_ue_map.h"

/// Rx antennas the spatial stage handles, the sensing library's chain ceiling
#define NR_AOA_MAX_ANT NR_SENSING_MAX_RX

/* Element spacing of the receive ULA, in wavelengths, and the per-chain phase
offsets removed before the steering vector is applied, in degrees.

The calibration is the blocking assumption of this whole module. The four X410 chains
share a sample clock, but the daughterboard LOs can settle at an arbitrary relative
phase on every retune, and until that offset is known and constant the steering vector
below describes an array that does not exist. Measure it from the direct path, which
is static and strong: arg(h_a[b_direct]) - arg(h_0[b_direct]) should be constant
across slots. If it is, it goes here. If it drifts within a run, this needs to become
a per-map estimate instead of a constant. */
#define NR_AOA_ELEMENT_SPACING 0.5
#define NR_AOA_CAL_DEG {0.0, 0.0, 0.0, 0.0}

/* Points of the spatial spectrum, uniform in sin(theta) over [-1, 1].

Uniform in sin(theta) and not in theta because the array response is a DFT over
sin(theta): the beamwidth is constant on that axis, so a uniform angle grid would
oversample near +-90 degrees and undersample nowhere useful. At Mr = 4 the beam is
1.0 wide in sin(theta), so 181 points oversample it ninety-fold; the cost is a few
thousand operations per cell and the spectrum comes out smooth enough to interpolate. */
#define NR_AOA_GRID 181

/* Local maxima examined before the resolution-cell merge picks the survivors. Only a
bound on the working buffer: a map with more maxima than this is noise, and the
strongest of them are kept. */
#define NR_AOA_MAX_CANDIDATES 256

/* When two local maxima are taken to be the same object: same delay bin to within
NR_AOA_SEP_BINS, and same angle to within NR_AOA_SEP_ANGLE_DEG.

Doppler is deliberately not part of the test, and the reason is geometric rather than
numerical. A delay bin is a bistatic ellipse and an angle is a ray from the array;
together they pin one point in space. Two cells agreeing on both are therefore looking
at one place, and one place holds one object, whatever the Doppler axis says about it.

The angle tolerance is well inside the 30 degree beamwidth of four antennas, so it can
never merge two directions the array could actually have told apart. */
#define NR_AOA_SEP_BINS 2
#define NR_AOA_SEP_ANGLE_DEG 10.0

/* Floor on the delay bins excluded above the direct path, used when the taper width
is unknown and as a lower bound on the width it implies.

The guard is normally computed from the grant, 2*idft_size/win_len bins; at a full
band grant that is about 2.5, so this floor binds only for the widest allocations.
Three bins is 7.3 m of excess range that can never be reported, which is the real
cost of the stage: a target closer than that to the direct path in bistatic range is
buried in its mainlobe and would be read as an artefact anyway. */
#define NR_AOA_LOS_GUARD_MIN 3

/* Spatial snapshot at one cell of the map: the slow-time transform evaluated at a
   single (bin, Doppler) point, once per Rx antenna, with the chain calibration
   removed.

   Each entry is one number carrying the whole coherent gain of the cell.

   slow  : [n_ant] conditioned slow-time samples, as nr_ue_sensing_range_doppler()
           exports them; every antenna must have been gathered over the same
           snapshots, which holds because they see the same reference symbols
   bin   : delay bin of the cell
   f_hz  : Doppler of the cell, on the grid the map reports
   x_re,
   x_im  : receive at least n_ant entries each */
void nr_ue_aoa_cell_vector(const nr_sensing_slowtime_t *slow,
                           int n_ant,
                           int bin,
                           double f_hz,
                           double *x_re,
                           double *x_im);

/* Angle of arrival from one spatial snapshot, by a periodogram over sin(theta).

   P(theta) = |a_R^H(theta) x|^2 with the ULA steering vector of 38.211 geometry,
   a_R(theta)[m] = exp(j*2*pi*(d/lambda)*m*sin(theta)), maximised over the grid and
   then refined by a parabolic fit through the peak and its two neighbours.

   quality_dB receives the peak over the mean of the spectrum. For a single plane
   wave that ratio is Mr, so 6 dB at four antennas; a cell holding two arrivals or
   just noise comes out well below it. Pass NULL to skip it.

   Returns the angle in degrees, 0 at the array normal, positive towards increasing
   element index. */
double nr_ue_aoa_estimate(const double *x_re, const double *x_im, int n_ant, double *quality_dB);

/* One cell to measure the angle at. bin is the delay bin (can be fractional, since the
   TDD detector refines it), f_hz the Doppler, snr_dB how far above the noise it is. */
typedef struct {
  double bin;
  double f_hz;
  double snr_dB;
} nr_aoa_cell_t;

/* Angles for the targets of one map.

   Where the targets come from:
     - cells != NULL (--sensing-tdd-detect on): the TDD detector's targets. Their
       Doppler is the true one, even beyond the edge of the map, because the detector
       already told the real peaks from the replicas.
     - cells == NULL (detector off): the AoA finds the peaks itself, with CFAR on the
       antenna-averaged map. Speeds are then only right inside the map's window.

   Then, for each cell, strongest first: skip it if it sits on the direct path, measure
   its angle from the complex samples of the antennas, and drop it if a stronger cell
   already reported the same place (same delay bin and same angle).

   map     : the antenna-averaged map, used for the peak search and the axes
   slow    : [n_ant] slow-time samples, one per Rx antenna, all from the same window
   n_ant   : 2 or more; below that there is no spatial information at all
   cells   : targets to measure, strongest first, or NULL to search the map
   n_cells : entries in cells
   out     : receives at most NR_SENSING_AOA_MAX entries

   Returns how many were written. */
int nr_ue_aoa_process(const nr_sensing_map_t *map,
                      const nr_sensing_slowtime_t *slow,
                      int n_ant,
                      const nr_aoa_cell_t *cells,
                      int n_cells,
                      nr_sensing_aoa_t out[NR_SENSING_AOA_MAX]);

#endif // __NR_UE_AOA__H__
