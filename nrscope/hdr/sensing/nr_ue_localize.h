#ifndef __NR_UE_LOCALIZE__
#define __NR_UE_LOCALIZE__

#include "nrscope/hdr/sensing/nr_ue_aoa.h"

/* Smallest cos(beta / 2) the bistatic velocity is computed for, beta the bistatic
angle at the target. The Doppler only sees the velocity along the bisector, scaled
by cos(beta / 2), so dividing it back out blows up as the target nears the baseline
between gNB and receiver (beta -> 180 deg, forward scatter), where the Doppler says
almost nothing about the speed. 0.17 is beta = 160 deg, a gain of about 6. */
#define NR_LOC_MIN_COS_HALF_BETA 0.17

typedef struct {
  float pos_x;
  float pos_y;
  /// distance from the receiver to the target, m
  float rho;
  int aoa_index;

  /* BISTATIC VELOCITY, see nr_ue_target_position(). Valid when vel_valid. */
  /// bistatic angle at the target, between the directions to the gNB and to the receiver, deg
  float beta_deg;
  /* Bearing of the inward bistatic bisector, the direction halfway between "towards
  the gNB" and "towards the receiver", in the same frame as the positions, deg */
  float bisector_deg;
  /* Velocity component along that bisector, m/s: positive moving towards the gNB and
  the receiver together, as the map's speed is positive approaching. The only part of
  the velocity one receiver's Doppler measures; the component across the bisector,
  along the ellipse of constant bistatic range, leaves the Doppler at zero. */
  float v_bisector_ms;
  /// false when the target is too close to the baseline, see NR_LOC_MIN_COS_HALF_BETA
  int vel_valid;
} nr_sensing_target_t;

/* Position, bistatic angle and bistatic velocity of each AoA cell.
   in        : the AoA cells; range_m is the excess path length over the direct path,
               speed_ms the map's speed (f * lambda / 2)
   pos_gnb,
   pos_ue    : gNB and receiver positions, m, in any planar frame
   boresight_deg : bearing of the receive array's normal in that frame; a cell's
               direction from the receiver is boresight_deg + angle_deg
   Returns the targets written to out, at most n_in: a cell with no geometric solution
   (non-positive excess range, or a ray that cannot meet the ellipse) is skipped. */
int nr_ue_target_position(const nr_sensing_aoa_t *in,
                          int n_in,
                          const float pos_gnb[2],
                          const float pos_ue[2],
                          float boresight_deg,   // global bearing of the array normal
                          nr_sensing_target_t *out);

#endif // __NR_UE_LOCALIZE__
