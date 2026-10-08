/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_localize.h"
#include "nrscope/hdr/sensing/nr_ue_aoa.h"

int nr_ue_target_position(const nr_sensing_aoa_t *in,
                          int n_in,
                          const float pos_gnb[2],
                          const float pos_ue[2],
                          float boresight_deg,
                          nr_sensing_target_t *out)
{

  /*
  Knowing the bistatic excess range delta_R and AoA for a potential target, we can compute its 2D position.
  Let R and T the positions of the UE and the gNB, "a" the vector R - T, and L = norm(a).
  Define D = L + delta_R, with delta_R the bistatic range computed from the periodogram.

  The position of target p = R + pho * u(theta), u(theta) = [cos(), sin()]^T
  Now pho = (D^2 - L^2) / (2 (D + a^T u(theta)))
  */

  const double ax = pos_ue[0] - pos_gnb[0];
  const double ay = pos_ue[1] - pos_gnb[1];
  const double L  = sqrt(ax * ax + ay * ay);

  int n_out = 0;
  for (int i = 0; i < n_in; i++) {
    const double dR = in[i].range_m;
    if (dR <= 0.0) continue;

    /* Array frame to global: angle_deg is measured off the array normal, the
    boresight says where that normal points. */
    const double th = (boresight_deg + in[i].angle_deg) * M_PI / 180.0;
    const double ux = cos(th);
    const double uy = sin(th);

    const double D   = L + dR;
    const double den = D + ax * ux + ay * uy;
    if (den <= 0.0) continue;

    const double rho = (D * D - L * L) / (2.0 * den);
    const double px  = pos_ue[0] + rho * ux;
    const double py  = pos_ue[1] + rho * uy;

    /*
    BISTATIC VELOCITY

    The bistatic range is R_b = |p - T| + |p - R|. Its rate, for a target moving with
    velocity v, is
        dR_b/dt = v . (e_T + e_R),   e_T = (p - T)/|p - T|,  e_R = (p - R)/|p - R|,
    the unit vectors from the gNB and from the receiver to the target. Their sum has
    length 2 cos(beta / 2), beta the bistatic angle between them, and points along the
    outward bisector b. So
        dR_b/dt = 2 cos(beta / 2) (v . b).
    The map's speed is f * lambda / 2 with f = -(dR_b/dt) / lambda (an echo whose path
    shortens advances in phase), i.e. speed = -(dR_b/dt) / 2, and therefore
        v . (-b) = speed / cos(beta / 2).
    -b is the inward bisector, so this is the velocity towards the two nodes together,
    positive approaching like the map's speed. With the gNB on the receiver, beta = 0
    and it is the speed itself, which is what the map's convention assumed all along.

    The component across the bisector, along the ellipse of constant R_b, does not move
    the Doppler at all and stays unknown from one receiver.
    */
    const double tx  = px - pos_gnb[0], ty = py - pos_gnb[1];
    const double nt  = sqrt(tx * tx + ty * ty);
    // e_R is u itself, the target lying along the ray from the receiver
    const double ex  = ux + (nt > 0.0 ? tx / nt : 0.0);
    const double ey  = uy + (nt > 0.0 ? ty / nt : 0.0);
    const double chb = sqrt(ex * ex + ey * ey) / 2.0; // |e_T + e_R| = 2 cos(beta / 2)

    nr_sensing_target_t o = {
        .aoa_index = i,
        .rho       = (float)rho,
        .pos_x     = (float)px,
        .pos_y     = (float)py,
        .beta_deg  = (float)(2.0 * acos(fmin(1.0, chb)) * 180.0 / M_PI),
    };
    if (nt > 0.0 && chb >= NR_LOC_MIN_COS_HALF_BETA) {
      o.bisector_deg  = (float)(atan2(-ey, -ex) * 180.0 / M_PI);
      o.v_bisector_ms = (float)(in[i].speed_ms / chb);
      o.vel_valid     = 1;
    }
    out[n_out++] = o;
  }
  return n_out;
}
