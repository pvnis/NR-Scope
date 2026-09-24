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

    out[n_out++] = (nr_sensing_target_t){
        .aoa_index = i,
        .rho = (float)rho,
        .pos_x = (float)(pos_ue[0] + rho * ux),
        .pos_y = (float)(pos_ue[1] + rho * uy),
    };
  }
  return n_out;
}
