#ifndef __NR_UE_LOCALIZE__
#define __NR_UE_LOCALIZE__

#include "nrscope/hdr/sensing/nr_ue_aoa.h"

typedef struct {
  float pos_x;
  float pos_y;
  float rho;
  int aoa_index;

} nr_sensing_target_t;

int nr_ue_target_position(const nr_sensing_aoa_t *in,
                          int n_in,
                          const float pos_gnb[2],
                          const float pos_ue[2],
                          float boresight_deg,   // global bearing of the array normal
                          nr_sensing_target_t *out);

#endif // __NR_UE_LOCALIZE__