#!/bin/bash
cd "$(dirname "$(readlink -f "$0")")/.."   # sensing/real
W="tracker.t_window_unix=[1795790234.976,1795790319.502]"
E="streams.dmrs.so_light_env={\"NO_TRS_OUT\":\"1\",\"MODE\":\"tab\",\"BA\":\"1\",\"BA_PH\":\"wide\",\"BA_G1\":\"1\"}"; D="streams.dmrs.est=\"/mnt/mydrive/field/real/lidar_carry3/est_dmrs_g1\""
for th in 15 12; do python3 real_fuse.py --config configs/lidar_carry3.json --name c3_g1_t$th --set "$W" "$E" "$D" streams.dmrs.omp_threshold_db=$th; done
echo ALLDONE
