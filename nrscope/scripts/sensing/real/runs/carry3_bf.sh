#!/bin/bash
# carry3: coherent 4-antenna beamforming detection (TBD cube, map_peaks) on the BA_G1 DMRS estimates, fused with the TRS
cd "$(dirname "$(readlink -f "$0")")/.."   # sensing/real
W="tracker.t_window_unix=[1795790234.976,1795790319.502]"; D="streams.dmrs.est=\"/mnt/mydrive/field/real/lidar_carry3/est_dmrs_g1\""
for wh in none map; do python3 tbd_track.py --config configs/lidar_carry3.json --caps c3_g1_t12 --name c3_bf_$wh --as-rx2 --set "$W" "$D" "tbd.output=\"map_peaks\"" "tbd.streams=[\"dmrs\"]" tbd.map_threshold=40 tbd.max_dets=4 "tbd.whiten=\"$wh\""; done
echo ALLDONE
