#!/bin/bash
# carry3: whitened 4-antenna beamforming with lower detection limits (find the walker as often as OMP, keep the sharper angle)
cd "$(dirname "$(readlink -f "$0")")/.."   # sensing/real
W="tracker.t_window_unix=[1795790234.976,1795790319.502]"; D="streams.dmrs.est=\"/mnt/mydrive/field/real/lidar_carry3/est_dmrs_g1\""
for th in 25 15; do python3 tbd_track.py --config configs/lidar_carry3.json --caps c3_g1_t12 --name c3_bfw_t$th --as-rx2 --set "$W" "$D" "tbd.output=\"map_peaks\"" "tbd.streams=[\"dmrs\"]" tbd.map_threshold=$th tbd.max_dets=6 "tbd.whiten=\"map\""; done
echo ALLDONE
