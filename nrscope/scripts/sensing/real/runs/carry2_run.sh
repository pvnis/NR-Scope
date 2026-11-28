#!/bin/bash
# carry2 (27 Nov): default chain + track-before-detect, scored on the labelled walker inside the FOV window
cd "$(dirname "$(readlink -f "$0")")/.."   # sensing/real
W="tracker.t_window_unix=[1795784394.942,1795784462.483]"
python3 real_new.py ~/GIOVANNI/real/cubotron_27_11/lidar_carry2 --set "$W"
python3 tbd_track.py --config configs/lidar_carry2.json --caps fused --name tbd_both --set "$W"
echo ALLDONE
