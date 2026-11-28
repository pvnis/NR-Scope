#!/bin/bash
# carry3 (27 Nov): wait for decode + grid, then the default chain and track-before-detect, scored on the labelled walker in the FOV window
for i in $(seq 900); do grep -q "^END" /tmp/carry3_dec.log && break; sleep 5; done
grep -q "grid=.*DONE" /tmp/carry3_dec.log || { echo "decode/grid not done: $(cat /tmp/carry3_dec.log)"; exit 1; }
cd "$(dirname "$(readlink -f "$0")")/.."   # sensing/real
W="tracker.t_window_unix=[1795790234.976,1795790319.502]"
python3 real_new.py ~/GIOVANNI/real/cubotron_27_11/lidar_carry3 --t-rec "2026-11-27 14:37:14.975" --iq-offset 0 --name c3_fused --set "$W"
python3 tbd_track.py --config configs/lidar_carry3.json --caps c3_fused --name c3_tbd --set "$W"
echo ALLDONE
