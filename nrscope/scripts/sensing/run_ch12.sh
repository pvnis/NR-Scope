#!/bin/bash
# 6f step 3: antennas 1 and 2 of three bursts with the blind/TRS/CSI-RS dumps (the CSI-RS span test needs 4 antennas).
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq
for b in b2 b9 b21; do /home/sens/field/replay4ch_blind.sh $R/$b 1 & /home/sens/field/replay4ch_blind.sh $R/$b 2; wait; echo "$(date +%T) $b" >> /home/sens/field/ana/ch12_replay.log; done
echo ALLDONE >> /home/sens/field/ana/ch12_replay.log
