#!/bin/bash
# 6f: replay every burst again with the blind DM-RS + TRS/CSI-RS raw dumps (antennas 0 and 3).
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq
for b in $(ls -d $R/b* | sort -t b -k2 -n); do
  /home/sens/field/replay4ch_blind.sh $b 0 & /home/sens/field/replay4ch_blind.sh $b 3; wait
  echo "$(date +%T) $b" >> /home/sens/field/ana/trs_replay.log
done
echo ALLDONE >> /home/sens/field/ana/trs_replay.log
