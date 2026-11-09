#!/bin/bash
# Replay every burst on all 4 antennas with the dumps of 2026-11-09 (blind DM-RS, TRS, window-shift GSW1), 2 at a time.
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq
for b in $(ls -d $R/b* | sed "s#.*/b##" | sort -n); do
  [ "$b" = 9 ] && continue
  /home/sens/field/replay4ch_blind.sh $R/b$b 0 & /home/sens/field/replay4ch_blind.sh $R/b$b 1; wait
  /home/sens/field/replay4ch_blind.sh $R/b$b 2 & /home/sens/field/replay4ch_blind.sh $R/b$b 3; wait
  echo "$(date +%T) b$b" >> /home/sens/field/ana/gsw_replay.log
done
echo ALLDONE >> /home/sens/field/ana/gsw_replay.log
