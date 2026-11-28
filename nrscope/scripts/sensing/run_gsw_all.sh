#!/bin/bash
# Replay every burst on all 4 antennas with the dumps of 2026-11-09 (blind DM-RS, TRS, window-shift GSW1), 2 at a time.
R=/home/sens/isac_data/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq
for b in $(ls -d $R/b* | sed "s#.*/b##" | sort -n); do
  [ "$b" = 9 ] && continue
  /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/replay4ch_blind.sh $R/b$b 0 & /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/replay4ch_blind.sh $R/b$b 1; wait
  /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/replay4ch_blind.sh $R/b$b 2 & /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/replay4ch_blind.sh $R/b$b 3; wait
  echo "$(date +%T) b$b" >> /home/sens/isac_data/ana/gsw_replay.log
done
echo ALLDONE >> /home/sens/isac_data/ana/gsw_replay.log
