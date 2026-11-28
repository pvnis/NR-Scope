#!/bin/bash
# 6e: replay every burst of the 5 Nov recording with the blind DM-RS dump (antennas 0 and 3), then the usage map.
R=/home/sens/isac_data/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq
for b in $(ls -d $R/b* | sort -t b -k2 -n); do
  [ -s $b/replay_blind_ch0/ue_dmrs.bin ] && grep -q DONE $b/replay_blind_ch0/nrscope.log && [ -s $b/replay_blind_ch3/ue_dmrs.bin ] && continue
  /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/replay4ch_blind.sh $b 0 & /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/replay4ch_blind.sh $b 3; wait
  echo "$(date +%T) $b" >> /home/sens/isac_data/ana/blind_replay.log
done
cd /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing && DIR=replay_blind_ch CH=0,3 nice python3 gs_usemap_v6e.py $R ana/usemap6e > ana/usemap6e.log 2>&1
echo ALLDONE >> /home/sens/isac_data/ana/blind_replay.log
