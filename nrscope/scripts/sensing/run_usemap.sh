#!/bin/bash
# Usage map of the current pipeline: light stream (TRS + CSI-RS + DM-RS, per-snapshot TRS reference) with and without
# the tablet, no injection; so_light.py writes every candidate symbol and its fate (SIDE), sensing_offline writes which
# symbols survive its alignment / random drop and which go into each map (NR_SENSING_TRACE).
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=/mnt/mydrive/field/usemap7; mkdir -p $O; cd /home/sens/field
one() {
  b=$1
  for m in tab notab; do
    SIDE=$O/b${b}_$m.side MINPRB=48 USE_DMRS=1 USE_CSI=1 MODE=$m nice python3 so_light.py $R/b$b $O/b${b}_$m.bin > $O/b${b}_$m.conv.log 2>&1
    NR_SENSING_TRACE=$O/b${b}_$m.trace nice $SO so/cfg_trs.yaml $O/b${b}_$m.bin --dump $O/b${b}_$m.map2d.csv > $O/b${b}_$m.run.log 2>&1
    rm -f $O/b${b}_$m.bin
  done
  echo "$(date +%T) b$b" >> $O/progress.log
}
for b in $(seq 0 25); do
  while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 4 ]; do sleep 5; done
  one $b &
  sleep 2
done
wait; echo ALLDONE >> $O/progress.log
