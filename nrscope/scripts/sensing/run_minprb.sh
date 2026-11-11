#!/bin/bash
# Narrow DM-RS symbols: leave out DM-RS symbols with fewer than 24 / 48 PRBs measured; light stream with and without
# the tablet, 3-target raw injection, random_drop off, all bursts.
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=/mnt/mydrive/field/minprb; mkdir -p $O; cd /home/sens/field; T="25,2.5,-20,15;60,-5,-22,-25;100,8,-25,30"
one() {
  b=$1
  for m in tab notab; do
    MINPRB=0,24,48 INJ="$T" USE_DMRS=1 USE_CSI=1 MODE=$m nice python3 so_light.py $R/b$b $O/b${b}_$m.bin > $O/b${b}_$m.log 2>&1
    for n in 0 24 48; do
      d=$O/m_b${b}_${m}_$n; mkdir -p $d; nice $SO so/cfg_trs.yaml $O/b${b}_$m.bin.m$n --dump $d/map2d.csv > $d/run.log 2>&1
      cp $O/b${b}_$m.bin.truth.csv $d/map2d.csv.truth.csv; rm -f $d/map2d.csv.snap.csv $O/b${b}_$m.bin.m$n
    done
  done
  REF_EP=snap INJ="$T" WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr_ep.py $R/b$b $O/b${b}_tabref.bin > $O/b${b}_tabref.log 2>&1
  d=$O/m_b${b}_tabref_0; mkdir -p $d; nice $SO so/cfg_sb8.yaml $O/b${b}_tabref.bin --dump $d/map2d.csv > $d/run.log 2>&1
  cp $O/b${b}_tabref.bin.truth.csv $d/map2d.csv.truth.csv; rm -f $d/map2d.csv.snap.csv $O/b${b}_tabref.bin
  echo "$(date +%T) b$b" >> $O/progress.log
}
for b in $(seq 0 25); do
  while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 4 ]; do sleep 5; done
  one $b & sleep 2
done
wait; echo ALLDONE >> $O/progress.log
