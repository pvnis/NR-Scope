#!/bin/bash
# Robustness: targets injected into the RAW signals (3 targets, mixed range/speed/level), every burst, every method.
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=/home/sens/field/so/robust; mkdir -p $O; cd /home/sens/field
T="25,2.5,-20,15;60,-5,-22,-25;100,8,-25,30"
one() {
  b=$1; B=$R/b$b
  INJ="$T" WSH=1 RUN=all VAR=raw DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr.py $B $O/b${b}_raw.bin > $O/b${b}_raw.log 2>&1
  INJ="$T" WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr.py $B $O/b${b}_trs.bin > $O/b${b}_trs.log 2>&1
  INJ="$T" USE_CSI=1 MODE=notab nice python3 so_light.py $B $O/b${b}_lnotab.bin > $O/b${b}_lnotab.log 2>&1
  INJ="$T" USE_CSI=1 MODE=tab nice python3 so_light.py $B $O/b${b}_ltab.bin > $O/b${b}_ltab.log 2>&1
  for v in raw trs lnotab ltab; do
    cfg=so/cfg_sb8.yaml; case $v in lnotab|ltab) cfg=so/cfg_trs.yaml;; esac
    d=$O/m_b${b}_$v; mkdir -p $d
    if [ -f $O/b${b}_$v.bin ]; then $SO $cfg $O/b${b}_$v.bin --dump $d/map2d.csv > $d/run.log 2>&1; cp $O/b${b}_$v.bin.truth.csv $d/map2d.csv.truth.csv 2>/dev/null; fi
    rm -f $O/b${b}_$v.bin
  done
  echo "$(date +%T) b$b" >> $O/progress.log
}
set -- $(ls -d $R/b* | sed "s#.*/b##" | grep -E "^[0-9]+$" | sort -n)
while [ $# -gt 0 ]; do one $1 & [ $# -gt 1 ] && one $2; wait; shift; [ $# -gt 0 ] && shift; done
echo ALLDONE >> $O/progress.log
