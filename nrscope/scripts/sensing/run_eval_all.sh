#!/bin/bash
# After run_gsw_all.sh: every burst through NR-Scope sensing_offline with injected targets, for each method.
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=/home/sens/field/so/all; mkdir -p $O; cd /home/sens/field
until grep -q ALLDONE ana/gsw_replay.log; do sleep 30; done
one() {   # burst
  b=$1; B=$R/b$b
  for v in raw proc trs; do WSH=1 RUN=all VAR=$v DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr.py $B $O/b${b}_$v.bin > $O/b${b}_$v.log 2>&1; done
  MODE=notab nice python3 so_light.py $B $O/b${b}_light.bin > $O/b${b}_light.log 2>&1
  USE_DMRS=0 nice python3 so_light.py $B $O/b${b}_trsonly.bin > $O/b${b}_trsonly.log 2>&1
  for v in raw proc trs light trsonly; do
    cfg=so/cfg_sb8.yaml; [ $v = light ] || [ $v = trsonly ] && cfg=so/cfg_trs.yaml
    for tg in slow fast; do
      if [ $tg = slow ]; then I1=range=40,speed=0.6; I2=range=80,speed=-0.8; else I1=range=40,speed=4; I2=range=80,speed=-3; fi
      for lv in -10 -20; do d=$O/m_b${b}_${v}_${tg}_$lv; mkdir -p $d; $SO $cfg $O/b${b}_$v.bin --dump $d/map2d.csv --inject $I1,level=$lv,angle=15 --inject $I2,level=$lv,angle=-20 > $d/run.log 2>&1; done
    done
    rm -f $O/b${b}_$v.bin
  done
  echo "$(date +%T) b$b" >> $O/progress.log
}
bs=$(ls -d $R/b* | sed "s#.*/b##" | sort -n); set -- $bs
while [ $# -gt 0 ]; do one $1 & [ $# -gt 1 ] && one $2; wait; shift; [ $# -gt 0 ] && shift; done
echo ALLDONE >> $O/progress.log
