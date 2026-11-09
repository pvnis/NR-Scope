#!/bin/bash
# Light-network streams on every burst (TRS + CSI-RS, TRS + CSI-RS + DM-RS without the tablet) through sensing_offline.
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=/home/sens/field/so/light; mkdir -p $O; cd /home/sens/field
one() {
  b=$1; B=$R/b$b
  USE_DMRS=0 USE_CSI=1 nice python3 so_light.py $B $O/b${b}_trscsi.bin > $O/b${b}_trscsi.log 2>&1
  USE_DMRS=1 USE_CSI=1 MODE=notab nice python3 so_light.py $B $O/b${b}_lightcsi.bin > $O/b${b}_lightcsi.log 2>&1
  for v in trscsi lightcsi; do
    for tg in slow fast faster; do
      case $tg in slow) I1=range=40,speed=0.6; I2=range=80,speed=-0.8;; fast) I1=range=40,speed=4; I2=range=80,speed=-3;; faster) I1=range=40,speed=7; I2=range=80,speed=-9;; esac
      for lv in -10 -20; do d=$O/m_b${b}_${v}_${tg}_$lv; mkdir -p $d; $SO so/cfg_trs.yaml $O/b${b}_$v.bin --dump $d/map2d.csv --inject $I1,level=$lv,angle=15 --inject $I2,level=$lv,angle=-20 > $d/run.log 2>&1; done
    done
    rm -f $O/b${b}_$v.bin
  done
  echo "$(date +%T) b$b" >> $O/progress.log
}
set -- $(ls -d $R/b* | sed "s#.*/b##" | grep -E "^[0-9]+$" | sort -n)
while [ $# -gt 0 ]; do one $1 & [ $# -gt 1 ] && one $2; wait; shift; [ $# -gt 0 ] && shift; done
echo ALLDONE >> $O/progress.log
