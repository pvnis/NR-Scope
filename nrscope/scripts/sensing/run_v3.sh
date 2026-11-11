#!/bin/bash
# After the AoA fixes (2026-11-10, OPEN_ISSUES 20): all phase ramps about DC, per-4-PRB-group inter-antenna calibration
# from the TRS (GCAL=1), tablet TRS reference per snapshot AND per precoder epoch.
# A) main table rows that the fixes can change (NR-Scope --inject after processing)
# B) three targets injected into the raw data, for the angles
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=/mnt/mydrive/field/main_v3; A=/mnt/mydrive/field/aoa_v3; mkdir -p $O $A; cd /home/sens/field; T="25,2.5,-20,15;60,-5,-22,-25;100,8,-25,30"
conv() {   # conv <variant> <burst dir> <out.bin> [INJ]
  v=$1; B=$2; out=$3; inj=$4
  case $v in
    trs)        GCAL=1 REF_EP=snap  INJ="$inj" WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr_ep.py $B $out ;;
    trsep)      GCAL=1 REF_EP=epoch INJ="$inj" WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr_ep.py $B $out ;;
    light)      GCAL=1 MINPRB=48 INJ="$inj" MODE=notab USE_CSI=0 nice python3 so_light.py $B $out ;;
    lightcsi)   GCAL=1 MINPRB=48 INJ="$inj" MODE=notab USE_CSI=1 nice python3 so_light.py $B $out ;;
    lighttab)   GCAL=1 MINPRB=48 INJ="$inj" MODE=tab USE_CSI=1 nice python3 so_light.py $B $out ;;
    lighttabep) GCAL=1 MINPRB=48 INJ="$inj" MODE=tab USE_CSI=1 TAB_EP=epoch nice python3 so_light.py $B $out ;;
  esac
}
cfg_of() { case $1 in light*) echo so/cfg_trs.yaml;; *) echo so/cfg_sb8.yaml;; esac; }
one() {
  b=$1; B=$R/b$b
  for v in trs trsep light lightcsi lighttab lighttabep; do
    conv $v $B $O/b${b}_$v.bin > $O/b${b}_$v.log 2>&1; cfg=$(cfg_of $v)
    for tg in slow fast faster; do
      case $tg in slow) I1=range=40,speed=0.6; I2=range=80,speed=-0.8;; fast) I1=range=40,speed=4; I2=range=80,speed=-3;; faster) I1=range=40,speed=7; I2=range=80,speed=-9;; esac
      for lv in -10 -20; do d=$O/m_b${b}_${v}_${tg}_$lv; mkdir -p $d; nice $SO $cfg $O/b${b}_$v.bin --dump $d/map2d.csv --inject $I1,level=$lv,angle=15 --inject $I2,level=$lv,angle=-20 > $d/run.log 2>&1; rm -f $d/map2d.csv.snap.csv; done
    done
    rm -f $O/b${b}_$v.bin
  done
  for v in trs trsep lightcsi lighttab lighttabep; do
    conv $v $B $A/b${b}_$v.bin "$T" > $A/b${b}_$v.log 2>&1
    d=$A/m_b${b}_$v; mkdir -p $d; nice $SO $(cfg_of $v) $A/b${b}_$v.bin --dump $d/map2d.csv > $d/run.log 2>&1
    cp $A/b${b}_$v.bin.truth.csv $d/map2d.csv.truth.csv; rm -f $d/map2d.csv.snap.csv $A/b${b}_$v.bin
  done
  echo "$(date +%T) b$b" >> $O/progress.log
}
for b in $(seq 0 25); do
  while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk '/MemAvailable/{print int($2/1e6)}' /proc/meminfo) -lt 4 ]; do sleep 5; done
  one $b & sleep 2
done
wait; echo ALLDONE >> $O/progress.log
