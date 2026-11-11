#!/bin/bash
# one burst, one reference variant, NR-Scope alignment on (sb8) and off (na): main slow rows + 3-target angles
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=${OUT:-/mnt/mydrive/field/hyb4}; mkdir -p $O; cd /home/sens/field; T="25,2.5,-20,15;60,-5,-22,-25;100,8,-25,30"; b=$1; B=$R/b$b; name=$2; shift 2
env "$@" GCAL=1 WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr_ep.py $B $O/b${b}_$name.bin > $O/b${b}_$name.log 2>&1
env "$@" INJ="$T" GCAL=1 WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr_ep.py $B $O/a${b}_$name.bin > $O/a${b}_$name.log 2>&1
for al in al na; do c=so/cfg_sb8.yaml; [ $al = na ] && c=so/cfg_sb8_noalign.yaml
  for tg in slow fast; do case $tg in slow) I1=range=40,speed=0.6; I2=range=80,speed=-0.8;; fast) I1=range=40,speed=4; I2=range=80,speed=-3;; esac
    for lv in -10 -20; do d=$O/m_b${b}_${name}${al}_${tg}_$lv; mkdir -p $d; nice $SO $c $O/b${b}_$name.bin --dump $d/map2d.csv --inject $I1,level=$lv,angle=15 --inject $I2,level=$lv,angle=-20 > $d/run.log 2>&1; rm -f $d/map2d.csv.snap.csv; done; done
  d=$O/p0_b${b}_${name}$al; mkdir -p $d; nice $SO $c $O/a${b}_$name.bin --dump $d/map2d.csv > $d/run.log 2>&1; cp $O/a${b}_$name.bin.truth.csv $d/map2d.csv.truth.csv; rm -f $d/map2d.csv.snap.csv
done
rm -f $O/b${b}_$name.bin $O/a${b}_$name.bin
echo "$(date +%T) DONE $b $name" >> $O/progress.log
