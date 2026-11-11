#!/bin/bash
# KPI gaps (OPEN_ISSUES 22): no injection at all -> detections per map/cell (Pfa upper bound: real movers count too),
# and maps with / without clutter removal on the same data -> clutter suppression ratio.
R=/home/sens/field/20261105_135044_sunrise_632_4ch_bursts_tabletdl_iq; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
O=/mnt/mydrive/field/kpi; mkdir -p $O; cd /home/sens/field
one() {
  b=$1; B=$R/b$b
  GCAL=1 REF_EP=snap WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr_ep.py $B $O/b${b}_trs.bin > $O/b${b}_trs.log 2>&1
  GCAL=1 MINPRB=48 MODE=tab USE_CSI=1 nice python3 so_light.py $B $O/b${b}_lighttab.bin > $O/b${b}_lighttab.log 2>&1
  GCAL=1 MINPRB=48 MODE=notab USE_CSI=1 nice python3 so_light.py $B $O/b${b}_lightcsi.bin > $O/b${b}_lightcsi.log 2>&1
  GCAL=1 REF_EP=hyb HYB_MAXEP=0.01 WSH=1 RUN=all VAR=trs DIR=replay_blind_ch CH=0,1,2,3 nice python3 gs_to_nrsr_ep.py $B $O/b${b}_trshyb.bin > $O/b${b}_trshyb.log 2>&1
  for v in trs lighttab lightcsi trshyb; do
    [ -f $O/b${b}_$v.bin ] || continue
    case $v in light*) c=so/cfg_trs;; trshyb) c=so/cfg_sb8_noalign;; *) c=so/cfg_sb8;; esac
    for cl in clut noclut; do cf=$c.yaml; [ $cl = noclut ] && cf=${c}_noclut.yaml
      d=$O/${cl}_b${b}_$v; mkdir -p $d; nice $SO $cf $O/b${b}_$v.bin --dump $d/map2d.csv > $d/run.log 2>&1; rm -f $d/map2d.csv.snap.csv; done
    rm -f $O/b${b}_$v.bin
  done
  echo "$(date +%T) b$b" >> $O/progress.log
}
for b in $(seq 0 25); do
  while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 4 ]; do sleep 5; done
  one $b & sleep 2
done
wait; echo ALLDONE >> $O/progress.log
