#!/bin/bash
# run_real_walk.sh <chunks dir> <out dir>: the sensing pipeline on a decoded recording (rec_pipe.sh decode): per chunk
# so_light (TRS + blind DM-RS of the other UEs, MODE=notab; PCI and TRS offset from the chunk) + sensing_offline
# (so/cfg_trs.yaml) -> <out>/b<k>.map2d.csv; 4 chunks at once. Angle calibration from the TRS by default (env SESSION=<session.json>; CAL=0: off). Progress: <out>/progress.log
CH=$1; O=$2; mkdir -p $O; L=$O/progress.log; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
cd /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing; echo "$(date +%T) start $CH" > $L
# angle calibration (default on): the X410 chain phase offsets from the TRS direct path of this recording
# (aoa_cal.py), written into a copy of CFG unless CFG already has an active aoa_cal_deg. Needs the recording
# session.json (gNB, lidar pose): env SESSION=<path>. Without it the angles are uncalibrated (warned).
CFG=${CFG:-so/cfg_trs.yaml}
if [ "${CAL:-1}" = 1 ] && ! grep -q "^  aoa_cal_deg:" $CFG; then
  if [ -n "${SESSION:-}" ] && [ -f "$SESSION" ]; then
    nice python3 /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/aoa_cal.py $CH $SESSION > $O/aoa_cal.log 2>&1; rc=$?
    line=$(grep "^aoa_cal_deg:" $O/aoa_cal.log)
    if [ $rc = 0 ] && [ -n "$line" ]; then
      awk -v L="  $line   # from the TRS direct path (aoa_cal.py)" '{print} /^sensing:/{print L}' $CFG > $O/cfg_cal.yaml; CFG=$O/cfg_cal.yaml
      echo "$(date +%T) angle calibration: $line" >> $L
    else echo "$(date +%T) WARNING: angle calibration failed or unstable (see aoa_cal.log): angles uncalibrated" >> $L; fi
  else echo "$(date +%T) WARNING: no SESSION=<session.json>: angles uncalibrated" >> $L; fi
fi
one() {
  b=$1; [ -s $O/b$b.map2d.csv ] && return
  SIDE=$O/b$b.side MINPRB=48 USE_DMRS=${USE_DMRS:-1} USE_CSI=0 MODE=${MODE:-notab} nice python3 so_light.py $CH/b$b $O/b$b.bin > $O/b$b.conv.log 2>&1
  NR_SENSING_TRACE=$O/b$b.trace nice $SO $CFG $O/b$b.bin --dump $O/b$b.map2d.csv > $O/b$b.run.log 2>&1
  rm -f $O/b$b.bin $O/b$b.map2d.csv.snap.csv; echo "$(date +%T) b$b" >> $L
}
for d in $(ls -d $CH/b* | sort -V); do b=${d##*/b}
  while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 5; done
  one $b & sleep 2
done
wait; echo "$(date +%T) ALLDONE" >> $L
