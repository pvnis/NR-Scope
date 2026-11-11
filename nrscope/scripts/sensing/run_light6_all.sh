#!/bin/bash
# merging ideas: one reference number per slot (TAB_EP/BLIND_EP=slot), blind record split into UE tracks (BLIND_TRACK)
cd /home/sens/field; export OUT=/mnt/mydrive/field/light6all
jobsq() { while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 4 ]; do sleep 5; done; }
for b in $(seq 0 25); do
  jobsq; ./run_light4.sh $b ltabslot MODE=tab USE_CSI=1 TAB_EP=slot BLIND_EP=slot TAB_ALIGN=0 & sleep 2
  jobsq; ./run_light4.sh $b ltabhybbs MODE=tab USE_CSI=1 TAB_EP=hyb BLIND_EP=slot HYB_MAXEP=0.01 TAB_ALIGN=0 & sleep 2
  jobsq; ./run_light4.sh $b ltabhybtrk MODE=tab USE_CSI=1 TAB_EP=hyb BLIND_EP=hyb BLIND_TRACK=1 HYB_MAXEP=0.01 TAB_ALIGN=0 & sleep 2
  jobsq; ./run_light4.sh $b lcsislot MODE=notab USE_CSI=1 BLIND_EP=slot TAB_ALIGN=0 & sleep 2
  jobsq; ./run_light4.sh $b lcsitrk MODE=notab USE_CSI=1 BLIND_EP=hyb BLIND_TRACK=1 HYB_MAXEP=0.01 TAB_ALIGN=0 & sleep 2
done
wait; echo ALLDONE >> $OUT/progress.log
