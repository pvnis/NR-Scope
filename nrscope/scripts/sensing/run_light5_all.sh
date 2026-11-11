#!/bin/bash
cd /home/sens/field; export OUT=/mnt/mydrive/field/light5all
jobsq() { while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 4 ]; do sleep 5; done; }
for b in $(seq 0 25); do
  jobsq; ./run_light4.sh $b ltabhyb MODE=tab USE_CSI=1 TAB_EP=hyb HYB_MAXEP=0.01 & sleep 2
  jobsq; ./run_light4.sh $b ltabhyb2 MODE=tab USE_CSI=1 TAB_EP=hyb BLIND_EP=hyb HYB_MAXEP=0.01 & sleep 2
  jobsq; ./run_light4.sh $b lcsihyb MODE=notab USE_CSI=1 BLIND_EP=hyb HYB_MAXEP=0.01 & sleep 2
  jobsq; ./run_light4.sh $b ltabhyb0 MODE=tab USE_CSI=1 TAB_EP=hyb HYB_MAXEP=0.01 TAB_ALIGN=0 & sleep 2
done
wait; echo ALLDONE >> $OUT/progress.log
