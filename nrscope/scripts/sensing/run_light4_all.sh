#!/bin/bash
cd /home/sens/field
jobsq() { while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 4 ]; do sleep 5; done; }
for b in $(seq 0 25); do
  jobsq; ./run_light4.sh $b ltab MODE=tab USE_CSI=1 & sleep 2
  jobsq; ./run_light4.sh $b ltabep MODE=tab USE_CSI=1 TAB_EP=epoch & sleep 2
  jobsq; ./run_light4.sh $b lcsi MODE=notab USE_CSI=1 & sleep 2
  jobsq; ./run_light4.sh $b lcsiep MODE=notab USE_CSI=1 BLIND_EP=epoch & sleep 2
  jobsq; ./run_light4.sh $b ltabep2 MODE=tab USE_CSI=1 TAB_EP=epoch BLIND_EP=epoch & sleep 2
done
wait; echo ALLDONE >> /mnt/mydrive/field/light4all/progress.log
