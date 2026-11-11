#!/bin/bash
cd /home/sens/field; B=${BURSTS:-20}
jobsq() { while [ $(jobs -r | wc -l) -ge 3 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 4 ]; do sleep 5; done; }
for b in $B; do
  jobsq; ./run_hyb4.sh $b snap REF_EP=snap & sleep 2
  jobsq; ./run_hyb4.sh $b hyb REF_EP=hyb & sleep 2
  jobsq; ./run_hyb4.sh $b hyb10 REF_EP=hyb HYB_MAXEP=0.01 & sleep 2
  jobsq; ./run_hyb4.sh $b hyb40 REF_EP=hyb HYB_MAXEP=0.04 & sleep 2
  jobsq; ./run_hyb4.sh $b epoch REF_EP=epoch & sleep 2
done
wait; echo ALLDONE >> ${OUT:-/mnt/mydrive/field/hyb4}/progress.log
