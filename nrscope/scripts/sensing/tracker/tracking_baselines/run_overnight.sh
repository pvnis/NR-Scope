#!/bin/bash
# The 5 other long runs of the 2-receiver scenario (density 3 and 6, seeds 1-3, minus dens6 seed1,
# done first), simulated one after another through the GUI, then every tracker on each
# (default settings, the tuning grids, the summary). Simulation and tracking never overlap
# (memory). Log: tracking_baselines/results/overnight.log
cd "$(dirname "$0")/.." || exit 1
LOG=tracking_baselines/results/overnight.log
echo "== $(date) simulating" >> $LOG
.venv/bin/python tools/batch.py bc_env_new_pos_2_receivers_long --only 3:1,3:2,6:2,3:3,6:3 >> $LOG 2>&1
for R in $(ls -d logs/*_bc_env_new_pos_2_receivers_long_dens*_seed* | xargs -n1 basename); do
  [ -f tracking_baselines/results/$R/summary.md ] && continue
  echo "== $(date) tracking $R" >> $LOG
  mkdir -p tracking_baselines/results/$R
  nice -n 19 env OMP_NUM_THREADS=1 .venv/bin/python tracking_baselines/harness.py $R > tracking_baselines/results/default_$R.log 2>&1
  nice -n 19 env OMP_NUM_THREADS=1 .venv/bin/python tracking_baselines/tune.py $R --jobs 4 > tracking_baselines/results/tune_$R.log 2>&1
  .venv/bin/python tracking_baselines/summary.py $R > /dev/null 2>> $LOG
done
echo "== $(date) done" >> $LOG
