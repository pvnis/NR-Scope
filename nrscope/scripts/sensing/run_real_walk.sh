#!/bin/bash
# run_real_walk.sh <chunks dir> <out dir>: the sensing pipeline on a decoded recording (rec_pipe.sh decode): per chunk
# so_light (TRS + blind DM-RS of the other UEs, MODE=notab; PCI and TRS offset from the chunk) + sensing_offline
# (so/cfg_trs.yaml) -> <out>/b<k>.map2d.csv; 4 chunks at once. Progress: <out>/progress.log
CH=$1; O=$2; mkdir -p $O; L=$O/progress.log; SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
cd /home/sens/field; echo "$(date +%T) start $CH" > $L
one() {
  b=$1; [ -s $O/b$b.map2d.csv ] && return
  SIDE=$O/b$b.side MINPRB=48 USE_DMRS=${USE_DMRS:-1} USE_CSI=0 MODE=${MODE:-notab} nice python3 so_light.py $CH/b$b $O/b$b.bin > $O/b$b.conv.log 2>&1
  NR_SENSING_TRACE=$O/b$b.trace nice $SO ${CFG:-so/cfg_trs.yaml} $O/b$b.bin --dump $O/b$b.map2d.csv > $O/b$b.run.log 2>&1
  rm -f $O/b$b.bin $O/b$b.map2d.csv.snap.csv; echo "$(date +%T) b$b" >> $L
}
for d in $(ls -d $CH/b* | sort -V); do b=${d##*/b}
  while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 5; done
  one $b & sleep 2
done
wait; echo "$(date +%T) ALLDONE" >> $L
