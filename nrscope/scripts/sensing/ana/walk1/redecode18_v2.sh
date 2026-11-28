#!/bin/bash
# redecode_v2.sh: 18 Nov decoded again with the joint DM-RS test (commit 54916ca) into a separate chunks folder; 4 at once
R=/mnt/mydrive/field/20261118_154424_sunrise_632_4ch_cont_tabletdl_iq; CH=/mnt/mydrive/field/20261118_154424_chunks_v2; L=$CH/progress.log
mkdir -p $CH; echo "$(date +%T) start" > $L
for k in $(seq 0 47); do t0=$(python3 -c "print($k * 2.5)")
  for c in 0 1 2 3; do
    while [ $(pgrep -xc nrscope) -ge 4 ]; do sleep 3; done
    PCI=632 /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/rec_direct/replay_chunk.sh $R $CH/b$k $c $t0 2.45 > /dev/null 2>&1 &
    sleep 1
  done
done
wait; echo "$(date +%T) ALLDONE" >> $L
