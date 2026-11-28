#!/bin/bash
# Resume the 18 Nov usage map after an interruption: replays only the incomplete chunk channels, then the rest.
REC=/home/sens/isac_data/20261118_154424_sunrise_632_4ch_cont_tabletdl_iq
CH=/mnt/mydrive/field/20261118_154424_chunks; O=/mnt/mydrive/field/usemap_20261118b; L=$O/progress.log
SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline
echo "$(date +%T) resume" >> $L
for k in $(seq 0 47); do
  t0=$(python3 -c "print($k * 2.5)")
  for c in 0 1 2 3; do
    d=$CH/b$k/replay_blind_ch$c
    if grep -q "^DONE" $d/nrscope.log 2>/dev/null && [ -s $d/ue_dmrs.bin ]; then continue; fi
    rm -rf $d
    while [ $(pgrep -xc nrscope) -ge 4 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 8 ]; do sleep 5; done
    /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/rec_direct/replay_chunk.sh $REC $CH/b$k $c $t0 2.45 > /dev/null 2>&1 & sleep 1
  done
done
wait; echo "$(date +%T) replay done" >> $L
cd /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing
one() { b=$1; for m in tab notab; do
    SIDE=$O/b${b}_$m.side MINPRB=48 USE_DMRS=1 USE_CSI=1 MODE=$m nice python3 so_light.py $CH/b$b $O/b${b}_$m.bin > $O/b${b}_$m.conv.log 2>&1
    NR_SENSING_TRACE=$O/b${b}_$m.trace nice $SO so/cfg_trs.yaml $O/b${b}_$m.bin --dump $O/b${b}_$m.map2d.csv > $O/b${b}_$m.run.log 2>&1
    rm -f $O/b${b}_$m.bin; done; echo "$(date +%T) b$b" >> $L; }
for b in $(seq 0 47); do
  while [ $(jobs -r | wc -l) -ge 4 ] || [ $(awk "/MemAvailable/{print int(\$2/1e6)}" /proc/meminfo) -lt 6 ]; do sleep 5; done
  one $b & sleep 2
done
wait
python3 gen_um7.py $O $O/json > $O/gen.log 2>&1
python3 build_um_page.py $CH $O $O/page > $O/page.log 2>&1
echo "$(date +%T) ALLDONE" >> $L
