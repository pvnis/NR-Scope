#!/bin/bash
# run_all.sh: after the recording ends, replay every burst x channel, then the multi-UE analysis on new + old recordings.
o=$(cat ~/field/.rec4_new); A=~/field/ana; mkdir -p $A
until [ -e $o/DONE ]; do sleep 10; done

echo "replay start $(date +%T)" > $A/run.log
for bd in $o/b*/; do for c in 0 1 2 3; do
  [ -s $bd/replay_ch$c/ue_dmrs.bin ] || bash ~/field/replay4ch.sh ${bd%/} $c >/dev/null 2>&1
done; echo "replayed $bd $(date +%T)" >> $A/run.log; done
cp ~/field/gs_all.py $A/
O=~/field/20261104_164101_sunrise_632_4ch_bursts_tabletdl_iq; T=~/field/20261105_132906_sunrise_632_4ch_bursts_tabletdl_iq
for bd in $T/b*/; do for c in 0 1 2 3; do [ -s $bd/replay_ch$c/ue_dmrs.bin ] || bash ~/field/replay4ch.sh ${bd%/} $c >/dev/null 2>&1; done; done
CH=0,3 python3 $A/gs_all.py $o $A/new03 > $A/new03.log 2>&1; echo "new03 done $(date +%T)" >> $A/run.log
CH=0,3 python3 $A/gs_all.py $T $A/test03 > $A/test03.log 2>&1
CH=0,3 python3 $A/gs_all.py $O $A/old03 > $A/old03.log 2>&1; echo "old03 done $(date +%T)" >> $A/run.log
python3 $A/gs_all.py $o $A/new > $A/new.log 2>&1; echo "new done $(date +%T)" >> $A/run.log
python3 $A/gs_all.py $O $A/old > $A/old.log 2>&1; echo "old done $(date +%T)" >> $A/run.log
echo ALLDONE >> $A/run.log
