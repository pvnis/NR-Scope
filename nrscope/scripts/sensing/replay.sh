#!/bin/bash
# replay.sh <recording dir> <name> [dump_rnti]: the recording through NR-Scope (file radio), lossless.
# Config: the recording's replay5 config with paths moved to <name>. Extra env (e.g. NRSCOPE_GS_RAW_SIZE) passes through.
R=$1; O=$R/$2; mkdir -p $O
sed "s#$R/replay5#$O#g" $R/replay5/config.yaml > $O/config.yaml
[ -n "$3" ] && sed -i "s#^  dump_rnti: .*#  dump_rnti: \"$3\"#" $O/config.yaml
rm -f $O/iq.fifo; mkfifo $O/iq.fifo
cd /home/sens/GIOVANNI/NR-Scope/build/nrscope/src
nice python3 -I /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sc16_to_fc32.py $R/iq.sc16 $O/iq.fifo &
stdbuf -oL -eL timeout 1800 ./nrscope $O/config.yaml > $O/nrscope.log 2>&1
rm -f $O/iq.fifo
echo DONE >> $O/nrscope.log
