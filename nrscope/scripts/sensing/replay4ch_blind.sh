#!/bin/bash
# replay4ch_blind.sh: as replay4ch.sh plus golden_sniffer.dump_blind (6e), output replay_blind_ch<c>.
# replay4ch.sh <4-channel recording dir> <channel 0..3>: one channel through NR-Scope (file radio, lossless),
# DM-RS of every UE dumped with both ports (GSD2), survey and grants. Output <rec>/replay_ch<c>/.
R=$1; C=$2; O=$R/replay_blind_ch$C; mkdir -p $O
TEMPLATE=/home/sens/field/20261103_155359_sunrise_632_tabletdl_iq/replay5
sed "s#$TEMPLATE#$O#g" $TEMPLATE/config.yaml > $O/config.yaml
sed -i "s#^  dump_rnti: .*#  dump_rnti: \"all\"#" $O/config.yaml
sed -i "/^  dump_rnti: /a\  dump_blind: true" $O/config.yaml
rm -f $O/iq.fifo; mkfifo $O/iq.fifo
cd /home/sens/GIOVANNI/NR-Scope/build/nrscope/src
nice python3 -I /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sc16_to_fc32.py $R/iq_ch$C.sc16 $O/iq.fifo &
NRSCOPE_GS_DUMP_P1001=1 stdbuf -oL -eL timeout 1200 ./nrscope $O/config.yaml > $O/nrscope.log 2>&1
rm -f $O/iq.fifo
echo DONE >> $O/nrscope.log
