#!/bin/bash
# replay_chunk.sh <4-ch recording dir> <out burst dir> <channel> <start s> <dur s>: one channel of one time slice
# through NR-Scope as replay4ch_blind.sh does (all RNTIs, blind dump, port 1001), output <out>/replay_blind_ch<c>.
# env PCI: the cell to decode (default the template's, 632)
R=$1; B=$2; C=$3; T0=$4; DUR=$5; O=$B/replay_blind_ch$C; mkdir -p $O
HERE=$(dirname "$(readlink -f "$0")")
sed "s#@OUT@#$O#g" $HERE/replay_template.yaml > $O/config.yaml   # the 3 Nov replay config, output paths filled in
sed -i "s#^  dump_rnti: .*#  dump_rnti: \"all\"#" $O/config.yaml
sed -i "/^  dump_rnti: /a\  dump_blind: true" $O/config.yaml
[ -n "${PCI:-}" ] && sed -i "s#^  pci: .*#  pci: $PCI#" $O/config.yaml
[ -n "${NIDS:-}" ] && sed -i "/^  dump_blind: true/a\  pdsch_dmrs_nids: [$NIDS]" $O/config.yaml   # env NIDS: extra PDSCH DM-RS scrambling IDs (blind records rnti 1, 2, ...)   # env PCI: the cell (default: the template, 632)
rm -f $O/iq.fifo; mkfifo $O/iq.fifo
cd /home/sens/GIOVANNI/NR-Scope/build/nrscope/src
nice python3 -I /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/rec_direct/sc16_slice.py $R/iq_ch$C.sc16 $O/iq.fifo $T0 $DUR &
NRSCOPE_GS_DUMP_P1001=1 stdbuf -oL -eL timeout 3600 ./nrscope $O/config.yaml > $O/nrscope.log 2>&1
rm -f $O/iq.fifo
echo DONE >> $O/nrscope.log
