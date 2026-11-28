#!/bin/bash
# try630.sh <IQ dir>: one slice of channel 0 through NR-Scope + GoldenSniffer (PCI env, default 630; T0, DUR), all RNTIs + blind dump, into /tmp/rtest/ch0
R=$1; T=/home/sens/isac_data/20261103_155359_sunrise_632_tabletdl_iq/replay5; O=/tmp/rtest/ch0
rm -rf /tmp/rtest; mkdir -p $O
sed "s#$T#$O#g" $T/config.yaml | sed "s#^  pci: .*#  pci: ${PCI:-630}#; s#^  dump_rnti: .*#  dump_rnti: \"all\"#" > $O/config.yaml
sed -i "/^  dump_rnti: /a\  dump_blind: true" $O/config.yaml
mkfifo $O/iq.fifo
cd /home/sens/GIOVANNI/NR-Scope/build/nrscope/src
python3 -I /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/rec_direct/sc16_slice.py $R/iq_ch0.sc16 $O/iq.fifo ${T0:-10} ${DUR:-2.45} &
timeout 600 ./nrscope $O/config.yaml > $O/nrscope.log 2>&1
rm -f $O/iq.fifo
echo "ue_dmrs.bin bytes: $(stat -c %s $O/ue_dmrs.bin)   grant lines: $(wc -l < $O/golden_sniffer_grants.csv)"
