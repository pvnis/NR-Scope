#!/bin/bash
# Figure 1 data: 18 Nov continuous capture, every chunk, the same processing apart from anchoring.
#   anch: so_light as in um_20261118b.sh (MODE=tab, MINPRB=48, USE_DMRS=1, USE_CSI=1): one gain per PRG and precoder epoch
#   raw : the same with PREC_IGNORE=1: one complex number per snapshot (common phase and scale only), no per-PRG gain
# Each through sensing_offline twice: cfg_trs.yaml (as the usage map: clutter removal on) and cfg_trs_noclut.yaml (clutter
# removal off: the window energy used as the direct-path reference). usage: fig1_run.sh <chunk numbers...>
CH=/mnt/mydrive/field/20261118_154424_chunks; O=/mnt/mydrive/field/ana_fig1; mkdir -p $O
SO=/home/sens/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline; F=/home/sens/isac_data/ana/fig1
cd /home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing
for b in "$@"; do
  for v in anch raw; do
    pi=0; [ $v = raw ] && pi=1
    PREC_IGNORE=$pi MINPRB=48 USE_DMRS=1 USE_CSI=1 MODE=tab nice python3 so_light.py $CH/b$b $O/b${b}_$v.bin > $O/b${b}_$v.conv.log 2>&1
    for c in clut noclut; do cf=$F/cfg_trs.yaml; [ $c = noclut ] && cf=$F/cfg_trs_noclut.yaml
      nice $SO $cf $O/b${b}_$v.bin --dump $O/b${b}_${v}_$c.map2d.csv > $O/b${b}_${v}_$c.run.log 2>&1; done
    rm -f $O/b${b}_$v.bin
  done
  echo "$(date +%T) b$b" >> $O/progress.log
done
