#!/bin/bash
# The twin grid masks from the 18 Nov continuous capture: every chunk, all records (b<k>.npz) and the tablet-only
# real / full-band variants (tablet 0x7688). usage: export18.sh <chunk>
k=$1; CH=/mnt/mydrive/field/20261118_154424_chunks; O=/mnt/mydrive/field/precsim_masks18; F=$HOME/isac_data/ana/masks18
python3 $F/masks_export.py $CH/b$k $O/b$k.npz > $O/b$k.log 2>&1 && for m in real full; do python3 $F/masks_tablet.py $O/b$k.npz $O/b${k}_tab_$m.npz $m 0x7688 >> $O/b$k.log 2>&1; done
