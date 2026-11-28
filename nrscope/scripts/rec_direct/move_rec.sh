#!/bin/bash
# Move a recording to /mnt/mydrive/field, verify (sizes + checksums), then replace it by a link.
set -u
S=$(readlink -m ${1:?usage: move_rec.sh <recording dir>}); D=/mnt/mydrive/field/$(basename $S); L=/home/sens/isac_data/ana/move_rec.log
echo "$(date +%T) start" > $L
[ -L $S ] && { echo "already a link" >> $L; exit 0; }
ionice -c3 nice rsync -a $S/ $D/ >> $L 2>&1 || { echo "rsync FAILED" >> $L; exit 1; }
echo "$(date +%T) copied, verifying" >> $L
a=$(cd $S && find . -type f | wc -l); b=$(cd $D && find . -type f | wc -l)
sa=$(du -sb $S | cut -f1); sd=$(du -sb $D | cut -f1)
echo "files $a / $b, bytes $sa / $sd" >> $L
diffs=$(ionice -c3 nice rsync -anc --out-format="%n" $S/ $D/ | grep -v "/$" | wc -l)
echo "$(date +%T) checksum differences: $diffs" >> $L
if [ "$a" = "$b" ] && [ "$sa" = "$sd" ] && [ "$diffs" = "0" ]; then
  rm -rf $S && ln -s $D $S && echo "$(date +%T) DONE: original removed, link $S -> $D" >> $L
else
  echo "MISMATCH: original kept" >> $L
fi
df -h / /mnt/mydrive >> $L
