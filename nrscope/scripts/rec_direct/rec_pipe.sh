#!/bin/bash
# rec_pipe.sh <store|decode|grid> <IQ capture dir> <lidar recording dir>: the processing of one recording (Real page buttons).
#   store   move the IQ to /mnt/mydrive/field (move_rec.sh: copy, verify, link) and copy the lidar folder next to it
#   decode  (env PCI, default 630) replay the IQ (already on the drive) through NR-Scope + GoldenSniffer: 2.45 s chunks every 2.5 s x 4 channels
#           (replay_chunk.sh, as the 18 Nov usage map), at most 4 at once; waits while ~/isac_data/.pause_decode exists
#           (rec_scene.sh sets it during a recording). Out: /mnt/mydrive/field/<IQ name>_chunks/b<k>/replay_blind_ch<c>
#   grid    the slot x PRB occupancy from the decoded records (gui/rec_grid.py) -> <lidar dir>/grid/
# Status: <lidar dir>/pipe_<step>.log (lines "PROGRESS done total", last line DONE or FAILED <why>); runs detached.
set -u
HERE=$(dirname "$(readlink -f "$0")")   # this script's folder: NR-Scope/nrscope/scripts/rec_direct
K=$1; IQ=$2; LD=$3; ST=$LD/pipe_$K.log; F=$HOME/isac_data
say() { echo "$(date +%T) $*" >> $ST; }
fail() { say "FAILED $*"; exit 1; }
: > $ST; say "start $K $IQ"
pgrep -x "rec_direct.*" >/dev/null && [ $K = store ] && fail "a recording is running: store after it"
case $K in
store)
  [ -L "$IQ" ] && { say "IQ already on the drive: $(readlink -f $IQ)"; } || {
    tot=$(du -sb "$IQ" | cut -f1); D=/mnt/mydrive/field/$(basename "$IQ")
    ( while sleep 5; do [ -d $D ] && say "PROGRESS $(du -sb $D 2>/dev/null | cut -f1) $tot"; [ -L "$IQ" ] && break; [ -e $ST.stop ] && break; done ) & mon=$!
    $HERE/move_rec.sh "$IQ"; touch $ST.stop; kill $mon 2>/dev/null; rm -f $ST.stop
    tail -3 $F/ana/move_rec.log >> $ST
    [ -L "$IQ" ] || fail "move_rec.sh did not complete (see ~/isac_data/ana/move_rec.log)"; }
  ( sudo -n fstrim / > /dev/null 2>&1 & )   # the NVMe copy is gone: trim, or the next recording writes too slowly
  D=$(readlink -f "$IQ")/lidar_$(basename "$LD" | sed 's/^lidar_//')
  ionice -c3 nice rsync -a --exclude 'pipe_*.log' "$LD/" "$D/" || fail "lidar copy to $D"
  say "lidar copied to $D"; say DONE ;;
decode)
  R=$(readlink -f "$IQ"); case $R in /mnt/mydrive/*) ;; *) fail "store the IQ on the drive first";; esac
  dur=$(python3 -c "
import re; m=open('$R/meta.txt').read(); n=int(re.search(r'samples_per_channel (\d+)',m).group(1)); r=float(re.search(r'rate_sps (\d+)',m).group(1)); print(n/r)")
  CH=/mnt/mydrive/field/$(basename "$R")_chunks; mkdir -p $CH
  export PCI=${PCI:-630} NIDS=${NIDS:-}; say "PCI $PCI, extra DM-RS IDs [${NIDS}]"
  jobs_=(); for k in $(python3 -c "import math; print(' '.join(str(k) for k in range(int(math.floor(($dur - 0.5) / 2.5)) + 1)))"); do
    t0=$(python3 -c "print($k * 2.5)"); d=$(python3 -c "print(round(min(2.45, $dur - $k * 2.5 - 0.02), 3))")
    for c in 0 1 2 3; do jobs_+=("$k $c $t0 $d"); done; done
  n=${#jobs_[@]}; say "PROGRESS 0 $n"; say "IQ $dur s -> $((n / 4)) chunks x 4 channels into $CH"
  for j in "${jobs_[@]}"; do set -- $j
    grep -q DONE $CH/b$1/replay_blind_ch$2/nrscope.log 2>/dev/null && continue    # already decoded (resume)
    while [ -e $F/.pause_decode ] || [ $(pgrep -xc nrscope) -ge 4 ] || [ $(awk '/MemAvailable/{print int($2/1e6)}' /proc/meminfo) -lt 8 ]; do sleep 5
      say "PROGRESS $(ls $CH/b*/replay_blind_ch*/nrscope.log 2>/dev/null | xargs -r grep -l DONE | wc -l) $n"; done
    $HERE/replay_chunk.sh "$R" $CH/b$1 $2 $3 $4 > /dev/null 2>&1 &
    sleep 1
  done
  wait; while [ $(pgrep -xc nrscope) -gt 0 ]; do sleep 5; done
  ok=$(ls $CH/b*/replay_blind_ch*/nrscope.log 2>/dev/null | xargs -r grep -l DONE | wc -l); say "PROGRESS $ok $n"
  [ $ok = $n ] || fail "$ok of $n chunk-channels decoded"
  echo "$CH" > $LD/chunks_dir; say DONE ;;
grid)
  CH=$(cat $LD/chunks_dir 2>/dev/null) || fail "decode first"
  $HOME/GIOVANNI/.venv/bin/python $HOME/GIOVANNI/gui/rec_grid.py "$CH" "$LD/grid" "$(readlink -f "$IQ")" >> $ST 2>&1 || fail "rec_grid.py"
  say DONE ;;
*) fail "unknown step $K";;
esac
