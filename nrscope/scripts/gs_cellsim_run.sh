#!/bin/bash
# gs_cellsim_run.sh <run_dir> <cell.ini> [gs_cellsim --key value ...]
#
# One GoldenSniffer evaluation on a synthetic cell, with no IQ on disk: gs_cellsim writes the
# fc32 into a FIFO that NR-Scope (file radio) reads while it is generated, then
# gs_compare.py --truth scores golden_sniffer.csv against the ground truth. Everything lands
# in <run_dir>: truth.csv (+ .events.csv), gs.csv, iq.fifo.yaml, gen.log, nrscope.log.
#   GS_CELLSIM   generator (default: build-sim/nrscope/src/tests/gs_cellsim next to this repo)
#   NRSCOPE_BIN  nrscope binary (default: build/nrscope/src/nrscope); a copy under another
#                name is safer when other scripts pkill nrscope
#   TASKSET      CPU list (default 0,3)
set -u
here=$(cd "$(dirname "$0")/../.." && pwd)
d=$1; ini=$(readlink -f "$2"); shift 2
G=${GS_CELLSIM:-$here/build-sim/nrscope/src/tests/gs_cellsim}
N=${NRSCOPE_BIN:-$here/build/nrscope/src/nrscope}
CPUS=${TASKSET:-0,3}
mkdir -p "$d" && cd "$d" && d=$(pwd)
rm -f iq.fifo iq.fifo.yaml; mkfifo iq.fifo
taskset -c $CPUS nice -n 19 "$G" -c "$ini" --out iq.fifo --truth truth.csv --gs_csv "$d/gs.csv" "$@" > gen.log 2>&1 &
gp=$!
for i in $(seq 100); do [ -s iq.fifo.yaml ] && break; sleep 0.1; done   # written before the FIFO opens
/usr/bin/time -f "nrscope wall %e s, cpu %U+%S s, maxrss %M kB" timeout 3600 taskset -c $CPUS nice -n 19 "$N" iq.fifo.yaml > nrscope.log 2>&1
wait $gp
rm -f iq.fifo
tail -1 gen.log
tail -1 nrscope.log
grep -E "Cell search: no cell|SIB 1 Decoded" nrscope.log | head -1
python3 "$here/nrscope/scripts/gs_compare.py" --truth truth.csv gs.csv --summary
