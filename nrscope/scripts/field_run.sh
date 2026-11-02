#!/bin/bash
# field_run.sh: one NR-Scope + GoldenSniffer run on a cell, everything kept in a run
# directory. Retries the X410 start (about one start in two fails on the management link).
#
# usage: field_run.sh --ssb HZ [--pci N] [--center HZ] [--gain DB] [--secs S] [--tag NAME]
#                     [--out DIR] [--set "key: value"]...   (extra golden_sniffer keys)
#   --ssb     SSB frequency (field_scan.py prints it)
#   --center  radio centre; default: the SSB frequency. After the first run, use the
#             rx_center_freq the "GoldenSniffer: carrier ..." line suggests.
#   --gain    rx gain (default 35); keep the "RX peak" of the [status] lines below ~-10 dBFS
#   --secs    run length (default 300)
# Output: <out>/<UTC date-time>_<tag>/{config.yaml, nrscope.log, golden_sniffer.csv,
#         msg4_bytes.log, summary.txt}
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
ssb="" pci="" center="" gain=35 secs=300 tag=cell out=$HOME/field sets=()
while [ $# -gt 0 ]; do
  case $1 in
    --ssb) ssb=$2; shift 2 ;;
    --pci) pci=$2; shift 2 ;;
    --center) center=$2; shift 2 ;;
    --gain) gain=$2; shift 2 ;;
    --secs) secs=$2; shift 2 ;;
    --tag) tag=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    --set) sets+=("$2"); shift 2 ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
done
[ -z "$ssb" ] && { sed -n 2,14p "$0"; exit 1; }
[ -z "$center" ] && center=$ssb
dir=$out/$(date -u +%Y%m%d_%H%M%S)_$tag
mkdir -p "$dir"
bin=$root/build/nrscope/src
cfg=$dir/config.yaml
python3 - "$root/nrscope/config/config_gs_field_x410.yaml" "$cfg" "$ssb" "$center" "$pci" "$gain" "$dir" "${sets[@]+"${sets[@]}"}" <<'EOF'
import re, sys
src, dst, ssb, center, pci, gain, d = sys.argv[1:8]
extra = sys.argv[8:]
s = open(src).read()
s = re.sub(r"(?m)^  ssb_freq: .*$", "  ssb_freq: %s" % ssb, s)
s = re.sub(r"(?m)^  rx_center_freq: .*$", "  rx_center_freq: %s" % center, s)
s = re.sub(r"(?m)^  rx_gain: .*$", "  rx_gain: %s" % gain, s)
s = re.sub(r"(?m)^  pci: .*\n", ("  pci: %s\n" % pci) if pci else "", s)
s = re.sub(r'(?m)^  csv: .*$', '  csv: "%s/golden_sniffer.csv"' % d, s)
i = s.index("golden_sniffer:\n")
for kv in extra:
    k = kv.split(":", 1)[0].strip()
    blk = s[i:]
    if re.search(r"(?m)^  %s:" % re.escape(k), blk):
        blk = re.sub(r"(?m)^  %s:.*$" % re.escape(k), "  " + kv.strip(), blk, count=1)
    else:
        blk = blk.replace("golden_sniffer:\n", "golden_sniffer:\n  %s\n" % kv.strip(), 1)
    s = s[:i] + blk
open(dst, "w").write(s)
EOF
echo "run directory: $dir"
cd "$bin" || exit 1
for try in 1 2 3 4 5; do
  sudo -n pkill -INT -x nrscope 2>/dev/null; sleep 2
  sudo -n rm -f msg4_bytes.log
  tmux kill-session -t field 2>/dev/null
  tmux new-session -d -s field "sudo stdbuf -oL -eL timeout --foreground -s INT $secs ./nrscope $cfg > $dir/nrscope.log 2>&1"
  for i in $(seq 40); do
    sleep 2
    grep -qE "GoldenSniffer ready|Cell search: no cell" "$dir/nrscope.log" && break
    pgrep -x nrscope >/dev/null || break
  done
  if grep -q "GoldenSniffer ready" "$dir/nrscope.log"; then
    echo "try $try: running ($(date -u +%T)); 'tmux attach -t field' to watch, Ctrl-C there to stop early"
    break
  fi
  if grep -q "Cell search: no cell" "$dir/nrscope.log"; then
    grep "Cell search: no cell" "$dir/nrscope.log"
    echo "no cell at this SSB frequency (wrong SSB/PCI, too weak, or gain): stopping"
    sudo -n pkill -INT -x nrscope; exit 2
  fi
  echo "try $try: did not start ($(grep -m1 -oE "rx xport timed out|AssertionError[^ ]*|Error[^\n]{0,60}" "$dir/nrscope.log" | head -1))"
  cp "$dir/nrscope.log" "$dir/nrscope_failed_start_$try.log"
done
grep -m1 "GoldenSniffer: carrier" "$dir/nrscope.log"
while pgrep -x nrscope >/dev/null; do sleep 5; done
cp -f msg4_bytes.log "$dir/" 2>/dev/null
{
  echo "ssb $ssb pci ${pci:-any} center $center gain $gain secs $secs"
  grep -m1 "GoldenSniffer: carrier" "$dir/nrscope.log"
  grep -oE "[0-9]+ slot\(s\) dropped so far" "$dir/nrscope.log" | tail -1
  grep -E "^\[status\]" "$dir/nrscope.log" | tail -1 | cut -c1-200
  grep -E "GoldenSniffer: [0-9]+ DCIs" "$dir/nrscope.log" | tail -1
  sed -n "/execution time/,/latency: slot -> tracked/p" "$dir/nrscope.log"
  echo "UEs (RNTI, DCIs, sizes) seen at least 10 times:"
  # the last report only (one is printed every report_s): the RNTI lines after the last
  # "GoldenSniffer: N DCIs" line
  last=$(grep -nE "GoldenSniffer: [0-9]+ DCIs" "$dir/nrscope.log" | tail -1 | cut -d: -f1)
  tail -n +"${last:-1}" "$dir/nrscope.log" | grep -E "^  RNTI 0x" | awk '$7+0 >= 10' | head -60
  echo "fallback states:"
  grep -cE "GoldenSniffer: UE .*-> CONVERGED" "$dir/nrscope.log" | sed "s/^/  converged: /"
  grep -cE "GoldenSniffer: UE .* discovered blind" "$dir/nrscope.log" | sed "s/^/  discovered: /"
  grep -c "hooray" "$dir/nrscope.log" | sed "s/^/  RRCSetup decoded (RACH): /"
} > "$dir/summary.txt"
cat "$dir/summary.txt"
