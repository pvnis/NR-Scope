#!/bin/bash
# field_record.sh: a short raw IQ recording (sc16) of a cell, for offline work later
# (replay through NR-Scope with the file radio, GoldenSniffer development on commercial
# signals). 122.88 Msps sc16 is 0.49 GB per second: keep recordings short.
#
# usage: field_record.sh --center HZ [--secs 4] [--gain 35] [--rate 122.88e6] [--tag NAME] [--out DIR]
# Output: <out>/<UTC date-time>_<tag>_iq/{iq.sc16, meta.txt}
set -u
center="" secs=4 gain=35 rate=122.88e6 tag=cell out=$HOME/field
args="type=x4xx,master_clock_rate=245760000,addr=192.168.20.2,mgmt_addr=192.168.1.190,use_dpdk=1234"
rx=/usr/local/lib/uhd/examples/rx_samples_to_file
while [ $# -gt 0 ]; do
  case $1 in
    --center) center=$2; shift 2 ;;
    --secs) secs=$2; shift 2 ;;
    --gain) gain=$2; shift 2 ;;
    --rate) rate=$2; shift 2 ;;
    --tag) tag=$2; shift 2 ;;
    --out) out=$2; shift 2 ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
done
[ -z "$center" ] && { sed -n 2,8p "$0"; exit 1; }
need_gb=$(python3 -c "print(int(${rate} * 4 * ${secs} / 1e9) + 2)")
free_gb=$(df -BG --output=avail "$out" 2>/dev/null | tail -1 | tr -dc 0-9)
[ -z "$free_gb" ] && { mkdir -p "$out"; free_gb=$(df -BG --output=avail "$out" | tail -1 | tr -dc 0-9); }
if [ "$free_gb" -lt "$need_gb" ]; then
  echo "not enough disk: $free_gb GB free, $need_gb needed"; exit 1
fi
pgrep -x nrscope >/dev/null && { echo "nrscope is running: stop it first (it holds the X410)"; exit 1; }
dir=$out/$(date -u +%Y%m%d_%H%M%S)_${tag}_iq
mkdir -p "$dir"
sudo -n "$rx" --args "$args" --rate "$rate" --freq "$center" --gain "$gain" --ant TX/RX0 --duration "$secs" \
  --type short --file "$dir/iq.sc16" --setup 0.5 --stats > "$dir/rx.log" 2>&1
st=$?
{
  echo "center_hz $center"
  echo "rate_sps $rate"
  echo "gain_db $gain"
  echo "secs $secs"
  echo "format sc16 interleaved I/Q, antenna TX/RX0, X410 reference about -4.02 ppm"
  echo "utc $(date -u +%FT%TZ)"
  echo "rx_samples_to_file exit $st"
  grep -iE "overflow|dropped|O+$|Received" "$dir/rx.log" | tail -3
} > "$dir/meta.txt"
ls -la "$dir"
cat "$dir/meta.txt"
