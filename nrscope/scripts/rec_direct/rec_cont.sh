#!/bin/bash
# rec_cont.sh: one continuous 4-channel IQ capture (rec_direct: gap-free, O_DIRECT to the NVMe) while the
# tablet downloads over its USB tether. Replaces the burst recorder rec4c.sh (rx_samples_to_file lost 0.37%).
#   SECS     capture length in s (default 100; keep <= ~120: the NVMe stays fast only up to ~310 GB)
#   TETHER_IF  tablet tether interface (default: autodetect the USB rndis/ncm interface)
#   NODL=1   record without generating tablet traffic
#   ANT      X410 port (default TX/RX0; the antennas may sit on RX1: check the level first)
#   GAIN     receive gain in dB (default 38)
# Output: ~/isac_data/<utc>_sunrise_632_4ch_cont_tabletdl_iq/{iq_ch<c>.sc16, meta.txt, gaps.csv, rec.log, dl_rate.log}
# Move to /mnt/mydrive afterwards with: move_rec.sh <dir>
set -u
SECS=${SECS:-100}; NODL=${NODL:-0}
HERE=$(dirname "$(readlink -f "$0")")   # this script's folder: NR-Scope/nrscope/scripts/rec_direct
REC=$HERE/rec_direct   # logs t0_wall_utc_s in meta.txt (lidar alignment without a fitted offset)
mkdir -p $HOME/isac_data   # data root (recordings, state files)
ARGS="type=x4xx,master_clock_rate=245760000,addr=192.168.20.2,mgmt_addr=192.168.1.190,use_dpdk=1234"
die() { echo "$*" >&2; [ -n "${out:-}" ] && echo "$*" > $out/FAILED; exit 1; }

[ -x $REC ] || die "rec_direct not built: $REC"
pgrep -x nrscope >/dev/null && die "nrscope is running"
pgrep -x "rec_direct.*" >/dev/null && die "rec_direct is already running"
need=$(python3 -c "print(int($SECS*4*122.88e6*4/1e9*1.05)+20)")
free=$(df -BG --output=avail $HOME/isac_data | tail -1 | tr -dc 0-9)
[ $free -ge $need ] || die "need ${need} GB free on $(df --output=target $HOME/isac_data | tail -1), have ${free} GB"

IF=${TETHER_IF:-}
if [ -z "$IF" ]; then
  for d in /sys/class/net/*; do
    drv=$(basename "$(readlink -f $d/device/driver 2>/dev/null)" 2>/dev/null)
    case "$drv" in rndis_host|cdc_ncm|cdc_ether) IF=$(basename $d);; esac
  done
fi
if [ $NODL = 0 ]; then
  [ -n "$IF" ] || die "no tablet tether interface found (plug the tablet, enable USB tethering, or set TETHER_IF / NODL=1)"
  ip -4 addr show dev $IF | grep -q inet || die "tether $IF has no IPv4 address yet"
fi

out=$HOME/isac_data/$(date -u +%Y%m%d_%H%M%S)_sunrise_632_4ch_cont_tabletdl_iq; mkdir -p $out
echo $out > $HOME/isac_data/.rec4_new
L=$out/dl_rate.log; : > $L
pids=()
if [ $NODL = 0 ]; then
  T=$((${SECS%.*}+60))
  for u in https://speedtest.init7.net/10GB.dd https://speedtest.init7.net/10GB.dd https://speedtest.init7.net/10GB.dd \
           https://nbg1-speed.hetzner.com/10GB.bin https://nbg1-speed.hetzner.com/10GB.bin https://ash-speed.hetzner.com/10GB.bin; do
    timeout $T curl -s --interface $IF -o /dev/null "$u" & pids+=($!)
  done
  ( p=$(cat /sys/class/net/$IF/statistics/rx_bytes)
    for i in $(seq $T); do sleep 1; n=$(cat /sys/class/net/$IF/statistics/rx_bytes)
      echo "$(date +%T.%3N) $(( (n-p)*8/1000000 )) Mbps" >> $L; p=$n; [ -e $out/STOP ] && break; done; echo DONE >> $L ) & pids+=($!)
  sleep 6
fi

echo "start $(date -u +%T.%3N) tether ${IF:-none}" > $out/rec.log
sudo -n $REC --dir $out --duration $SECS --args "$ARGS" --rate 122.88e6 --freq 3750e6 --gain ${GAIN:-38} \
  --ant ${ANT:-TX/RX0} --channels 0,1,2,3 --setup 1 >> $out/rec.log 2>&1
rc=$?
touch $out/STOP; for p in "${pids[@]}"; do kill $p 2>/dev/null; done; wait 2>/dev/null

{ echo "utc_end $(date -u +%FT%TZ)"; echo "tether ${IF:-none}"; echo "tablet_traffic $([ $NODL = 0 ] && echo yes || echo no)"
  [ -s $L ] && awk '/Mbps/{s+=$2;n++} END{if(n) printf "dl_mean_mbps %.0f\n", s/n}' $L; } >> $out/meta.txt
g=$(sed -n 's/^gaps \([0-9]*\)/\1/p' $out/meta.txt)
if [ $rc -ne 0 ]; then echo "rec_direct rc $rc" > $out/FAILED
elif [ "${g:-1}" != 0 ]; then echo "gaps $g (zero-filled, see gaps.csv)" > $out/WARN; echo DONE > $out/DONE
else echo DONE > $out/DONE; fi
tail -3 $out/rec.log
