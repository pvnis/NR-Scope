#!/bin/bash
# tablet_dl.sh: heavy downlink traffic through the tablet's USB tether (6 parallel downloads), so the tablet is a busy UE
# on the cell while a capture runs; the received rate goes to a log once per second. Used by rec_cont.sh (tablet on this
# host) and by rec_scene.sh over ssh (DL_HOST: the tablet on another receiver's host).
# usage: tablet_dl.sh check                             prints the tether interface, or the reason and exit 1
#        tablet_dl.sh run <seconds> <log> [<stopfile>]  until <seconds> have passed or <stopfile> exists
# env TETHER_IF: the interface (default: the USB rndis/ncm/ether interface)
set -u
IF=${TETHER_IF:-}
if [ -z "$IF" ]; then
  for d in /sys/class/net/*; do
    drv=$(basename "$(readlink -f $d/device/driver 2>/dev/null)" 2>/dev/null)
    case "$drv" in rndis_host|cdc_ncm|cdc_ether) IF=$(basename $d);; esac
  done
fi
case ${1:-} in
check)
  [ -n "$IF" ] || { echo "no tablet tether interface found on $(hostname) (plug the tablet, enable USB tethering, or set TETHER_IF / NODL=1)"; exit 1; }
  ip -4 addr show dev $IF | grep -q inet || { echo "tether $IF on $(hostname) has no IPv4 address yet"; exit 1; }
  echo $IF ;;
run)
  T=${2:?seconds}; L=${3:?log}; STOPF=${4:-}; : > $L
  [ -n "$IF" ] || { echo "no tether interface" >> $L; exit 1; }
  pids=()
  for u in https://speedtest.init7.net/10GB.dd https://speedtest.init7.net/10GB.dd https://speedtest.init7.net/10GB.dd \
           https://nbg1-speed.hetzner.com/10GB.bin https://nbg1-speed.hetzner.com/10GB.bin https://ash-speed.hetzner.com/10GB.bin; do
    timeout $T curl -s --interface $IF -o /dev/null "$u" & pids+=($!)
  done
  trap 'kill "${pids[@]}" 2>/dev/null; echo DONE >> $L; exit 0' TERM INT
  p=$(cat /sys/class/net/$IF/statistics/rx_bytes)
  for i in $(seq $T); do sleep 1; n=$(cat /sys/class/net/$IF/statistics/rx_bytes)
    echo "$(date +%T.%3N) $(( (n-p)*8/1000000 )) Mbps" >> $L; p=$n; [ -n "$STOPF" ] && [ -e "$STOPF" ] && break; done
  kill "${pids[@]}" 2>/dev/null; wait 2>/dev/null; echo DONE >> $L ;;
*) echo "usage: tablet_dl.sh check | run <seconds> <log> [<stopfile>]"; exit 2 ;;
esac
