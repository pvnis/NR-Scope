#!/bin/bash
# rec4c.sh: 4-channel IQ bursts (one streamer, sample-aligned, into RAM, then moved to disk) while the tablet
# downloads. Guards: waits for free RAM before each burst, watchdog on the radio, discards incomplete bursts,
# flushes the disk after each move. Keeps going until NB good bursts or MAXTRY attempts.
set -u
NB=${NB:-8}; SECS=${SECS:-2.45}; MAXTRY=${MAXTRY:-$((NB*2))}
IF=$(ip route show default | sed -n "s/.* dev \([^ ]*\).*/\1/p" | head -1)
out=$HOME/field/$(date -u +%Y%m%d_%H%M%S)_sunrise_632_4ch_bursts_tabletdl_iq; mkdir -p $out; L=$out/dl_rate.log; : > $L
echo $out > $HOME/field/.rec4_new
shm=/dev/shm/rec4; sudo -n rm -rf $shm
pgrep -x nrscope >/dev/null && { echo "nrscope running" > $out/FAILED; exit 1; }
need=$(( 4*4*${SECS%.*}*123 + 4*4*123 ))   # MB for 4 channels x SECS (rounded up) + margin
need=6200
T=$((MAXTRY*40+60))
for u in https://speedtest.init7.net/10GB.dd https://speedtest.init7.net/10GB.dd https://speedtest.init7.net/10GB.dd \
         https://nbg1-speed.hetzner.com/10GB.bin https://nbg1-speed.hetzner.com/10GB.bin https://ash-speed.hetzner.com/10GB.bin; do
  timeout $T curl -s -o /dev/null "$u" &
done
( p=$(cat /sys/class/net/$IF/statistics/rx_bytes); for i in $(seq $T); do sleep 1; n=$(cat /sys/class/net/$IF/statistics/rx_bytes); echo "$(date +%T.%3N) $(( (n-p)*8/1000000 )) Mbps" >> $L; p=$n; [ -e $out/STOP ] && break; done; echo DONE >> $L ) &
sleep 6
good=0; try=0; nbad=0; exp=$(python3 -c "print(int($SECS*122.88e6*0.95))")
while [ $good -lt $NB ] && [ $try -lt $MAXTRY ]; do
  try=$((try+1))
  for w in $(seq 60); do av=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo); [ $av -ge $need ] && break; sleep 1; done
  bd=$out/b$good; mkdir -p $bd $shm; sudo -n chmod 777 $shm
  echo "try $try burst $good start $(date -u +%T.%3N) memavail ${av}MB" >> $out/bursts.log
  cd $shm && sudo -n timeout -k 5 25 /usr/local/lib/uhd/examples/rx_samples_to_file --args "type=x4xx,master_clock_rate=245760000,addr=192.168.20.2,mgmt_addr=192.168.1.190,use_dpdk=1234" \
    --rate 122.88e6 --freq 3750e6 --gain 38 --ant TX/RX0 --channels "0,1,2,3" --duration $SECS --type short --file iq.sc16 --setup 0.5 --stats > $bd/rx.log 2>&1
  rc=$?
  n=$(sed -n "s/^Received \([0-9]*\) samples in \([0-9.]*\) seconds.*/\1/p" $bd/rx.log | head -1); s=$(sed -n "s/^Received \([0-9]*\) samples in \([0-9.]*\) seconds.*/\2/p" $bd/rx.log | head -1)
  sz=$(stat -c %s $shm/iq_ch3.sc16 2>/dev/null || echo 0)
  if [ "${n:-0}" -ge $exp ] && python3 -c "import sys; sys.exit(0 if float('${s:-99}')<$SECS+0.3 else 1)" && [ $sz -ge $((exp*4)) ]; then
    for c in 0 1 2 3; do sudo -n mv $shm/iq_ch$c.sc16 $bd/ && sudo -n chown sens:sens $bd/iq_ch$c.sc16; done
    echo "  good: rc $rc received $n in $s s" >> $out/bursts.log; good=$((good+1)); nbad=0
  else
    echo "  BAD: rc $rc received ${n:-?} in ${s:-?} s, file $sz bytes; discarded" >> $out/bursts.log; mv $bd/rx.log $out/bad_try$try.log; rmdir $bd; nbad=$((nbad+1))
    if [ $rc -ge 124 ] || [ $nbad -ge 2 ]; then
      sudo -n pkill -9 -f "[r]x_samples_to_file"; sudo -n rm -f /dev/hugepages/rtemap_*
      timeout 120 ssh -o BatchMode=yes root@192.168.1.190 "systemctl restart usrp-hwd" && sleep 20
      echo "  restarted usrp-hwd on the X410" >> $out/bursts.log; nbad=0
    fi
  fi
  sudo -n rm -rf $shm; sync
  sleep 2
done
touch $out/STOP; pkill -f "[c]url -s -o /dev/null"
{ echo "center_hz 3750e6"; echo "rate_sps 122.88e6"; echo "gain_db 38"; echo "bursts $good good of $try tries x $SECS s"
  echo "format sc16 interleaved I/Q, 4 channels (A:0 A:1 B:0 B:1, antenna TX/RX0 each), one streamer per burst (sample-aligned), files b<n>/iq_ch<c>.sc16"
  echo "utc $(date -u +%FT%TZ)"; cat $out/bursts.log; } > $out/meta.txt
echo DONE > $out/DONE
