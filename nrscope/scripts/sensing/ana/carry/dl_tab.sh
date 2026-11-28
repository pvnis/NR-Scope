#!/system/bin/sh
# On the tablet (adb push to /data/local/tmp, run detached): 4 continuous downloads for DUR s plus a logger that writes,
# every second, the time and every rmnet interface's received bytes to /data/local/tmp/dl.log.
DUR=${1:-2700}
L=/data/local/tmp/dl.log; : > $L
for u in https://speedtest.init7.net/10GB.dd https://nbg1-speed.hetzner.com/10GB.bin https://speedtest.init7.net/10GB.dd https://ash-speed.hetzner.com/10GB.bin; do
  ( end=$(( $(date +%s) + DUR )); while [ $(date +%s) -lt $end ]; do curl -s -o /dev/null --max-time 300 "$u"; done ) &
done
( end=$(( $(date +%s) + DUR )); while [ $(date +%s) -lt $end ]; do echo "$(date +%s) $(grep rmnet /proc/net/dev | tr -s ' ' | tr '\n' '|')" >> $L; sleep 1; done ) &
wait
echo "END $(date +%s)" >> $L
