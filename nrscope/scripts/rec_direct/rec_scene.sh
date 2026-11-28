#!/bin/bash
# rec_scene.sh: one recording of the scene, run by hand on sens8. Prints START, then GO when the lidar is writing frames AND
# the 4-channel IQ is streaming (rec_cont.sh, re-started up to 3 times: the X410 start fails about one time in two), then
# STOP when the IQ capture ends. Ctrl-C stops everything it started.
# usage: rec_scene.sh <name> [seconds]     env as rec_cont.sh: NODL=1 (no tablet download), ANT (default RX1), GAIN (default 45)
#        EXPECT_PCI (default 630: the cell the tablet must be on, checked over adb)
#        SESSION (default: the newest test in ~/GIOVANNI/real)  -> lidar frames in ~/GIOVANNI/real/<SESSION>/lidar_<name>
set -u
HERE=$(dirname "$(readlink -f "$0")")   # this script's folder: NR-Scope/nrscope/scripts/rec_direct
NAME=${1:?usage: rec_scene.sh <name> [seconds]}; SECS=${2:-120}; export ANT=${ANT:-RX1} GAIN=${GAIN:-45} NODL=${NODL:-0} SECS
G=$HOME/GIOVANNI; SESSION=${SESSION:-$(cd $G/real && ls -t */session.json | head -1 | cut -d/ -f1)}
LD=$G/real/$SESSION/lidar_$NAME; LP=""; RP=""
say() { echo; echo "==================== $* ===================="; echo; }
stop_all() { [ -n "$RP" ] && kill $RP 2>/dev/null; [ -n "$LP" ] && kill $LP 2>/dev/null
  for p in $(pgrep -x "rec_direct.*"); do sudo -n kill -INT $p; done; say "STOPPED (Ctrl-C)"; exit 1; }
trap stop_all INT TERM
trap "rm -f $HOME/isac_data/.pause_decode" EXIT   # the decode resumes when this recording ends
[ -e "$LD" ] && { echo "$LD exists: choose another name"; exit 1; }
# a decode running (rec_pipe.sh): pause it (no new chunk starts) and wait for its NR-Scope runs to end; it resumes at exit
if pgrep -f "[r]ec_pipe.sh decode" >/dev/null; then
  touch $HOME/isac_data/.pause_decode; echo "decode running: pausing it, waiting for its NR-Scope runs to finish (a few minutes)..."
  for i in $(seq 600); do pgrep -x nrscope >/dev/null || break; sleep 1; done
  pgrep -x nrscope >/dev/null && { rm -f $HOME/isac_data/.pause_decode; echo "NR-Scope still running after 10 min: not recording"; exit 1; }
  echo "decode paused"
fi
pgrep -x "rec_direct.*" >/dev/null && { echo "rec_direct already running"; exit 1; }
pgrep -x nrscope >/dev/null && { echo "nrscope is running"; exit 1; }
tablet_cell() {   # the tablet's serving cell (adb, read-only): "NR PCI 630" / "LTE PCI 64" / "unknown"
  adb shell dumpsys telephony.registry 2>/dev/null | grep -o "CellInfo\(Nr\|Lte\):{ *mRegistered=YES[^}]*mPci *= *[0-9]*" | head -1 |
    sed -E 's/CellInfo(Nr|Lte).*mPci *= *([0-9]+)/\1 PCI \2/; s/^Nr/NR/; s/^Lte/LTE/' | grep . || echo "unknown"
}
TAB=$(tablet_cell); echo "tablet serving cell: $TAB"
MODE=$(adb shell cmd phone get-allowed-network-types-for-users 2>/dev/null | tr -d "\r" | head -1); echo "tablet allowed network types: ${MODE:-unknown}"
if [ $NODL = 0 ] && [ "$MODE" != "NR" ]; then
  read -r -p "the tablet is not in NR-only mode (${MODE:-unknown}): it may fall back to LTE; record anyway? [y/N] " ans; [ "$ans" = y ] || { echo "not recording"; exit 1; }
fi
if [ $NODL = 0 ] && [ "$TAB" != "NR PCI ${EXPECT_PCI:-630}" ]; then
  read -r -p "the tablet is not on NR PCI ${EXPECT_PCI:-630} ($TAB): record anyway? [y/N] " ans; [ "$ans" = y ] || { echo "not recording"; exit 1; }
fi
say "START  $NAME: ${SECS} s, antenna $ANT, gain $GAIN, tablet download $([ $NODL = 0 ] && echo on || echo off), session $SESSION"

echo "lidar: starting (frames into $LD)"
(cd $G && exec .venv/bin/python tools/lidar_record.py $LD $((SECS + 15)) > /tmp/rec_scene_lidar.log 2>&1) & LP=$!
for i in $(seq 30); do sleep 1; ls $LD/*.ply >/dev/null 2>&1 && break; done
ls $LD/*.ply >/dev/null 2>&1 || { echo "lidar: NO FRAMES after 30 s (see /tmp/rec_scene_lidar.log)"; kill $LP 2>/dev/null; say "STOP: lidar failed, nothing recorded"; exit 1; }
echo "lidar: recording"

ok=0
for try in 1 2 3; do
  prev=$(cat $HOME/isac_data/.rec4_new 2>/dev/null)
  echo "radio: starting the IQ capture (try $try)"
  (exec $HERE/rec_cont.sh > /tmp/rec_scene_iq.log 2>&1) & RP=$!
  D=""; for i in $(seq 30); do sleep 1; n=$(cat $HOME/isac_data/.rec4_new 2>/dev/null); [ "$n" != "$prev" ] && { D=$n; break; }; done
  s0=0; for i in $(seq 45); do sleep 1
    kill -0 $RP 2>/dev/null || break
    s=$(stat -c %s $D/iq_ch0.sc16 2>/dev/null || echo 0); [ $s -gt 0 ] && [ $s -gt $s0 ] && [ $s0 -gt 0 ] && { ok=1; break; }; s0=$s
  done
  [ $ok = 1 ] && break
  echo "radio: did not start (try $try): $(tail -2 $D/rec.log 2>/dev/null | tr '\n' ' ')"
  kill $RP 2>/dev/null; for p in $(pgrep -x "rec_direct.*"); do sudo -n kill -INT $p; done; wait $RP 2>/dev/null; RP=""
  sleep 5
done
[ $ok = 1 ] || { kill $LP 2>/dev/null; say "STOP: the radio did not start in 3 tries, nothing recorded"; exit 1; }

say "GO  --  start moving now  (${SECS} s)"
t0=$(date +%s)
while kill -0 $RP 2>/dev/null; do sleep 1; e=$(( $(date +%s) - t0 ))
  [ $((e % 10)) = 0 ] && [ $e -lt $SECS ] && echo "  ${e} s / ${SECS} s"
done
say "STOP  --  stop moving"
grep -E "^gaps |^overflows " $D/meta.txt 2>/dev/null; ls $D | grep -E "^(DONE|WARN|FAILED)$"
echo "IQ: $D"
echo "lidar: $LD (keeps writing a little longer, then stops by itself)"
