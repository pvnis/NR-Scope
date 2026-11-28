#!/bin/bash
# rec_scene.sh: one recording of the scene, run by hand on sens8 (receiver 1 + lidar). With RX2 set (default sens9) a second
# receiver records at the same time: its capture runs detached on that host (rec_try.sh: 3 tries, survives a dropped ssh link),
# ssh only starts, watches and stops it; the recordings are aligned afterwards by their wall-clock starts (t0_wall_utc_s).
# Prints START, then GO when the lidar is writing frames AND every receiver is streaming, then STOP after the walk time.
# Both captures last SECS + MARGIN, so a receiver that starts later than the other does not cut the walk. Ctrl-C stops everything.
# usage: rec_scene.sh <name> [walk seconds, default 90]
#   env as rec_cont.sh: NODL=1 (no tablet download), ANT (default RX1), GAIN (default 39)
#   RX2        second receiver's host (default sens9; RX2= records with sens8 alone); ANT2 / GAIN2 (default ANT / GAIN)
#   MARGIN     extra capture seconds on each receiver (default 30; walk + margin <= ~150 s: the NVMe slows past ~310 GB)
#   NOLIDAR=1  no lidar (radio-only test)
#   EXPECT_PCI (default 630: the cell the tablet must be on, checked over adb)
#   SESSION    (default: the newest test in ~/GIOVANNI/real)  -> lidar frames in ~/GIOVANNI/real/<SESSION>/lidar_<name>
# Out: the IQ folders (~/isac_data/<utc>_..._iq on each host), listed in <lidar folder>/receivers.txt
set -u
HERE=$(dirname "$(readlink -f "$0")")   # this script's folder: NR-Scope/nrscope/scripts/rec_direct (same path on every host)
NAME=${1:?usage: rec_scene.sh <name> [seconds]}; SECS=${2:-90}; MARGIN=${MARGIN:-30}; CAP=$((SECS + MARGIN))
export ANT=${ANT:-RX1} GAIN=${GAIN:-39} NODL=${NODL:-0}; export SECS=$CAP   # rec_cont.sh records CAP seconds
RX2=${RX2-sens9}; ANT2=${ANT2:-$ANT}; GAIN2=${GAIN2:-$GAIN}; NOLIDAR=${NOLIDAR:-0}; REC_CONT=${REC_CONT:-$HERE/rec_cont.sh}
G=$HOME/GIOVANNI; SESSION=${SESSION:-$(cd $G/real && ls -t */session.json | head -1 | cut -d/ -f1)}
LD=$G/real/$SESSION/lidar_$NAME; LP=""; RP=""; R2PID=""; D2=""
SSH2="ssh -o BatchMode=yes -o ConnectTimeout=5 $RX2"
say() { echo; echo "==================== $* ===================="; echo; }
rx2() { timeout 15 $SSH2 "$@" 2>/dev/null; }
stop_rx2() { [ -n "$RX2" ] && [ -n "$R2PID" ] && rx2 "kill $R2PID 2>/dev/null; for p in \$(pgrep -x rec_direct); do sudo -n kill -INT \$p; done"; }
stop_all() { [ -n "$RP" ] && kill $RP 2>/dev/null; [ -n "$LP" ] && kill $LP 2>/dev/null
  for p in $(pgrep -x "rec_direct.*"); do sudo -n kill -INT $p; done; stop_rx2; say "STOPPED (Ctrl-C)"; exit 1; }
trap stop_all INT TERM
trap "rm -f $HOME/isac_data/.pause_decode" EXIT   # the decode resumes when this recording ends
[ -e "$LD" ] && { echo "$LD exists: choose another name"; exit 1; }
# receiver 2: reachable, idle, and both clocks on time (the alignment uses the wall clocks)
clock_off() { chronyc tracking 2>/dev/null | awk '/System time/{print $4}'; }
if [ -n "$RX2" ]; then
  [ "$(rx2 hostname)" = "$RX2" ] || { echo "receiver 2 ($RX2) unreachable over ssh: fix the link, or RX2= to record with sens8 alone"; exit 1; }
  rx2 "pgrep -x rec_direct || pgrep -x nrscope" >/dev/null && { echo "receiver 2 ($RX2) is already recording or decoding"; exit 1; }
  h1=$(git -C $HERE rev-parse HEAD); h2=$(rx2 "git -C $HERE rev-parse HEAD")
  [ "$h1" = "$h2" ] || { echo "receiver 2 ($RX2) runs NR-Scope ${h2:0:7}, this host ${h1:0:7}: git pull there (or here) first"; exit 1; }
  o1=$(clock_off); o2=$(rx2 "chronyc tracking | awk '/System time/{print \$4}'")
  python3 -c "import sys; sys.exit(0 if abs(float('${o1:-1}')) < 0.05 and abs(float('${o2:-1}')) < 0.05 else 1)" ||
    { echo "clock offset too large or unknown (sens8 ${o1:-?} s, $RX2 ${o2:-?} s; need < 0.05 s): check chrony on both"; exit 1; }
  echo "receiver 2: $RX2 reachable and idle; clocks: $(hostname) ${o1} s, $RX2 ${o2} s from NTP"
fi
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
if [ $NODL = 0 ]; then
  TAB=$(tablet_cell); echo "tablet serving cell: $TAB"
  MODE=$(adb shell cmd phone get-allowed-network-types-for-users 2>/dev/null | tr -d "\r" | head -1); echo "tablet allowed network types: ${MODE:-unknown}"
  if [ "$MODE" != "NR" ]; then
    read -r -p "the tablet is not in NR-only mode (${MODE:-unknown}): it may fall back to LTE; record anyway? [y/N] " ans; [ "$ans" = y ] || { echo "not recording"; exit 1; }
  fi
  if [ "$TAB" != "NR PCI ${EXPECT_PCI:-630}" ]; then
    read -r -p "the tablet is not on NR PCI ${EXPECT_PCI:-630} ($TAB): record anyway? [y/N] " ans; [ "$ans" = y ] || { echo "not recording"; exit 1; }
  fi
fi
WALK=$((CAP - MARGIN))
say "START  $NAME: walk ${WALK} s (captures ${CAP} s), receivers $(hostname)${RX2:+ + $RX2}, antenna $ANT, gain $GAIN, tablet download $([ $NODL = 0 ] && echo on || echo off), session $SESSION"

if [ $NOLIDAR = 0 ]; then
  echo "lidar: starting (frames into $LD)"
  (cd $G && exec .venv/bin/python tools/lidar_record.py $LD $((CAP + 15)) > /tmp/rec_scene_lidar.log 2>&1) & LP=$!
  for i in $(seq 30); do sleep 1; ls $LD/*.ply >/dev/null 2>&1 && break; done
  ls $LD/*.ply >/dev/null 2>&1 || { echo "lidar: NO FRAMES after 30 s (see /tmp/rec_scene_lidar.log)"; kill $LP 2>/dev/null; say "STOP: lidar failed, nothing recorded"; exit 1; }
  echo "lidar: recording"
else mkdir -p $LD; echo "lidar: off (NOLIDAR=1)"; fi

# receiver 2: start its capture detached on its host, then receiver 1 here; both start in parallel
if [ -n "$RX2" ]; then
  prev2=$(rx2 "cat ~/isac_data/.rec4_new")
  R2PID=$(rx2 "GAIN=$GAIN2 ANT=$ANT2 ${REC_CONT2:+REC_CONT=$REC_CONT2} setsid nohup $HERE/rec_try.sh $CAP 1 3 $NAME > /dev/null 2>&1 < /dev/null & echo \$!")
  [ -n "$R2PID" ] && echo "radio 2: starting on $RX2 (rec_try pid $R2PID, log ~/isac_data/rec_try_$NAME.log there)" || echo "radio 2: COULD NOT START on $RX2"
fi
ok=0
for try in 1 2 3; do
  prev=$(cat $HOME/isac_data/.rec4_new 2>/dev/null)
  echo "radio 1: starting the IQ capture (try $try)"
  (exec $REC_CONT > /tmp/rec_scene_iq.log 2>&1) & RP=$!
  D=""; for i in $(seq 30); do sleep 1; n=$(cat $HOME/isac_data/.rec4_new 2>/dev/null); [ "$n" != "$prev" ] && { D=$n; break; }; done
  s0=0; for i in $(seq 45); do sleep 1
    kill -0 $RP 2>/dev/null || break
    s=$(stat -c %s $D/iq_ch0.sc16 2>/dev/null || echo 0); [ $s -gt 0 ] && [ $s -gt $s0 ] && [ $s0 -gt 0 ] && { ok=1; break; }; s0=$s
  done
  [ $ok = 1 ] && break
  echo "radio 1: did not start (try $try): $(tail -2 $D/rec.log 2>/dev/null | tr '\n' ' ')"
  kill $RP 2>/dev/null; for p in $(pgrep -x "rec_direct.*"); do sudo -n kill -INT $p; done; wait $RP 2>/dev/null; RP=""
  sleep 5
done
[ $ok = 1 ] || { kill $LP 2>/dev/null; stop_rx2; say "STOP: radio 1 did not start in 3 tries, nothing recorded"; exit 1; }
t1=$(date +%s); echo "radio 1: streaming ($D)"
# wait for receiver 2 to stream (its file growing); give up when waiting longer would cut the walk by more than 30 s,
# and go on with receiver 1 alone
t2=""
if [ -n "$RX2" ] && [ -n "$R2PID" ]; then
  s0=0; while [ $(date +%s) -lt $((t1 + MARGIN + 30)) ]; do
    st=$(rx2 "d=\$(cat ~/isac_data/.rec4_new); [ \"\$d\" != \"$prev2\" ] && echo \$d \$(stat -c %s \$d/iq_ch0.sc16 2>/dev/null || echo 0); grep -q 'RESULT FAILED' ~/isac_data/rec_try_$NAME.log 2>/dev/null && echo FAILED")
    case "$st" in *FAILED*) echo "radio 2: FAILED on $RX2 (3 tries)"; break;; esac
    set -- $st; s=${2:-0}; D2=${1:-}
    [ "$s" -gt 0 ] && [ "$s0" -gt 0 ] && [ "$s" -gt "$s0" ] && { t2=$(date +%s); echo "radio 2: streaming on $RX2 ($D2)"; break; }
    s0=$s; sleep 2
  done
  [ -z "$t2" ] && { say "WARNING: receiver 2 ($RX2) is NOT recording: this run has receiver 1 only (Ctrl-C to abort)"; D2=""; }
fi
# the walk must end before the first capture does: captures end CAP s after their start
tend=$((t1 + CAP)); [ -n "$t2" ] && [ $((t2 + CAP)) -lt $tend ] && tend=$((t2 + CAP))
now=$(date +%s); W=$WALK; [ $((now + W)) -gt $tend ] && W=$((tend - now))
if [ $W -lt 20 ]; then kill $RP $LP 2>/dev/null; for p in $(pgrep -x "rec_direct.*"); do sudo -n kill -INT $p; done; stop_rx2
  say "STOP: only ${W} s of common recording left (the receivers started too far apart): run again"; exit 1; fi
{ echo "rx1 $(hostname) $D"; [ -n "$D2" ] && echo "rx2 $RX2 $D2"; echo "walk_start_unix $now"; echo "walk_s $W"; } > $LD/receivers.txt
say "GO  --  start moving now  (${W} s)"
t0=$(date +%s)
while :; do sleep 1; e=$(( $(date +%s) - t0 )); [ $e -ge $W ] && break
  [ $((e % 10)) = 0 ] && echo "  ${e} s / ${W} s"
done
say "STOP  --  stop moving"
echo "captures finishing (their margin)..."; wait $RP 2>/dev/null
grep -E "^gaps |^overflows " $D/meta.txt 2>/dev/null; ls $D | grep -E "^(DONE|WARN|FAILED)$"
echo "IQ 1: $(hostname) $D"
if [ -n "$D2" ]; then
  for i in $(seq 60); do rx2 "grep -q '^RESULT' ~/isac_data/rec_try_$NAME.log" && break; sleep 2; done
  echo "IQ 2: $RX2 $D2  ($(rx2 "grep '^RESULT' ~/isac_data/rec_try_$NAME.log; grep -E '^gaps |^overflows ' $D2/meta.txt" | tr '\n' ' '))"
fi
[ $NOLIDAR = 0 ] && echo "lidar: $LD (keeps writing a little longer, then stops by itself)"
exit 0
