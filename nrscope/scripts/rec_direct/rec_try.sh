#!/bin/bash
# on sens8: rec_cont.sh with retries, and a guard for rec_direct hanging at exit (after its final "done:" line,
# i.e. with every sample written): 20 s grace, then SIGINT, then SIGKILL by PID. usage: rec_try.sh <SECS> <NODL> <tries> <label>
SECS=$1; NODL=$2; N=$3; LAB=$4; LOG=$HOME/isac_data/rec_try_$LAB.log; : > $LOG
for i in $(seq $N); do
  old=$(cat $HOME/isac_data/.rec4_new 2>/dev/null)
  echo "$(date -u +%T) try $i: SECS=$SECS NODL=$NODL" >> $LOG
  SECS=$SECS NODL=$NODL "$(dirname "$(readlink -f "$0")")"/rec_cont.sh >> $LOG 2>&1 & rc_pid=$!
  sleep 5; d=$(cat $HOME/isac_data/.rec4_new 2>/dev/null)
  end=$(( $(date +%s) + SECS + 120 ))
  while kill -0 $rc_pid 2>/dev/null; do
    if grep -q "^done:" $d/rec.log 2>/dev/null; then
      sleep 20; p=$(pgrep -x "rec_direct.*")
      if [ -n "$p" ]; then echo "$(date -u +%T) rec_direct hung after done: SIGINT/SIGKILL $p" >> $LOG
        sudo -n kill -INT $p; sleep 5; pgrep -x "rec_direct.*" > /dev/null && sudo -n kill -9 $(pgrep -x "rec_direct.*"); fi
      break
    fi
    if [ $(date +%s) -gt $end ]; then echo "$(date -u +%T) timeout without done: killing" >> $LOG
      p=$(pgrep -x "rec_direct.*"); [ -n "$p" ] && sudo -n kill -9 $p; break; fi
    sleep 2
  done
  wait $rc_pid 2>/dev/null
  if [ -e "$d/DONE" ]; then echo "$(date -u +%T) OK $d $(cat $d/WARN 2>/dev/null)" >> $LOG; echo "RESULT OK $d" >> $LOG; exit 0; fi
  echo "$(date -u +%T) try $i failed: $(cat $d/FAILED 2>/dev/null) | $(tail -2 $d/rec.log 2>/dev/null | tr '\n' ' ')" >> $LOG
  p=$(pgrep -x "rec_direct.*"); [ -n "$p" ] && sudo -n kill -9 $p; sleep 10
done
echo "RESULT FAILED" >> $LOG; exit 1
