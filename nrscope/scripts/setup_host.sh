#!/bin/bash
# setup_host.sh: make this machine (sens8, sens9, ...) a field host identical to the others, from the repos alone.
# Idempotent: run after cloning and after every pull that touches it. Code lives only in the repos; data in ~/isac_data.
#   - system Python packages the chain uses (apt, list below); ~/GIOVANNI/.venv = the system Python (no pip)
#   - ~/isac_data (data root: recordings, decode outputs, the tracker's runs in paper_prec/logs)
#   - sensing/tracker/logs -> ~/isac_data/paper_prec/logs (the tracker reads and writes its runs there)
#   - builds the recorder (rec_direct/rec_direct) and checks NR-Scope's build is at the checked-out commit
# Host-specific settings stay out of the repos: /root/.config/uhd.conf (dpdk_mac), netplan, the receiver pose.
set -u
S=$(dirname "$(readlink -f "$0")"); NR=$(cd "$S/../.." && pwd); G=$HOME/GIOVANNI; D=$HOME/isac_data
ok() { echo "OK    $*"; }; fix() { echo "FIX   $*"; }; fail() { echo "FAIL  $*"; }
PKGS="python3-numpy python3-scipy python3-matplotlib python3-sympy python3-munkres python3-gps python3-lz4 python3-lxml
      python3-tornado python3-bs4 python3-yaml adb android-udev-rules"
miss=$(for p in $PKGS; do dpkg-query -W -f='${Status}' $p 2>/dev/null | grep -q "ok installed" || echo $p; done | tr "\n" " ")
if [ -n "$miss" ]; then sudo -n env DEBIAN_FRONTEND=noninteractive apt-get install -y $miss >/dev/null 2>&1 && fix "installed $miss" || fail "apt install $miss"; else ok "python packages"; fi
if [ -d $G ]; then
  [ -x $G/.venv/bin/python ] && ok "GIOVANNI venv" || { python3 -m venv --without-pip --system-site-packages $G/.venv && fix "created $G/.venv"; }
else fail "~/GIOVANNI not cloned (git clone git@gitlab.com:pvnis/passive-isac-twin.git ~/GIOVANNI)"; fi
mkdir -p $D/paper_prec/logs $D/ana && ok "data root $D"
L=$S/sensing/tracker/logs
if [ -L $L ] && [ "$(readlink $L)" = "$D/paper_prec/logs" ]; then ok "tracker logs link"
elif [ -e $L ] && [ ! -L $L ]; then fail "$L is a real folder: move it to $D/paper_prec/logs first"
else ln -sfn $D/paper_prec/logs $L && fix "tracker logs -> $D/paper_prec/logs"; fi
R=$S/rec_direct
if [ ! -x $R/rec_direct ] || [ $R/rec_direct.cpp -nt $R/rec_direct ]; then
  g++ -O2 -std=c++17 -o $R/rec_direct $R/rec_direct.cpp -luhd -lpthread && fix "built rec_direct" || fail "rec_direct build"
else ok "rec_direct built"; fi
# NR-Scope binaries: make -n lists what a build would redo (compiles or links) = the binary is older than its sources
for t in nrscope sensing_offline; do
  if [ ! -d $NR/build ]; then fail "NR-Scope not configured ($NR/build)"; break; fi
  n=$(cd $NR/build && make -n $t 2>/dev/null | grep -cE "Building CXX|Linking CXX")
  [ "$n" = 0 ] && ok "$t build up to date" || fail "$t older than its sources ($n steps): cd $NR/build && make -j\$(nproc) $t"
done
echo "host $(hostname): NR-Scope $(git -C $NR rev-parse --short HEAD) $(git -C $NR status -s | grep -q . && echo '(uncommitted changes)'), GIOVANNI $(git -C $G rev-parse --short HEAD 2>/dev/null)"
