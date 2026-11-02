#!/bin/bash
# field_prep.sh: before a field session with NR-Scope on sens6 + X410 (run as the sens
# user; uses sudo). Re-applies what does not survive a reboot or a power cycle of the
# X410, and checks the rest. Prints OK/FIX/FAIL per item; nothing here reboots anything.
#
#  - X410 management IPv4 192.168.1.190 on its RJ45 eth0 (DHCP-only there: gone after a
#    power cycle; UHD's DPDK discovery needs it). Re-added over the IPv6 link-local
#    address, which needs no configuration.
#  - rcuc/3..13 at SCHED_FIFO 95: NR-Scope's busy threads otherwise starve the RCU
#    callbacks of their CPU (RCU stalls; once it took sens6 off the network).
#  - the USB wifi interrupt on CPU 0, away from the workers.
#  - UHD's DPDK I/O thread on an isolated CPU (/root/.config/uhd.conf dpdk_lcore=2).
#  - CPU governor, FFTW wisdom, disk space, the X410 links.
X410_IP=192.168.1.190
X410_LL=fe80::280:2fff:fe36:ecdb
MGMT_IF=enp128s31f6
ok()   { echo "OK    $*"; }
fix()  { echo "FIX   $*"; }
fail() { echo "FAIL  $*"; FAILED=1; }
FAILED=0
sudo -n true 2>/dev/null || { echo "needs passwordless sudo"; exit 1; }

# X410 management
if ping -c1 -W1 $X410_IP >/dev/null 2>&1; then
  ok "X410 management $X410_IP"
elif timeout 15 ssh -o BatchMode=yes -o StrictHostKeyChecking=no -o ConnectTimeout=5 root@$X410_LL%$MGMT_IF \
     "ip addr add $X410_IP/24 dev eth0" 2>/dev/null && ping -c1 -W2 $X410_IP >/dev/null 2>&1; then
  fix "X410 management: added $X410_IP to its eth0 (over IPv6 link-local)"
else
  fail "X410 management unreachable (RJ45 cable from sens6 $MGMT_IF to the X410 eth0? X410 booted?)"
fi
carrier=$(cat /sys/class/net/$MGMT_IF/carrier_changes 2>/dev/null)
echo "      ($MGMT_IF carrier changes since boot: $carrier; a growing number means a bad RJ45 contact)"
if ip link show enp2s0f1np1 2>/dev/null | grep -q "state UP"; then ok "X410 SFP data link enp2s0f1np1 up"; else fail "X410 SFP data link enp2s0f1np1 down"; fi

# RCU callback threads
n=0
for p in $(pgrep -x "rcuc/([3-9]|1[0-3])"); do
  sudo -n chrt -f -p 95 $p && n=$((n + 1))
done
[ $n -ge 11 ] && ok "rcuc/3..13 at FIFO 95 ($n threads)" || fail "rcuc priorities: only $n of 11 set"

# wifi USB interrupt
irq=$(grep -m1 xhci_hcd /proc/interrupts | cut -d: -f1 | tr -d " ")
if [ -n "$irq" ]; then
  echo 0 | sudo -n tee /proc/irq/$irq/smp_affinity_list >/dev/null 2>&1
  [ "$(cat /proc/irq/$irq/effective_affinity_list)" = "0" ] && ok "xhci IRQ $irq on CPU 0" || fail "xhci IRQ $irq not on CPU 0"
fi

# UHD DPDK lcore
if sudo -n grep -qE "^dpdk_lcore=2" /root/.config/uhd.conf && sudo -n grep -qE "^dpdk_corelist=1,2" /root/.config/uhd.conf; then
  ok "UHD DPDK I/O on isolated CPU 2"
else
  fail "/root/.config/uhd.conf: dpdk_corelist must be 1,2 and dpdk_lcore 2 (CPUs 1-2 are the isolated ones)"
fi

# governor
gov=$(cat /sys/devices/system/cpu/cpu7/cpufreq/scaling_governor 2>/dev/null)
if [ "$gov" = performance ]; then ok "CPU governor performance"; else
  for g in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do echo performance | sudo -n tee $g >/dev/null; done
  fix "CPU governor set to performance (was $gov)"
fi

# FFTW wisdom (saved by NR-Scope at exit: the first run of the day measures the plans)
w=$(sudo -n stat -c %s /root/.srsran_fftwisdom 2>/dev/null || echo 0)
[ "$w" -gt 200 ] && ok "FFTW wisdom $w bytes" || echo "NOTE  no FFTW wisdom yet: the first run measures the FFT plans (some dropped slots at start)"

# disk
free_gb=$(df -BG --output=avail / | tail -1 | tr -dc 0-9)
[ "$free_gb" -ge 3 ] && ok "disk: $free_gb GB free" || fail "disk: only $free_gb GB free (raw IQ records need ~0.5 GB/s)"

# nothing left running on the X410
pgrep -x nrscope >/dev/null && fail "an nrscope is running (pid $(pgrep -x nrscope | tr "\n" " "))" || ok "no nrscope running"
[ $FAILED = 0 ] && echo "ready" || echo "NOT ready"
exit $FAILED
