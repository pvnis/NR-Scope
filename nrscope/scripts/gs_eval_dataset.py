#!/usr/bin/env python3
"""Run NR-Scope + GoldenSniffer on GoldenSniffer dataset captures and score each against
its gNB log (gs_compare.py).

For every configuration: the SSB is located in the capture (find_ssb.py), a config is
written from a template, the sc16 capture is converted to fc32 on the fly through a FIFO
(no fc32 copy on disk), NR-Scope replays it, and gs_compare.py scores the DCIs.

  mode always:   blind search in every slot (GoldenSniffer alone)
  mode fallback: discovery for --startup-s, then the learned UEs are decoded by the
                 tracked decode (the path used live once UEs are known)

usage: gs_eval_dataset.py --conf 3,17 [--mode always|fallback|both] [--cpus 4-13]
       [--captures DIR] [--template YAML] [--nrscope BIN] [--out DIR] [--startup-s 2]
"""
import argparse
import glob
import os
import re
import subprocess
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", ".."))
ap = argparse.ArgumentParser()
ap.add_argument("--conf", required=True, help="comma separated configuration ids")
ap.add_argument("--mode", default="both", choices=["always", "fallback", "both"])
ap.add_argument("--captures", default=os.path.expanduser("~/GIOVANNI/5G-GoldenSniffer/captures"))
ap.add_argument("--template", default=os.path.join(root, "nrscope", "config", "config_gs_file.yaml"))
ap.add_argument("--nrscope", default=os.path.join(root, "build", "nrscope", "src", "nrscope"))
ap.add_argument("--find-ssb", default=os.path.expanduser("~/GIOVANNI/tools/gs/find_ssb.py"))
ap.add_argument("--out", default="/tmp/gs_eval")
ap.add_argument("--cpus", default="4-13")
ap.add_argument("--startup-s", type=float, default=2.0)
ap.add_argument("--timeout", type=int, default=600)
ap.add_argument("--extra", default="", help="extra golden_sniffer keys, 'k: v; k: v'")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)


def set_key(text, key, value):
    """replaces the first 'key: ...' line (keeping its indentation)"""
    pat = re.compile(r"^(\s*)" + re.escape(key) + r":.*$", re.M)
    if not pat.search(text):
        raise SystemExit("template has no key %s" % key)
    return pat.sub(lambda m: "%s%s: %s" % (m.group(1), key, value), text, count=1)


def gs_block(text, keys):
    """sets keys in the golden_sniffer: block, appending the missing ones"""
    i = text.index("golden_sniffer:")
    head, block = text[:i], text[i:]
    for k, v in keys.items():
        pat = re.compile(r"^(\s+)" + re.escape(k) + r":.*$", re.M)
        if pat.search(block):
            block = pat.sub(lambda m: "%s%s: %s" % (m.group(1), k, v), block, count=1)
        else:
            block = block.replace("golden_sniffer:\n", "golden_sniffer:\n  %s: %s\n" % (k, v), 1)
    return head + block


template = open(a.template).read()
rows = []
for cid in [int(x) for x in a.conf.split(",") if x]:
    files = glob.glob(os.path.join(a.captures, "iq_*_conf%d_0.sc16" % cid))
    log = os.path.join(a.captures, "gnb_conf%d.log" % cid)
    if not files or not os.path.exists(log):
        print("conf%d: capture or gNB log missing" % cid)
        continue
    cap = files[0]
    f0 = float(re.search(r"iq_([\d.]+)M_", os.path.basename(cap)).group(1)) * 1e6
    scs = 15 if f0 < 3e9 else 30
    band = 2 if f0 < 3e9 else 78
    out = subprocess.run(["python3", "-I", a.find_ssb, cap, str(f0), "23.04e6", str(scs), "40"],
                         capture_output=True, text=True).stdout.splitlines()
    m = re.match(r"ssb_freq ([\d.]+) MHz", out[0]) if out else None
    if not m:
        print("conf%d: no SSB found" % cid)
        continue
    ssb = int(round(float(m.group(1)) * 1e6))
    for mode in (["always", "fallback"] if a.mode == "both" else [a.mode]):
        tag = "conf%d_%s" % (cid, mode)
        fifo = os.path.join(a.out, tag + ".fifo")
        csv = os.path.join(a.out, tag + ".csv")
        yml = os.path.join(a.out, tag + ".yaml")
        lg = os.path.join(a.out, tag + ".log")
        for p in (fifo, csv):
            if os.path.exists(p):
                os.remove(p)
        os.mkfifo(fifo)
        t = template
        t = set_key(t, "rf_args", '"rx_file=%s,base_srate=23.04e6"' % fifo)
        t = set_key(t, "ssb_freq", str(ssb))
        t = set_key(t, "rx_center_freq", str(int(round(f0))))
        t = set_key(t, "band_list", "[%d]" % band)
        t = set_key(t, "scs_index", "0" if scs == 15 else "1")
        keys = {"mode": mode, "csv": '"%s"' % csv}
        if mode == "fallback":
            keys.update({"startup_s": a.startup_s, "scan_period_s": 30, "quiet_s": 1})
        for kv in [x for x in a.extra.split(";") if x.strip()]:
            k, v = kv.split(":", 1)
            keys[k.strip()] = v.strip()
        t = gs_block(t, keys)
        open(yml, "w").write(t)
        conv = subprocess.Popen(["taskset", "-c", "0,3", "nice", "python3", "-I",
                                 os.path.join(here, "sc16_to_fc32.py"), cap, fifo])
        t0 = time.time()
        with open(lg, "w") as fo:
            try:
                subprocess.run(["taskset", "-c", a.cpus, "nice", "stdbuf", "-oL", "-eL", a.nrscope, yml],
                               stdout=fo, stderr=subprocess.STDOUT, timeout=a.timeout)
            except subprocess.TimeoutExpired:
                fo.write("\n[gs_eval] timeout\n")
        conv.kill()
        os.remove(fifo)
        dt = time.time() - t0
        sc = subprocess.run(["python3", os.path.join(here, "gs_compare.py"), log, csv, "--scs-khz", str(scs),
                             "--wrap-by-sfn"], capture_output=True, text=True).stdout if os.path.exists(csv) else ""
        open(os.path.join(a.out, tag + ".score"), "w").write(sc)
        tot = re.findall(r"all\s+sent\s+(\d+)\s+found\s+(\d+).*?extra \(never sent\) (\d+)\s+duplicates (\d+)", sc)
        sent = sum(int(x[0]) for x in tot)
        found = sum(int(x[1]) for x in tot)
        extra = sum(int(x[2]) for x in tot)
        dup = sum(int(x[3]) for x in tot)
        other = re.search(r"never scheduled: (\{.*\})", sc)
        nother = sum(int(v) for v in re.findall(r": (\d+)", other.group(1))) if other else 0
        states = re.findall(r"GoldenSniffer: UE 0x[0-9a-f]+ (DISCOVERED -> CONVERGED|re-converged|STALE -> CONVERGED)",
                            open(lg).read())
        rows.append((tag, ssb / 1e6, sent, found, extra, dup, nother, len(states), dt))
        print("%-18s SSB %9.3f  sent %6d  found %6d (%5.1f%%)  extra %d  dup %d  other-RNTI DCIs %d  converged %d  "
              "(%.0f s)" % (tag, ssb / 1e6, sent, found, 100.0 * found / max(sent, 1), extra, dup, nother,
                            len(states), dt), flush=True)
print("\n%-18s %8s %8s %7s %6s %5s %6s" % ("run", "sent", "found", "%", "extra", "dup", "other"))
for r in rows:
    print("%-18s %8d %8d %6.1f%% %6d %5d %6d" % (r[0], r[2], r[3], 100.0 * r[3] / max(r[2], 1), r[4], r[5], r[6]))
