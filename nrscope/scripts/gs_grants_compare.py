#!/usr/bin/env python3
"""Score GoldenSniffer's DL grants (golden_sniffer_grants.csv: the PRBs of tracked UEs' DCI
1_1, from the FDRA found on air) against a srsRAN/OCUDU gNB log ("Slot decisions", DL
grants with rb=[a..b), BWP-relative).

usage: gs_grants_compare.py gnb.log golden_sniffer_grants.csv [--scs-khz 30] [--rnti 0x4611]
Per RNTI: grants reported, how many the gNB sent in that slot to that RNTI, how many have
exactly its PRBs, and the mean overlap (Jaccard) of the others.
"""
import argparse
import csv
import re
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("log")
ap.add_argument("grants")
ap.add_argument("--scs-khz", type=int, default=30)
ap.add_argument("--rnti", default="")
a = ap.parse_args()
spf = 10 * a.scs_khz // 15

sched_re = re.compile(r"\[\s*(\d+)\.(\d+)\]\s+Slot decisions")
dl_re = re.compile(r"DL: ue=\d+ (?:c-|tc-)?rnti=(0x[0-9a-fA-F]+).*?rb=\[(\d+)\.\.(\d+)\)")
sent = defaultdict(list)  # (sfn, slot, rnti) -> [(a, b)]
for line in open(a.log, errors="replace"):
    m = sched_re.search(line)
    if not m:
        continue
    sfn, slot = int(m.group(1)), int(m.group(2))
    for seg in re.split(r"(?=\b(?:DL|UL): ue=)", line):
        g = dl_re.match(seg.strip())
        if g:
            sent[(sfn, slot, g.group(1).lower())].append((int(g.group(2)), int(g.group(3))))

rows = list(csv.DictReader(open(a.grants)))
want = set(x.lower() for x in a.rnti.split(",") if x)
stat = defaultdict(lambda: {"n": 0, "sent": 0, "exact": 0, "jac": []})
for r in rows:
    rnti = r["rnti"].lower()
    if want and rnti not in want:
        continue
    st = stat[rnti]
    st["n"] += 1
    key = (int(r["sfn"]), int(r["slot"]), rnti)
    if key not in sent:
        continue
    st["sent"] += 1
    lo, n = int(r["rb_start"]), int(r["nof_rb"])
    mine = set(range(lo, lo + n))
    best = 0.0
    for (ga, gb) in sent[key]:
        theirs = set(range(ga, gb))
        j = len(mine & theirs) / max(len(mine | theirs), 1)
        best = max(best, j)
    if best == 1.0:
        st["exact"] += 1
    else:
        st["jac"].append(best)
for rnti, st in sorted(stat.items()):
    others = sum(st["jac"]) / len(st["jac"]) if st["jac"] else float("nan")
    print("RNTI %s  grants %6d  in the gNB log %6d  exact PRBs %6d (%5.1f%%)  others: mean overlap %.2f" %
          (rnti, st["n"], st["sent"], st["exact"], 100.0 * st["exact"] / max(st["sent"], 1), others))
