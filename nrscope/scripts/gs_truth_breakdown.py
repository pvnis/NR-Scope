#!/usr/bin/env python3
"""Where GoldenSniffer misses, against gs_cellsim's ground truth.

usage: gs_truth_breakdown.py truth.csv golden_sniffer.csv [--rnti 0x1234 [--window 0.5]]

A DCI counts as found when the CSV has one of the same RNTI, slot, direction (first payload
bit) and size. Without --rnti: found / sent by AL, format, search space (uss/css) and DCIs
in the slot, the per-UE spread, and which GoldenSniffer path found them (search: the
whole-slot blind search; tracked: the direct decode of CONVERGED UEs). With --rnti: that
UE over time, per window, with the path of each hit and the AL of each miss.
"""
import argparse
import csv
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("truth")
ap.add_argument("csv")
ap.add_argument("--rnti", default="")
ap.add_argument("--window", type=float, default=0.5)
a = ap.parse_args()
truth = [r for r in csv.DictReader(open(a.truth)) if r["rnti"] != "0xffff"]
gs = list(csv.DictReader(open(a.csv)))
spf = 20 if any(int(r["slot"]) >= 10 for r in csv.DictReader(open(a.truth))) else 10


def unwrap(rows):
    prev, out = None, []
    for r in rows:
        f = int(r["sfn"])
        if prev is not None:
            base = prev - prev % 1024
            f = min((base + d + f for d in (-1024, 0, 1024)), key=lambda v: abs(v - prev))
        prev = f
        out.append(f * spf + int(r["slot"]))
    return out


got = defaultdict(list)
gsl = unwrap(gs)
for s, r in zip(gsl, gs):
    got[(s, r["rnti"].lower())].append(r)
lo, hi = (min(gsl), max(gsl)) if gsl else (0, -1)


def match(s, r):
    want = 1 if r["dir"] == "DL" else 0
    for g in got.get((s, r["rnti"].lower()), []):
        if g["payload_hex"] and (int(g["payload_hex"][0], 16) >> 3) == want and int(g["size"]) == int(r["size"]):
            return g
    return None


ts = unwrap(truth)
if a.rnti:
    tot, ok, why = Counter(), Counter(), defaultdict(Counter)
    for s, r in zip(ts, truth):
        if r["rnti"].lower() != a.rnti.lower():
            continue
        b = int(float(r["t_s"]) / a.window)
        g = match(s, r)
        tot[b] += 1
        ok[b] += g is not None
        why[b][g["source"] if g else "miss AL" + r["al"]] += 1
    for b in sorted(tot):
        print("%6.2f-%6.2f s: %4d/%4d (%3.0f%%) %s" % (b * a.window, (b + 1) * a.window, ok[b], tot[b],
                                                    100.0 * ok[b] / tot[b], dict(why[b])))
    raise SystemExit
load = Counter(ts)
by, src, ue = defaultdict(Counter), Counter(), Counter()
for s, r in zip(ts, truth):
    if not lo <= s <= hi:
        continue
    g = match(s, r)
    for k in ("al", "format", "ss"):
        by[k][(r[k], g is not None)] += 1
    by["dcis/slot"][(min(load[s], 8), g is not None)] += 1
    ue[(r["rnti"], g is not None)] += 1
    if g:
        src[g["source"]] += 1
for k in ("al", "format", "ss", "dcis/slot"):
    vals = sorted({x[0] for x in by[k]}, key=lambda v: (len(str(v)), str(v)))
    print("%-9s " % k + "  ".join("%s: %d/%d (%.0f%%)" % (v, by[k][(v, True)], by[k][(v, True)] + by[k][(v, False)],
          100.0 * by[k][(v, True)] / max(1, by[k][(v, True)] + by[k][(v, False)])) for v in vals))
rates = sorted((100.0 * ue[(u, True)] / max(1, ue[(u, True)] + ue[(u, False)]), u) for u in {x[0] for x in ue})
if rates:
    print("per UE: min %.0f%% (%s), median %.0f%%, max %.0f%%" % (rates[0][0], rates[0][1], rates[len(rates) // 2][0],
                                                               rates[-1][0]))
print("found by:", dict(src))
