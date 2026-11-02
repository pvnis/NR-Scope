#!/usr/bin/env python3
"""Score GoldenSniffer's DCIs (golden_sniffer.csv) against a srsRAN/OCUDU gNB log.

Every "PDCCH: rnti=... format=... al=..." line of the gNB log is a DCI that was sent.
The CSV's SFNs are aligned to the log's (the log spans several SFN wraps: the wrap with
the most matching (SFN, slot, RNTI) is taken), then, over the slots the CSV covers:
  found    = DCIs sent and reported (same slot and RNTI)
  missed   = sent, not reported
  extra    = reported, never sent (wrong RNTI or slot: false detections)
per RNTI and format; plus whether the AL and the DCI direction (first payload bit, the
format identifier of 0_0/1_0) agree.

The log's "PDCCH: rnti=... format=... al=..." lines are used, or, if the PHY does not log
them (OCUDU at PHY info), the MAC "Slot decisions" grants (DL/UL per RNTI, no AL).

usage: gs_compare.py gnb.log golden_sniffer.csv [--scs-khz 15] [--rnti 0x4601,0x4602]
"""
import argparse
import csv
import re
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("log")
ap.add_argument("csv")
ap.add_argument("--scs-khz", type=int, default=15)
ap.add_argument("--rnti", default="", help="comma separated; default: every C-RNTI in the CSV")
a = ap.parse_args()
spf = 10 * a.scs_khz // 15

line_re = re.compile(r"\[\s*(\d+)\.(\d+)\]\s+PDCCH:\s+rnti=(0x[0-9a-fA-F]+).*?format=(\S+).*?al=(\d+)")
def unwrapper():
    """SFN -> absolute frame, nearest to the previous one (rows can come out of order)"""
    st = {"prev": None}

    def f(sfn):
        if st["prev"] is None:
            st["prev"] = sfn
            return sfn
        base = st["prev"] - st["prev"] % 1024
        best = min((base + d + sfn for d in (-1024, 0, 1024)), key=lambda v: abs(v - st["prev"]))
        st["prev"] = best
        return best
    return f


# OCUDU/srsRAN MAC "Slot decisions" lines, when the PHY does not log PDCCH: one grant
# per "DL: ue=N c-rnti=..." / "UL: ue=N rnti=..." entry; format by direction, AL unknown.
sched_re = re.compile(r"\[\s*(\d+)\.(\d+)\]\s+Slot decisions")
grant_re = re.compile(r"(DL|UL): ue=\d+ (?:c-|tc-)?rnti=(0x[0-9a-fA-F]+)")

sent = []  # (abs_slot, rnti, fmt, al)
uw = unwrapper()
uw2 = unwrapper()
sched = []
for line in open(a.log, errors="replace"):
    if "PDCCH:" in line:
        m = line_re.search(line)
        if not m:
            continue
        sfn, slot = int(m.group(1)), int(m.group(2))
        sent.append((uw(sfn) * spf + slot, m.group(3).lower(), m.group(4), int(m.group(5))))
    elif "Slot decisions" in line:
        m = sched_re.search(line)
        if not m:
            continue
        sfn, slot = int(m.group(1)), int(m.group(2))
        base = uw2(sfn) * spf + slot
        # one segment per grant; an UL grant is logged at its PUSCH slot, its DCI is k2 earlier
        for seg in re.split(r"(?=\b(?:DL|UL): ue=)", line):
            m2 = grant_re.match(seg.strip())
            if not m2:
                continue
            d, r = m2.group(1), m2.group(2)
            k2 = re.search(r"\bk2=(\d+)", seg)
            at = base - (int(k2.group(1)) if d == "UL" and k2 else 0)
            sched.append((at, r.lower(), "1_x" if d == "DL" else "0_x", -1))
if not sent:
    print("no PDCCH lines in the log: using the MAC slot decisions (AL not logged)")
    sent = sched

rows = list(csv.DictReader(open(a.csv)))
got = []
uw = unwrapper()
for r in rows:
    sfn, slot = int(r["sfn"]), int(r["slot"])
    bit0 = int(r["payload_hex"][0], 16) >> 3 if r["payload_hex"] else -1
    got.append((uw(sfn) * spf + slot, r["rnti"].lower(), int(r["al"]), bit0, r["kind"]))

rntis = set(x.lower() for x in a.rnti.split(",") if x) or {g[1] for g in got if g[4] == "C"}
sent_set = defaultdict(list)
for s in sent:
    sent_set[(s[0], s[1])].append(s)

# SFN wrap alignment
best = (-1, 0)
max_epoch = sent[-1][0] // (1024 * spf) + 1 if sent else 1
for e in range(-2, max_epoch + 2):
    sh = e * 1024 * spf
    n = sum(1 for g in got if (g[0] + sh, g[1]) in sent_set)
    if n > best[0]:
        best = (n, sh)
shift = best[1]
got = [(g[0] + shift,) + g[1:] for g in got]
lo, hi = min(g[0] for g in got), max(g[0] for g in got)
print("aligned on SFN wrap %d; CSV covers %.2f s" % (shift // (1024 * spf), (hi - lo + 1) / (spf * 100.0)))

got_set = defaultdict(list)
for g in got:
    got_set[(g[0], g[1])].append(g)
for rnti in sorted(rntis):
    tot, fnd, al_ok, dir_ok = Counter(), Counter(), Counter(), Counter()
    for s in sent:
        if s[1] != rnti or not lo <= s[0] <= hi:
            continue
        tot[s[2]] += 1
        cands = got_set.get((s[0], rnti), [])
        if cands:
            want = 1 if s[2].startswith("1") else 0
            same_dir = [g for g in cands if g[3] == want]
            g = same_dir[0] if same_dir else cands[0]
            fnd[s[2]] += 1
            al_ok[s[2]] += s[3] < 0 or g[2] == s[3]
            dir_ok[s[2]] += g[3] == want
    extra = sum(1 for g in got if g[1] == rnti and (g[0], rnti) not in sent_set)
    print("RNTI %s" % rnti)
    for f in sorted(tot):
        print("  format %-4s sent %6d  found %6d (%5.1f%%)  AL agrees %d  direction agrees %d" %
              (f, tot[f], fnd[f], 100.0 * fnd[f] / max(tot[f], 1), al_ok[f], dir_ok[f]))
    T, F = sum(tot.values()), sum(fnd.values())
    print("  all         sent %6d  found %6d (%5.1f%%)  extra (never sent) %d" % (T, F, 100.0 * F / max(T, 1), extra))
other = Counter(g[1] for g in got if g[1] not in rntis and g[4] == "C")
if other:
    print("C-RNTIs reported that the log never scheduled:", dict(other))
