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

SFN wraps (10.24 s) are resolved from the wall clocks when both sides have them (the
log's line timestamps, UTC, and the CSV's timestamp column): SFN continuity alone puts
everything after a gap or a resync longer than 5.12 s in the wrong wrap.
--wrap-by-sfn falls back to continuity.

usage: gs_compare.py gnb.log golden_sniffer.csv [--scs-khz 15] [--rnti 0x4601,0x4602]
"""
import argparse
import csv
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone

ap = argparse.ArgumentParser()
ap.add_argument("log")
ap.add_argument("csv")
ap.add_argument("--scs-khz", type=int, default=15)
ap.add_argument("--rnti", default="", help="comma separated; default: every C-RNTI in the CSV")
ap.add_argument("--wrap-by-sfn", action="store_true", help="resolve SFN wraps by continuity, not by the clocks")
a = ap.parse_args()
spf = 10 * a.scs_khz // 15

line_re = re.compile(r"\[\s*(\d+)\.(\d+)\]\s+PDCCH:\s+rnti=(0x[0-9a-fA-F]+).*?format=(\S+).*?al=(\d+)")
ts_re = re.compile(r"^(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d+)")


def log_time(line):
    """the gNB log line's timestamp (UTC) in s, or None"""
    m = ts_re.match(line)
    if not m:
        return None
    return datetime.fromisoformat(m.group(1)).replace(tzinfo=timezone.utc).timestamp()


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

# entries: [frame, slot offset within the frame's numbering, rnti, fmt, al] plus (time, sfn)
sent, sent_tf = [], []
sched, sched_tf = [], []
uw = unwrapper()
uw2 = unwrapper()
for line in open(a.log, errors="replace"):
    if "PDCCH:" in line:
        m = line_re.search(line)
        if not m:
            continue
        sfn, slot = int(m.group(1)), int(m.group(2))
        sent.append((uw(sfn) * spf + slot, m.group(3).lower(), m.group(4), int(m.group(5))))
        sent_tf.append((log_time(line), sfn, None))  # continuity frame filled in below
    elif "Slot decisions" in line:
        m = sched_re.search(line)
        if not m:
            continue
        sfn, slot = int(m.group(1)), int(m.group(2))
        fr = uw2(sfn)
        base = fr * spf + slot
        t = log_time(line)
        # one segment per grant; an UL grant is logged at its PUSCH slot, its DCI is k2 earlier
        for seg in re.split(r"(?=\b(?:DL|UL): ue=)", line):
            m2 = grant_re.match(seg.strip())
            if not m2:
                continue
            d, r = m2.group(1), m2.group(2)
            k2 = re.search(r"\bk2=(\d+)", seg)
            at = base - (int(k2.group(1)) if d == "UL" and k2 else 0)
            sched.append((at, r.lower(), "1_x" if d == "DL" else "0_x", -1))
            sched_tf.append((t, sfn, fr))
if not sent:
    print("no PDCCH lines in the log: using the MAC slot decisions (AL not logged)")
    sent, sent_tf = sched, sched_tf
else:
    # recompute the continuity frames of the PDCCH lines for the clock correction below
    uw = unwrapper()
    sent_tf = [(t, sfn, uw(sfn)) for t, sfn, _ in sent_tf]

rows = list(csv.DictReader(open(a.csv)))
got, got_tf = [], []
uw = unwrapper()
for r in rows:
    sfn, slot = int(r["sfn"]), int(r["slot"])
    bit0 = int(r["payload_hex"][0], 16) >> 3 if r["payload_hex"] else -1
    fr = uw(sfn)
    got.append((fr * spf + slot, r["rnti"].lower(), int(r["al"]), bit0, r["kind"]))
    got_tf.append((float(r["timestamp"]) if r.get("timestamp") else None, sfn, fr))


def clock_frames(entries, tf, c):
    """moves each entry to the wrap its timestamp says (whole wraps only)"""
    out = []
    for e, (t, sfn, fr) in zip(entries, tf):
        want = sfn + 1024 * round((t * 100.0 - c - sfn) / 1024.0)
        out.append((e[0] + (want - fr) * spf,) + e[1:])
    return out


by_clock = (not a.wrap_by_sfn and sent and got and all(x[0] is not None for x in sent_tf) and
            all(x[0] is not None for x in got_tf))
if by_clock:
    # SFN phase against the clock, from the log: circular median of (100 t - sfn) mod 1024,
    # so the timestamps may be off by up to +-5 s before a wrap flips
    ph = [(t * 100.0 - sfn) % 1024 for t, sfn, _ in sent_tf]
    ref = sorted(ph)[len(ph) // 2]
    c = ref + sorted(((p - ref + 512) % 1024) - 512 for p in ph)[len(ph) // 2]
    sent = clock_frames(sent, sent_tf, c)
    got = clock_frames(got, got_tf, c)
    print("SFN wraps from the clocks")

rntis = set(x.lower() for x in a.rnti.split(",") if x) or {g[1] for g in got if g[4] == "C"}
sent_set = defaultdict(list)
for s in sent:
    sent_set[(s[0], s[1])].append(s)

# SFN wrap alignment (by the clocks: absorbs a clock offset or a lag of whole wraps)
best = (-1, 0)
lo_s = min(s[0] for s in sent) // (1024 * spf) if sent else 0
hi_s = max(s[0] for s in sent) // (1024 * spf) + 1 if sent else 1
lo_g = min(g[0] for g in got) // (1024 * spf) if got else 0
for e in range(lo_s - lo_g - 3, hi_s - lo_g + 3):
    sh = e * 1024 * spf
    n = sum(1 for g in got if (g[0] + sh, g[1]) in sent_set)
    if n > best[0]:
        best = (n, sh)
shift = best[1]
got = [(g[0] + shift,) + g[1:] for g in got]
lo, hi = min(g[0] for g in got), max(g[0] for g in got)
print("aligned with a shift of %d wraps; CSV covers %.2f s" % (shift // (1024 * spf), (hi - lo + 1) / (spf * 100.0)))

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
    dup = sum(n - 1 for n in Counter((g[0], g[1], g[3]) for g in got if g[1] == rnti).values() if n > 1)
    print("RNTI %s" % rnti)
    for f in sorted(tot):
        print("  format %-4s sent %6d  found %6d (%5.1f%%)  AL agrees %d  direction agrees %d" %
              (f, tot[f], fnd[f], 100.0 * fnd[f] / max(tot[f], 1), al_ok[f], dir_ok[f]))
    T, F = sum(tot.values()), sum(fnd.values())
    print("  all         sent %6d  found %6d (%5.1f%%)  extra (never sent) %d  duplicates %d" %
          (T, F, 100.0 * F / max(T, 1), extra, dup))
other = Counter(g[1] for g in got if g[1] not in rntis and g[4] == "C")
if other:
    print("C-RNTIs reported that the log never scheduled:", dict(other))
