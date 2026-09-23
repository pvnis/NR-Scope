#!/usr/bin/env python3
"""Cross-validate NR-Scope's decoded DCIs against the gNB's own scheduler log.

Both sides stamp every transmission with the same SFN.slot, which gives an exact
join key -- but SFN only counts 0..1023, so it repeats every 10.24 s. Any capture
longer than that aliases: grants minutes apart collapse onto the same key and get
"matched" against each other, producing garbage. So the wall clock is used to
pick the right SFN cycle, and SFN.slot for the precise phase within it. Clock
accuracy only has to be better than +/- 5 s for this to work.

Assumes both machines keep UTC (use --skew if they do not).

Usage:
    python3 validate_dci.py gnb.log logs_nrscope.csv [--rnti 0x4601]
"""

import argparse
import csv
import datetime as dt
import re
import sys
from collections import Counter, defaultdict

SLOTS_PER_FRAME = 20          # 30 kHz SCS
SLOT_SECONDS = 0.0005         # 0.5 ms
CYCLE_SLOTS = 1024 * SLOTS_PER_FRAME
CYCLE_SECONDS = CYCLE_SLOTS * SLOT_SECONDS      # 10.24 s

GNB_PDCCH_RE = re.compile(
    r"\[\s*\d+\.\d+\]\s+PDCCH:\s+rnti=0x(?P<rnti>[0-9a-fA-F]+)\s+"
    r"ss_id=(?P<ss>\d+)\s+format=(?P<fmt>\S+)"
)

GNB_RE = re.compile(
    r"^(?P<time>\d{4}-\d{2}-\d{2}T[\d:.]+).*?"
    r"\[\s*(?P<sfn>\d+)\.(?P<slot>\d+)\]\s+(?P<ch>PDSCH|PUSCH):\s+"
    r"rnti=0x(?P<rnti>[0-9a-fA-F]+).*?"
    r"h_id=(?P<harq>\d+).*?"
    r"prb=\[(?P<prb0>\d+),\s*(?P<prb1>\d+)\)\s+"
    r"symb=\[(?P<sym0>\d+),\s*(?P<sym1>\d+)\)\s+"
    r"mod=(?P<mod>\S+)\s+rv=(?P<rv>\d+).*?tbs=(?P<tbs>\d+)"
)

FIELDS = ["prb_start", "prb_len", "sym_start", "sym_len", "rv", "harq"]


def unwrap(epoch, sfn, slot, ref):
    """Absolute slot counter. SFN.slot fixes the phase; the clock picks the cycle."""
    phase = (int(sfn) * SLOTS_PER_FRAME + int(slot)) % CYCLE_SLOTS
    approx = (epoch - ref) / SLOT_SECONDS
    cycle = round((approx - phase) / CYCLE_SLOTS)
    return cycle * CYCLE_SLOTS + phase


def parse_gnb(path, want):
    recs = []
    # The gNB logs a grant's PDCCH just before its PDSCH/PUSCH, so the most
    # recent PDCCH for that RNTI and direction tells us which search space the
    # grant used: format *_1 is the UE-specific (dedicated) space, *_0 the
    # common one. NR-Scope only decodes the dedicated space, so the two must be
    # counted separately or its recall looks far worse than it is.
    last_pdcch = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            p = GNB_PDCCH_RE.search(line)
            if p:
                fmt = p.group("fmt")
                last_pdcch[(int(p.group("rnti"), 16),
                            "DL" if fmt.startswith("1") else "UL")] = fmt
                continue
            m = GNB_RE.match(line)
            if not m:
                continue
            rnti = int(m.group("rnti"), 16)
            if want is not None and rnti != want:
                continue
            t = dt.datetime.fromisoformat(m.group("time")).replace(tzinfo=dt.timezone.utc)
            direction = "DL" if m.group("ch") == "PDSCH" else "UL"
            fmt = last_pdcch.get((rnti, direction), "?")
            recs.append({
                "fmt": fmt,
                "dedicated": fmt.endswith("_1"),
                "epoch": t.timestamp(), "sfn": int(m.group("sfn")), "slot": int(m.group("slot")),
                "rnti": rnti, "dir": "DL" if m.group("ch") == "PDSCH" else "UL",
                "prb_start": int(m.group("prb0")),
                "prb_len": int(m.group("prb1")) - int(m.group("prb0")),
                "sym_start": int(m.group("sym0")),
                "sym_len": int(m.group("sym1")) - int(m.group("sym0")),
                "mod": m.group("mod").upper(), "rv": int(m.group("rv")),
                "tbs_bits": int(m.group("tbs")) * 8,   # gNB logs BYTES, NR-Scope logs BITS
                "harq": int(m.group("harq")),
            })
    return recs


def parse_nrscope(path, want):
    recs = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if not r.get("rnti") or r["rnti"] == "rnti":
                continue
            rnti = int(r["rnti"])
            if want is not None and rnti != want:
                continue
            # NR-Scope logs the PDCCH slot; k is the offset to the data slot.
            total = int(r["system_frame_index"]) * SLOTS_PER_FRAME + int(r["slot_index"]) + int(r["k"])
            recs.append({
                "epoch": float(r["timestamp"]),
                "sfn": (total // SLOTS_PER_FRAME) % 1024, "slot": total % SLOTS_PER_FRAME,
                "rnti": rnti, "dir": "DL" if r["dci_format"].startswith("1") else "UL",
                "prb_start": int(r["frequency_start"]), "prb_len": int(r["frequency_length"]),
                "sym_start": int(r["time_start"]), "sym_len": int(r["time_length"]),
                "mod": r["modulation"].upper(), "rv": int(r["redundancy_version"]),
                "tbs_bits": int(r["transport_block_size"]), "harq": int(r["harq_id"]),
            })
    return recs


def estimate_skew(gnb, nrs):
    """Recover the clock offset between the two machines.

    SFN.slot is generated by the radio, not either host clock, so it is the same
    on both sides regardless of how far the machines' clocks have drifted. Match
    on SFN.slot plus identical allocation fields -- a coincidence across eight
    fields is implausible -- then the spread of the time differences tells us
    whether we have found a genuine constant offset.

    Returns (skew, n_matches, spread) or None.
    """
    G = {}
    for r in gnb:
        G.setdefault((r["sfn"], r["slot"], r["dir"], r["rnti"]), []).append(r)
    diffs = []
    for r in nrs:
        for c in G.get((r["sfn"], r["slot"], r["dir"], r["rnti"]), []):
            if all(c[f] == r[f] for f in FIELDS):
                diffs.append(c["epoch"] - r["epoch"])
                break
    if len(diffs) < 5:
        return None
    diffs.sort()
    med = diffs[len(diffs) // 2]
    # Ignore outliers from SFN aliasing: keep differences within one SFN cycle.
    near = [d for d in diffs if abs(d - med) < CYCLE_SECONDS / 2]
    if len(near) < 5:
        return None
    return med, len(near), max(near) - min(near)


def clusters(recs, gap=60):
    """Split records into separate capture runs (NR-Scope appends to its CSV)."""
    if not recs:
        return []
    recs = sorted(recs, key=lambda r: r["epoch"])
    out, cur = [], [recs[0]]
    for r in recs[1:]:
        if r["epoch"] - cur[-1]["epoch"] > gap:
            out.append(cur)
            cur = [r]
        else:
            cur.append(r)
    out.append(cur)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gnb_log")
    ap.add_argument("nrscope_csv")
    ap.add_argument("--rnti", help="hex (0x4601) or decimal; default: all")
    ap.add_argument("--skew", type=float, default=0.0,
                    help="seconds to add to NR-Scope timestamps if the machines' clocks differ")
    ap.add_argument("--show", type=int, default=5)
    a = ap.parse_args()

    want = int(a.rnti, 0) if a.rnti else None
    for path, what in ((a.gnb_log, "gNB log"), (a.nrscope_csv, "NR-Scope CSV")):
        try:
            open(path).close()
        except OSError as e:
            sys.exit(f"Cannot read the {what}: {e.strerror}: {path}\n"
                     f"(NR-Scope only writes logs_nrscope.csv once it has decoded a DCI, "
                     f"so a run that never locked onto the cell leaves no file.)")
    gnb = parse_gnb(a.gnb_log, None)
    nrs = parse_nrscope(a.nrscope_csv, None)
    for r in nrs:
        r["epoch"] += a.skew

    if want is None:
        # Pick the C-RNTI both sides saw most. One-row RNTIs in the CSV are
        # blind-decode false positives, so the dominant one is the real UE.
        gc = Counter(r["rnti"] for r in gnb)
        nc = Counter(r["rnti"] for r in nrs)
        common = set(gc) & set(nc)
        if not common:
            sys.exit("No RNTI appears in both logs. Are they from the same run?")
        want = max(common, key=lambda x: min(gc[x], nc[x]))
        others = [f"0x{r:04x}" for r in sorted(set(gc) | set(nc), key=lambda x: -nc.get(x, 0))
                  if r != want and nc.get(r, 0) > 1]
        if others:
            print(f"also present: {' '.join(others)}  (--rnti to select)")

    gnb = [r for r in gnb if r["rnti"] == want]
    nrs = [r for r in nrs if r["rnti"] == want]
    if not gnb or not nrs:
        sys.exit(f"RNTI 0x{want:04x} not present in both logs.")

    # The two machines need not agree on the time of day. SFN.slot comes from the
    # radio, so it lets us measure the offset and correct for it rather than
    # wrongly declaring a simultaneous capture non-overlapping.
    if not a.skew:
        est = estimate_skew(gnb, nrs)
        if est:
            skew, n, spread = est
            if abs(skew) > 1.0:
                print(f"clock skew detected: the gNB host is {abs(skew):.1f} s "
                      f"{'behind' if skew < 0 else 'ahead of'} this one "
                      f"({n} SFN-matched grants, spread {spread*1000:.0f} ms). Correcting.\n")
                for r in nrs:
                    r["epoch"] += skew

    if not gnb:
        sys.exit("No PDSCH/PUSCH lines parsed from the gNB log.")
    if not nrs:
        sys.exit("No DCI rows parsed from the NR-Scope CSV.")

    fmt = lambda e: dt.datetime.fromtimestamp(e, dt.timezone.utc).strftime("%H:%M:%S")

    runs = clusters(nrs)
    if len(runs) > 1:
        print(f"NOTE: the CSV holds {len(runs)} separate runs (NR-Scope appends). "
              f"Delete logs_nrscope.csv between runs.")
        for c in runs:
            print(f"      {fmt(c[0]['epoch'])} -> {fmt(c[-1]['epoch'])}  ({len(c)} rows)")
        print()

    g_lo, g_hi = min(r["epoch"] for r in gnb), max(r["epoch"] for r in gnb)
    # keep the NR-Scope run that overlaps the gNB log; if none does, the closest
    # one, so the error message names the run the user most likely meant.
    def overlap_score(c):
        n = sum(1 for r in c if g_lo - 1 <= r["epoch"] <= g_hi + 1)
        gap = max(g_lo - c[-1]["epoch"], c[0]["epoch"] - g_hi, 0)
        return (n, -gap)

    nrs = max(runs, key=overlap_score)

    lo = max(g_lo, min(r["epoch"] for r in nrs))
    hi = min(g_hi, max(r["epoch"] for r in nrs))
    if hi <= lo:
        print(f"gNB log      : {fmt(g_lo)} -> {fmt(g_hi)}")
        print(f"NR-Scope run : {fmt(nrs[0]['epoch'])} -> {fmt(nrs[-1]['epoch'])}")
        print()
        sys.exit("The two captures DO NOT OVERLAP in wall-clock time -- nothing to compare.\n"
                 "Record the gNB log so that it spans the NR-Scope capture, then retry.")

    print(f"rnti 0x{want:04x}   compared over {fmt(lo)}-{fmt(hi)} ({hi-lo:.0f} s, "
          f"where both were running)")
    print()

    ref = lo
    g = {}
    for r in gnb:
        if lo <= r["epoch"] <= hi:
            g.setdefault((unwrap(r["epoch"], r["sfn"], r["slot"], ref), r["rnti"], r["dir"]), r)
    n = {}
    for r in nrs:
        if lo <= r["epoch"] <= hi:
            n.setdefault((unwrap(r["epoch"], r["sfn"], r["slot"], ref), r["rnti"], r["dir"]), r)

    both = sorted(set(g) & set(n))
    ded = [k for k in g if g[k]["dedicated"]]
    fb = [k for k in g if not g[k]["dedicated"]]
    ded_hit = sum(1 for k in ded if k in n)
    fb_hit = sum(1 for k in fb if k in n)
    false_pos = len(set(n) - set(g))

    print(f"{'':<22}{'gNB':>6}{'matched':>10}")
    print(f"  {'dedicated (1_1/0_1)':<20}{len(ded):>6}{ded_hit:>10}"
          f"{'   <- what NR-Scope decodes' if ded else ''}")
    print(f"  {'fallback  (1_0/0_0)':<20}{len(fb):>6}{fb_hit:>10}"
          f"{'   not implemented' if fb else ''}")
    print()

    if not both:
        sys.exit("Nothing matched. If the machines' clocks differ by more than ~5 s, pass --skew.")

    ALLOC = ["prb_start", "prb_len", "sym_start", "sym_len"]
    alloc_ok = sum(1 for k in both if all(g[k][f] == n[k][f] for f in ALLOC))
    hq_ok = sum(1 for k in both if g[k]["harq"] == n[k]["harq"] and g[k]["rv"] == n[k]["rv"])
    tbs_ok = sum(1 for k in both if g[k]["tbs_bits"] == n[k]["tbs_bits"])
    mod_ok = sum(1 for k in both if g[k]["mod"] == n[k]["mod"])
    pct = lambda x: f"{x:>3}/{len(both):<3} {100*x/len(both):>6.1f}%"

    print(f"  PRB + symbol allocation   {pct(alloc_ok)}   <- DMRS positions depend on this")
    print(f"  HARQ id / RV              {pct(hq_ok)}")
    print(f"  modulation                {pct(mod_ok)}")
    print(f"  transport block size      {pct(tbs_ok)}")
    print(f"  false positives           {false_pos:>3}")
    print()

    pairs = Counter((g[k]["mod"], n[k]["mod"]) for k in both if g[k]["mod"] != n[k]["mod"])
    if any(x == "256QAM" for x, _ in pairs):
        print(f"  note: {sum(pairs.values())} grants use mcs-Table qam256 but NR-Scope assumed "
              f"64QAM\n        (also skews their TBS) -- set rrc_recfg_config in config.yaml")
        print()

    bad = Counter()
    ex = defaultdict(list)
    for k in both:
        for f in ALLOC:
            if g[k][f] != n[k][f]:
                bad[f] += 1
                if len(ex[f]) < a.show:
                    ex[f].append((k, g[k][f], n[k][f]))
    if bad:
        print("  ALLOCATION ERRORS (these would corrupt DMRS extraction):")
        for f, c in bad.most_common():
            print(f"    {f:<10} {c}")
            for k, gv, nv in ex[f]:
                print(f"      SFN {g[k]['sfn']}.{g[k]['slot']} {k[2]}: gNB={gv} nrscope={nv}")
        print()

    verdict = ("allocations are exact -- safe to use for DMRS"
               if alloc_ok == len(both) else
               f"{len(both)-alloc_ok} of {len(both)} allocations are WRONG -- not safe for DMRS")
    print(f"VERDICT: {verdict}")


if __name__ == "__main__":
    main()
