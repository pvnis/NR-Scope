#!/usr/bin/env python3
"""Benchmark NR-Scope's DCI decoding against the gNB's own log of what it sent.

The gNB (OCUDU / srsRAN Project) logs, at info level, every PDCCH and PDSCH it
transmits:

    [PHY] [I] [  894.15] PDCCH: rnti=0x4602 ss_id=2 format=1_1 cce=2 al=2 ...
    [PHY] [I] [  894.15] PDSCH: rnti=0x4602 h_id=1 k1=4 prb=[0, 2) symb=[1, 14) mod=64QAM rv=0 bg=BG1 tbs=372 ...
    [SCHED] [I] [  894.15] Slot decisions ... DL: ue=1 c-rnti=0x4602 h_id=1 ss_id=2 rb=[0..2) ... tbs=372 ri=2 ...

Those are the ground truth. NR-Scope's recording_mode writes one row per DCI it
decoded to DCIs/dci_<run>_pci<N>.csv. This script joins the two per RNTI and per
slot and reports:

  * detection: of the downlink DCIs the gNB sent to a UE while NR-Scope was
    tracking it, how many NR-Scope found (recall), and how many NR-Scope rows
    have no gNB counterpart (false positives);
  * accuracy: for every matched DCI, whether each decoded field agrees with the
    gNB: aggregation level, CCE, HARQ id, RV, PRBs, symbols, modulation, TBS and
    number of layers. The PDSCH allocation is what the sensing path needs to
    regenerate DM-RS, so a mismatch here matters even when detection is perfect;
  * a timeline that puts the gNB's RRC events (rrcReconfiguration, rrcRelease)
    next to where detection or accuracy changes.

The gNB log is used only to grade the sniffer. Nothing from it feeds back into
NR-Scope.

Aligning the two logs
---------------------
The sniffer and gNB machines have unsynchronised clocks, so rows are joined on
SFN.slot, never on wall-clock time. SFN wraps every 10.24 s, so each side's
wall clock is used only to unwrap its own SFN into a continuous slot count. The
offset between the two clocks is then estimated from the data: every pair of
rows with the same RNTI and SFN.slot votes for (t_nrscope - t_gnb), and true
pairs agree to within milliseconds while chance pairs scatter in steps of
10.24 s. The estimate and its support are printed; --clock-offset overrides it.

usage:
  dci_vs_gnb.py GNB_LOG NRSCOPE_DCI_CSV [--msg4 NRSCOPE_MSG4_CSV] [--rnti 0x4602]
                [--expiry 5] [--pdcch NRSCOPE_PDCCH_CSV] [--out-csv joined.csv] [--json summary.json]

Why a DCI was missed
--------------------
With --pdcch, the PDCCH candidate CSV that NR-Scope writes when
log_config.record_pdcch_candidates is on is looked up at the exact position
(slot, aggregation level, CCE) where the gNB sent each DCI. Each miss then falls
in one class:

  slot not searched  NR-Scope evaluated no candidate for this RNTI in that slot
                     (slot dropped, skipped as uplink, or RNTI not yet known)
  not a candidate    the slot was searched, but not at the gNB's position:
                     NR-Scope's search space or candidate hashing differs
  DM-RS below thr.   searched there, but the PDCCH DM-RS energy or correlation
                     was under the threshold: the PDCCH was not where or as
                     NR-Scope expected it (timing, frequency, fade)
  CRC failed         the DM-RS matched but the decoded bits failed the CRC:
                     wrong DCI size, scrambling, or a poor channel estimate

The same measurements for the detected DCIs are printed as the baseline.
"""
import argparse
import csv
import datetime
import json
import re
import statistics
import sys
from collections import Counter, defaultdict

SLOTS_PER_FRAME = 20  # 30 kHz SCS; set with --scs-khz
SFN_PERIOD = 1024
DL_FORMATS = ("1_0", "1_1", "1_2")

RE_TS = r"^(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d+)"
RE_SLOT = r"\[\s*(\d+)\.(\d+)\]"
RE_PDCCH = re.compile(RE_TS + r" \[PHY\s*\] \[\w\] " + RE_SLOT +
                      r" PDCCH: rnti=0x([0-9a-fA-F]+) ss_id=(\d+) format=(\d_\d) cce=(\d+) al=(\d+)")
RE_PDSCH = re.compile(RE_TS + r" \[PHY\s*\] \[\w\] " + RE_SLOT +
                      r" PDSCH: rnti=0x([0-9a-fA-F]+) h_id=(\d+) k1=(\d+) prb=\[(\d+), (\d+)\) "
                      r"symb=\[(\d+), (\d+)\) mod=(\w+) rv=(\d+) bg=\w+ tbs=(\d+)")
RE_SCHED = re.compile(RE_TS + r" \[SCHED\s*\] \[\w\] " + RE_SLOT + r" Slot decisions")
RE_SCHED_DL = re.compile(r"DL: ue=\d+ c-rnti=0x([0-9a-fA-F]+) h_id=(\d+) ss_id=(\d+) rb=\[(\d+)\.\.(\d+)\) "
                         r"k1=\d+ newtx=(\w+) rv=(\d+) tbs=(\d+) ri=(\d+)")
RE_RRC = re.compile(RE_TS + r" \[RRC\s*\] \[\w\] ue=\d+ c-rnti=0x([0-9a-fA-F]+): (\w+ \w+) (\w+)")


def gnb_time(s):
    return datetime.datetime.fromisoformat(s).timestamp()


def parse_gnb(path):
    """PDCCH, PDSCH and SCHED DL decisions keyed by (rnti, sfn, slot), plus RRC events."""
    pdcch, pdsch, sched, rrc = [], [], [], []
    with open(path, errors="replace") as f:
        for line in f:
            if "] PDCCH: " in line:
                m = RE_PDCCH.match(line)
                if m:
                    pdcch.append(dict(t=gnb_time(m[1]), sfn=int(m[2]), slot=int(m[3]), rnti=int(m[4], 16),
                                      ss_id=int(m[5]), format=m[6], cce=int(m[7]), al=int(m[8])))
            elif "] PDSCH: " in line:
                m = RE_PDSCH.match(line)
                if m:
                    pdsch.append(dict(t=gnb_time(m[1]), sfn=int(m[2]), slot=int(m[3]), rnti=int(m[4], 16),
                                      h_id=int(m[5]), prb=(int(m[7]), int(m[8])), symb=(int(m[9]), int(m[10])),
                                      mod=m[11], rv=int(m[12]), tbs_bits=8 * int(m[13])))
            elif "Slot decisions" in line:
                m = RE_SCHED.match(line)
                if m:
                    for d in RE_SCHED_DL.finditer(line):
                        sched.append(dict(t=gnb_time(m[1]), sfn=int(m[2]), slot=int(m[3]), rnti=int(d[1], 16),
                                          h_id=int(d[2]), ri=int(d[9])))
            elif "[RRC" in line:
                m = RE_RRC.match(line)
                if m:
                    rrc.append(dict(t=gnb_time(m[1]), rnti=int(m[2], 16), channel=m[3], msg=m[4]))
    return pdcch, pdsch, sched, rrc


def parse_nrscope(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["direction"] != "DL":
                continue
            prbs = parse_prbs(r["prbs"])
            ts, tl = int(r["time_start"] or 0), int(r["time_length"] or 0)
            rows.append(dict(t=float(r["timestamp"]), sfn=int(r["sfn"]), slot=int(r["slot"]), rnti=int(r["rnti"]),
                             format=r["dci_format"], ss=r["ss_type"], al=int(r["aggregation_level"]),
                             cce=int(r["cce"]), h_id=int(r["harq_id"]), rv=int(r["dci_rv"]), k0=int(r["k"] or 0),
                             prb=prbs, symb=(ts, ts + tl), mod=r["modulation"], tbs_bits=int(r["tbs"] or 0),
                             layers=int(r["nof_layers"] or 0), mcs_table=r["mcs_table"], dci_mcs=r["dci_mcs"]))
    return rows


def parse_prbs(s):
    """'0-1' -> (0, 2); contiguous allocations only, anything else is kept as the raw string."""
    m = re.fullmatch(r"(\d+)-(\d+)", s.strip())
    if m:
        return (int(m[1]), int(m[2]) + 1)
    m = re.fullmatch(r"(\d+)", s.strip())
    if m:
        return (int(m[1]), int(m[1]) + 1)
    return s


STAGE_RANK = {"no_measure": 0, "epre": 1, "corr": 2, "crc_fail": 3, "crc_ok": 4}


def parse_pdcch(path, keep):
    """Candidate rows whose (rnti, sfn, slot) is in keep; the file can be millions of rows."""
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (int(r["rnti"]), int(r["sfn"]), int(r["slot"]))
            if key not in keep:
                continue
            rows.append(dict(t=float(r["timestamp"]), rnti=key[0], sfn=key[1], slot=key[2], ca=r["ca_variant"] == "1",
                             al=int(r["aggregation_level"]), cce=int(r["cce"]), epre=float(r["epre_dBfs"]),
                             corr=float(r["norm_corr"]), sync_us=float(r["sync_error_us"]), stage=r["stage"]))
    return rows


def diagnose(pd_at_slot, pd_at_pos, rnti, g):
    """Best candidate measurement at the gNB's position, and the class of the outcome."""
    here = pd_at_pos.get((rnti, g["abs"], g["al"], g["cce"]))
    if here:
        best = max(here, key=lambda c: (STAGE_RANK[c["stage"]], c["corr"]))
        cls = {"crc_ok": "decoded", "crc_fail": "CRC failed"}.get(best["stage"], "DM-RS below thr.")
        return cls, best
    if pd_at_slot.get((rnti, g["abs"])):
        return "not a candidate", None
    return "slot not searched", None


def parse_msg4(path):
    with open(path) as f:
        return [dict(t=float(r["timestamp"]), sfn=int(r["sfn"]), slot=int(r["slot"]), rnti=int(r["c_rnti"] or r["tc_rnti"]))
            for r in csv.DictReader(f)]


class Clock:
    """Turns (wall time on the sniffer's clock, SFN, slot) into one continuous slot index.

    Built on the sniffer's first row, so both sides unwrap against the same
    reference once the gNB's times are shifted onto the sniffer's clock. Any row
    within +-5.12 s of its true air time unwraps correctly."""

    def __init__(self, t0, sfn0, slot0):
        self.ref = t0 - self.in_period(sfn0, slot0) / slot_rate()

    @staticmethod
    def in_period(sfn, slot):
        return sfn * SLOTS_PER_FRAME + slot

    def abs_slot(self, t, sfn, slot):
        s = self.in_period(sfn, slot)
        period = SFN_PERIOD * SLOTS_PER_FRAME
        k = round(((t - self.ref) - s / slot_rate()) / (period / slot_rate()))
        return k * period + s

    def time_of(self, abs_slot):
        return self.ref + abs_slot / slot_rate()


def slot_rate():
    return SLOTS_PER_FRAME * 100.0  # slots per second


def estimate_offset(nr_rows, gnb_rows):
    """Clock offset t_nrscope - t_gnb, voted by rows sharing RNTI and SFN.slot."""
    by_key = defaultdict(list)
    for g in gnb_rows:
        by_key[(g["rnti"], g["sfn"], g["slot"])].append(g["t"])
    votes = [n["t"] - tg for n in nr_rows for tg in by_key.get((n["rnti"], n["sfn"], n["slot"]), ())]
    if not votes:
        return None, 0, 0
    bins = Counter(round(v / 0.05) for v in votes)
    best, support = bins.most_common(1)[0]
    near = [v for v in votes if abs(v - best * 0.05) <= 0.1]
    return statistics.median(near), len(near), len(votes)


def fmt_pct(a, b):
    return f"{100.0 * a / b:5.1f}%" if b else "   n/a"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("gnb_log")
    ap.add_argument("dci_csv")
    ap.add_argument("--msg4", help="NR-Scope msg4 CSV of the same run; marks when each RNTI became known")
    ap.add_argument("--rnti", action="append", help="only these C-RNTIs (hex or decimal); default: every RNTI NR-Scope decoded")
    ap.add_argument("--expiry", type=float, default=5.0,
                    help="the run's rnti_expiry_s: NR-Scope stops searching an RNTI this many seconds after "
                         "its last DCI (default 5)")
    ap.add_argument("--scs-khz", type=int, default=30)
    ap.add_argument("--clock-offset", type=float, help="t_nrscope - t_gnb in seconds, instead of estimating it")
    ap.add_argument("--examples", type=int, default=5, help="mismatch examples printed per field")
    ap.add_argument("--pdcch", help="NR-Scope PDCCH candidate CSV of the same run (record_pdcch_candidates): "
                                    "says why each missed DCI was missed")
    ap.add_argument("--out-csv", help="write one row per gNB DL DCI with its NR-Scope match")
    ap.add_argument("--json", help="write the summary as JSON, for comparing runs")
    a = ap.parse_args()

    global SLOTS_PER_FRAME
    SLOTS_PER_FRAME = 10 * a.scs_khz // 15

    pdcch, pdsch, sched, rrc = parse_gnb(a.gnb_log)
    nr = parse_nrscope(a.dci_csv)
    msg4 = parse_msg4(a.msg4) if a.msg4 else []
    if not nr:
        sys.exit("no downlink DCIs in the NR-Scope CSV")
    if not pdcch:
        sys.exit("no PDCCH lines in the gNB log (PHY must log at info level)")

    rntis = sorted({n["rnti"] for n in nr})
    if a.rnti:
        rntis = [int(x, 0) for x in a.rnti]

    # --- clock alignment -------------------------------------------------
    if a.clock_offset is not None:
        offset, support, nvotes = a.clock_offset, None, None
    else:
        offset, support, nvotes = estimate_offset(nr, [p for p in pdcch if p["format"] in DL_FORMATS])
        if offset is None:
            sys.exit("no RNTI/SFN.slot in common between the two logs: different sessions, or a different cell?")
    nr.sort(key=lambda r: r["t"])
    clk = Clock(nr[0]["t"], nr[0]["sfn"], nr[0]["slot"])
    for rows in (pdcch, pdsch, sched):
        for g in rows:
            g["abs"] = clk.abs_slot(g["t"] + offset, g["sfn"], g["slot"])
    for n in nr:
        n["abs"] = clk.abs_slot(n["t"], n["sfn"], n["slot"])
    for m in msg4:
        m["abs"] = clk.abs_slot(m["t"], m["sfn"], m["slot"])

    def rel(abs_slot, origin):
        return (abs_slot - origin) / slot_rate()

    print(f"gNB log      : {a.gnb_log}")
    print(f"NR-Scope DCIs: {a.dci_csv} ({len(nr)} DL rows)")
    if support is None:
        print(f"clock offset : {offset:+.3f} s (given)")
    else:
        print(f"clock offset : {offset:+.3f} s (t_nrscope - t_gnb), supported by {support} of {nvotes} candidate pairs")
        if support < 0.5 * len(nr):
            print("  WARNING: weak support; check both logs are from the same session")

    pd_at_slot, pd_at_pos = defaultdict(list), defaultdict(list)
    if a.pdcch:
        keep = {(p["rnti"], p["sfn"], p["slot"]) for p in pdcch if p["format"] in DL_FORMATS}
        for c in parse_pdcch(a.pdcch, keep):
            c["abs"] = clk.abs_slot(c["t"], c["sfn"], c["slot"])
            pd_at_slot[(c["rnti"], c["abs"])].append(c)
            pd_at_pos[(c["rnti"], c["abs"], c["al"], c["cce"])].append(c)
        print(f"PDCCH cands. : {a.pdcch} ({sum(len(v) for v in pd_at_slot.values())} rows in slots the gNB used)")

    pdsch_at = {(p["rnti"], p["abs"]): p for p in pdsch}
    ri_at = {(s["rnti"], s["abs"], s["h_id"]): s["ri"] for s in sched}
    summary = {"gnb_log": a.gnb_log, "dci_csv": a.dci_csv, "clock_offset_s": offset, "rntis": {}}
    joined_rows = []

    for rnti in rntis:
        g_dl = sorted((p for p in pdcch if p["rnti"] == rnti and p["format"] in DL_FORMATS), key=lambda p: p["abs"])
        n_dl = sorted((n for n in nr if n["rnti"] == rnti), key=lambda n: n["abs"])
        if not g_dl:
            print(f"\n=== C-RNTI 0x{rnti:04x}: no DL PDCCH in the gNB log")
            continue

        # When NR-Scope knew the RNTI: its Msg4, else its first DCI. Everything
        # the gNB sent before (the RACH, Msg4 itself) is out of scope.
        known_rows = [m for m in msg4 if m["rnti"] == rnti]
        known = [m["abs"] for m in known_rows]
        start = min(known) + 1 if known else (n_dl[0]["abs"] if n_dl else None)
        if start is None:
            print(f"\n=== C-RNTI 0x{rnti:04x}: NR-Scope never learned this RNTI")
            continue
        last_found = n_dl[-1]["abs"] if n_dl else start
        expiry_at = last_found + int(a.expiry * slot_rate())
        origin = start

        nr_at = defaultdict(list)
        for n in n_dl:
            nr_at[n["abs"]].append(n)

        in_scope = [g for g in g_dl if g["abs"] >= start]
        searched = [g for g in in_scope if g["abs"] <= expiry_at]
        after_expiry = [g for g in in_scope if g["abs"] > expiry_at]

        # --- detection -----------------------------------------------------
        matched, missed = [], []
        used = set()
        for g in searched:
            cands = [n for n in nr_at.get(g["abs"], []) if id(n) not in used]
            if cands:
                n = cands[0]
                used.add(id(n))
                matched.append((g, n))
            else:
                missed.append(g)
        false_pos = [n for n in n_dl if id(n) not in used and n["abs"] >= start]

        print(f"\n=== C-RNTI 0x{rnti:04x}")
        first = min(known_rows, key=lambda m: m["abs"]) if known_rows else n_dl[0]
        print(f"tracked from : SFN {first['sfn']}.{first['slot']} "
              f"({'NR-Scope Msg4' if known_rows else 'first NR-Scope DCI, no msg4 CSV'}); times below are relative to it")
        print(f"               last NR-Scope DCI at +{rel(last_found, origin):.2f} s;"
              f" gNB's last DL DCI at +{rel(g_dl[-1]['abs'], origin):.2f} s")
        print(f"gNB DL DCIs  : {len(in_scope)} after the RNTI was known; {len(searched)} while NR-Scope searched it,"
              f" {len(after_expiry)} after it expired (not searched)")
        print(f"detected     : {len(matched)} / {len(searched)} = {fmt_pct(len(matched), len(searched))} recall")
        print(f"false pos.   : {len(false_pos)} NR-Scope DCIs with no gNB PDCCH for this RNTI in that slot")

        by_kind = defaultdict(lambda: [0, 0])
        for g in searched:
            by_kind[(g["format"], g["ss_id"])][1] += 1
        for g, _ in matched:
            by_kind[(g["format"], g["ss_id"])][0] += 1
        for (fmt, ss), (hit, tot) in sorted(by_kind.items()):
            print(f"  format {fmt} ss_id {ss}: {hit:4d} / {tot:4d} = {fmt_pct(hit, tot)}")
        # Per PDCCH position: a location NR-Scope never decodes points at its view
        # of the CORESET or the search-space hashing rather than at the radio.
        by_loc = defaultdict(lambda: [0, 0])
        for g in searched:
            by_loc[(g["al"], g["cce"])][1] += 1
        for g, _ in matched:
            by_loc[(g["al"], g["cce"])][0] += 1
        print("  by AL/CCE (detected/sent): " + "  ".join(
            f"{al}/{cce}:{hit}/{tot}" + ("!" if hit == 0 else "") for (al, cce), (hit, tot) in sorted(by_loc.items())))

        if missed and not a.pdcch:
            print(f"missed ({len(missed)}), first {min(len(missed), 10)}:")
            for g in missed[:10]:
                print(f"  +{rel(g['abs'], origin):7.3f} s  SFN {g['sfn']:4d}.{g['slot']:<2d} format {g['format']} "
                      f"ss_id {g['ss_id']} al {g['al']} cce {g['cce']}")

        miss_classes = Counter()
        if a.pdcch:
            # Baseline: the winning candidate of every detected DCI
            base = [diagnose(pd_at_slot, pd_at_pos, rnti, g)[1] for g, _ in matched]
            base = [b for b in base if b]
            if base:
                q = lambda xs, f: statistics.quantiles(xs, n=20)[0 if f == 5 else 18] if len(xs) >= 2 else xs[0]
                corr, epre, sync = [b["corr"] for b in base], [b["epre"] for b in base], [b["sync_us"] for b in base]
                print(f"detected DCIs at the gNB's position ({len(base)}): norm_corr median {statistics.median(corr):.3f}"
                      f" (5th pct {q(corr, 5):.3f}), EPRE median {statistics.median(epre):+.1f} dBfs"
                      f" (5th pct {q(epre, 5):+.1f}), sync error median {statistics.median(sync):+.3f} us")
            if missed:
                print(f"missed ({len(missed)}), with what NR-Scope measured at the gNB's position:")
            for g in missed:
                cls, best = diagnose(pd_at_slot, pd_at_pos, rnti, g)
                miss_classes[cls] += 1
                meas = (f"stage {best['stage']:8s} norm_corr {best['corr']:.3f} EPRE {best['epre']:+6.1f} dBfs "
                        f"sync {best['sync_us']:+.3f} us" if best else
                        f"{len(pd_at_slot.get((rnti, g['abs']), []))} other candidates evaluated in that slot")
                print(f"  +{rel(g['abs'], origin):7.3f} s  SFN {g['sfn']:4d}.{g['slot']:<2d} {g['format']} ss_id {g['ss_id']}"
                      f" al {g['al']} cce {g['cce']:2d}: {cls:17s} {meas}")
            if miss_classes:
                print("  misses by class: " + ", ".join(f"{k} {v}" for k, v in miss_classes.most_common()))

        # --- accuracy ------------------------------------------------------
        fields = ["format", "al", "cce", "h_id", "rv", "prb", "symb", "mod", "tbs_bits", "layers"]
        agree = {f: [0, 0] for f in fields}
        examples = defaultdict(list)
        first_bad = {}
        for g, n in matched:
            p = pdsch_at.get((rnti, n["abs"] + n["k0"]))
            truth = {"format": g["format"], "al": g["al"], "cce": g["cce"]}
            if p:
                truth.update(h_id=p["h_id"], rv=p["rv"], prb=p["prb"], symb=p["symb"], mod=p["mod"], tbs_bits=p["tbs_bits"])
                ri = ri_at.get((rnti, p["abs"], p["h_id"]))
                if ri is not None:
                    truth["layers"] = ri
            got = {f: n[f] for f in fields}
            for f in fields:
                if f not in truth:
                    continue
                agree[f][1] += 1
                if got[f] == truth[f]:
                    agree[f][0] += 1
                else:
                    first_bad.setdefault(f, n["abs"])
                    if len(examples[f]) < a.examples:
                        examples[f].append((rel(n["abs"], origin), n["sfn"], n["slot"], got[f], truth[f]))
            joined_rows.append(dict(rnti=f"0x{rnti:04x}", sfn=g["sfn"], slot=g["slot"], t_rel=round(rel(g["abs"], origin), 4),
                                    detected=1, **{f"gnb_{f}": truth.get(f, "") for f in fields},
                                    **{f"nr_{f}": got[f] for f in fields}))
        for g in missed:
            joined_rows.append(dict(rnti=f"0x{rnti:04x}", sfn=g["sfn"], slot=g["slot"], t_rel=round(rel(g["abs"], origin), 4),
                                    detected=0, gnb_format=g["format"], gnb_al=g["al"], gnb_cce=g["cce"]))

        print("field agreement on detected DCIs (NR-Scope vs gNB):")
        for f in fields:
            ok, tot = agree[f]
            line = f"  {f:9s}: {ok:4d} / {tot:4d} = {fmt_pct(ok, tot)}"
            if f in first_bad:
                line += f"   first mismatch at +{rel(first_bad[f], origin):.3f} s"
            print(line)
            for ex in examples.get(f, []):
                print(f"      +{ex[0]:7.3f} s SFN {ex[1]}.{ex[2]}: NR-Scope {ex[3]}  gNB {ex[4]}")

        # --- timeline: RRC events next to detection, per second ------------
        print("gNB RRC events (moved onto NR-Scope's clock; ~ms accuracy, they carry no SFN):")
        for e in sorted((e for e in rrc if e["rnti"] == rnti), key=lambda e: e["t"]):
            t_ev = e["t"] + offset
            print(f"  {rel((t_ev - clk.ref) * slot_rate(), origin):+8.3f} s  {e['channel']} {e['msg']}")
        buckets = defaultdict(lambda: [0, 0, 0])
        for g in searched:
            buckets[int(rel(g["abs"], origin))][0] += 1
        for g, _ in matched:
            buckets[int(rel(g["abs"], origin))][1] += 1
        for g in after_expiry:
            buckets[int(rel(g["abs"], origin))][2] += 1
        if buckets:
            print("per second since the RNTI was known: gNB DL DCIs / detected  (x = sent after NR-Scope's expiry)")
            secs = range(0, max(buckets) + 1)
            print("  " + " ".join(
                f"{s}:{buckets[s][1]}/{buckets[s][0]}" + (f"x{buckets[s][2]}" if buckets[s][2] else "")
                for s in secs if any(buckets[s])))

        summary["rntis"][f"0x{rnti:04x}"] = dict(
            gnb_dl_dcis_searched=len(searched), detected=len(matched), missed=len(missed),
            after_expiry=len(after_expiry), false_positives=len(false_pos),
            miss_classes=dict(miss_classes),
            recall=len(matched) / len(searched) if searched else None,
            by_format={f"{k[0]}/ss{k[1]}": v for k, v in by_kind.items()},
            by_al_cce={f"al{k[0]}/cce{k[1]}": v for k, v in by_loc.items()},
            field_agreement={f: (agree[f][0] / agree[f][1] if agree[f][1] else None) for f in fields})

    if a.out_csv and joined_rows:
        keys = list(dict.fromkeys(k for r in joined_rows for k in r))
        with open(a.out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(sorted(joined_rows, key=lambda r: (r["rnti"], r["t_rel"])))
        print(f"\njoined rows written to {a.out_csv}")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"summary written to {a.json}")


if __name__ == "__main__":
    main()
