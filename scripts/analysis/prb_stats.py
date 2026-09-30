#!/usr/bin/env python3
"""How much of the carrier the gNB gave each downlink PDSCH, from an NR-Scope run.

Every row of DCIs/dci_<run>_pci<N>.csv is one DCI NR-Scope decoded. For a
downlink DCI, nof_prb is the number of PRBs of the PDSCH it schedules and prbs
the range (e.g. 0-75). NR-Scope gets them from the DCI's frequency-domain
resource assignment (type 1: start and length coded as one RIV, TS 38.214
5.1.2.2.2), decoded with srsRAN. The gNB benchmark (scripts/benchmark) checks
these columns against the gNB's own log.

usage: prb_stats.py [DCI_CSV] [--carrier-prb 273] [--min-rate 1000]
  DCI_CSV       default: the newest file in DCIs/
  --min-rate    only count seconds with at least this many DL DCIs, i.e. while
                traffic is running (0 counts everything)
"""
import argparse
import csv
import glob
import os
import statistics
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dci_csv", nargs="?")
    ap.add_argument("--carrier-prb", type=int, default=273, help="carrier bandwidth in PRBs (default 273 = 100 MHz at 30 kHz)")
    ap.add_argument("--min-rate", type=int, default=1000, help="count only seconds with at least this many DL DCIs (default 1000)")
    a = ap.parse_args()

    path = a.dci_csv or max(glob.glob(os.path.join(ROOT, "DCIs", "dci_*.csv")), key=os.path.getmtime)
    rows = [r for r in csv.DictReader(open(path)) if r["direction"] == "DL" and r["nof_prb"]]
    if not rows:
        raise SystemExit(f"{path}: no downlink grants")

    # Keep the seconds where traffic was running
    t0 = float(rows[0]["timestamp"])
    per_second = Counter(int(float(r["timestamp"]) - t0) for r in rows)
    busy = [r for r in rows if per_second[int(float(r["timestamp"]) - t0)] >= a.min_rate]
    if not busy:
        raise SystemExit(f"{path}: no second with {a.min_rate}+ DL DCIs; try --min-rate 0")

    prb = [int(r["nof_prb"]) for r in busy]
    q = statistics.quantiles(prb, n=20) if len(prb) > 1 else [prb[0]] * 19
    full = a.carrier_prb
    seconds = sum(1 for s, n in per_second.items() if n >= a.min_rate)

    print(f"{os.path.basename(path)}: {len(busy)} DL grants over {seconds} s with >= {a.min_rate} DL DCIs/s")
    print(f"  PRBs per PDSCH: median {statistics.median(prb):.0f} of {full} ({100 * statistics.median(prb) / full:.0f}% of the carrier),"
          f" mean {statistics.mean(prb):.0f}, 5th pct {q[0]:.0f}, 95th pct {q[18]:.0f}, max {max(prb)}")
    print(f"  full band ({full} PRBs): {100 * sum(p == full for p in prb) / len(prb):.1f}% of grants;"
          f" >= 90% of it: {100 * sum(p >= 0.9 * full for p in prb) / len(prb):.1f}%")
    print(f"  most common sizes: " + ", ".join(f"{p} PRB ({100 * n / len(prb):.0f}%)" for p, n in Counter(prb).most_common(5)))
    starts = Counter(r["prbs"].split("-")[0].split(",")[0] for r in busy if r["prbs"])
    print(f"  allocation starts at PRB: " + ", ".join(f"{s} ({100 * n / len(busy):.0f}%)" for s, n in starts.most_common(3)))
    print(f"  layers: " + ", ".join(f"{k} ({100 * n / len(busy):.0f}%)" for k, n in Counter(r["nof_layers"] for r in busy).most_common()) +
          f"; MCS median {statistics.median(int(r['dci_mcs']) for r in busy):.0f}; PDSCH symbols: " +
          ", ".join(f"{k} ({100 * n / len(busy):.0f}%)" for k, n in Counter(r["time_length"] for r in busy).most_common(2)))
    # DM-RS pattern: configuration type (comb), CDM groups without data, symbols
    def pct(counter, fmt):
        return ", ".join(f"{fmt(k)} ({100 * v / len(busy):.1f}%)" for k, v in counter.most_common())
    comb = {"1": "type 1 = comb-2 (every other subcarrier, 6 RE/PRB)",
            "2": "type 2 = 2 adjacent subcarriers per CDM group (4 RE/PRB)"}
    print("  DM-RS type: " + pct(Counter(r["dmrs_type"] for r in busy), lambda k: comb.get(k, f"type {k}")))
    cdm = {"1": "1 (other comb carries data)", "2": "2 (other comb left empty)", "3": "3 (type 2 only)"}
    print("  CDM groups without data: " + pct(Counter(r["nof_dmrs_cdm_groups"] for r in busy), lambda k: cdm.get(k, k)))
    if busy[0].get("dmrs_symbols") is not None:
        print("  DM-RS symbols: " + pct(Counter(r["dmrs_symbols"] for r in busy if r["dmrs_symbols"]), lambda k: "{" + k + "}"))
    print("  DM-RS additional positions: " + pct(Counter(r["dmrs_add_pos"] for r in busy), lambda k: f"pos{k}") +
          "; length: " + pct(Counter(r["dmrs_len"] for r in busy), str))

    # What that means for sensing: the DM-RS of a grant only sees the grant's bandwidth
    scs_hz = 30e3
    bw = statistics.median(prb) * 12 * scs_hz
    print(f"  median grant bandwidth {bw / 1e6:.1f} MHz -> range resolution c/(2B) about {3e8 / (2 * bw):.1f} m"
          f" (full carrier: {3e8 / (2 * full * 12 * scs_hz):.1f} m)")


if __name__ == "__main__":
    main()
