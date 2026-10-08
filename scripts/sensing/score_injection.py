#!/usr/bin/env python3
# SPDX-License-Identifier: LicenseRef-CSSL-1.0
"""Score the maps of a sensing_offline --inject run against the injected truth.

    score_injection.py <map2d.csv> [--truth <map2d.csv.truth.csv>] [--per-map]

sensing_offline writes the maps (<dump>) and, with --inject, the trajectory of every
injected target, one row per slot (<dump>.truth.csv). For each map this script takes
the truth over that map's window -- the window ends at the map's frame.slot and lasts
t_span_s -- and asks whether a detection sits on it:

  detections  the TDD detector's accepted markers (verdict 0) when the map has any
              markers, otherwise the AoA rows; ranges are taken as excess range over
              the direct path, the axis --inject's range= is given on
  hit         a detection within the range tolerance of the target's mean range over
              the window (half the range it covered in the window plus 1.5 bins) and
              within the speed tolerance of its mean speed (half the speed it covered
              plus 2 Doppler resolution cells, lambda / (2 t_span) each)
  other       an accepted detection that matches no injected target: real movers,
              clutter residue, or the copies the injection itself makes (a weak copy
              of every real echo, at its range + the target's and its speed + the
              target's). Those within NEAR_LOS_M of the direct path are counted apart,
              since that is where the residue lives.

Prints the detection rate per target and stream, the mean error of the hits, and the
count of other detections per map.
"""
import argparse
import bisect
import csv
import math
import sys
from collections import defaultdict

C = 299792458.0
NEAR_LOS_M = 3.0
SLOTS_PER_HYPERFRAME = 1024 * 20  # 30 kHz


def read_truth(path):
    """{target: sorted list of (slot_abs, active, range_m, speed_ms)}, and (frame, slot) -> [slot_abs]."""
    tr = defaultdict(list)
    where = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            sa = int(r["slot_abs"])
            tr[int(r["target"])].append((sa, int(r["active"]), float(r["range_m"]), float(r["speed_ms"])))
            if r["target"] == "0":
                where[(int(r["frame"]), int(r["slot"]))].append(sa)
    for v in tr.values():
        v.sort()
    return tr, where


def read_maps(path):
    """One dict per map row of nr_ue_sensing_dump_map()'s format."""
    maps = []
    with open(path) as f:
        for row in csv.reader(l for l in f if not l.startswith("#")):
            nb, nf, nm, na = int(row[12]), int(row[13]), int(row[14]), int(row[15])
            mpb, los = float(row[9]), float(row[16])
            off = 21 + nb * nf
            markers = [tuple(float(x) for x in row[off + 4 * i: off + 4 * i + 4]) for i in range(nm)]
            off += 4 * nm
            aoa = [tuple(float(x) for x in row[off + 5 * i: off + 5 * i + 5]) for i in range(na)]
            if markers:
                # marker range is from the window start; excess range is from the direct path
                det = [(m[0] - los * mpb, m[1], m[2]) for m in markers if int(m[3]) == 0]
            else:
                det = [(a[0], a[1], a[3]) for a in aoa]
            maps.append(dict(frame=int(row[0]), slot=int(row[1]), aarx=int(row[2]), layer=int(row[4]),
                             t_span=float(row[8]), mpb=mpb, carrier=float(row[11]), det=det))
    return maps


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("maps")
    ap.add_argument("--truth", help="default <maps>.truth.csv")
    ap.add_argument("--per-map", action="store_true", help="print one line per map")
    args = ap.parse_args()

    truth, where = read_truth(args.truth or args.maps + ".truth.csv")
    maps = read_maps(args.maps)
    if not maps:
        sys.exit("no maps in " + args.maps)
    slot_s = 0.5e-3  # 30 kHz

    stats = defaultdict(lambda: dict(n=0, hit=0, dr=0.0, dv=0.0))
    others = defaultdict(lambda: dict(n=0, far=0, near=0))
    prev_sa = None
    for m in maps:
        # the map's slot_abs: the candidate nearest the previous map, so runs past the
        # 10.24 s SFN wrap stay unambiguous (maps are dumped in about time order)
        cand = where.get((m["frame"], m["slot"]))
        if not cand:
            continue
        sa_end = cand[0] if prev_sa is None else min(cand, key=lambda c: abs(c - prev_sa))
        prev_sa = sa_end
        sa_start = sa_end - int(round(m["t_span"] / slot_s))
        lam = C / m["carrier"]
        dv_cell = lam / (2.0 * m["t_span"])
        stream = ("rxavg" if m["aarx"] < 0 else "rx%d" % m["aarx"]) + " L%d" % m["layer"]
        used = set()
        line = []
        for tgt, rows in truth.items():
            i0 = bisect.bisect_left(rows, (sa_start,))
            i1 = bisect.bisect_right(rows, (sa_end, 2))
            win = [r for r in rows[i0:i1] if r[1]]
            # a target present for under half the window is not scored on this map
            if len(win) < 0.5 * max(1, i1 - i0) or not win:
                continue
            rs = [r[2] for r in win]
            vs = [r[3] for r in win]
            r_mean, v_mean = sum(rs) / len(rs), sum(vs) / len(vs)
            tol_r = (max(rs) - min(rs)) / 2 + 1.5 * m["mpb"]
            tol_v = (max(vs) - min(vs)) / 2 + 2.0 * dv_cell
            best = None
            for j, d in enumerate(m["det"]):
                er, ev = d[0] - r_mean, d[1] - v_mean
                if abs(er) <= tol_r and abs(ev) <= tol_v:
                    score = (er / tol_r) ** 2 + (ev / tol_v) ** 2
                    if best is None or score < best[0]:
                        best = (score, j, er, ev)
            st = stats[(tgt, stream)]
            st["n"] += 1
            if best:
                st["hit"] += 1
                st["dr"] += abs(best[2])
                st["dv"] += abs(best[3])
                used.add(best[1])
                line.append("T%d hit (%+.1f m, %+.2f m/s)" % (tgt, best[2], best[3]))
            else:
                line.append("T%d MISS (at %.1f m, %+.2f m/s)" % (tgt, r_mean, v_mean))
        ot = others[stream]
        ot["n"] += 1
        for j, d in enumerate(m["det"]):
            if j in used:
                continue
            ot["near" if abs(d[0]) < NEAR_LOS_M else "far"] += 1
        if args.per_map:
            rest = ["%.1f m %+.2f m/s %.0f dB" % d for j, d in enumerate(m["det"]) if j not in used]
            print("%4d.%-2d %-9s %s | other: %s" % (m["frame"], m["slot"], stream, "; ".join(line), ", ".join(rest)))

    if args.per_map:
        print()
    print("target  stream      maps  detected  mean |range err|  mean |speed err|")
    for (tgt, stream), st in sorted(stats.items()):
        h = st["hit"]
        print("T%-6d %-10s %5d  %4d %4.0f%%  %13s  %16s" % (
            tgt, stream, st["n"], h, 100.0 * h / st["n"] if st["n"] else 0,
            "%.2f m" % (st["dr"] / h) if h else "-", "%.3f m/s" % (st["dv"] / h) if h else "-"))
    print("\nother accepted detections per map (not an injected target)")
    for stream, ot in sorted(others.items()):
        print("  %-10s far from the direct path %.2f, within %.0f m of it %.2f  (%d maps)" % (
            stream, ot["far"] / ot["n"], NEAR_LOS_M, ot["near"] / ot["n"], ot["n"]))


if __name__ == "__main__":
    main()
