#!/usr/bin/env python3
"""Spread of a target over the Doppler axis of its own range bin.

    streak.py LOGDIR [EXCESS_M SPEED_MS]     (without a fixed target: the ground-truth car echo)

For each map: the target's peak over the median of its range row (how far it stands above its
own smear), the Doppler cells within 10 dB of it, and its row's median over the map's median."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze  # noqa: E402

d = sys.argv[1]
maps = analyze.read_maps(os.path.join(d, "map2d.csv"))
fixed = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else None
if not fixed:
    truth = analyze.read_truth(os.path.join(d, "ground_truth.csv"))
    times = np.array(sorted(truth))
pk, wide, rowup, hit, cnt = [], [], [], 0, 0
for m in maps:
    if fixed:
        e, v = fixed
    else:
        gt = analyze.truth_at(truth, times, m["t"], m["lam"])
        if "excess" not in gt or abs(gt["speed"]) < 1:
            continue
        e, v = gt["excess"], gt["speed"]
    b = int(round(m["bin_los"] + e / m["m_per_bin"]))
    if b >= m["n_bins"]:
        continue
    row = m["power"][b]
    db = 10 * np.log10(row / np.median(row))
    k = int(np.argmax(db))
    pk.append(db[k]); wide.append(np.sum(db > db[k] - 10))
    rowup.append(10 * np.log10(np.median(row) / np.median(m["power"])))
    cnt += 1
    hit += abs(m["speeds"][k] - v) <= 2 * abs(m["speed_step"]) * 2
print(f"{os.path.basename(d.rstrip('/')):18s} {cnt:3d} maps: peak over its row {np.median(pk):5.1f} dB, "
      f"cells within 10 dB {np.median(wide):4.0f}, row over map median {np.median(rowup):5.1f} dB, "
      f"row peak at the right speed {hit}/{cnt}")
