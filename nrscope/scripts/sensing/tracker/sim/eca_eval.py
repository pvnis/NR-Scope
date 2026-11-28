#!/usr/bin/env python3
"""Compare clutter settings on the synthetic-target runs.

    eca_eval.py LOGDIR TARGET_EXCESS_M TARGET_SPEED_MS [CLUTTER_EXCESS_M]

Per run: maps where NR-Scope's detector reports the target (within 2 range bins, 0.7 m/s), false
targets per map, the map's own peak at the target, and with CLUTTER_EXCESS_M the residue the
stepping reflector leaves in its row away from zero speed (median, dB over the map's median)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze  # noqa: E402

d, te, tv = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
ce = float(sys.argv[4]) if len(sys.argv) > 4 else None
maps = analyze.read_maps(os.path.join(d, "map2d.csv"))
found = mpk = 0
false, resid = [], []
for m in maps:
    det = [(r[0] - m["bin_los"] * m["m_per_bin"], r[1]) for r in m["markers"] if int(r[3]) == 0]
    ok = [x for x in det if abs(x[0] - te) <= 2 * m["m_per_bin"] and abs(x[1] - tv) <= 0.7]
    found += bool(ok)
    false.append(len(det) - (1 if ok else 0))
    b = int(round(m["bin_los"] + te / m["m_per_bin"]))
    k = int(np.argmax(m["power"][b]))
    mpk += abs(m["speeds"][k] - tv) <= 0.7
    if ce is not None:
        bc = int(round(m["bin_los"] + ce / m["m_per_bin"]))
        away = np.abs(m["speeds"]) > 1.0
        resid.append(10 * np.log10(np.median(m["power"][bc, away]) / np.median(m["power"])))
line = (f"{os.path.basename(d.rstrip('/')):14s} {len(maps):3d} maps: detector found the target {found}/{len(maps)}, "
        f"false targets/map {np.mean(false):.2f}, map peak at the target {mpk}/{len(maps)}")
if resid:
    line += f", reflector residue {np.median(resid):5.1f} dB"
print(line)
