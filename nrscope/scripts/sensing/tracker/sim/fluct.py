#!/usr/bin/env python3
"""Clutter / amplitude-fluctuation analysis of NR-Scope's raw snapshots (map2d.csv.snap.csv).

    fluct.py LOGDIR

For one window of slow-time delay profiles h[m, b] (snapshot m, delay bin b):
  static        mean over snapshots, per bin
  residual A    h - mean                          (what mean clutter removal leaves)
  residual B    h - g_m * mean, g_m = <h_m, mean> / |mean|^2
                                                  (a common complex gain per snapshot removed)
and the car echo the ground truth expects at its bin, at the window's time.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze  # noqa: E402

d = sys.argv[1]
rows = [l.rstrip().split(",") for l in open(os.path.join(d, "map2d.csv.snap.csv")) if not l.startswith("#")]
t = np.array([int(r[1]) for r in rows]) / float(rows[0][9])
kfirst = np.array([int(r[4]) for r in rows]); npil = np.array([int(r[5]) for r in rows])
nb = int(rows[0][7]); fc = float(rows[0][10]); lam = 299792458 / fc
mpb = 299792458 / (int(rows[0][6]) * int(rows[0][3]) * float(rows[0][8]))
H = np.array([[float(r[11 + 2 * b]) + 1j * float(r[12 + 2 * b]) for b in range(nb)] for r in rows])
print(f"{len(H)} snapshots, {t[0]:.3f}-{t[-1]:.3f} s, {nb} bins of {mpb:.2f} m; grants: "
      + ", ".join(f"k_first {k}/{n} pilots x{c}" for (k, n), c in zip(*np.unique(np.c_[kfirst, npil], axis=0, return_counts=True))))

mean = H.mean(0)
los = int(np.argmax(abs(mean)))
g = (H @ mean.conj()) / np.vdot(mean, mean)
resA = H - mean
resB = H - np.outer(g, mean)
P = lambda x: 10 * np.log10(np.mean(abs(x) ** 2, axis=0) + 1e-30)
ref = P(H)[los]
print(f"direct path at bin {los}: amplitude std {abs(H[:, los]).std() / abs(H[:, los]).mean():.3f}, "
      f"phase std {np.std(np.angle(H[:, los] * mean[los].conj())):.3f} rad")
print(f"common gain g_m: |g| std {abs(g).std():.3f}, phase std {np.std(np.angle(g)):.3f} rad")

truth = analyze.read_truth(os.path.join(d, "ground_truth.csv"))
times = np.array(sorted(truth))
gt = analyze.truth_at(truth, times, t.mean(), lam)
bcar = los + gt["excess"] / mpb if "excess" in gt else None
print(f"truth at {t.mean():.2f} s: car echo {gt.get('excess', float('nan')):.1f} m excess (bin {bcar:.1f}), "
      f"{gt.get('gain_db', float('nan')):.1f} dB re direct, {gt.get('speed', float('nan')):.2f} m/s")
print("\n bin  excess   static  resid(mean)  resid(gain)   [dB re direct path power]")
for b in range(nb):
    mark = " <- direct" if b == los else (" <- car" if bcar is not None and abs(b - bcar) < 0.5 else "")
    if b < los + 6 or mark or b % 5 == 0:
        print(f" {b:3d} {(b - los) * mpb:7.1f}  {P(mean[None])[b] - ref:7.1f}  {P(resA)[b] - ref:9.1f}  {P(resB)[b] - ref:11.1f}{mark}")

# Doppler: coherent integration gain at the car's Doppler vs the residual floor
fd = gt.get("speed", 0) * 2 / lam
if bcar is not None:
    bc = int(round(bcar))
    for name, R in (("mean removed", resA), ("gain removed", resB)):
        spec = lambda f: abs(np.sum(R[:, bc] * np.exp(-2j * np.pi * f * (t - t.mean())))) ** 2 / len(t)
        fgrid = np.linspace(-500, 500, 1001)
        s = np.array([spec(f) for f in fgrid])
        print(f"{name}: car bin {bc}, power at the car's Doppler {fd:.0f} Hz {10 * np.log10(spec(fd) + 1e-30) - ref:.1f} dB, "
              f"median over Doppler {10 * np.log10(np.median(s) + 1e-30) - ref:.1f} dB, expected echo after integration "
              f"{gt['gain_db'] + 10 * np.log10(len(t)):.1f} dB (all re direct path per snapshot)")
