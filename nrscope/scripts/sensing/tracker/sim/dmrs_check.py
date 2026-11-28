#!/usr/bin/env python3
"""Step 2: recompute NR-Scope's DM-RS snapshots from the exact IQ it received (bridge --tap).

    dmrs_check.py LOGDIR [N_ID]

For each snapshot in map2d.csv.snap.csv that falls inside tap.npz: demodulate that OFDM symbol,
least-squares estimate on the DM-RS (type 1, port 1000, nSCID 0, TS 38.211 7.4.1.1), IDFT to a
delay profile like NR-Scope's, and compare with NR-Scope's own profile.
"""
import os
import sys

import numpy as np

d = sys.argv[1]
NID = int(sys.argv[2]) if len(sys.argv) > 2 else 1
N, SLOT, CP0, CP = 4096, 61440, 352, 288
NSC = 3276  # 273 PRB


def gold(cinit, n):
    x1 = np.zeros(n + 1600 + 31, np.uint8); x2 = np.zeros_like(x1)
    x1[0] = 1
    for i in range(31):
        x2[i] = (cinit >> i) & 1
    for i in range(n + 1600):
        x1[i + 31] = x1[i + 3] ^ x1[i]
        x2[i + 31] = x2[i + 3] ^ x2[i + 2] ^ x2[i + 1] ^ x2[i]
    return (x1[1600:1600 + n] ^ x2[1600:1600 + n]).astype(float)


def dmrs(slot, l, n_pilots):
    cinit = ((1 << 17) * (14 * slot + l + 1) * (2 * NID + 1) + 2 * NID) % (1 << 31)
    c = gold(cinit, 2 * n_pilots)
    return ((1 - 2 * c[0::2]) + 1j * (1 - 2 * c[1::2])) / np.sqrt(2)


def sym_start(l):  # first sample of symbol l's useful part, from the slot start
    return CP0 + l * (N + CP) + (0 if l == 0 else (CP0 - CP) * 0)


tap = np.load(os.path.join(d, "tap.npz"))
x, n0 = tap["x"], int(tap["n0"][0])
rows = [r.rstrip().split(",") for r in open(os.path.join(d, "map2d.csv.snap.csv")) if not r.startswith("#")]
out = []
for r in rows:
    ts, kf, npil, idft, nb = int(r[1]), int(r[4]), int(r[5]), int(r[6]), int(r[7])
    h_nrs = np.array([float(r[11 + 2 * b]) + 1j * float(r[12 + 2 * b]) for b in range(nb)])
    slot_abs, within = divmod(ts, SLOT)
    starts = np.array([CP0 + l * (N + CP) for l in range(14)])  # useful-part start of each symbol
    l = int(np.argmin(abs(starts - within))) if np.min(abs(starts - within)) < CP + 1 else None
    a = slot_abs * SLOT + (starts[l] if l is not None else within) - n0
    if l is None or a < 0 or a + N > len(x):
        continue
    X = np.fft.fftshift(np.fft.fft(x[a:a + N])) / np.sqrt(N)
    grid = X[N // 2 - NSC // 2:N // 2 + NSC // 2]          # carrier subcarriers 0..3275
    k = np.arange(kf, kf + 2 * npil, 2)                     # even comb (delta 0)
    m = k // 2                                              # sequence index from point A
    seq = dmrs(slot_abs % 20, l, m.max() + 1)[m]
    H = grid[k] * np.conj(seq) * 2                          # / (1/sqrt2)^2 amplitude
    prof = np.fft.ifft(H, idft)[:nb]
    if not np.any(prof) or not np.any(h_nrs):
        continue
    c = abs(np.vdot(prof, h_nrs)) / np.linalg.norm(prof) / np.linalg.norm(h_nrs)
    pk = abs(prof) ** 2
    out.append((within - starts[l], l, c, 10 * np.log10(pk[0] / pk[1:].sum())))
o = np.array(out)
print(f"{len(o)} snapshots compared (symbol offset of NR-Scope's t_sample vs symbol start: {sorted(set(o[:, 0].astype(int)))})")
print(f"symbols: {sorted(set(o[:, 1].astype(int)))}")
print(f"my estimate: direct-path bin over all other bins {np.median(o[:, 3]):.1f} dB (a clean flat channel: large)")
print(f"correlation my profile vs NR-Scope's: median {np.median(o[:, 2]):.3f}, 10% {np.percentile(o[:, 2], 10):.3f}, 90% {np.percentile(o[:, 2], 90):.3f}")
