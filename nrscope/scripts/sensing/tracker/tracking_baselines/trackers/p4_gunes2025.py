"""P4: Gunes & Caliskan, "A New Plot Fusion Approach for Multi-static Multi-Target Tracking",
IEEE RadarConf 2025: plot fusion, then a GNN tracker.
Fusion (their ToA + azimuth variant): every combination of one detection per receiver (rx1's scan
with rx2's nearest scan, within max_dt) is localised by least squares on the bistatic ranges and
angles of both (their eq. 10; theirs fmincon, here scipy least_squares from the mean of the two
single-receiver positions); combinations taken greedily by lowest normalised residual, each
detection used once (deghosting, after Malanowski), and only if the residual passes a chi^2
test (assumed). The fused plots go to a textbook GNN tracker on positions (track.run_xy_tracker:
KF, Hungarian, 3 hits to confirm, 3 misses to delete; theirs is not described further).
One receiver: the plots are its own localised detections (no fusion). Their abstract names a
multi-Bernoulli filter, their body and results a GNN: we follow the body."""
import numpy as np
from scipy.optimize import least_squares

from . import common as C


def fuse(m1, g1, m2, g2, R1, R2):
    def res(p):
        x = np.r_[p, 0.0, 0.0]
        h1, h2 = C.track.bistatic_h(x, *g1)[0], C.track.bistatic_h(x, *g2)[0]
        r = [(m1["e"] - h1[0]) / np.sqrt(R1[0, 0]), (m2["e"] - h2[0]) / np.sqrt(R2[0, 0])]
        for m, h, R in ((m1, h1, R1), (m2, h2, R2)):
            if m["theta"] is not None:
                r.append(((m["theta"] - h[2] + np.pi) % (2 * np.pi) - np.pi) / np.sqrt(R[2, 2]))
        return np.array(r)
    p0 = (np.asarray(m1["xy"]) + np.asarray(m2["xy"])) / 2
    s = least_squares(res, p0)
    return s.x, float(s.fun @ s.fun), len(s.fun)


def run(frames, max_dt=0.08, gate_p=0.999):
    rxs = {f[5] for f in frames}
    plots = []
    if len(rxs) == 1:
        plots = [(t, [m["xy"] for m in meas if m["xy"] is not None]) for t, meas, rb, sc, geo, rx in frames]
    else:
        by = {r: [f for f in frames if f[5] == r] for r in rxs}
        r1, r2 = sorted(rxs)
        t2 = np.array([f[0] for f in by[r2]])
        for f in frames:
            t, meas, rb, sc, geo, rx = f
            if rx != r1:
                plots.append(None)   # rx2's own scans: fused into rx1's; reports the latest tracks
                continue
            j = int(np.argmin(np.abs(t2 - t)))
            fused = []
            if abs(t2[j] - t) <= max_dt:
                g = by[r2][j]
                A = [m for m in meas if m["xy"] is not None]
                B = [m for m in g[1] if m["xy"] is not None]
                Ra, Rb = C.Rz(rb, sc), C.Rz(g[2], g[3])
                cand = []
                for a, ma in enumerate(A):
                    for b, mb in enumerate(B):
                        if np.hypot(*(np.asarray(ma["xy"]) - mb["xy"])) > 20:
                            continue
                        p, c2, k = fuse(ma, geo, mb, g[4], Ra, Rb)
                        if c2 < C.gate(k - 2, gate_p):
                            cand.append((c2, a, b, p))
                ua, ub = set(), set()
                for c2, a, b, p in sorted(cand, key=lambda c: c[0]):
                    if a not in ua and b not in ub:
                        ua.add(a), ub.add(b)
                        fused.append(p)
            plots.append((t, fused))
    res = iter(C.track.run_xy_tracker([p for p in plots if p is not None]))
    out, last = [], []
    for p in plots:
        if p is not None:
            last = next(res)
        out.append(last)
    return out
