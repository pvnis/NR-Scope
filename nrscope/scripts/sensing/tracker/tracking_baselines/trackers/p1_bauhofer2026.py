"""P1: the tracker of Bauhofer et al., "Experimental Demonstration of Multi-Target Tracking in
Integrated Sensing and Communication", EuCNC/6G Summit 2026: a GM-PHD filter (Vo & Ma 2006) in the
range / range-rate domain of one sensor, no angle, as in their Table II:
state [e, e_dot], CV with Q = diag(0.01, 0.01) per second, H: (e, v) with v = -e_dot / 2 (our
bistatic speed), p_S = 0.93 and p_D = 0.8 per scan, no spawning, merging radius 0.4 (Euclidean),
at most 40 components, extraction of weights > 0.5.
Adaptations (theirs reaches 15-60 m only): the birth intensity, their 1e-6 N([20, 0.5],
diag(100, 100)), tiled every 20 m over our 0-100 m; clutter density as for every tracker here.
Their tracker has no angle, hence no position: for scoring, each estimate takes the angle of the
nearest measurement of the scan (within 3 sigma), and is localised on its ellipse; none, no
position. Single-sensor: for the fused setup each receiver keeps its own filter (no fusion)."""
import numpy as np

from . import common as C

H = np.array([[1.0, 0.0], [0.0, -0.5]])


def run(frames, q=0.01, p_s=0.93, p_d=0.8, clutter=2.0, w_birth=1e-6, merge_r=0.4, jmax=40):
    filt = {}     # receiver -> components [w, x, P, label]
    out, nid = [], [1]
    births = [np.array([e, 0.5]) for e in range(10, 100, 20)]
    for t, meas, rb, sc, geo, rx in frames:
        f = filt.setdefault(rx, dict(comps=[], t=None))
        dt = 0.0 if f["t"] is None else t - f["t"]
        f["t"] = t
        F = np.array([[1.0, dt], [0.0, 1.0]])
        Q = np.diag([q, q]) * dt
        comps = [[w * p_s, F @ x, F @ P @ F.T + Q, l] for w, x, P, l in f["comps"]]
        for b in births:
            comps.append([w_birth, b.copy(), np.diag([100.0, 100.0]), nid[0]])
            nid[0] += 1
        R = C.Rz(rb, sc)[:2, :2]
        new = [[w * (1 - p_d), x, P, l] for w, x, P, l in comps]
        for m in meas:
            z = np.array([m["e"], m["v"]])
            upd = []
            for w, x, P, l in comps:
                S = H @ P @ H.T + R
                nu = z - H @ x
                K = P @ H.T @ np.linalg.inv(S)
                upd.append([w * p_d * C.gauss(nu, S), x + K @ nu, (np.eye(2) - K @ H) @ P, l])
            tot = C.kappa(clutter, 2) + sum(u[0] for u in upd)
            new += [[u[0] / tot, *u[1:]] for u in upd]
        new = sorted([c for c in new if c[0] > 1e-12], key=lambda c: -c[0])
        merged = []
        while new:
            c0 = new[0]
            near = [c for c in new if np.linalg.norm(c[1] - c0[1]) <= merge_r]
            ids = {id(c) for c in near}
            new = [c for c in new if id(c) not in ids]
            x, P = C.moment_match([c[0] for c in near], [c[1] for c in near], [c[2] for c in near])
            merged.append([sum(c[0] for c in near), x, P, c0[3]])
        f["comps"] = merged[:jmax]
        est, seen = [], set()
        Rv = C.Rz(rb, sc)
        T, Rr, L, bore = geo
        for c in f["comps"]:
            w, x, P, l = c
            if w <= 0.5:
                continue
            if l in seen:   # a label shared by two targets: a new one
                c[3] = l = nid[0]
                nid[0] += 1
            seen.add(l)
            e, v = x[0], -x[1] / 2
            best, bd = None, 9.0
            for m in meas:
                if m["theta"] is None:
                    continue
                d = (m["e"] - e) ** 2 / Rv[0, 0] + (m["v"] - v) ** 2 / Rv[1, 1]
                if d < bd:
                    best, bd = m, d
            if best is None:
                continue
            p = C.track.localize(e, np.rad2deg(best["theta"]), T, Rr, np.rad2deg(bore))
            if p is not None:
                est.append((l, p[0], p[1], 0.0, 0.0))
        out.append(est)
    return out


SINGLE_SENSOR = True
