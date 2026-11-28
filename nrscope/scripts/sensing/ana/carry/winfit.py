#!/usr/bin/env python3
"""Sliding-window trajectory fit (offline test): in each window, RANSAC + robust least squares of a constant-velocity
trajectory to all detections (bistatic excess range, speed, angle) of the TRS and DM-RS streams of one run; the fitted
position at the window centre is scored against the lidar truth, with a time-shifted chance level.
Run from paper_prec (its python): winfit.py RUN W0 W1 [win_s step_s]"""
import csv, json, sys, numpy as np
sys.path[:0] = ["sim", "tracking_baselines"]
import track as TK
from scipy.optimize import least_squares
RUN, W0, W1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
WIN = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0; STEP = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
KMAX = int(sys.argv[6]) if len(sys.argv) > 6 else 5
Z = 1.0; rng = np.random.default_rng(0)
st = json.load(open(RUN + "/settings.json")); t0 = float(st["t0_unix"])
G = np.array(st["gnb"], float); R = np.array(st["rx"], float); bore = float(st["boresight_deg"]); L0 = np.linalg.norm(G - R)
def dets(mf, sig_a):
    unix = [float(q["unix_s"]) for q in csv.DictReader(open(RUN + "/" + mf + ".time.csv"))]; out = []
    for l, u in zip([l for l in open(RUN + "/" + mf) if not l.startswith("#")], unix):
        v = l.strip().split(","); nb, nf, nm, na = int(v[12]), int(v[13]), int(v[14]), int(v[15]); t = u - float(v[8]) / 2
        if not W0 - WIN <= t <= W1 + WIN: continue
        off = 21 + nb * nf; A = np.array(v[off + nm * 4: off + nm * 4 + na * 5], float).reshape(-1, 5)
        for a in A:
            if abs(a[1]) >= 0.3 and a[0] >= 4: out.append((t, a[0], a[1], a[2], sig_a))
    return out
D = np.array(dets("map2d.csv", 6.0) + dets("map2d_rx2.csv", 10.0))   # columns: t, excess, speed, angle, sigma_angle
print("detections kept (|v| >= 0.3 m/s, range >= 4 m): %d" % len(D))
def pred(p, t, tc):   # p = x0, y0, vx, vy at the window centre tc
    x = np.c_[p[0] + p[2] * (t - tc), p[1] + p[3] * (t - tc), np.full(len(t), Z)]; V = np.c_[np.full(len(t), p[2]), np.full(len(t), p[3]), np.zeros(len(t))]
    dg, dr = x - G, x - R; ng, nr = np.linalg.norm(dg, axis=1), np.linalg.norm(dr, axis=1)
    e = ng + nr - L0; s = -0.5 * np.sum(V * (dg / ng[:, None] + dr / nr[:, None]), 1)
    a = (np.degrees(np.arctan2(x[:, 1] - R[1], x[:, 0] - R[0])) - bore + 180) % 360 - 180
    return e, s, a
def resid(p, d, tc):
    e, s, a = pred(p, d[:, 0], tc)
    return np.r_[(d[:, 1] - e) / 1.0, (d[:, 2] - s) / 0.2, ((d[:, 3] - a + 180) % 360 - 180) / d[:, 4]]
def fit_one(d, pos, ok, tc):
    best, bn = None, 0
    for _ in range(400):
        i, j = rng.choice(ok, 2, replace=False)
        if abs(d[i, 0] - d[j, 0]) < 0.2: continue
        v = (pos[j] - pos[i]) / (d[j, 0] - d[i, 0])
        if np.hypot(*v) > 3.0: continue   # people: < 3 m/s
        p = np.r_[pos[i] + v * (tc - d[i, 0]), v]
        e, s, a = pred(p, d[:, 0], tc)
        inl = (np.abs(d[:, 1] - e) <= 3) & (np.abs(d[:, 2] - s) <= 0.5) & (np.abs((d[:, 3] - a + 180) % 360 - 180) <= 25)
        if inl.sum() > bn: bn, best = int(inl.sum()), (p, inl)
    if best is None or bn < 4: return None, None, bn
    r = least_squares(resid, best[0], args=(d[best[1]], tc), loss="soft_l1")
    return r.x, best[1], bn
def fit(tc, K=KMAX):   # sequential RANSAC: up to K trajectories, each with its own detections
    d = D[np.abs(D[:, 0] - tc) <= WIN / 2]; out = []
    if len(d) < 4: return out
    pos = []
    for row in d:
        p = TK.localize(row[1], row[3], G, R, bore); pos.append(None if p is None else np.asarray(p)[:2])
    pos = np.array(pos, dtype=object); keep = np.ones(len(d), bool)
    for _ in range(K):
        ok = [i for i in np.flatnonzero(keep) if pos[i] is not None]
        if len(ok) < 4: break
        p, inl, n = fit_one(d, pos, ok, tc)
        if p is None: break
        idx = np.flatnonzero(keep); out.append((p, n)); keep[inl & keep] = False
    return out
# truth (receiver-1 file: the lidar walker position)
gt = [x for x in csv.DictReader(open(RUN + "/ground_truth.csv")) if x["via_car"] == "1" and x["car_x"] != "nan"]
tt = np.array([float(x["t_s"]) for x in gt]) + t0; P = np.array([[float(x["car_x"]), float(x["car_y"])] for x in gt])
res = []
for tc in np.arange(W0 + WIN / 2, W1 - WIN / 2 + 1e-9, STEP):
    res.append((tc, fit(tc)))
def score(shift, rad):
    hit = tot = 0
    for tc, F in res:
        ts = W0 + ((tc + shift - W0) % (W1 - W0)); k = np.argmin(abs(tt - ts))
        if abs(tt[k] - ts) > 0.5: continue
        tot += 1; hit += any(np.hypot(*(p[:2] - P[k])) <= rad for p, n in F)
    return 100 * hit / max(tot, 1), tot
for rad in (3, 6):
    tr, n = score(0, rad); ch = [score(s, rad)[0] for s in (12, 24, 36, 48, 60, 72)]
    print("window %.1f s: fitted position within %d m of the walker in %4.1f %% of %d windows (chance %4.1f +- %3.1f %%)" % (WIN, rad, tr, n, np.mean(ch), np.std(ch)))
err = [min(np.hypot(*(p[:2] - P[np.argmin(abs(tt - tc))])) for p, n in F) for tc, F in res if F]
print("windows with trajectories %d of %d (median %d per window); nearest trajectory to the walker: median %.1f m (p25 %.1f, p75 %.1f)" % (len(err), len(res), int(np.median([len(F) for _, F in res if F])), np.median(err), np.percentile(err, 25), np.percentile(err, 75)))
