"""P2: the tracking core of BSense (Meng et al., "Needle in a Haystack: Tracking UAVs from Massive
Noise in Real-World 5G-A Base Station Data", MobiSys 2026), without its learned parts.
Object level (their Sec. 5): every localised object (our measurement group = their DBSCAN cluster)
gets two confidences, C <- 0.9 (C_prev + I):
  spatial: I = 1 when Hungarian matching to the previous scan's objects (gate g_s) finds it;
  velocity: I = 1 when the displacement velocity's predicted bistatic speed agrees with the
  measured one (|dv| < 2 m/s, ratio < 2, same sign; their radial check made bistatic).
Objects with C_s < 1.7 and C_v < 1.1 are dropped. Trajectory level: IMM (CV + coordinated turn;
theirs CV + CTRV with UKFs, here EKFs on the same bistatic model), Hungarian association
(chi^2 gate), confirmed after L = 6 consecutive associations, deleted after max_miss misses in
the field of view (their deletion rule is not given).
Not reproduced: the per-cube noise fingerprint (needs their 5 signal metrics) and TrajFormer
(needs labelled training data). g_s: theirs 15 m at 1 km range; a parameter here."""
import numpy as np
from scipy.optimize import linear_sum_assignment

from . import common as C

PI = np.array([[0.95, 0.05], [0.05, 0.95]])   # IMM model switching (assumed: not given)


def ct_predict(x, P, dt, q, q_w, turn):
    """5-state [px, py, vx, vy, w]; turn False: CV (w held at 0)."""
    xx = x.copy()
    if not turn:
        xx[4] = 0.0
    f = lambda s: C.track.ct_step(s, dt)
    F = np.empty((5, 5))
    for k in range(5):
        d = np.zeros(5)
        d[k] = 1e-5
        F[:, k] = (f(xx + d) - f(xx - d)) / 2e-5
    if not turn:
        F[:, 4] = 0.0
        F[4, 4] = 0.0
    G = np.zeros((5, 5))
    _, Q4 = C.predict(np.zeros(4), np.zeros((4, 4)), dt, q)
    G[:4, :4] = Q4
    G[4, 4] = q_w * dt if turn else 1e-6
    return f(xx), F @ P @ F.T + G


def run(frames, q=0.35, q_w=0.1, g_s=6.0, L=6, max_miss=3, gate_p=0.999):
    prev = {}      # receiver -> (t, [(xy, Cs, Cv)])
    tracks, out, nid = [], [], 1
    for t, meas, rb, sc, geo, rx in frames:
        R = C.Rz(rb, sc)
        # object level: confidences against this receiver's previous scan
        objs = [m for m in meas if m["xy"] is not None]
        pt, po = prev.get(rx, (None, []))
        cs, cv = np.zeros(len(objs)), np.zeros(len(objs))
        if objs and po:
            d = np.array([[np.hypot(*(o["xy"] - p[0])) for p in po] for o in objs])
            for i, j in zip(*linear_sum_assignment(d)):
                if d[i, j] <= g_s:
                    cs[i] = 0.9 * (po[j][1] + 1)
                    vel = (objs[i]["xy"] - po[j][0]) / max(t - pt, 1e-3)
                    _, H, _ = C.track.bistatic_h(np.r_[objs[i]["xy"], 0, 0], *geo)
                    vp, vm = float(H[1, 2:] @ vel), objs[i]["v"]
                    ok = abs(vp - vm) < 2 and np.sign(vp) == np.sign(vm) and max(abs(vp), abs(vm)) < 2 * max(min(abs(vp), abs(vm)), 1e-3)
                    cv[i] = 0.9 * (po[j][2] + ok)
                else:
                    cs[i], cv[i] = 0.0, 0.0
        prev[rx] = (t, [(o["xy"], a, b) for o, a, b in zip(objs, cs, cv)])
        keep = [o for o, a, b in zip(objs, cs, cv) if not (a < 1.7 and b < 1.1)]
        # trajectory level: IMM tracks
        for tr in tracks:
            dt = t - tr["t"]
            tr["t"] = t
            mu = tr["mu"]
            cbar = PI.T @ mu
            mix = (PI * mu[:, None]) / np.maximum(cbar[None, :], 1e-12)
            new = []
            for j in range(2):
                x0, P0 = C.moment_match(mix[:, j], tr["xs"], tr["Ps"])
                new.append(ct_predict(x0, P0, dt, q, q_w, turn=j == 1))
            tr["xs"], tr["Ps"] = [n[0] for n in new], [n[1] for n in new]
            tr["mu"] = cbar
        zs = [C.zvec(m) for m in keep]
        cost = np.full((len(tracks), len(keep)), 1e9)
        cache = {}
        for i, tr in enumerate(tracks):
            x, P = C.moment_match(tr["mu"], tr["xs"], tr["Ps"])
            for j, z in enumerate(zs):
                nu, S, H = C.innov(x[:4], P[:4, :4], z, geo, R)
                d2 = float(nu @ np.linalg.solve(S, nu))
                if d2 < C.gate(len(z), gate_p):
                    cost[i, j], cache[i, j] = d2, True
        pairs = [(i, j) for i, j in zip(*linear_sum_assignment(cost)) if (i, j) in cache] if cache else []
        hit, used = {i for i, _ in pairs}, {j for _, j in pairs}
        for i, j in pairs:
            tr = tracks[i]
            lik = []
            for k in range(2):
                x, P = tr["xs"][k], tr["Ps"][k]
                nu, S, H4 = C.innov(x[:4], P[:4, :4], zs[j], geo, R)
                H = np.zeros((len(zs[j]), 5))
                H[:, :4] = H4
                K = P @ H.T @ np.linalg.inv(S)
                tr["xs"][k], tr["Ps"][k] = x + K @ nu, (np.eye(5) - K @ H) @ P
                lik.append(C.gauss(nu, S))
            mu = tr["mu"] * np.array(lik)
            tr["mu"] = mu / mu.sum() if mu.sum() > 0 else np.array([0.5, 0.5])
            tr["run"] += 1
            tr["miss"] = 0
            tr["conf"] |= tr["run"] >= L
        for i, tr in enumerate(tracks):
            x, _ = C.moment_match(tr["mu"], tr["xs"], tr["Ps"])
            if i not in hit and C.in_fov(x, geo):
                tr["miss"] += 1
                tr["run"] = 0
        tracks = [tr for tr in tracks if tr["miss"] < max_miss and (tr["conf"] or tr["miss"] == 0)]
        for j, m in enumerate(keep):
            if j not in used and not any((i, j) in cache for i in range(len(tracks))):
                x0, P0 = C.birth(m, geo)
                x5, P5 = np.r_[x0, 0.0], np.pad(P0, ((0, 1), (0, 1)))
                P5[4, 4] = 0.3 ** 2
                tracks.append(dict(id=nid, xs=[x5.copy(), x5.copy()], Ps=[P5.copy(), P5.copy()], mu=np.array([0.5, 0.5]),
                                   t=t, run=1, miss=0, conf=False))
                nid += 1
        res = []
        for tr in tracks:
            if tr["conf"]:
                x, _ = C.moment_match(tr["mu"], tr["xs"], tr["Ps"])
                res.append((tr["id"], *map(float, x[:4])))
        out.append(res)
    return out
