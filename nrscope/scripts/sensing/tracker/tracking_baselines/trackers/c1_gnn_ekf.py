"""C1: the textbook tracker. EKF per track on the bistatic measurements, global nearest neighbour
association (chi^2 gate on the Mahalanobis distance, Hungarian assignment), M-of-N track logic:
confirmed after `confirm` hits, deleted after `max_miss` consecutive misses in the receiver's field
of view (a target out of a receiver's view is not missed by it). Unassigned localised
measurements start tentative tracks; a tentative track missed once is dropped.
Refs: Blackman & Popoli 1999 (GNN, M/N logic); Bar-Shalom, Willett & Tian 2011.
Also the tracking core of BSense (Meng et al. 2026) and of Gunes & Caliskan 2025."""
import numpy as np
from scipy.optimize import linear_sum_assignment

from . import common as C


def run(frames, q=0.35, confirm=3, max_miss=3, gate_p=0.999):
    tracks, out, nid = [], [], 1
    for t, meas, rb, sc, geo, rx in frames:
        R = C.Rz(rb, sc)
        for tr in tracks:
            tr["x"], tr["P"] = C.predict(tr["x"], tr["P"], t - tr["t"], q)
            tr["t"] = t
        zs = [C.zvec(m) for m in meas]
        cost = np.full((len(tracks), len(meas)), 1e9)
        cache = {}
        for i, tr in enumerate(tracks):
            for j, z in enumerate(zs):
                nu, S, H = C.innov(tr["x"], tr["P"], z, geo, R)
                d2 = float(nu @ np.linalg.solve(S, nu))
                if d2 < C.gate(len(z), gate_p):
                    cost[i, j], cache[i, j] = d2, (nu, S, H)
        pairs = []
        if cache:
            pairs = [(i, j) for i, j in zip(*linear_sum_assignment(cost)) if (i, j) in cache]
        hit, used = {i for i, _ in pairs}, {j for _, j in pairs}
        for i, j in pairs:
            tr = tracks[i]
            tr["x"], tr["P"] = C.update(tr["x"], tr["P"], *cache[i, j])
            tr["hits"] += 1
            tr["miss"] = 0
            tr["conf"] |= tr["hits"] >= confirm
        for i, tr in enumerate(tracks):
            if i not in hit and C.in_fov(tr["x"], geo):
                tr["miss"] += 1
        tracks = [tr for tr in tracks if tr["miss"] < max_miss and (tr["conf"] or tr["miss"] == 0)]
        gated = used | {j for _, j in cache}
        for j, m in enumerate(meas):
            if j not in gated and m["xy"] is not None:
                x0, P0 = C.birth(m, geo)
                tracks.append(dict(id=nid, x=x0, P=P0, t=t, hits=1, miss=0, conf=False))
                nid += 1
        out.append([(tr["id"], *map(float, tr["x"])) for tr in tracks if tr["conf"]])
    return out
