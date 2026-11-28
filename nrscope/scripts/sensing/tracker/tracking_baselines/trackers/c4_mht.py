"""C4: multiple hypothesis tracking (Reid 1979), hypothesis-oriented, with the M best assignments
of each hypothesis by Murty's algorithm (Cox & Hingorani 1996) and track scores as log-likelihood
ratios (Blackman 2004; Blackman & Popoli 1999). Per scan, every kept global hypothesis spawns its
k best assignments of the scan's measurements to: its tracks (score + log(pd g / kappa)), new
tentative tracks (score log(mu_n / clutter)), with misses scoring log(1 - pd) (pd inside the
field of view, 0 outside). The n_hyp best global hypotheses survive (normalised). Tracks are
confirmed when their score exceeds s_conf (SPRT threshold) and deleted when it falls s_del below
its peak. Reported: the confirmed tracks of the best hypothesis."""
import numpy as np
from scipy.optimize import linear_sum_assignment

from . import common as C

BIG = 1e6


def kbest(cost, k):
    """Murty's k best assignments of the rows of a rectangular cost matrix (rows <= columns):
    [(total, cols)]."""
    if cost.shape[0] == 0:
        return [(0.0, np.zeros(0, int))]

    def solve(c):
        r, cl = linear_sum_assignment(c)
        tot = c[r, cl].sum()
        return (tot, cl) if tot < BIG / 2 else None
    first = solve(cost)
    if first is None:
        return []
    out, queue = [], [(first[0], first[1], cost)]
    while queue and len(out) < k:
        queue.sort(key=lambda q: q[0])
        tot, cl, c = queue.pop(0)
        out.append((tot, cl))
        c = c.copy()
        for i in range(len(cl)):   # partition: fix rows < i, forbid (i, cl[i])
            ci = c.copy()
            ci[i, cl[i]] = BIG
            s = solve(ci)
            if s is not None:
                queue.append((s[0], s[1], ci))
            keep = c[i, cl[i]]
            c[i, :] = BIG
            c[:, cl[i]] = BIG
            c[i, cl[i]] = keep
    return out


def run(frames, q=0.35, p_d=0.9, clutter=2.0, mu_n=0.05, n_hyp=10, k=5, s_conf=4.6, s_del=6.9, gate_p=0.999):
    hyps = [(0.0, [])]    # (log weight, tracks); track: dict(id, x, P, s, peak, conf)
    out, nid, t_prev = [], [1], None
    for t, meas, rb, sc, geo, rx in frames:
        R = C.Rz(rb, sc)
        dt = 0.0 if t_prev is None else t - t_prev
        t_prev = t
        zs = [C.zvec(m) for m in meas]
        M = len(meas)
        new_hyps = []
        for lw, tracks in hyps:
            pred = []
            for tr in tracks:
                x, P = C.predict(tr["x"], tr["P"], dt, q)
                pred.append(dict(tr, x=x, P=P))
            n = len(pred)
            # rows: tracks; columns: measurements (M), then one "missed" column per track (n).
            # A measurement no track takes is a false alarm or a new tentative track (below).
            cost = np.full((n, M + n), BIG)
            cache = {}
            for i, tr in enumerate(pred):
                pd = p_d if C.in_fov(tr["x"], geo) else 0.0
                cost[i, M + i] = -np.log(max(1 - pd, 1e-9))
                if pd == 0:
                    continue
                for j, z in enumerate(zs):
                    nu, S, H = C.innov(tr["x"], tr["P"], z, geo, R)
                    if float(nu @ np.linalg.solve(S, nu)) < C.gate(len(z), gate_p):
                        cache[i, j] = (nu, S, H)
                        cost[i, j] = -np.log(max(pd * C.gauss(nu, S) / C.kappa(clutter, len(z)), 1e-300))
            for tot, cl in kbest(cost, k):
                tr_new = []
                for i, tr in enumerate(pred):
                    j = cl[i]
                    tr = dict(tr)
                    if j < M:
                        tr["x"], tr["P"] = C.update(tr["x"], tr["P"], *cache[i, j])
                        tr["s"] += -cost[i, j]
                    else:
                        tr["s"] += -cost[i, M + i]
                    tr["peak"] = max(tr["peak"], tr["s"])
                    tr["conf"] = tr["conf"] or tr["s"] >= s_conf
                    if tr["s"] >= tr["peak"] - s_del:
                        tr_new.append(tr)
                used = {cl[i] for i in range(n) if cl[i] < M}
                for j, m in enumerate(meas):
                    if j not in used and m["xy"] is not None:
                        x0, P0 = C.birth(m, geo)
                        s0 = np.log(mu_n / clutter)
                        tr_new.append(dict(id=nid[0], x=x0, P=P0, s=s0, peak=s0, conf=False))
                        nid[0] += 1
                new_hyps.append((lw - tot, tr_new))
        new_hyps.sort(key=lambda h: -h[0])
        hyps = new_hyps[:n_hyp]
        top = hyps[0][0]
        hyps = [(lw - top, trs) for lw, trs in hyps]
        out.append([(tr["id"], *map(float, tr["x"])) for tr in hyps[0][1] if tr["conf"]])
    return out
