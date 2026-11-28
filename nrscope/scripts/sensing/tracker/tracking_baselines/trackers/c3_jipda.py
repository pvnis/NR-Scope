"""C3: JIPDA, joint integrated probabilistic data association (Musicki & Evans 2004): JPDA
(Fortmann, Bar-Shalom & Scheffe 1983) with a track existence probability (IPDA, Musicki et al.
1994), for an unknown number of targets. Per scan: gating, clusters of tracks sharing
measurements, exact enumeration of the joint association events in each cluster (each track at
most one measurement, each measurement at most one track), the marginal association
probabilities, then per track a moment-matched mixture (missed + each measurement) and the
existence update. Unassigned localised measurements (in no gate) start tracks with r_birth.
Reported when r >= r_confirm; deleted when r < r_delete. p_d inside the field of view, 0 outside."""
import itertools

import numpy as np

from . import common as C


def run(frames, q=0.35, p_s=0.9, p_d=0.9, clutter=2.0, r_birth=0.1, r_confirm=0.9, r_delete=0.01, gate_p=0.999):
    tracks, out, nid = [], [], 1
    for t, meas, rb, sc, geo, rx in frames:
        R = C.Rz(rb, sc)
        for tr in tracks:
            dt = t - tr["t"]
            tr["r"] *= p_s ** dt
            tr["x"], tr["P"] = C.predict(tr["x"], tr["P"], dt, q)
            tr["t"] = t
        zs = [C.zvec(m) for m in meas]
        pd = [p_d if C.in_fov(tr["x"], geo) else 0.0 for tr in tracks]
        lr = {}    # (i, j): r pd g / kappa, the measurement's likelihood ratio
        cache = {}
        for i, tr in enumerate(tracks):
            if pd[i] == 0:
                continue
            for j, z in enumerate(zs):
                nu, S, H = C.innov(tr["x"], tr["P"], z, geo, R)
                if float(nu @ np.linalg.solve(S, nu)) < C.gate(len(z), gate_p):
                    cache[i, j] = (nu, S, H)
                    lr[i, j] = tr["r"] * pd[i] * C.gauss(nu, S) / C.kappa(clutter, len(z))
        beta = joint_marginals(len(tracks), [tr["r"] * d for tr, d in zip(tracks, pd)], lr)
        for i, tr in enumerate(tracks):
            b0, bj = beta[i]
            miss = b0 * tr["r"] * (1 - pd[i]) / max(1 - tr["r"] * pd[i], 1e-12)   # exists, not detected
            ws, xs, Ps = [miss], [tr["x"]], [tr["P"]]
            for j, b in bj.items():
                xu, Pu = C.update(tr["x"], tr["P"], *cache[i, j])
                ws.append(b), xs.append(xu), Ps.append(Pu)
            tr["r"] = min(sum(ws), 1.0)
            if tr["r"] > 0:
                tr["x"], tr["P"] = C.moment_match(ws, xs, Ps)
            tr["conf"] |= tr["r"] >= r_confirm
        tracks = [tr for tr in tracks if tr["r"] >= r_delete]
        gated = {j for _, j in cache}
        for j, m in enumerate(meas):
            if j not in gated and m["xy"] is not None:
                x0, P0 = C.birth(m, geo)
                tracks.append(dict(id=nid, x=x0, P=P0, t=t, r=r_birth, conf=False))
                nid += 1
        out.append([(tr["id"], *map(float, tr["x"])) for tr in tracks if tr["conf"] and tr["r"] >= r_confirm / 2])
    return out


def joint_marginals(n, rpd, lr, max_events=200000):
    """Exact JPDA marginals per cluster. Event weight: prod over tracks of (lr[i, j] if assigned j,
    else 1 - r pd). Returns per track (beta_0, {j: beta_j})."""
    res = [(1.0, {}) for _ in range(n)]
    adj = {}
    for i, j in lr:
        adj.setdefault(i, set()).add(j)
    seen = set()
    for i0 in adj:
        if i0 in seen:
            continue
        ti, mj, stack = set(), set(), [i0]   # the cluster: tracks linked through shared measurements
        while stack:
            i = stack.pop()
            if i in ti:
                continue
            ti.add(i)
            for j in adj[i]:
                mj.add(j)
                stack += [k for k in adj if j in adj[k] and k not in ti]
        seen |= ti
        ti = sorted(ti)
        opts = [[None] + sorted(adj[i]) for i in ti]
        tot, b0, bj = 0.0, np.zeros(len(ti)), [dict() for _ in ti]
        count = 0
        for ev in itertools.product(*opts):
            js = [j for j in ev if j is not None]
            if len(js) != len(set(js)):
                continue
            w = 1.0
            for k, (i, j) in enumerate(zip(ti, ev)):
                w *= (1 - rpd[i]) if j is None else lr[i, j]
            tot += w
            for k, j in enumerate(ev):
                if j is None:
                    b0[k] += w
                else:
                    bj[k][j] = bj[k].get(j, 0.0) + w
            count += 1
            if count > max_events:
                break
        for k, i in enumerate(ti):
            res[i] = (b0[k] / tot, {j: v / tot for j, v in bj[k].items()})
    return res
