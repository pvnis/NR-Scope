"""C6: extended-target GGIW-PHD filter (Granstrom & Orguner, IEEE TSP 2012 "A PHD filter for
tracking multiple extended targets using random matrices"; GGIW: Granstrom & Orguner 2012,
"Estimation and maintenance of measurement rates"; survey: Granstrom, Baum & Reuter, JAIF 2017).
Each component: weight, Poisson measurement rate ~ Gamma(a, b), kinematics N(x, P), extent
X ~ IW(v, V) in the plane. Per scan: prediction (rate and extent forgetting, eta and tau),
distance partitioning of the localised measurements (single linkage at several thresholds, plus
all-singletons), and the ET-PHD update over every partition and cell:
  missed:     w (1 - p_d + p_d (b/(b+1))^a)
  cell W:     w p_d Gamma-Poisson(|W|) N(z_W; h, S_W) / (prod kappa) / d_W
  partition:  weight prod_W d_W,  d_W = [|W| == 1] + sum_j (cell likelihoods).
Kinematics: EKF on the stacked bistatic measurements of the cell, each with the extent's spread
(z X, z = 1/4, through the position Jacobian) and the limb spread in speed as in our tracker;
extent: Feldmann et al. 2011 update with the cell's localised positions. Births: one per cell
of the coarsest partition. Estimates: weight > 0.5; labels by tagging (Panta et al. 2009).
Measurements without a position (no angle) are not used (they cannot be partitioned)."""
import math

import numpy as np

from . import common as C

Z = 0.25
X0 = 0.3 ** 2      # birth extent (a person)
LIMB = 1.0         # speed spread of a person's limbs (m/s)


def partitions(xy, dists):
    """Single-linkage partitions of the points at each distance threshold (distinct ones)."""
    n = len(xy)
    out = {tuple(range(n)): [[i] for i in range(n)]}
    for d in dists:
        parent = list(range(n))

        def root(i):
            while parent[i] != i:
                i = parent[i]
            return i
        for i in range(n):
            for j in range(i + 1, n):
                if np.hypot(*(xy[i] - xy[j])) <= d:
                    parent[root(j)] = root(i)
        cells = {}
        for i in range(n):
            cells.setdefault(root(i), []).append(i)
        out[tuple(root(i) for i in range(n))] = list(cells.values())
    uniq = {}
    for cells in out.values():
        uniq[tuple(sorted(tuple(c) for c in cells))] = cells
    return list(uniq.values())


def cell_update(c, ms, geo, R, p_d, clutter):
    """Likelihood (times prod kappa^-1) and updated component for cell ms."""
    w, a, b, x, P, v, V, l = c
    Xh = V / max(v - 6, 1.0)
    zs, Hs, Rs, nus = [], [], [], []
    for m in ms:
        z = C.zvec(m)
        zh, H4, Rh = C.track.bistatic_h(x, *geo)
        k = len(z)
        nu = z - zh[:k]
        if k == 3:
            nu[2] = (nu[2] + np.pi) % (2 * np.pi) - np.pi
        Rk = R[:k, :k] + Rh[:k, :k] + Z * H4[:k, :2] @ Xh @ H4[:k, :2].T
        Rk[1, 1] += LIMB ** 2 if abs(x[2:4] @ x[2:4]) < 2.5 ** 2 else 0.0
        nus.append(nu), Hs.append(H4[:k]), Rs.append(Rk)
        zs.append(k)
    H = np.vstack(Hs)
    nu = np.concatenate(nus)
    Rb = np.zeros((len(nu), len(nu)))
    o = 0
    for Rk in Rs:
        Rb[o:o + len(Rk), o:o + len(Rk)] = Rk
        o += len(Rk)
    S = H @ P @ H.T + Rb
    d2 = float(nu @ np.linalg.solve(S, nu))
    if d2 > C.gate(len(nu)):
        return 0.0, None
    g = math.exp(-0.5 * d2) / math.sqrt(np.linalg.det(2 * np.pi * S))
    n = len(ms)
    gp = math.exp(math.lgamma(a + n) - math.lgamma(a) + a * math.log(b) - (a + n) * math.log(b + 1) - math.lgamma(n + 1))
    kap = np.prod([C.kappa(clutter, k) for k in zs])
    L = p_d * gp * g / kap
    K = P @ H.T @ np.linalg.inv(S)
    xu, Pu = x + K @ nu, (np.eye(4) - K @ H) @ P
    pts = np.array([m["xy"] for m in ms])
    Vu = V
    if n >= 1:   # Feldmann-style extent update from the positions' spread about the predicted centre
        zbar = pts.mean(0)
        eps = zbar - x[:2]
        Ns = np.outer(eps, eps) * (Xh.trace() / max(Xh.trace() + P[:2, :2].trace(), 1e-9))
        Zs = (pts - zbar).T @ (pts - zbar)
        Vu = V + Ns + Zs
    return L, [None, a + n, b + 1, xu, Pu, v + n, Vu, l]


def run(frames, q=0.35, p_s=0.9, p_d=0.9, clutter=2.0, w_birth=0.05, eta=1.25, tau=2.0,
        dists=(1.0, 2.0, 3.0, 5.0), prune=1e-5, U=4.0, jmax=60):
    comps, out, nid, t_prev = [], [], 1, None   # [w, a, b, x, P, v, V, label]
    for t, meas, rb, sc, geo, rx in frames:
        R = C.Rz(rb, sc)
        dt = 0.0 if t_prev is None else t - t_prev
        t_prev = t
        f = math.exp(-dt / tau)
        for c in comps:
            c[0] *= p_s ** dt
            c[1], c[2] = c[1] / eta ** dt, c[2] / eta ** dt
            c[3], c[4] = C.predict(c[3], c[4], dt, q)
            v0 = c[5]
            c[5] = 6 + f * (v0 - 6)
            c[6] = c[6] * (c[5] - 6) / max(v0 - 6, 1e-9)
        ms = [m for m in meas if m["xy"] is not None]
        xy = [np.asarray(m["xy"]) for m in ms]
        pd = [p_d if C.in_fov(c[3], geo) else 0.0 for c in comps]
        new = [[c[0] * (1 - d + d * (c[2] / (c[2] + 1)) ** c[1]), *c[1:]] for c, d in zip(comps, pd)]
        if ms:
            parts = partitions(xy, dists)
            cache = {}
            dW = {}
            for cells in parts:
                for cell in cells:
                    key = tuple(cell)
                    if key in dW:
                        continue
                    ups = []
                    for c, d in zip(comps, pd):
                        if d == 0:
                            continue
                        L, cu = cell_update(c, [ms[i] for i in cell], geo, R, d, clutter)
                        if cu is not None and L > 0:
                            cu[0] = c[0] * L
                            ups.append(cu)
                    cache[key] = ups
                    dW[key] = (len(cell) == 1) + sum(u[0] for u in ups)
            wp = np.array([np.prod([dW[tuple(cl)] for cl in cells]) for cells in parts])
            wp = wp / wp.sum() if wp.sum() > 0 else np.ones(len(parts)) / len(parts)
            for cells, wpi in zip(parts, wp):
                for cell in cells:
                    key = tuple(cell)
                    for u in cache[key]:
                        new.append([wpi * u[0] / dW[key], *u[1:]])
            coarse = min(parts, key=len)
            for cell in coarse:
                m = ms[cell[0]]
                x0, P0 = C.birth(m, geo)
                x0[:2] = np.mean([xy[i] for i in cell], 0)
                new.append([w_birth, 2.0, 1.0, x0, P0, 10.0, np.eye(2) * X0 * 4, nid])
                nid += 1
        comps = merge([c for c in new if c[0] > prune], U, jmax)
        est, seen = [], set()
        for c in sorted(comps, key=lambda c: -c[0]):
            if c[0] > 0.5:
                if c[7] in seen:   # a label shared by two targets: a new one
                    c[7] = nid
                    nid += 1
                seen.add(c[7])
                est.append((c[7], *map(float, c[3])))
        out.append(est)
    return out


def merge(comps, U, jmax):
    comps = sorted(comps, key=lambda c: -c[0])
    out = []
    while comps:
        c0 = comps[0]
        Pi = np.linalg.inv(c0[4])
        near = [c for c in comps if float((c[3] - c0[3]) @ Pi @ (c[3] - c0[3])) <= U]
        ids = {id(c) for c in near}
        comps = [c for c in comps if id(c) not in ids]
        ws = np.array([c[0] for c in near])
        w = ws.sum()
        x, P = C.moment_match(ws, [c[3] for c in near], [c[4] for c in near])
        a = float(ws @ [c[1] for c in near] / w)
        b = float(ws @ [c[2] for c in near] / w)
        v = float(ws @ [c[5] for c in near] / w)
        V = sum(wi * c[6] for wi, c in zip(ws, near)) / w
        out.append([w, a, b, x, P, v, V, c0[7]])
    return out[:jmax]
