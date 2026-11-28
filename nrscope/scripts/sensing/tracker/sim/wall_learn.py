"""Learn the reflecting walls from the tracks themselves, without any map (the tracker's
"ghost walls: learned" option).

A target's echo that bounces once more on a wall on its way to NR-Scope (target -> wall -> NR-Scope)
is tracked as a ghost: a second track, farther in bistatic range, that moves with the first.
For every pair of tracks (P, ghost X) in a map, if X were P's out-bounce echo, the bounce point
W lies on NR-Scope's ray towards X at the distance where gNB -> P -> W -> NR-Scope has X's path
length, and the wall through W reflects P -> W into W -> NR-Scope (normal = bisector). The pair is
kept only if X's bistatic range rate is the one P's velocity gives along gNB -> P -> W. Every
kept pair votes for a wall line (normal angle, distance from the origin); a line voted for in at
least MIN_SECONDS different seconds, by at least MIN_PAIRS pairs that each voted for it at least PAIR_MIN times, is a wall, as long as the bounce points that voted for it
(plus EXT_M at each end). Coincidences (two people walking together) do not pile up on one line.
"""
import numpy as np

MIN_LIFE = 10          # maps a track must live to take part
RATE_TOL_MS = 0.5      # ghost's range rate vs the one predicted from its parent (m/s)
MIN_SEP_M = 5.0        # bounce point at least this far from the parent
DPHI = np.deg2rad(3)   # vote cell: wall angle
DRHO = 1.0             # vote cell: wall distance (m)
MIN_SECONDS = 3        # distinct seconds of votes for a wall
MIN_PAIRS = 2          # different pairs of tracks per wall (a wall reflects everyone who passes; one pair can be a coincidence)
PAIR_MIN = 5           # votes one pair of tracks must give the same wall (a reflection persists, a coincidence does not)
EXT_M = 5.0            # wall beyond the bounce points seen


def _u(v):
    return v / max(np.linalg.norm(v), 1e-9)


def learn(frames, gnb, rx):
    """frames: [(t, [(id, x, y, vx, vy, ...)])] confirmed tracks per map; gnb, rx: positions.
    Returns walls [dict(a, b, votes, seconds)] (a, b: the segment's ends, local metres)."""
    T, R = np.asarray(gnb, float)[:2], np.asarray(rx, float)[:2]
    life = {}
    for _, tr in frames:
        for q in tr:
            life[q[0]] = life.get(q[0], 0) + 1
    votes = []   # (t, phi, rho, W, pair)
    for t, tr in frames:
        tr = [q for q in tr if life[q[0]] >= MIN_LIFE]
        for p in tr:
            P, vp = np.array(p[1:3], float), np.array(p[3:5], float)
            Lp = np.linalg.norm(P - T) + np.linalg.norm(P - R)
            for g in tr:
                if g[0] == p[0]:
                    continue
                X, vg = np.array(g[1:3], float), np.array(g[3:5], float)
                L = np.linalg.norm(X - T) + np.linalg.norm(X - R)
                if L <= Lp + 1:
                    continue
                d, D, q = _u(X - R), L - np.linalg.norm(P - T), P - R
                den = 2 * (D - d @ q)
                if abs(den) < 1e-6:
                    continue
                s = (D * D - q @ q) / den   # |P - (R + s d)| = D - s
                if not 1 < s < D:
                    continue
                W = R + s * d
                if np.linalg.norm(W - P) < MIN_SEP_M:
                    continue
                rate = vg @ (_u(X - T) + _u(X - R))
                if abs(rate) < 0.3 or abs(vp @ (_u(P - T) + _u(P - W)) - rate) > RATE_TOL_MS:
                    continue
                n = _u(_u(P - W) + _u(R - W))
                phi, rho = np.arctan2(n[1], n[0]), n @ W
                if rho < 0:
                    phi, rho = phi + np.pi, -rho
                votes.append((t, phi % (2 * np.pi), rho, W, (p[0], g[0])))
    if not votes:
        return []
    ts = np.array([v[0] for v in votes])
    ph = np.array([v[1] for v in votes])
    rh = np.array([v[2] for v in votes])
    Ws = np.array([v[3] for v in votes])
    pair = [v[4] for v in votes]
    cells = {}
    for i, k in enumerate(zip(np.round(ph / DPHI).astype(int), np.round(rh / DRHO).astype(int))):
        cells.setdefault(k, []).append(i)
    walls, used = [], set()
    near = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    for k in sorted(cells, key=lambda k: -len(cells[k])):
        if k in used:
            continue
        idx = [i for a, b in near for i in cells.get((k[0] + a, k[1] + b), [])]
        n_pair = {}
        for i in idx:
            n_pair[pair[i]] = n_pair.get(pair[i], 0) + 1
        idx = [i for i in idx if n_pair[pair[i]] >= PAIR_MIN]
        if len({pair[i] for i in idx}) < MIN_PAIRS:
            continue
        secs = len(set(np.round(ts[idx]).astype(int)))
        if secs < MIN_SECONDS:
            continue
        used.update((k[0] + a, k[1] + b) for a, b in near)
        p, r = np.median(ph[idx]), np.median(rh[idx])
        n = np.array([np.cos(p), np.sin(p)])
        along = np.array([-n[1], n[0]])
        s = Ws[idx] @ along
        walls.append(dict(a=(r * n + (s.min() - EXT_M) * along).round(2).tolist(),
                          b=(r * n + (s.max() + EXT_M) * along).round(2).tolist(), votes=len(idx), seconds=secs))
    return walls


def as_walls(walls):
    """learn()'s walls in the tracker's form (a, b, base z, top z), or None."""
    if not walls:
        return None
    A, B = np.array([w["a"] for w in walls]), np.array([w["b"] for w in walls])
    return A, B, np.full(len(A), -1e9), np.full(len(A), 1e9)
