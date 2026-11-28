#!/usr/bin/env python3
"""Scenario engine: many people and vehicles moving in a scene, generated once from a spec.

    scenario.py generate SPEC.json OUT.npz      trajectories of every agent
    scenario.py preview SPEC.json               a decimated JSON of them (what the GUI animates)

The spec (the GUI's Scenarios page writes it):
    {"scene": "simple_street_canyon", "seed": 1, "duration_s": 120, "dt": 0.05,
     "groups": [
       {"kind": "pedestrians", "rect": [x0, y0, x1, y1], "n": 10, "speed": [0.8, 1.5],
        "stop_prob": 0.3, "stop_s": [1, 6], "turn_rate": 0.05,
        "group_frac": 0.4, "group_size": [2, 4], "split_rate": 0.01, "rcs": 0.5},
       {"kind": "lane", "points": [[x, y], ...], "rate_per_min": 6, "speed": [6, 12],
        "mix": {"bike": 0.1, "car": 0.6, "van": 0.2, "truck": 0.1}},
       {"kind": "walkway", "points": [[x, y], ...], "rate_per_min": 6, "width_m": 3,
        "both_ways": 0.5, "speed": [1.0, 1.5], "group_frac": 0.3, "group_size": [2, 3],
        "min_spacing_m": 0, "rcs": 0.5}]}

People (one zone per "pedestrians" group) walk between random waypoints inside their zone
with the social-force model (Helbing & Molnar 1995): each is pulled towards its goal at its
own preferred speed and pushed away from the others and from the buildings. At a waypoint a
person stops with probability stop_prob for a random time, then picks a new one; turn_rate
(per second) is the rate of changing one's mind mid-way (a turn). A fraction group_frac
walks in groups: followers keep a place beside their leader, stop when it stops, and leave
the group at split_rate per second. min_spacing_m > 0 keeps people who are not in the same
group at least that far apart (a stronger push from 2 m before it, and a step that would
bring two closer than it is not taken); 0 leaves only the usual ~1 m push.

Walkers (one "walkway" group per walkway: a polyline of width width_m, like a lane for
people) arrive at random (Poisson, rate_per_min), a share both_ways of them at the far end,
walk it at their own speed, each at a random offset across its width with a slight sway,
and leave at the other end. A share group_frac arrives as a group of group_size walking side
by side. Walkers keep out of each other's way (min_spacing_m, at least 0.8 m): they step aside
within the width, and slow down behind a slower one going their way. A walkway ignores the
obstacles: draw it on open ground.

Vehicles (one "lane" group per lane, a polyline driven from its first point to its last)
arrive at random (Poisson, rate_per_min) with a class drawn from mix and a desired speed
drawn from speed, and follow the vehicle ahead with the Intelligent Driver Model (Treiber,
Hennecke & Helbing 2000): they keep a safe time gap, brake behind slower vehicles and leave
at the lane's end. A vehicle enters only when there is room.

Every agent has a class (CLASSES: size and radar cross-section), a birth and a death time
(a zone's people live the whole run, a walkway's come and go), and its position and velocity every dt. Agents are numbered
from 1 in the order of their birth, which is also their number in the ground truth.
"""
import json
import sys

import numpy as np

# length, width, height (m) of the scatterer body; RCS (m^2) at 3.5 GHz, order of magnitude
CLASSES = {
    "person": dict(size=(0.5, 0.5, 1.8), rcs=0.5),
    "bike": dict(size=(1.8, 0.6, 1.5), rcs=2.0),
    "car": dict(size=(4.4, 1.8, 1.4), rcs=10.0),
    "van": dict(size=(5.5, 2.0, 2.2), rcs=20.0),
    "truck": dict(size=(9.0, 2.5, 3.5), rcs=50.0),
}
# IDM: max acceleration, comfortable braking (m/s^2), time gap (s), standstill gap (m)
IDM = {"bike": (1.0, 2.0, 1.2, 1.5), "car": (1.5, 2.0, 1.5, 2.0), "van": (1.2, 2.0, 1.6, 2.0),
       "truck": (0.8, 1.5, 1.8, 3.0)}

PED_DEFAULTS = dict(n=10, speed=[0.8, 1.5], stop_prob=0.3, stop_s=[1.0, 6.0], turn_rate=0.05,
                    group_frac=0.4, group_size=[2, 4], split_rate=0.01, min_spacing_m=0.0, rcs=0.5)
WALK_DEFAULTS = dict(rate_per_min=6.0, width_m=3.0, both_ways=0.5, speed=[1.0, 1.5],
                     group_frac=0.3, group_size=[2, 3], min_spacing_m=0.0, rcs=0.5)
LANE_DEFAULTS = dict(rate_per_min=6.0, speed=[6.0, 12.0], mix={"bike": 0.1, "car": 0.6, "van": 0.2, "truck": 0.1})


# ---- geometry -------------------------------------------------------------------------------

def inside(p, poly):
    """p inside the polygon (even-odd rule)."""
    x, y, c = p[0], p[1], False
    for i in range(len(poly)):
        (x1, y1), (x2, y2) = poly[i], poly[i - 1]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            c = not c
    return c


def seg_dist(p, a, b):
    """Distance from p to the segment ab and the unit vector from it to p."""
    ab, ap = b - a, p - a
    t = np.clip(ap @ ab / max(ab @ ab, 1e-12), 0, 1)
    d = p - (a + t * ab)
    n = np.linalg.norm(d)
    return n, d / n if n > 1e-9 else np.zeros(2)


def crosses(a, b, poly):
    """The segment ab crosses or enters the polygon."""
    if inside(b, poly):
        return True
    for i in range(len(poly)):
        c, d = np.asarray(poly[i]), np.asarray(poly[i - 1])
        o = lambda p, q, r: np.sign((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
        if o(a, b, c) != o(a, b, d) and o(c, d, a) != o(c, d, b):
            return True
    return False


class Polyline:
    def __init__(self, pts):
        self.p = np.asarray(pts, float)[:, :2]
        seg = np.diff(self.p, axis=0)
        self.len = np.linalg.norm(seg, axis=1)
        self.cum = np.r_[0, np.cumsum(self.len)]
        self.length = float(self.cum[-1])
        self.dir = seg / np.maximum(self.len, 1e-9)[:, None]

    def at(self, s):
        """Point and unit direction at arc length s."""
        i = int(np.clip(np.searchsorted(self.cum, s, side="right") - 1, 0, len(self.len) - 1))
        return self.p[i] + self.dir[i] * (s - self.cum[i]), self.dir[i]


# ---- people ---------------------------------------------------------------------------------

def simulate_people(g, buildings, rng, steps, dt):
    """One pedestrian zone: returns [(x, y, vx, vy) arrays of length steps] per person.
    buildings: footprints people never enter and walk around (buildings, and a real site's
    planters, benches, walls, tree trunks...: generate's obstacles)."""
    g = dict(PED_DEFAULTS, **g)
    x0, y0, x1, y1 = g["rect"]
    x0, x1, y0, y1 = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
    n = int(g["n"])
    if n <= 0:
        return []
    polys = [q for q in (np.asarray(b["footprint"], float) for b in buildings)   # the ones near the zone
             if q[:, 0].max() > x0 - 5 and q[:, 0].min() < x1 + 5 and q[:, 1].max() > y0 - 5 and q[:, 1].min() < y1 + 5]
    box = np.array([[q[:, 0].min(), q[:, 1].min(), q[:, 0].max(), q[:, 1].max()] for q in polys]).reshape(-1, 4)
    # every edge of every polygon, for the push (vectorised: a site has hundreds of small obstacles)
    ea = np.concatenate([np.roll(q, 1, axis=0) for q in polys]) if polys else np.zeros((0, 2))
    eb = np.concatenate(polys) if polys else np.zeros((0, 2))

    def free(p, m=0.0):   # outside every footprint (bounding boxes first: most are far)
        near = np.flatnonzero((box[:, 0] - m <= p[0]) & (p[0] <= box[:, 2] + m) & (box[:, 1] - m <= p[1]) & (p[1] <= box[:, 3] + m))
        return not any(inside(p, polys[i]) for i in near)

    def push(P):
        """The obstacles' push on every person: (n, 2), from every edge closer than 2 m."""
        if not len(ea):
            return np.zeros_like(P)
        ab = eb - ea                                                     # (E, 2)
        ap = P[:, None, :] - ea[None, :, :]                              # (n, E, 2)
        t = np.clip((ap * ab).sum(2) / np.maximum((ab * ab).sum(1), 1e-12), 0, 1)
        d = ap - t[:, :, None] * ab[None]                                # from the edge to the person
        dist = np.linalg.norm(d, axis=2)
        w = np.where(dist < 2.0, 5.0 * np.exp((0.3 - dist) / 0.2), 0.0) / np.maximum(dist, 1e-9)
        return (w[:, :, None] * d).sum(1)

    def waypoint(frm):
        for _ in range(60):
            w = np.array([rng.uniform(x0, x1), rng.uniform(y0, y1)])
            if free(w, 1.0) and not any(crosses(frm, w, polys[i]) for i in seg_near(frm, w)):   # 1 m clear: within the obstacles' push one never gets there
                return w
        for _ in range(60):   # walled in (pushed into a pocket): any free point, the push walks them out
            w = np.array([rng.uniform(x0, x1), rng.uniform(y0, y1)])
            if free(w, 0.4):
                return w
        return frm.copy()

    def seg_near(a, b):   # the obstacles whose box the segment ab's box overlaps
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        return np.flatnonzero((box[:, 0] <= hi[0]) & (box[:, 2] >= lo[0]) & (box[:, 1] <= hi[1]) & (box[:, 3] >= lo[1]))

    sp_min = float(g.get("min_spacing_m") or 0.0)
    pos, tries = [], 0
    while len(pos) < n:
        p = np.array([rng.uniform(x0, x1), rng.uniform(y0, y1)])
        tries += 1
        if tries > 20000 or (free(p, 0.4) and all(np.linalg.norm(p - q) > max(0.8, sp_min) for q in pos)):
            pos.append(p)
    pos = np.array(pos)
    v0 = rng.uniform(*g["speed"], n)
    # groups: leaders followed by 1..size-1 followers, side by side
    leader = -np.ones(n, int)
    offset = np.zeros((n, 2))
    i, n_grouped = 0, int(round(g["group_frac"] * n))
    while i < n_grouped - 1:
        size = int(rng.integers(g["group_size"][0], g["group_size"][1] + 1))
        size = min(size, n_grouped - i)
        if size < 2:
            break
        for k in range(1, size):
            leader[i + k] = i
            offset[i + k] = [-0.3 * k, 0.7 * (k if k % 2 else -k) / 2 + (0.35 if k % 2 else -0.35)]
            pos[i + k] = pos[i] + offset[i + k] * 0.5
            v0[i + k] = v0[i]
        i += size
    goal = np.array([waypoint(p) for p in pos])
    stop_until = np.full(n, -1.0)
    stuck = np.zeros(n)
    vel = np.zeros((n, 2))
    out = np.zeros((n, steps, 4))
    tau, A, B, R = 0.5, 2.0, 0.3, 0.6
    for k in range(steps):
        t = k * dt
        des = np.zeros((n, 2))
        for j in range(n):
            if leader[j] >= 0 and rng.random() < g["split_rate"] * dt:   # leaves its group
                leader[j] = -1
                goal[j] = waypoint(pos[j])
            if leader[j] >= 0:
                L = leader[j]
                h = vel[L] / np.linalg.norm(vel[L]) if np.linalg.norm(vel[L]) > 0.2 else np.array([1.0, 0.0])
                tgt = pos[L] + offset[j, 0] * h + offset[j, 1] * np.array([-h[1], h[0]])
                d = tgt - pos[j]
                dn = np.linalg.norm(d)
                sp = min(1.5 * v0[j], v0[j] * dn / 0.7 + np.linalg.norm(vel[L]))
                des[j] = d / dn * sp if dn > 0.05 else vel[L]
                continue
            if t < stop_until[j]:
                continue
            d = goal[j] - pos[j]
            dn = np.linalg.norm(d)
            if dn < 0.5:
                if rng.random() < g["stop_prob"]:
                    stop_until[j] = t + rng.uniform(*g["stop_s"])
                goal[j] = waypoint(pos[j])
                continue
            if rng.random() < g["turn_rate"] * dt:
                goal[j] = waypoint(pos[j])
            des[j] = d / dn * v0[j] * min(1.0, dn / 1.0 + 0.3)
        f = (des - vel) / tau
        # people push each other away (less within a group)
        diff = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(diff, axis=2) + np.eye(n) * 1e6
        same = (leader[:, None] == np.arange(n)[None, :]) | (leader[None, :] == np.arange(n)[:, None]) | \
               ((leader[:, None] == leader[None, :]) & (leader[:, None] >= 0))
        r = np.where(same, 0.45, max(R, sp_min))
        mag = A * np.exp((r - dist) / np.where(same, B, B if sp_min <= 0 else max(B, 0.7)))
        f += (mag[:, :, None] * diff / dist[:, :, None]).sum(axis=1)
        # and the buildings and obstacles push them away
        f += push(pos)
        vel = vel + f * dt
        sp = np.linalg.norm(vel, axis=1)
        cap = 1.6 * v0
        vel *= np.minimum(1, cap / np.maximum(sp, 1e-9))[:, None]
        new = pos + vel * dt
        # stay out of buildings and inside the zone (a margin of 2 m)
        for j in range(n):
            if not free(new[j]):
                new[j], vel[j] = pos[j], vel[j] * 0
                stuck[j] += dt
            elif des[j] @ des[j] > 0.09 and vel[j] @ vel[j] < 0.01:   # wants to walk but held still (pull and push cancel)
                stuck[j] += dt
            else:
                stuck[j] = 0.0
            if stuck[j] > 3.0 and leader[j] < 0:   # stuck against an obstacle: somewhere else
                goal[j], stuck[j] = waypoint(pos[j]), 0.0
        new[:, 0] = np.clip(new[:, 0], x0 - 2, x1 + 2)
        new[:, 1] = np.clip(new[:, 1], y0 - 2, y1 + 2)
        if sp_min > 0:   # no step that brings two people of different groups closer than sp_min
            for j in range(n):
                dn = np.linalg.norm(new[j] - new, axis=1)
                do = np.linalg.norm(pos[j] - pos, axis=1)
                bad = (dn < sp_min) & (dn < do) & ~same[j]
                bad[j] = False
                if bad.any():   # slide: drop the part of the step towards them, else stay
                    v = vel[j].copy()
                    for q in np.flatnonzero(bad):
                        u = (pos[q] - pos[j]) / max(do[q], 1e-9)
                        v -= max(v @ u, 0.0) * u
                    trial = pos[j] + v * dt
                    dt_ = np.linalg.norm(trial - new, axis=1)
                    ok = free(trial) and not ((dt_ < sp_min) & (dt_ < do) & ~same[j] & (np.arange(n) != j)).any()
                    new[j], vel[j] = (trial, v) if ok else (pos[j], vel[j] * 0)
                    stuck[j] = 0.0 if ok and v @ v > 0.04 else stuck[j] + dt
        pos = new
        out[:, k, :2], out[:, k, 2:] = pos, vel
    return [out[j] for j in range(n)]


# ---- vehicles -------------------------------------------------------------------------------

def simulate_lane(g, rng, steps, dt):
    """One lane: returns [(class, k_birth, array of (x, y, vx, vy))] per vehicle."""
    g = dict(LANE_DEFAULTS, **g)
    lane = Polyline(g["points"])
    mix = {c: float(w) for c, w in g["mix"].items() if c in IDM and float(w) > 0} or {"car": 1.0}
    names, w = list(mix), np.array(list(mix.values()))
    rate = float(g["rate_per_min"]) / 60
    arrivals, t = [], rng.exponential(1 / rate) if rate > 0 else np.inf
    t = min(t, 1.0)   # the first one comes soon
    while t < steps * dt:
        arrivals.append(t)
        t += rng.exponential(1 / rate)
    queue = []   # (class, v0) waiting for room at the lane's start
    cars = []    # dict(cls, L, v0, s, v, k0, rec)
    done = []
    ai = 0
    for k in range(steps):
        t = k * dt
        while ai < len(arrivals) and arrivals[ai] <= t:
            c = names[rng.choice(len(names), p=w / w.sum())]
            queue.append((c, rng.uniform(*g["speed"]) * (0.6 if c == "bike" else 1.0)))
            ai += 1
        if queue:
            c, v0 = queue[0]
            L = CLASSES[c]["size"][0]
            last = cars[-1] if cars else None
            room = last is None or last["s"] - (last["L"] + L) / 2 > IDM[c][3] + 0.5 * v0 * IDM[c][2]
            if room:
                v = v0 if last is None else min(v0, last["v"])
                cars.append(dict(cls=c, L=L, v0=v0, s=0.0, v=v, k0=k, rec=[]))
                queue.pop(0)
        # IDM, front to back
        for i, c in enumerate(cars):
            a, b, T, s0 = IDM[c["cls"]]
            acc = a * (1 - (c["v"] / c["v0"]) ** 4)
            if i > 0:
                ld = cars[i - 1]
                gap = max(ld["s"] - c["s"] - (ld["L"] + c["L"]) / 2, 0.1)
                ss = s0 + max(0.0, c["v"] * T + c["v"] * (c["v"] - ld["v"]) / (2 * np.sqrt(a * b)))
                acc -= a * (ss / gap) ** 2
            c["acc"] = max(acc, -9.0)
        for c in cars:
            c["v"] = max(0.0, c["v"] + c["acc"] * dt)
            c["s"] += c["v"] * dt
            p, h = lane.at(c["s"])
            c["rec"].append((p[0], p[1], c["v"] * h[0], c["v"] * h[1]))
        while cars and cars[0]["s"] > lane.length:
            done.append(cars.pop(0))
    return [(c["cls"], c["k0"], np.array(c["rec"])) for c in done + cars if c["rec"]]


# ---- walkways --------------------------------------------------------------------------------

def simulate_walkway(g, rng, steps, dt):
    """One walkway: returns [(k_birth, array of (x, y, vx, vy))] per walker."""
    g = dict(WALK_DEFAULTS, **g)
    path = Polyline(g["points"])
    half = max(float(g["width_m"]), 0.0) / 2
    sp_min = max(float(g.get("min_spacing_m") or 0.0), 0.8)
    rate = float(g["rate_per_min"]) / 60
    arrivals, t = [], min(rng.exponential(1 / rate), 1.0) if rate > 0 else np.inf
    while t < steps * dt:
        arrivals.append(t)
        t += rng.exponential(1 / rate)
    walkers, done, ai = [], [], 0   # dict(fwd, v0, s, lat, home, phase, k0, rec)
    for k in range(steps):
        t = k * dt
        while ai < len(arrivals) and arrivals[ai] <= t:
            fwd = rng.random() >= g["both_ways"]
            size = int(rng.integers(g["group_size"][0], g["group_size"][1] + 1)) if rng.random() < g["group_frac"] else 1
            v0 = rng.uniform(*g["speed"])
            lat0 = rng.uniform(-half, half)
            for m in range(max(size, 1)):   # a group: side by side, 0.7 m apart, the same pace
                lat = float(np.clip(lat0 + 0.7 * (m - (size - 1) / 2), -half, half))
                walkers.append(dict(fwd=fwd, v0=v0 * rng.uniform(0.97, 1.03), s=-0.3 * m, lat=lat, home=lat,
                                    phase=rng.uniform(0, 2 * np.pi), k0=k, rec=[], pos=None))
            ai += 1
        # positions now, to keep out of each other's way
        for w in walkers:
            sp = w["s"] if w["fwd"] else path.length - w["s"]
            p, h = path.at(min(max(sp, 0.0), path.length))
            w["h"], w["n"] = (h if w["fwd"] else -h), np.array([-h[1], h[0]])   # n: the side lat is measured to, either way
            w["pos"] = p + np.array([-h[1], h[0]]) * w["lat"]
        for w in walkers:
            v, push = w["v0"], 0.0
            for o in walkers:
                if o is w:
                    continue
                d = o["pos"] - w["pos"]
                along, across = d @ w["h"], d @ w["n"]
                ahead = sp_min + (1.0 if o["fwd"] == w["fwd"] else 5.0)   # one coming the other way: step aside early
                if 0 < along < ahead and abs(across) < sp_min:
                    if o["fwd"] == w["fwd"]:   # someone slower ahead, going my way: keep my distance
                        v = min(v, o["v0"] * max(0.0, (along - 0.5 * sp_min) / sp_min))
                    push -= np.sign(across or 1.0) * (sp_min - abs(across))   # and step aside
            w["lat"] += float(np.clip(push, -0.6, 0.6)) * dt + (w["home"] - w["lat"]) * 0.2 * dt
            w["lat"] = float(np.clip(w["lat"], -half, half))
            w["v"] = v
        for w in walkers:
            w["s"] += w["v"] * dt
            sway = 0.04 * np.sin(2 * np.pi * 0.9 * t + w["phase"])   # the body sways at the step rate
            sp = w["s"] if w["fwd"] else path.length - w["s"]
            if 0.0 <= w["s"]:   # a group's followers enter a step behind: born at their first step on it
                if not w["rec"]:
                    w["k0"] = k
                p, h = path.at(min(sp, path.length))
                nrm = np.array([-h[1], h[0]])
                xy = p + nrm * (w["lat"] + sway)
                vel = (h if w["fwd"] else -h) * w["v"]
                w["rec"].append((xy[0], xy[1], vel[0], vel[1]))
        for w in [w for w in walkers if w["s"] > path.length]:
            walkers.remove(w)
            done.append(w)
    return [(w["k0"], np.array(w["rec"])) for w in done + walkers if w["rec"]]


# ---- the whole scenario ---------------------------------------------------------------------

def generate(spec, buildings=(), obstacles=()):
    """Every agent of the spec (people avoid the buildings and the obstacles, footprints both): dict(dt, duration_s, agents=[dict(id, cls, rcs, size, k0, xy)]),
    xy an array (n, 4) of x, y, vx, vy from step k0 on."""
    dt = float(spec.get("dt", 0.05))
    steps = int(round(float(spec.get("duration_s", 120)) / dt))
    rng = np.random.default_rng(int(spec.get("seed", 1)))
    agents = []
    for g in spec.get("groups", []):
        if g.get("kind") == "pedestrians":
            for tr in simulate_people(g, list(buildings) + list(obstacles), rng, steps, dt):
                agents.append(dict(cls="person", rcs=float(g.get("rcs") or CLASSES["person"]["rcs"]), k0=0, xy=tr))
        elif g.get("kind") == "walkway" and len(g.get("points", [])) >= 2:
            for k0, tr in simulate_walkway(g, rng, steps, dt):
                agents.append(dict(cls="person", rcs=float(g.get("rcs") or CLASSES["person"]["rcs"]), k0=k0, xy=tr))
        elif g.get("kind") == "lane" and len(g.get("points", [])) >= 2:
            for c, k0, tr in simulate_lane(g, rng, steps, dt):
                agents.append(dict(cls=c, rcs=CLASSES[c]["rcs"], k0=k0, xy=tr))
    agents.sort(key=lambda a: a["k0"])
    for i, a in enumerate(agents, 1):
        a["id"], a["size"] = i, CLASSES[a["cls"]]["size"]
    return dict(dt=dt, duration_s=steps * dt, agents=agents)


def save(sc, path):
    a = sc["agents"]
    np.savez_compressed(path, dt=sc["dt"], duration_s=sc["duration_s"],
                        cls=np.array([x["cls"] for x in a]), rcs=np.array([x["rcs"] for x in a]),
                        size=np.array([x["size"] for x in a]).reshape(-1, 3),
                        k0=np.array([x["k0"] for x in a], int), n=np.array([len(x["xy"]) for x in a], int),
                        xy=np.concatenate([x["xy"] for x in a]) if a else np.zeros((0, 4)))


class Agents:
    """The saved trajectories, played back: at(t) gives each agent's (pos, vel) or None."""

    def __init__(self, path):
        z = np.load(path)
        self.dt = float(z["dt"])
        self.cls = [str(c) for c in z["cls"]]
        self.rcs = z["rcs"].astype(float)
        self.size = z["size"]
        self.k0 = z["k0"]
        cut = np.r_[0, np.cumsum(z["n"])]
        self.xy = [z["xy"][cut[i]:cut[i + 1]] for i in range(len(self.cls))]
        self.steps = int(round(float(z["duration_s"]) / self.dt))   # a person whose record ends before has left

    def __len__(self):
        return len(self.cls)

    def at(self, t):
        out = []
        for k0, xy, c in zip(self.k0, self.xy, self.cls):
            u = t / self.dt - k0
            if u < 0 or (u > len(xy) - 1 and not (c == "person" and k0 + len(xy) >= self.steps)):
                out.append(None)     # not born yet, or gone
                continue
            u = min(u, len(xy) - 1)  # people stay at their last place
            i = min(int(u), len(xy) - 2) if len(xy) > 1 else 0
            f = u - i
            s = xy[i] if len(xy) == 1 else (1 - f) * xy[i] + f * xy[i + 1]
            out.append((np.array([s[0], s[1], 0.0]), np.array([s[2], s[3], 0.0])))
        return out


def preview(sc, every_s=0.25):
    """Positions every every_s, for animation: dict(dt, frames=n, agents=[dict(id, cls, k0, pts)])."""
    step = max(1, int(round(every_s / sc["dt"])))
    return dict(dt=sc["dt"] * step, duration_s=sc["duration_s"],
                agents=[dict(id=a["id"], cls=a["cls"], k0=int(np.ceil(a["k0"] / step)),
                             pts=np.round(a["xy"][(-a["k0"]) % step::step, :2], 2).tolist()) for a in sc["agents"]])


def main():
    what, spec = sys.argv[1], json.load(open(sys.argv[2]))
    bld = json.load(open(sys.argv[4])).get("buildings", []) if len(sys.argv) > 4 else []
    sc = generate(spec, bld)
    if what == "generate":
        save(sc, sys.argv[3])
        print(f"{len(sc['agents'])} agents over {sc['duration_s']:g} s -> {sys.argv[3]}")
    else:
        json.dump(preview(sc), sys.stdout)


if __name__ == "__main__":
    main()
