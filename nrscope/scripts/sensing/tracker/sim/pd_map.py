"""Where each receiver can see a target: line of sight from the gNB and from the receiver,
inside the array's +-FOV_DEG, on a grid over the site (the digital twin). The tracker's
detection probability P_D(x) per receiver: a track behind a building is not expected to be
detected, so its misses there are not evidence against it. Also where wall ghosts of the
targets on a scenario's paths would appear (ghosts), for placing the receivers."""
import os
import numpy as np

FOV_DEG = 60.0     # the array's useful half-width (NR-Scope's mirror and grating limit)
STEP_M = 1.0       # grid step
SPAN_M = 80.0      # the grid reaches this far beyond the gNB and the receiver
TARGET_Z = 1.0     # a person's chest, a car's body


def rays(scene, o, d, maxt, test=True):
    """Mitsuba rays from o along d (unit, (n, 3)): blocked before maxt (test) or the hit distance."""
    import mitsuba as mi
    f = lambda a: mi.Float(np.ascontiguousarray(a, dtype=np.float32))
    ray = mi.Ray3f(mi.Point3f(*(f(o[:, k]) for k in range(3))), mi.Vector3f(*(f(d[:, k]) for k in range(3))))
    ray.maxt = f(maxt)
    if test:
        return np.array(scene.mi_scene.ray_test(ray)).astype(bool).ravel()
    si = scene.mi_scene.ray_intersect(ray)
    return np.where(np.array(si.is_valid()).astype(bool).ravel(), np.array(si.t).ravel(), np.inf)


def ground(scene, xy):
    """The lowest surface under each (x, y), as SionnaChannel._ground (0 where none)."""
    bb = scene.mi_scene.bbox()
    lo, hi = float(bb.min[2]) - 1, float(bb.max[2]) + 1
    xy = np.atleast_2d(xy)
    hit = rays(scene, np.c_[xy, np.full(len(xy), lo)], np.tile([0.0, 0.0, 1.0], (len(xy), 1)), np.full(len(xy), hi - lo), test=False)
    return np.where(np.isfinite(hit), lo + hit, 0.0)


def los_grid(scene, gnb, rx, bore_deg):
    """(x0, y0, step, seen[ny, nx]) for one receiver; scene: a loaded Sionna scene; gnb, rx
    and the targets at their heights above the ground (as the simulator places them)."""
    T, R = (np.r_[p[:2], ground(scene, p[:2])[0] + p[2]] for p in (np.asarray(gnb, float), np.asarray(rx, float)))
    lo, hi = np.minimum(T, R)[:2] - SPAN_M, np.maximum(T, R)[:2] + SPAN_M
    xs, ys = np.arange(lo[0], hi[0], STEP_M), np.arange(lo[1], hi[1], STEP_M)
    X, Y = np.meshgrid(xs, ys)
    xy = np.c_[X.ravel(), Y.ravel()]
    P = np.c_[xy, ground(scene, xy) + TARGET_Z]
    seen = np.ones(len(P), bool)
    for end in (T, R):
        d = end - P
        n = np.linalg.norm(d, axis=1)
        seen &= ~rays(scene, P, d / n[:, None], n - 0.05)
    a = np.degrees(np.arctan2(P[:, 1] - R[1], P[:, 0] - R[0])) - bore_deg
    seen &= np.abs((a + 180) % 360 - 180) <= FOV_DEG
    return float(xs[0]), float(ys[0]), STEP_M, seen.reshape(X.shape)


def site_grids(logdir, scene_name, geos):
    """One grid per receiver geo (gnb, rx, bore_deg), cached in the run's folder; None
    without a site scene."""
    if not str(scene_name).startswith("site:"):
        return None
    path = os.path.join(logdir, "pd_grids.npz")
    key = np.array([v for g in geos for v in (*g[0], *g[1], g[2])], float)
    if os.path.exists(path):
        z = np.load(path)
        if len(z["key"]) == len(key) and np.allclose(z["key"], key):
            return [(*z[f"o{k}"], STEP_M, z[f"g{k}"]) for k in range(len(geos))]
    from sionna_channel import load_base
    scene = load_base(scene_name)
    grids = [los_grid(scene, *g) for g in geos]
    np.savez(path, key=key, **{f"o{k}": np.array(g[:2]) for k, g in enumerate(grids)},
             **{f"g{k}": g[3] for k, g in enumerate(grids)})
    return grids


def seen_at(grid, x, y):
    """The grid's value at (x, y); outside it, not seen."""
    x0, y0, step, g = grid
    i, j = int(round((y - y0) / step)), int(round((x - x0) / step))
    return bool(0 <= i < g.shape[0] and 0 <= j < g.shape[1] and g[i, j])


def ghosts(scene_name, gnb, rxs, paths, step=1.0):
    """Where the single-bounce wall ghosts of targets on the paths appear, per receiver, before
    any run: each path point the gNB lights is a target; its echoes that also bounce on a wall
    (track.ghost_echoes, the image method; wall heights and blocked bounces ignored) are
    localised as the tracker would (track.localize). rxs: [(rx, bore_deg)]; paths: polylines
    [[x, y], ...]. Returns per receiver [[target x, y, ghost x, y, removed]]: removed when the
    ghost lands inside a building (no track starts there) or another receiver sees its place
    directly (it misses there, which kills the track)."""
    import analyze
    import track
    from sionna_channel import load_base
    from scenario import inside
    scene, walls, foot = load_base(scene_name), analyze.site_walls(scene_name), analyze.site_footprints(scene_name)
    grids = [los_grid(scene, gnb, rx, b) for rx, b in rxs]
    pts = []
    for line in paths:
        p = np.asarray(line, float)
        for a, b in zip(p[:-1], p[1:]):
            n = max(1, int(np.linalg.norm(b - a) / step))
            pts += [a + (b - a) * k / n for k in range(n)]
    if not pts or walls is None:
        return [[] for _ in rxs]
    P = np.array(pts)
    P3 = np.c_[P, ground(scene, P) + TARGET_Z]
    Tg = np.r_[gnb[:2], ground(scene, np.asarray(gnb[:2], float))[0] + gnb[2]]
    d = Tg - P3
    n = np.linalg.norm(d, axis=1)
    lit = ~rays(scene, P3, d / n[:, None], n - 0.05)
    out = []
    for k, (rx, b) in enumerate(rxs):
        T, R = track.as3(gnb), track.as3(rx)
        res = []
        for x, y in P[lit]:
            for e, _, th in track.ghost_echoes([x, y, 0.0, 0.0], T, R, float(np.linalg.norm(R - T)), np.deg2rad(b), walls):
                if abs(np.degrees(th)) > FOV_DEG:
                    continue
                g = track.localize(e, np.degrees(th), gnb, rx, b)
                if g is not None:
                    res.append([round(float(x), 1), round(float(y), 1), round(float(g[0]), 1), round(float(g[1]), 1),
                                any(inside(g[:2], f) for f in foot) or any(seen_at(grids[j], g[0], g[1]) for j in range(len(rxs)) if j != k)])
        out.append(res)
    return out


if __name__ == "__main__":   # pd_map.py ghosts spec.json out.json (the GUI's Ghosts button)
    import json
    import sys
    spec = json.load(open(sys.argv[2]))
    rxs = [(spec["rx"], float(spec["boresight_deg"]))] + ([(spec["rx2"], float(spec["boresight2_deg"]))] if spec.get("rx2") else [])
    paths = [g["points"] for g in spec.get("groups", []) if g.get("points")]
    json.dump(dict(ghosts=ghosts(spec["scene"], spec["gnb"], rxs, paths)), open(sys.argv[3], "w"))
