"""Line-of-sight map for the tracker's detection probability, from the site file alone (no ray tracer): terrain
(height grid), buildings, walls and other footprint objects (prisms between base and top), trees (cylinders of radius
r between base and top). A grid point at target_z_m above the local ground is seen when the straight lines to the gNB
and to NR-Scope (each at its height above its local ground) are blocked by none of them, and it lies inside the
array's +-fov_deg. Written as the tracker's own cache (pd_grids.npz: key, o0, g0) in a run folder, keyed by the run's
gNB / NR-Scope / boresight. los.grid_file: use a given map (checked against the key) instead.

  python3 los.py [--config ...] RUN_DIR [--compare OTHER.npz]"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C


def site(cfg, origin):
    """everything in the run frame (the importer's east/north metres from the origin)"""
    sys.path.insert(0, os.path.join(cfg["paths"]["paper"], "sim")); import import_capture
    sd = os.path.join(cfg["paths"]["paper"], "sim", "sites", cfg["recording"]["site"])
    meta = json.load(open(os.path.join(sd, "site.json"))); g = json.load(open(os.path.join(sd, "geometry.json")))
    sh = np.asarray(import_capture.enu(np.array([meta["lat"], meta["lon"]]), np.asarray(origin, float)), float)
    prisms = [(np.asarray(o["footprint"], float) + sh, float(o["base"]), float(o["top"]))
              for o in g["buildings"] + [o for o in g["objects"] if "footprint" in o] if len(o["footprint"]) >= 3]
    trees = np.array([[o["x"] + sh[0], o["y"] + sh[1], o["r"], o["base"], o["top"]] for o in g["objects"] if o.get("kind") == "tree"], float).reshape(-1, 5)
    t = g["terrain"]; ticks = np.asarray(t["x"], float); Z = np.asarray(t["z"], float)

    def ground_fn(Zg):   # bilinear on the height grid (site frame ticks, both axes)
        def gnd(x, y):
            st = ticks[1] - ticks[0]
            fx = np.clip((np.asarray(x) - sh[0] - ticks[0]) / st, 0, len(ticks) - 1.001); fy = np.clip((np.asarray(y) - sh[1] - ticks[0]) / st, 0, len(ticks) - 1.001)
            i, j = fy.astype(int), fx.astype(int); a, b = fy - i, fx - j
            return Zg[i, j] * (1 - a) * (1 - b) + Zg[i, j + 1] * (1 - a) * b + Zg[i + 1, j] * a * (1 - b) + Zg[i + 1, j + 1] * a * b
        return gnd
    cands = [ground_fn(Z), ground_fn(Z.T)]   # orientation not stored: keep the one that puts the building bases on the ground
    err = [float(np.mean([abs(f(*p.mean(0)) - b) for p, b, _ in prisms])) for f in cands]
    return prisms, trees, cands[int(np.argmin(err))], min(err)


def blocked(P, end, prisms, trees, ground, n_terrain=40):
    """P (n, 3), end (3,): is the segment P -> end blocked?"""
    out = np.zeros(len(P), bool); d = end[None, :2] - P[:, :2]; dz = end[2] - P[:, 2]
    for poly, base, top in prisms:
        a, e = poly, np.roll(poly, -1, axis=0) - poly
        for k in range(len(a)):
            den = d[:, 0] * e[k, 1] - d[:, 1] * e[k, 0]; ok = np.abs(den) > 1e-9; dd = np.where(ok, den, 1.0); w = a[k] - P[:, :2]
            t = (w[:, 0] * e[k, 1] - w[:, 1] * e[k, 0]) / dd; u = (w[:, 0] * d[:, 1] - w[:, 1] * d[:, 0]) / dd
            z = P[:, 2] + t * dz
            out |= ok & (t > 1e-6) & (t < 1 - 1e-6) & (u >= 0) & (u <= 1) & (z < top) & (z > base)
    L2 = np.maximum((d ** 2).sum(1), 1e-9)
    for x, y, r, base, top in trees:
        t = np.clip(((x - P[:, 0]) * d[:, 0] + (y - P[:, 1]) * d[:, 1]) / L2, 0, 1)
        dist = np.hypot(P[:, 0] + t * d[:, 0] - x, P[:, 1] + t * d[:, 1] - y); z = P[:, 2] + t * dz
        out |= (dist < r) & (t > 1e-3) & (t < 1 - 1e-3) & (z > base) & (z < top)
    for t in np.linspace(0.02, 0.98, n_terrain):
        out |= P[:, 2] + t * dz < ground(P[:, 0] + t * d[:, 0], P[:, 1] + t * d[:, 1])
    return out


def write(cfg, run_dir):
    s = json.load(open(os.path.join(run_dir, "settings.json"))); L = cfg["los"]
    key = np.array([*s["gnb"], *s["rx"], float(s["boresight_deg"])], float)
    if L.get("grid_file"):
        z = np.load(L["grid_file"])
        if not np.allclose(z["key"], key, atol=1e-6):
            raise ValueError("line-of-sight map %s was made for another geometry: key %s vs %s" % (L["grid_file"], np.round(z["key"], 3), np.round(key, 3)))
        np.savez(os.path.join(run_dir, "pd_grids.npz"), key=key, o0=z["o0"], g0=z["g0"]); return z["g0"]
    origin = json.load(open(os.path.join(run_dir, "geometry_origin.json")))["origin"]
    prisms, trees, ground, _ = site(cfg, origin)
    T0, R0 = np.asarray(s["gnb"], float), np.asarray(s["rx"], float); bore = float(s["boresight_deg"])
    T = np.r_[T0[:2], ground(T0[0], T0[1]) + T0[2]]; R = np.r_[R0[:2], ground(R0[0], R0[1]) + R0[2]]
    lo, hi = np.minimum(T0, R0)[:2] - L["span_m"], np.maximum(T0, R0)[:2] + L["span_m"]
    xs, ys = np.arange(lo[0], hi[0], L["step_m"]), np.arange(lo[1], hi[1], L["step_m"]); X, Y = np.meshgrid(xs, ys)
    P = np.c_[X.ravel(), Y.ravel(), ground(X.ravel(), Y.ravel()) + L["target_z_m"]]
    seen = ~blocked(P, T, prisms, trees, ground) & ~blocked(P, R, prisms, trees, ground)
    a = np.degrees(np.arctan2(P[:, 1] - R[1], P[:, 0] - R[0])) - bore
    seen &= np.abs((a + 180) % 360 - 180) <= L["fov_deg"]
    g = seen.reshape(X.shape)
    np.savez(os.path.join(run_dir, "pd_grids.npz"), key=key, o0=np.array([xs[0], ys[0]]), g0=g)
    return g


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("run_dir"); ap.add_argument("--compare"); ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[])
    a = ap.parse_args(); c = C.load(a.config, a.set)
    gerr = site(c, json.load(open(os.path.join(a.run_dir, "geometry_origin.json")))["origin"])[3]
    g = write(c, a.run_dir); print("%s: %d x %d, seen %.1f %% (terrain fit to the building bases: %.2f m)" % (a.run_dir, *g.shape, 100 * g.mean(), gerr))
    if a.compare:
        z = np.load(a.compare); h = z["g0"]
        if h.shape == g.shape:
            print("vs %s: agree on %.1f %% of cells; seen %.1f %% here, %.1f %% there; seen only here %.1f %%, only there %.1f %%"
                  % (a.compare, 100 * (g == h).mean(), 100 * g.mean(), 100 * h.mean(), 100 * (g & ~h).mean(), 100 * (h & ~g).mean()))
