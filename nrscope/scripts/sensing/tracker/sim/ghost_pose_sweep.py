"""How exactly the wall map must be registered to the receiver for the tracker's multipath
ghost prediction (ghost_map=1) to help: re-tracks a run's maps with the site's walls moved
rigidly about NR-Scope by a heading error (deg) and a position error (m), and prints the
scores next to no walls (ghost_map=0), walls learned from the run (2) and the exact walls.

The error is applied to the walls, not to the receiver, so localisation itself stays exact
and only the ghost prediction sees the error (what a wrong lidar/OSM fit does). Each
nonzero error is averaged over its signs (heading) or four directions (position).

Usage: python3 ghost_pose_sweep.py <run dir> [--deg 0,1,2,5,10] [--m 0,0.5,1,2,5] [--json out.json]
"""
import argparse
import json
import sys

import numpy as np

import analyze

KEYS = ("false_per_map", "missed_per_map", "gospa_m", "false_tracks", "n_tracks")


def moved(walls, rx, deg, dxy):
    """walls rotated by deg about rx (x-y) and shifted by dxy."""
    a, b, z0, z1 = walls
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    rot = np.array([[c, -s], [s, c]])
    f = lambda p: (p - rx) @ rot.T + rx + dxy
    return f(a), f(b), z0, z1


def score(run, ghost_map, walls=None):
    exact = analyze.site_walls
    if walls is not None:
        analyze.site_walls = lambda scene: walls
    try:
        s = analyze.summary(run, ekf={"ghost_map": ghost_map})
    finally:
        analyze.site_walls = exact
    p = s.get("position") or {}
    out = {k: float(p[k]) for k in KEYS if k in p}
    for c in (s.get("per_class") or {}).get("position") or []:
        out[c["cls"] + "_coverage"] = float(c["coverage"])
    return out


def mean(rows):
    return {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--deg", default="0,1,2,5,10")
    ap.add_argument("--m", default="0,0.5,1,2,5")
    ap.add_argument("--json")
    a = ap.parse_args()
    st = json.load(open(f"{a.run}/settings.json"))
    scene = st.get("scene") or (st.get("scenario") or {}).get("scene")
    walls = analyze.site_walls(scene)
    if walls is None:
        sys.exit(f"no site walls for scene {scene!r}")
    rx = np.asarray(st["rx"][:2], float)

    res = [dict(case="no walls", **score(a.run, 0)), dict(case="learned walls", **score(a.run, 2))]
    for deg in map(float, a.deg.split(",")):
        for m in map(float, a.m.split(",")):
            signs = [1, -1] if deg else [0]
            dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)] if m else [(0, 0)]
            rows = [score(a.run, 1, moved(walls, rx, sg * deg, m * np.array(d, float)))
                    for sg in signs for d in dirs]
            res.append(dict(case=f"map {deg:g} deg {m:g} m", deg=deg, m=m, **mean(rows)))
            print(f"  {res[-1]['case']} done", file=sys.stderr)

    cols = [k for k in res[0] if k != "case"]
    print(f"{'case':22s}" + "".join(f"{k:>18s}" for k in cols))
    for r in res:
        print(f"{r['case']:22s}" + "".join(f"{r.get(k, float('nan')):18.3f}" for k in cols))
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
