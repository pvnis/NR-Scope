#!/usr/bin/env python3
"""Scene helpers for the GUI.

  scene_tools.py geometry OUT.json
      Building footprints of the street canyon (2D convex hulls, metres) and the floor extent.

  scene_tools.py render OUT.png --gnb x,y,z --rx x,y,z --car x,y,z [--car ...] [--paths] [--view top|street]
      3D picture of the scene rendered by Sionna, with the gNB, NR-Scope and the cars (one --car each);
      --paths also draws the ray-traced propagation paths from the gNB to NR-Scope;
      --cam x,y,z --look x,y,z place the camera (the GUI passes its 3D view's camera).

  scene_tools.py paths OUT.json --gnb x,y,z --rx x,y,z --car x,y,z [--car ...]
      The ray-traced paths as polylines (gNB, interaction points, NR-Scope) with their gain
      in dB and the car each touched (0 = none), for the GUI's interactive 3D view.
"""
import argparse
import json
import os
import sys

import mitsuba as mi
import numpy as np
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sionna_channel import build_scene, pattern_gain, propagation  # noqa: E402

import sionna.rt as rt  # noqa: E402


def geometry(out, scene="simple_street_canyon"):
    """Footprints of the buildings (2D convex hulls of each connected piece, with its height
    above the ground) and the ground's extent. The ground: objects named floor/ground/terrain/
    plane, or flat ones (under 0.5 m of height); its height is their median vertex height."""
    name0 = scene
    scene, cars = build_scene((0, 0, 10), (1, 1, 2), n_cars=0, scene=scene)
    shapes, ground, every = [], [], []
    for name, obj in scene.objects.items():
        v = np.array(obj.mi_mesh.vertex_positions_buffer()).reshape(-1, 3)
        every.append(v)
        if any(k in name.lower() for k in ("floor", "ground", "terrain", "plane")) or np.ptp(v[:, 2]) < 0.5:
            ground.append(v)
            continue
        # One object can hold several separate buildings: split it into connected pieces.
        faces = np.array(obj.mi_mesh.faces_buffer()).reshape(-1, 3)
        parent = list(range(len(v)))
        def root(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for a, b, c in faces:
            parent[root(b)] = root(a)
            parent[root(c)] = root(a)
        groups = {}
        for i in range(len(v)):
            groups.setdefault(root(i), []).append(i)
        for idx in groups.values():
            pts = v[idx]
            if len(pts) < 3 or np.ptp(pts[:, 0]) < 0.5 or np.ptp(pts[:, 1]) < 0.5 or np.ptp(pts[:, 2]) < 1.5:
                continue
            hull = ConvexHull(pts[:, :2])
            shapes.append(dict(name=name, height=float(pts[:, 2].max()), base=float(pts[:, 2].min()),
                               footprint=[[round(float(x), 2), round(float(y), 2)] for x, y in pts[hull.vertices, :2]]))
    g = np.concatenate(ground) if ground else np.concatenate(every)
    floor = [float(g[:, 0].min()), float(g[:, 1].min()), float(g[:, 0].max()), float(g[:, 1].max())]
    json.dump(dict(scene=name0, floor=floor, ground_z=float(np.median(g[:, 2])) if ground else 0.0, buildings=shapes),
              open(out, "w"))


def place_cars(gnb, rx, car_pos, scene="simple_street_canyon"):
    scene, cars = build_scene(gnb, rx, n_cars=len(car_pos), scene=scene)
    for car, p in zip(cars, car_pos):
        car.position = mi.Point3f(*map(float, p))
    return scene, cars


def paths(out, gnb, rx, car_pos, scene="simple_street_canyon"):
    scene, cars = place_cars(gnb, rx, car_pos, scene)
    p = rt.PathSolver()(scene, max_depth=2, diffuse_reflection=True, samples_per_src=200000)
    valid = np.array(p.valid).reshape(-1).astype(bool)
    vert = np.array(p.vertices).reshape(p.vertices.shape[0], -1, 3)
    objs = np.array(p.objects).reshape(p.objects.shape[0], -1)
    a = np.array(p.a[0]).reshape(-1) + 1j * np.array(p.a[1]).reshape(-1)
    ids = {c.object_id: k for k, c in enumerate(cars, 1)}
    lines = []
    for i in np.flatnonzero(valid):
        hops = [vert[d, i].tolist() for d in range(vert.shape[0]) if objs[d, i] != 0xFFFFFFFF]
        car = next((ids[o] for o in objs[:, i] if o in ids), 0)
        lines.append(dict(pts=[list(map(float, gnb))] + hops + [list(map(float, rx))],
                          gain_db=float(20 * np.log10(abs(a[i]) + 1e-30)), car=car))
    json.dump(dict(paths=lines), open(out, "w"))


def rays(out, gnb, rx, agents_file, scene="simple_street_canyon", prop=None, keep=400):
    """The scene's ray-traced paths gNB -> NR-Scope (Sionna, no targets) and every agent's echo
    gNB -> agent -> NR-Scope as the chain makes it: a point at the agent's centre, its gain from
    the bistatic radar equation, hidden when a building blocks either leg. gnb, rx, the agents:
    heights above the ground under them (as the Scenarios page gives them). prop: the
    propagation settings (sionna_channel.propagation). Only the `keep` strongest paths are
    written (with the total count and power)."""
    from scenario import CLASSES
    prop = propagation(prop)
    sc, _ = build_scene((0, 0, 1e3), (1, 1, 1e3), n_cars=0, scene=scene, prop=prop)
    top = float(np.array(sc.mi_scene.bbox().max)[2]) + 10.0
    bot = float(np.array(sc.mi_scene.bbox().min)[2]) - 10.0
    f = lambda a: mi.Float(np.ascontiguousarray(a, dtype=np.float32))
    def cast(o, d, maxt, test=True):
        ray = mi.Ray3f(mi.Point3f(f(o[:, 0]), f(o[:, 1]), f(o[:, 2])), mi.Vector3f(f(d[:, 0]), f(d[:, 1]), f(d[:, 2])))
        ray.maxt = f(maxt)
        if test:
            return np.array(sc.mi_scene.ray_test(ray)).astype(bool).ravel()
        si = sc.mi_scene.ray_intersect(ray)
        return np.where(np.array(si.is_valid()).astype(bool).ravel(), np.array(si.t).ravel(), np.inf)
    def ground(xy):   # the lowest surface under each point, as the chain (sionna_channel._ground)
        xy = np.atleast_2d(np.asarray(xy, float))
        h = cast(np.c_[xy, np.full(len(xy), bot)], np.tile([0.0, 0.0, 1.0], (len(xy), 1)), np.full(len(xy), top - bot), test=False)
        return np.where(np.isfinite(h), bot + h, 0.0)
    ag = json.load(open(agents_file)) if agents_file else []
    g, r = np.array(gnb, float), np.array(rx, float)
    g[2] += ground(g[:2])[0]; r[2] += ground(r[:2])[0]
    sc.transmitters["gnb"].position = mi.Point3f(*map(float, g))
    sc.receivers["nrscope"].position = mi.Point3f(*map(float, r))
    dif = bool(prop["diffraction"])
    p = rt.PathSolver()(sc, max_depth=int(prop["max_depth"]), diffuse_reflection=True, samples_per_src=200000,
                        diffraction=dif, edge_diffraction=dif)
    valid = np.array(p.valid).reshape(-1).astype(bool)
    vert = np.array(p.vertices).reshape(p.vertices.shape[0], -1, 3)
    objs = np.array(p.objects).reshape(p.objects.shape[0], -1)
    a = np.array(p.a[0]).reshape(-1) + 1j * np.array(p.a[1]).reshape(-1)
    inter = np.array(p.interactions).reshape(p.interactions.shape[0], -1) if hasattr(p, "interactions") else None
    KIND = {1: "specular", 2: "diffuse", 4: "refraction", 8: "diffraction"}
    static = []
    idx = np.flatnonzero(valid)
    total_db = float(10 * np.log10((np.abs(a[idx]) ** 2).sum() + 1e-30))
    for i in idx[np.argsort(-np.abs(a[idx]))][:keep]:
        hops = [vert[d, i].tolist() for d in range(vert.shape[0]) if objs[d, i] != 0xFFFFFFFF]
        kinds = sorted({KIND.get(int(inter[d, i]), "?") for d in range(inter.shape[0]) if inter[d, i]}) if inter is not None else []
        static.append(dict(pts=[g.tolist()] + hops + [r.tolist()], gain_db=float(20 * np.log10(abs(a[i]) + 1e-30)),
                           bounces=len(hops), kinds=kinds))
    lam = 3e8 / float(np.array(sc.frequency).ravel()[0])
    echoes = []
    if ag:
        xy = np.array([[x["x"], x["y"]] for x in ag], float)
        h = np.array([1.0 if x["cls"] == "person" else 0.5 * CLASSES.get(x["cls"], CLASSES["car"])["size"][2] for x in ag])
        c = np.c_[xy, ground(xy) + h]
        seen = np.ones(len(ag), bool)
        for end in (g, r):
            d = end - c; n = np.linalg.norm(d, axis=1)
            seen &= ~cast(c, d / n[:, None], n - 0.05)
        for k, x in enumerate(ag):
            lg, lr = np.linalg.norm(c[k] - g), np.linalg.norm(c[k] - r)
            rcs = CLASSES.get(x["cls"], CLASSES["car"])["rcs"]
            amp = lam * np.sqrt(rcs) / ((4 * np.pi) ** 1.5 * lg * lr)
            amp *= np.sqrt(pattern_gain(prop["pattern"], c[k] - g)[0] * pattern_gain(prop["pattern"], c[k] - r)[0])
            echoes.append(dict(id=x.get("id"), cls=x["cls"], pts=[g.tolist(), c[k].tolist(), r.tolist()],
                               gain_db=float(20 * np.log10(amp)), seen=bool(seen[k]),
                               excess_m=float(lg + lr - np.linalg.norm(g - r))))
    json.dump(dict(paths=static, n_paths=int(len(idx)), total_db=total_db, prop=prop, echoes=echoes,
                   gnb=g.tolist(), rx=r.tolist()), open(out, "w"))


def render(out, gnb, rx, car_pos, with_paths, view, cam_pos=None, look_at=None, scene="simple_street_canyon"):
    scene, cars = place_cars(gnb, rx, car_pos, scene)
    paths = None
    if with_paths:
        paths = rt.PathSolver()(scene, max_depth=2, diffuse_reflection=True, samples_per_src=200000)
    centre = (np.array(gnb) + np.array(rx)) / 2
    if cam_pos is not None:
        cam = rt.Camera(position=mi.Point3f(*map(float, cam_pos)), look_at=mi.Point3f(*map(float, look_at)))
    elif view == "top":
        cam = rt.Camera(position=mi.Point3f(float(centre[0]), float(centre[1]) - 1, 220.0),
                        look_at=mi.Point3f(float(centre[0]), float(centre[1]), 0.0))
    else:
        # down the street (it runs along x), from behind the gNB and above the roofs
        cam = rt.Camera(position=mi.Point3f(float(centre[0]) - 110, 0.0, 75.0),
                        look_at=mi.Point3f(float(centre[0]), 0.0, 0.0))
    scene.render_to_file(camera=cam, filename=out, paths=paths, resolution=(960, 640), num_samples=128)


def mesh(out, scene="simple_street_canyon"):
    """The scene's real triangles for the GUI's 3D view: per object, flat vertex (x,y,z, cm
    precision) and triangle index lists, and whether it is ground (same rule as geometry)."""
    scene, cars = build_scene((0, 0, 10), (1, 1, 2), n_cars=0, scene=scene)
    objs = []
    for name, obj in scene.objects.items():
        m = obj.mi_mesh
        v = np.array(m.vertex_positions_buffer()).reshape(-1, 3)
        f = np.array(m.faces_buffer()).astype(int).ravel()
        ground = any(k in name.lower() for k in ("floor", "ground", "terrain", "plane")) or np.ptp(v[:, 2]) < 0.5
        mat = getattr(obj.radio_material, "name", "")
        objs.append(dict(name=name, ground=bool(ground), material=mat, v=np.round(v, 2).ravel().tolist(), f=f.tolist()))
    json.dump(dict(objects=objs), open(out, "w"), separators=(",", ":"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["geometry", "mesh", "render", "paths", "rays"])
    ap.add_argument("--agents", help="rays: a JSON list of agents [{id, cls, x, y}]")
    ap.add_argument("--propagation", default=None, help="rays: propagation settings (JSON, sionna_channel.PROPAGATION)")
    ap.add_argument("out")
    xyz = lambda v: tuple(float(c) for c in v.split(","))
    ap.add_argument("--gnb", type=xyz, default=(-40, 0, 8))
    ap.add_argument("--rx", type=xyz, default=(10, -6, 2))
    ap.add_argument("--car", type=xyz, action="append", help="a car's position; repeat for more cars")
    ap.add_argument("--paths", action="store_true")
    ap.add_argument("--view", choices=["top", "street"], default="street")
    ap.add_argument("--cam", type=xyz, default=None, help="camera position (with --look)")
    ap.add_argument("--look", type=xyz, default=(0, 0, 0), help="point the camera looks at")
    ap.add_argument("--scene", default="simple_street_canyon", help="a Sionna scene name or a Mitsuba XML file")
    a = ap.parse_args()
    if a.what == "geometry":
        geometry(a.out, a.scene)
    elif a.what == "mesh":
        mesh(a.out, a.scene)
    elif a.what == "rays":
        rays(a.out, a.gnb, a.rx, a.agents, a.scene, a.propagation)
    elif a.what == "paths":
        paths(a.out, a.gnb, a.rx, a.car if a.car is not None else [(-20, 3, 0)], a.scene)
    else:
        render(a.out, a.gnb, a.rx, a.car if a.car is not None else [(-20, 3, 0)], a.paths, a.view, a.cam, a.look, a.scene)


if __name__ == "__main__":
    main()
