#!/usr/bin/env python3
"""A real site as a Sionna scene, from OpenStreetMap (buildings, roads) and, in Switzerland,
swisstopo (terrain, aerial photo). No Blender.

  osm_site.py NAME LAT LON [--radius 200] [--level-m 3.5]

writes sim/sites/NAME/:
  scene.xml, meshes/*.ply   the Mitsuba scene Sionna loads (scene "site:NAME")
  geometry.json             building footprints, floor, terrain grid for the GUI
  roads.json                roads, paths and crossings (local metres, with ground height)
  site.json                 origin (lat, lon, ground altitude), settings; edit "heights"
                            ({"BC": 21.0}) or "materials" ({"BC": "glass"}) and re-run to override OSM
                            (keys: a name, or an OSM id like "w123456" / "n789")
  aerial.jpg                swisstopo SWISSIMAGE of the area (with its bounds in site.json)
  osm.json                  the raw OSM answer, kept so a re-run does not download again

Local frame: x east, y north, metres from (LAT, LON); z up, 0 = the ground at the origin.
Re-running an existing NAME rebuilds it in place: the origin stays (so radios and scenarios keep
their place), OSM is downloaded again (your latest edits; Overpass lags a few minutes), and the
edits in site.json are kept. --new-origin moves the origin to LAT LON; --cached reuses osm.json.

From OSM: buildings and building:part (height / min_height, materials from building:material,
building:facade:material or material), trees (natural=tree, tree_row: height, diameter_crown),
artworks (tourism=artwork: sculptures), walls (barrier=wall), street lamps, roads, and any
feature type listed in site.json "features" (e.g. {"tags": {"man_made": "planter"}, "height": 0.8,
"width": 1.5, "material": "concrete"}): add a line there and rebuild, no code change.
Windows: site.json "facades": {"BC": {"window_w": 1.6, "window_h": 2.0, "sill_m": 0.9,
"spacing_m": 2.4, "level_m": 3.5}} cuts each wall of BC into glass panes and wall (FACADE).
Trees OSM lacks go in site.json "extra_trees": [{"name": "diagonal alley", "points": [[x, y], ...], "spacing_m": 8,
"offsets_m": [-4, 4], "height_m": 9, "crown_m": 5, "clear_of": [[[x, y], ...]], "clear_m": 3}]: a row every spacing_m
along the polyline, one per side offset, none within clear_m of the clear_of polylines (roads, paths); or "at": [[x, y], ...] for single trees.
Shapes OSM lacks go in site.json "custom": {"name": "BC roof overhang", "footprint": [[x, y], ...],
"bottom_m": 16.5, "top_m": 17.5, "material": "concrete"} (metres above the ground under it).
OSM buildings that are wrong go in site.json "remove": ["w1567595005", ...] (id or name); a
roof one walks under is then a custom shape.
Building heights: OSM "height", else "building:levels" x level_m, else 10 m. Each building
stands from below the lowest ground under it (or from min_height) to the median ground under it
plus its height.
"""
import argparse
import io
import json
import math
import os
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
SITES = os.path.join(HERE, "sites")
UA = {"User-Agent": "giovanni-isac-sim/0.1 (research; sionna scenes)"}
R_EARTH = 6378137.0
OVERPASS = ["https://overpass-api.de/api/interpreter", "https://overpass.private.coffee/api/interpreter",
            "https://maps.mail.ru/osm/tools/overpass/api/interpreter"]
TERRAIN_STEP_M = 10.0
DEFAULT_HEIGHT_M = 10.0
GROUND_MATERIAL = "medium_dry_ground"
BUILDING_MATERIAL = "concrete"
# OSM material words -> Sionna ITU materials (and our "foliage")
MATERIALS = {"concrete": "concrete", "cement_block": "concrete", "plaster": "concrete", "reinforced_concrete": "concrete",
             "brick": "brick", "glass": "glass", "mirror": "glass", "metal": "metal", "steel": "metal",
             "aluminium": "metal", "aluminum": "metal", "copper": "metal", "bronze": "metal", "iron": "metal",
             "metal_plates": "metal", "wood": "wood", "timber_framing": "wood", "stone": "marble",
             "sandstone": "marble", "granite": "marble", "limestone": "marble", "marble": "marble",
             "plasterboard": "plasterboard"}
ITU = {"concrete", "brick", "plasterboard", "wood", "glass", "ceiling_board", "chipboard", "plywood", "marble",
       "floorboard", "metal", "very_dry_ground", "medium_dry_ground", "wet_ground"}
# foliage is not an ITU material: a lossy, rough low-permittivity slab, about 1 dB per metre of
# crown at 3.5 GHz (ITU-R P.833 order of magnitude). Rough; change it in site.json "foliage".
FOLIAGE = dict(relative_permittivity=2.0, conductivity=0.001, thickness=2.0, scattering_coefficient=0.4)
# slab thickness per material (m): what a wave crossing a wall or a pane goes through
THICKNESS = dict(concrete=0.3, brick=0.25, marble=0.3, plasterboard=0.015, wood=0.05, glass=0.02, metal=0.005,
                 medium_dry_ground=1.0)
# window pattern on a facade (site.json "facades": {"BC": {...}}); missing keys from here
FACADE = dict(window_w=1.6, window_h=2.0, sill_m=0.9, spacing_m=2.4, level_m=3.5, glass="glass")
DEFAULTS = dict(tree_height_m=8.0, tree_crown_m=5.0, tree_row_spacing_m=8.0, artwork_height_m=3.0,
                artwork_width_m=1.5, wall_height_m=1.5, lamp_height_m=5.0)


def get(url, data=None, timeout=90):
    req = urllib.request.Request(url, data=data, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


class Frame:
    """WGS84 <-> local metres around (lat0, lon0) (equirectangular, exact to cm over a km)."""

    def __init__(self, lat0, lon0):
        self.lat0, self.lon0 = lat0, lon0
        self.kx = R_EARTH * math.cos(math.radians(lat0)) * math.pi / 180
        self.ky = R_EARTH * math.pi / 180

    def xy(self, lat, lon):
        return (lon - self.lon0) * self.kx, (lat - self.lat0) * self.ky

    def latlon(self, x, y):
        return self.lat0 + y / self.ky, self.lon0 + x / self.kx


def lv95(lat, lon):
    """WGS84 -> Swiss LV95 (swisstopo's approximate formulas, about 1 m)."""
    p, l = (lat * 3600 - 169028.66) / 10000, (lon * 3600 - 26782.5) / 10000
    e = 2600072.37 + 211455.93 * l - 10938.51 * l * p - 0.36 * l * p * p - 44.54 * l ** 3
    n = 1200147.07 + 308807.95 * p + 3745.25 * l * l + 76.63 * p * p - 194.56 * l * l * p + 119.79 * p ** 3
    return e, n


def in_switzerland(lat, lon):
    return 45.8 < lat < 47.9 and 5.9 < lon < 10.6


def overpass(lat, lon, radius, rules=()):
    a = f"(around:{radius},{lat},{lon})"
    q = (f'[out:json][timeout:90];(way["building"]{a};way["building:part"]{a};'
         f'node["natural"="tree"]{a};way["natural"="tree_row"]{a};node["tourism"="artwork"]{a};'
         f'way["tourism"="artwork"]{a};way["barrier"~"^(wall|retaining_wall)$"]{a};node["highway"="street_lamp"]{a};'
         + "".join(f'nwr{"".join(f"[{json.dumps(k)}={json.dumps(v)}]" for k, v in r["tags"].items())}{a};' for r in rules)
         + f'way["highway"](around:{radius + 50},{lat},{lon}););out body;>;out skel qt;')
    data = urllib.parse.urlencode({"data": q}).encode()
    err = None
    for attempt in range(2):   # the public servers are often busy: ask them all at once, keep the first answer
        pool = ThreadPoolExecutor(len(OVERPASS))
        jobs = {pool.submit(lambda u: json.loads(get(u, data, timeout=60)), u): u for u in OVERPASS}
        try:
            for j in as_completed(jobs):
                try:
                    return j.result()
                except Exception as e:
                    print(f"  {jobs[j]}: {e}", flush=True)
                    err = e
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
    raise RuntimeError(f"OpenStreetMap (Overpass) unreachable: {err}")


def terrain(frame, half):
    """Ground altitude on a grid (swisstopo DTM, one profile request per row), or None."""
    if not in_switzerland(frame.lat0, frame.lon0):
        return None
    xs = np.arange(-half, half + 1e-6, TERRAIN_STEP_M)
    z = np.empty((len(xs), len(xs)))
    for i, y in enumerate(xs):
        a, b = lv95(*frame.latlon(xs[0], y)), lv95(*frame.latlon(xs[-1], y))
        geom = json.dumps({"type": "LineString", "coordinates": [a, b]})
        url = ("https://api3.geo.admin.ch/rest/services/profile.json?sr=2056&nb_points=%d&geom=%s"
               % (len(xs), urllib.parse.quote(geom)))
        prof = json.loads(get(url))
        alt = np.array([p["alts"]["DTM2"] if "DTM2" in p["alts"] else p["alts"]["COMB"] for p in prof], float)
        d = np.array([p["dist"] for p in prof], float)
        z[i] = np.interp(np.linspace(0, d[-1], len(xs)), d, alt)
    return xs, z   # z[iy, ix]


def aerial(frame, half, out, zoom=19):
    """swisstopo SWISSIMAGE tiles (web mercator) over the square, stitched; returns the image's
    bounds in local metres [x0, y0, x1, y1]."""
    def tile(lat, lon):
        n = 2 ** zoom
        return ((lon + 180) / 360 * n,
                (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n)
    def latlon(tx, ty):
        n = 2 ** zoom
        return math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * ty / n)))), tx / n * 360 - 180
    from PIL import Image
    tx0, ty0 = tile(*frame.latlon(-half, half))
    tx1, ty1 = tile(*frame.latlon(half, -half))
    tx0, ty0, tx1, ty1 = int(tx0), int(ty0), int(tx1) + 1, int(ty1) + 1
    img = Image.new("RGB", ((tx1 - tx0) * 256, (ty1 - ty0) * 256))
    for x in range(tx0, tx1):
        for y in range(ty0, ty1):
            url = f"https://wmts.geo.admin.ch/1.0.0/ch.swisstopo.swissimage/default/current/3857/{zoom}/{x}/{y}.jpeg"
            img.paste(Image.open(io.BytesIO(get(url))), ((x - tx0) * 256, (y - ty0) * 256))
    img.save(out, quality=88)
    la0, lo0 = latlon(tx0, ty1)
    la1, lo1 = latlon(tx1, ty0)
    (x0, y0), (x1, y1) = frame.xy(la0, lo0), frame.xy(la1, lo1)
    return [x0, y0, x1, y1]


def triangulate(poly):
    """Ear clipping of a simple counterclockwise polygon [(x, y)]: triangles as index triples."""
    idx = list(range(len(poly)))
    P = np.asarray(poly, float)
    cross = lambda a, b, c: (P[b, 0] - P[a, 0]) * (P[c, 1] - P[a, 1]) - (P[b, 1] - P[a, 1]) * (P[c, 0] - P[a, 0])
    def inside(p, a, b, c):
        return cross(a, b, p) >= 0 and cross(b, c, p) >= 0 and cross(c, a, p) >= 0
    tris, guard = [], 0
    while len(idx) > 3 and guard < 10000:
        guard += 1
        n = len(idx)
        for k in range(n):
            a, b, c = idx[k - 1], idx[k], idx[(k + 1) % n]
            if cross(a, b, c) <= 1e-12:
                continue
            if any(inside(p, a, b, c) for p in idx if p not in (a, b, c)):
                continue
            tris.append((a, b, c))
            idx.pop(k)
            break
        else:   # degenerate leftovers: fan them
            break
    tris += [(idx[0], idx[i], idx[i + 1]) for i in range(1, len(idx) - 1)]
    return tris


def write_ply(path, v, f):
    v = np.asarray(v, "<f4")
    with open(path, "wb") as fh:
        fh.write((f"ply\nformat binary_little_endian 1.0\nelement vertex {len(v)}\n"
                  "property float x\nproperty float y\nproperty float z\n"
                  f"element face {len(f)}\nproperty list uchar int vertex_indices\nend_header\n").encode())
        fh.write(v.tobytes())
        for t in f:
            fh.write(np.uint8(3).tobytes() + np.asarray(t, "<i4").tobytes())


def num(v, default=None):
    """A number from an OSM value ("12", "12 m", "3,5"), else default."""
    try:
        return float(str(v).lower().replace("m", "").replace(",", ".").strip())
    except (TypeError, ValueError):
        return default


def material_of(tags, default):
    for k in ("building:facade:material", "building:material", "material"):
        m = MATERIALS.get(str(tags.get(k, "")).lower())
        if m:
            return m
    return default


def prism(P, base, top):
    """A closed vertical prism over the counterclockwise polygon P: vertices, triangles."""
    m = len(P)
    v = np.r_[np.c_[P, np.full(m, base)], np.c_[P, np.full(m, top)]]
    f = []
    for i in range(m):
        j = (i + 1) % m
        f += [(i, j, m + j), (i, m + j, m + i)]
    f += [(m + a, m + b, m + c) for a, b, c in triangulate(P)]
    f += [(c, b, a) for a, b, c in triangulate(P)]
    return v, f


def ccw(P):
    P = np.asarray(P, float)
    return P[::-1] if np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]) < 0 else P


def circle(x, y, r, n=8):
    a = np.arange(n) * 2 * np.pi / n
    return np.c_[x + r * np.cos(a), y + r * np.sin(a)]


def band(L, w):
    """A polyline L widened to a w-wide strip (a wall), as a polygon."""
    L = np.asarray(L, float)
    d = np.gradient(L, axis=0)
    nrm = np.c_[-d[:, 1], d[:, 0]] / np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-9)[:, None] * w / 2
    return ccw(np.r_[L + nrm, (L - nrm)[::-1]])


def tree(x, y, z0, h, crown):
    """A trunk (prism) and a crown (8-sided double cone between 35% of h and h)."""
    tv, tf = prism(circle(x, y, 0.2, 6), z0 - 0.3, z0 + 0.45 * h)
    r, zb, zt = crown / 2, z0 + 0.35 * h, z0 + h
    ring = np.c_[circle(x, y, r), np.full(8, (zb + zt) / 2)]
    cv = np.r_[ring, [[x, y, zb], [x, y, zt]]]
    cf = [(i, 9, (i + 1) % 8) for i in range(8)] + [((i + 1) % 8, 8, i) for i in range(8)]
    return (tv, tf), (cv, cf)


def facade(P, base, top, zg, spec, wall):
    """A building with windows: walls cut in a grid of panes (spec["glass"]) and wall (`wall`),
    roof and floor of `wall`. Windows window_w x window_h every spacing_m along each wall
    (centred), sill_m above each floor, floors every level_m from the ground zg. Returns
    {material: (vertices, triangles)}."""
    sp = dict(FACADE, **spec)
    out = {}
    def put(mat, quad):
        v, f = out.setdefault(mat, ([], []))
        o = len(v); v.extend(quad); f += [(o, o + 1, o + 2), (o, o + 2, o + 3)]
    rows, k = [], 0
    while True:   # window rows: (bottom, top), the top row at least 0.3 m under the roof
        z0 = zg + k * sp["level_m"] + sp["sill_m"]
        if z0 + sp["window_h"] > top - 0.3:
            break
        if z0 >= base:
            rows.append((z0, z0 + sp["window_h"]))
        k += 1
    for i in range(len(P)):
        a, b = np.asarray(P[i], float), np.asarray(P[(i + 1) % len(P)], float)
        L = float(np.hypot(*(b - a)))
        n = int(L // sp["spacing_m"]) if rows else 0
        u0 = (L - n * sp["spacing_m"]) / 2
        cols = [(u0 + c * sp["spacing_m"] + (sp["spacing_m"] - sp["window_w"]) / 2,) for c in range(n)]
        cols = [(c[0], c[0] + sp["window_w"]) for c in cols]
        ub = sorted({0.0, L, *[x for c in cols for x in c]})
        zb = sorted({base, top, *[z for r in rows for z in r]})
        e = (b - a) / max(L, 1e-9)
        for ua, ubb in zip(ub, ub[1:]):
            in_col = any(c0 <= ua + 1e-6 and ubb <= c1 + 1e-6 for c0, c1 in cols)
            pa, pb = a + e * ua, a + e * ubb
            for za, zt in zip(zb, zb[1:]):
                glass = in_col and any(r0 <= za + 1e-6 and zt <= r1 + 1e-6 for r0, r1 in rows)
                put(sp["glass"] if glass else wall, [[*pa, za], [*pb, za], [*pb, zt], [*pa, zt]])
    v, f = out.setdefault(wall, ([], []))
    m = len(P)
    for z, flip in ((top, False), (base, True)):
        o = len(v); v.extend([[*p, z] for p in P])
        f += [(o + a_, o + c_, o + b_) if flip else (o + a_, o + b_, o + c_) for a_, b_, c_ in triangulate(P)]
    return {k: (np.asarray(vv, float), ff) for k, (vv, ff) in out.items()}


def inside(pt, P):
    x, y = pt; c = False
    for i in range(len(P)):
        (x1, y1), (x2, y2) = P[i - 1], P[i]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            c = not c
    return c


def height_of(tags, level_m):
    for k in ("height", "building:height"):
        try:
            return float(str(tags[k]).replace("m", "").strip())
        except (KeyError, ValueError):
            pass
    try:
        return float(tags["building:levels"]) * level_m + float(tags.get("roof:levels", 0)) * level_m
    except (KeyError, ValueError):
        return DEFAULT_HEIGHT_M


def build(name, lat, lon, radius=None, level_m=3.5, new_origin=False, cached=False):
    out = os.path.join(SITES, name)
    os.makedirs(os.path.join(out, "meshes"), exist_ok=True)
    cfg_path = os.path.join(out, "site.json")
    cfg = json.load(open(cfg_path)) if os.path.exists(cfg_path) else {}
    if cfg.get("lat") is not None and not new_origin:   # rebuild in place: same origin, same local metres
        if abs(cfg["lat"] - lat) > 1e-7 or abs(cfg["lon"] - lon) > 1e-7:
            print(f"existing site: keeping its origin {cfg['lat']:.6f}, {cfg['lon']:.6f}")
        lat, lon = cfg["lat"], cfg["lon"]
    radius = float(radius or cfg.get("radius") or 200.0)   # a rebuild keeps its radius unless given
    cfg.update(name=name, lat=lat, lon=lon, radius=radius, level_m=level_m)
    for k in ("heights", "materials", "facades"):
        cfg.setdefault(k, {})
    cfg["defaults"] = dict(DEFAULTS, **cfg.get("defaults", {}))
    cfg["foliage"] = dict(FOLIAGE, **cfg.get("foliage", {}))
    # your own OSM feature types: every object with these tags becomes a solid (point: a hexagon
    # "width" wide; closed way: its outline; open way: a "width"-wide strip), "height" high
    cfg.setdefault("features", [dict(tags={"man_made": "planter"}, height=0.8, width=1.5, material="concrete"),
                                dict(tags={"amenity": "bench"}, height=0.5, width=0.6, material="wood")])
    cfg.setdefault("custom", [])
    cfg.setdefault("remove", [])
    rule_of = lambda t: next((r for r in cfg["features"] if all(t.get(k) == v for k, v in r["tags"].items())), None)
    D = cfg["defaults"]
    fr = Frame(lat, lon)
    half = radius + 30

    raw_path = os.path.join(out, "osm.json")
    if cached and os.path.exists(raw_path) and json.load(open(raw_path)).get("_q") == [lat, lon, radius]:
        osm = json.load(open(raw_path))
    else:
        print("downloading OpenStreetMap...")
        try:
            osm = overpass(lat, lon, radius, cfg["features"])
            osm["_q"] = [lat, lon, radius]
            json.dump(osm, open(raw_path, "w"))
        except RuntimeError as e:
            if not os.path.exists(raw_path):
                raise
            print(f"{e}\nusing the last download (osm.json); your latest OSM edits are not in it", flush=True)
            osm = json.load(open(raw_path))
    nodes = {e["id"]: fr.xy(e["lat"], e["lon"]) for e in osm["elements"] if e["type"] == "node"}
    ways = [e for e in osm["elements"] if e["type"] == "way"]
    points = [e for e in osm["elements"] if e["type"] == "node" and e.get("tags")]

    print("terrain...")
    ter = terrain(fr, half)
    if ter is None:
        alt0, ground = 0.0, (lambda x, y: np.zeros(np.shape(x)))
        gx = np.array([-half, half]); gz = np.zeros((2, 2))
    else:
        gx, alt = ter
        f = RegularGridInterpolator((gx, gx), alt, bounds_error=False, fill_value=None)
        alt0 = float(f([[0.0, 0.0]])[0])
        gz = alt - alt0
        g2 = RegularGridInterpolator((gx, gx), gz, bounds_error=False, fill_value=None)
        ground = lambda x, y: g2(np.c_[np.ravel(y), np.ravel(x)]).reshape(np.shape(x))
    cfg["ground_altitude_m"] = alt0

    for fn in os.listdir(os.path.join(out, "meshes")):   # objects deleted from OSM must go
        os.remove(os.path.join(out, "meshes", fn))
    # meshes per material (Sionna merges them by material anyway): {material: [(v, f), ...]}
    parts = {}
    add = lambda mat, vf: parts.setdefault(mat, []).append(vf)
    X, Y = np.meshgrid(gx, gx)
    n = len(gx)
    add(GROUND_MATERIAL, (np.c_[X.ravel(), Y.ravel(), gz.ravel()],
                          [(i * n + j, i * n + j + 1, (i + 1) * n + j + 1) for i in range(n - 1) for j in range(n - 1)] +
                          [(i * n + j, (i + 1) * n + j + 1, (i + 1) * n + j) for i in range(n - 1) for j in range(n - 1)]))
    near = lambda P: np.hypot(*np.asarray(P, float).reshape(-1, 2).mean(0)) <= radius
    key = lambda e: f"{e['type'][0]}{e['id']}"
    def pick(table, e, label, default):   # site.json override by OSM id, then name, else OSM
        return table.get(key(e), table.get(label, default))
    poly = lambda w: [nodes[i] for i in w["nodes"][:-1] if i in nodes] if w["nodes"][0] == w["nodes"][-1] else None

    foot, objs = [], []
    part_ways = [w for w in ways if "building:part" in w.get("tags", {}) and poly(w) and len(poly(w)) >= 3]
    part_polys = [ccw(poly(w)) for w in part_ways]
    solids = []   # (way, polygon, kind) extruded the same way
    for w in ways:
        t = w.get("tags", {})
        if "building" in t and poly(w) and len(poly(w)) >= 3:
            solids.append((w, ccw(poly(w)), "building"))
    solids += [(w, P, "part") for w, P in zip(part_ways, part_polys)]
    for w, P, kind in solids:
        t = w["tags"]
        label = t.get("name") or t.get("ref") or key(w)
        if not near(P) or key(w) in cfg["remove"] or label in cfg["remove"]:
            continue
        h = float(pick(cfg["heights"], w, label, height_of(t, level_m)))
        min_h = num(t.get("min_height"), (num(t.get("building:min_level"), 0) or 0) * level_m)
        mat = pick(cfg["materials"], w, label, material_of(t, BUILDING_MATERIAL))
        gpts = ground(P[:, 0], P[:, 1])
        base = float(np.median(gpts)) + min_h if min_h else float(gpts.min()) - 0.5
        top = float(np.median(gpts)) + h
        # an outline whose 3D shape is described by building:part pieces is drawn by those pieces
        has_parts = kind == "building" and any(inside(Q.mean(0), P) for Q in part_polys)
        fac = cfg["facades"].get(key(w), cfg["facades"].get(label))
        if not has_parts and top > base:
            if fac:
                for m_, vf in facade(P, base, top, float(np.median(gpts)), fac, mat).items():
                    add(m_, vf)
            else:
                add(mat, prism(P, base, top))
        if kind == "building":
            foot.append(dict(name=label, osm=key(w), height=h, base=base, top=top, material=mat,
                             footprint=P.round(2).tolist(), levels=t.get("building:levels"), osm_height=t.get("height")))
        else:
            objs.append(dict(kind="part", name=label, osm=key(w), base=base, top=top, material=mat,
                             footprint=P.round(2).tolist()))

    def put_tree(e, x, y, t, label):
        h = float(pick(cfg["heights"], e, label, num(t.get("height"), D["tree_height_m"])))
        c = num(t.get("diameter_crown"), D["tree_crown_m"])
        z0 = float(ground(np.array([x]), np.array([y]))[0])
        trunk, crown = tree(x, y, z0, h, c)
        add("wood", trunk); add("foliage", crown)
        objs.append(dict(kind="tree", name=t.get("species") or t.get("genus") or "tree", osm=key(e),
                         x=round(x, 2), y=round(y, 2), r=c / 2, base=z0, top=z0 + h, material="foliage"))

    def put_rule(e, r, t, label, P):
        h = float(pick(cfg["heights"], e, label, num(t.get("height"), r.get("height", 1.0))))
        mat = pick(cfg["materials"], e, label, material_of(t, r.get("material", "concrete")))
        gp = ground(P[:, 0], P[:, 1]); z0 = float(np.median(gp))
        add(mat, prism(P, float(gp.min()) - 0.3, z0 + h))
        kind = "/".join(f"{k}={v}" for k, v in r["tags"].items())
        objs.append(dict(kind=kind, name=label, osm=key(e), footprint=P.round(2).tolist(), base=z0, top=z0 + h, material=mat))

    for e in points:
        t = e["tags"]; x, y = nodes[e["id"]]
        if np.hypot(x, y) > radius:
            continue
        label = t.get("name") or key(e)
        z0 = float(ground(np.array([x]), np.array([y]))[0])
        if rule_of(t):
            put_rule(e, rule_of(t), t, label, ccw(circle(x, y, num(t.get("width"), rule_of(t).get("width", 1.0)) / 2, 6)))
        elif t.get("natural") == "tree":
            put_tree(e, x, y, t, label)
        elif t.get("tourism") == "artwork":
            h = float(pick(cfg["heights"], e, label, num(t.get("height"), D["artwork_height_m"])))
            wd = num(t.get("width"), D["artwork_width_m"])
            mat = pick(cfg["materials"], e, label, material_of(t, "concrete"))
            P = ccw(circle(x, y, wd / 2, 6))
            add(mat, prism(P, z0 - 0.3, z0 + h))
            objs.append(dict(kind="artwork", name=label, osm=key(e), footprint=P.round(2).tolist(),
                             base=z0, top=z0 + h, material=mat))
        elif t.get("highway") == "street_lamp":
            h = num(t.get("height"), D["lamp_height_m"])
            add("metal", prism(ccw(circle(x, y, 0.08, 6)), z0 - 0.3, z0 + h))
            objs.append(dict(kind="lamp", name="lamp", osm=key(e), x=round(x, 2), y=round(y, 2), r=0.3,
                             base=z0, top=z0 + h, material="metal"))

    for w in ways:
        t = w.get("tags", {})
        L = np.array([nodes[i] for i in w["nodes"] if i in nodes])
        if len(L) < 2 or not near(L):
            continue
        label = t.get("name") or key(w)
        if rule_of(t) and "building" not in t:
            r = rule_of(t)
            put_rule(w, r, t, label, ccw(poly(w)) if poly(w) and len(poly(w)) >= 3 else band(L, num(t.get("width"), r.get("width", 0.5))))
        elif t.get("natural") == "tree_row":
            seg = np.hypot(*np.diff(L, axis=0).T); s_ = np.r_[0, np.cumsum(seg)]
            for d in np.arange(0, s_[-1] + 1e-6, D["tree_row_spacing_m"]):
                put_tree(w, float(np.interp(d, s_, L[:, 0])), float(np.interp(d, s_, L[:, 1])), t, label)
        elif t.get("tourism") == "artwork" or t.get("barrier") in ("wall", "retaining_wall"):
            art = t.get("tourism") == "artwork"
            h = float(pick(cfg["heights"], w, label, num(t.get("height"), D["artwork_height_m"] if art else D["wall_height_m"])))
            mat = pick(cfg["materials"], w, label, material_of(t, "concrete"))
            P = ccw(poly(w)) if poly(w) and len(poly(w)) >= 3 else band(L, num(t.get("width"), 0.3))
            gp = ground(P[:, 0], P[:, 1]); z0 = float(np.median(gp))
            add(mat, prism(P, float(gp.min()) - 0.3, z0 + h))
            objs.append(dict(kind="artwork" if art else "wall", name=label, osm=key(w), footprint=P.round(2).tolist(),
                             base=z0, top=z0 + h, material=mat))

    # your own trees (site.json "extra_trees"): rows along a polyline (local metres), one per side offset
    def _dist_to(p, Q):
        best = 1e9
        for a, b in zip(Q[:-1], Q[1:]):
            d = b - a; t = min(max(float((p - a) @ d / max(d @ d, 1e-9)), 0.0), 1.0); best = min(best, float(np.hypot(*(p - (a + t * d)))))
        return best
    for r in cfg.get("extra_trees", []):
        h, c = float(r.get("height_m", D["tree_height_m"])), float(r.get("crown_m", D["tree_crown_m"]))
        for x, y in r.get("at", []):   # single trees at given positions
            z0 = float(ground(np.array([x]), np.array([y]))[0])
            trunk, crown = tree(x, y, z0, h, c)
            add("wood", trunk); add("foliage", crown)
            objs.append(dict(kind="tree", name=r.get("name", "extra tree"), osm="extra", x=round(x, 2), y=round(y, 2),
                             r=c / 2, base=z0, top=z0 + h, material="foliage"))
        if len(r.get("points", [])) < 2:
            continue
        L = np.array(r["points"], float); seg = np.hypot(*np.diff(L, axis=0).T); s_ = np.r_[0, np.cumsum(seg)]
        for d in np.arange(float(r.get("start_m", 0.0)), s_[-1] + 1e-6, float(r.get("spacing_m", D["tree_row_spacing_m"]))):
            k = min(int(np.searchsorted(s_, d, side="right")) - 1, len(seg) - 1)
            u = (L[k + 1] - L[k]) / max(seg[k], 1e-9); nrm = np.array([-u[1], u[0]])
            p0 = np.array([np.interp(d, s_, L[:, 0]), np.interp(d, s_, L[:, 1])])
            for off in r.get("offsets_m", [0.0]):
                x, y = (p0 + float(off) * nrm).tolist()
                if any(_dist_to(np.array([x, y]), np.array(q, float)) < float(r.get("clear_m", 3.0)) for q in r.get("clear_of", [])):
                    continue   # keep roads and paths free (clear_of: polylines)
                z0 = float(ground(np.array([x]), np.array([y]))[0])
                trunk, crown = tree(x, y, z0, h, c)
                add("wood", trunk); add("foliage", crown)
                objs.append(dict(kind="tree", name=r.get("name", "extra tree"), osm="extra", x=round(x, 2), y=round(y, 2),
                                 r=c / 2, base=z0, top=z0 + h, material="foliage"))

    # your own shapes (site.json "custom"): a footprint in local metres (x east, y north, as the
    # GUI's cursor shows), from bottom_m to top_m above the ground under it. A roof overhang is a
    # thin slab high up; a missing kiosk is a box from 0.
    for c in cfg["custom"]:
        P = ccw(c["footprint"])
        gp = ground(P[:, 0], P[:, 1]); z0 = float(np.median(gp))
        bot = z0 + c.get("bottom_m", 0.0) if c.get("bottom_m", 0.0) > 0 else float(gp.min()) - 0.3
        add(c.get("material", "concrete"), prism(P, bot, z0 + c["top_m"]))
        objs.append(dict(kind="custom", name=c.get("name", "custom"), footprint=P.round(2).tolist(),
                         clear_m=float(c.get("bottom_m", 0.0)),   # its underside above the ground
                         base=bot, top=z0 + c["top_m"], material=c.get("material", "concrete")))

    shapes = []
    for mat, lst in parts.items():
        V, F, o = [], [], 0
        for v, f in lst:
            V.append(np.asarray(v)); F += [tuple(np.asarray(t) + o) for t in f]; o += len(v)
        fname = "ground" if mat == GROUND_MATERIAL else f"m_{mat}"
        write_ply(os.path.join(out, "meshes", fname + ".ply"), np.concatenate(V), F)
        shapes.append((fname, mat))

    fo = cfg["foliage"]
    def bsdf(m):
        if m in ITU:
            return (f'\t<bsdf type="itu-radio-material" id="{m}">\n\t\t<string name="type" value="{m}"/>\n'
                    f'\t\t<float name="thickness" value="{THICKNESS.get(m, 0.1)}"/>\n\t</bsdf>')
        return (f'\t<bsdf type="radio-material" id="{m}">\n' +
                "".join(f'\t\t<float name="{k}" value="{v}"/>\n' for k, v in fo.items()) + '\t</bsdf>')
    xml = ['<scene version="2.1.0">'] + [bsdf(m) for m in sorted({m for _, m in shapes})]
    xml += [f'\t<shape type="ply" id="mesh-{s}">\n\t\t<string name="filename" value="meshes/{s}.ply"/>\n'
            f'\t\t<boolean name="face_normals" value="true"/>\n\t\t<ref id="{m}" name="bsdf"/>\n\t</shape>' for s, m in shapes]
    xml.append("</scene>")
    open(os.path.join(out, "scene.xml"), "w").write("\n".join(xml) + "\n")

    roads = []
    for w in ways:
        t = w.get("tags", {})
        if "highway" not in t:
            continue
        P = np.array([nodes[i] for i in w["nodes"] if i in nodes])
        if len(P) < 2:
            continue
        z = ground(P[:, 0], P[:, 1])
        roads.append(dict(id=w["id"], highway=t["highway"], name=t.get("name"), lanes=t.get("lanes"),
                          oneway=t.get("oneway"), footway=t.get("footway"), width=t.get("width"),
                          pts=np.c_[P, z].round(2).tolist()))
    json.dump(dict(roads=roads), open(os.path.join(out, "roads.json"), "w"))

    if cached and os.path.exists(os.path.join(out, "aerial.jpg")) and cfg.get("aerial_bounds"):
        pass
    else:
        try:
            print("aerial photo...")
            cfg["aerial_bounds"] = aerial(fr, half, os.path.join(out, "aerial.jpg")) if in_switzerland(lat, lon) else None
        except Exception as e:   # the scene works without it
            print("no aerial photo:", e)
            cfg["aerial_bounds"] = None
    json.dump(dict(scene=f"site:{name}", floor=[-half, -half, half, half], ground_z=0.0, buildings=foot, objects=objs,
                   terrain=dict(x=gx.tolist(), z=np.round(gz, 2).tolist())),
              open(os.path.join(out, "geometry.json"), "w"))
    json.dump(cfg, open(cfg_path, "w"), indent=1)
    cnt = {}
    for o in objs:
        cnt[o["kind"]] = cnt.get(o["kind"], 0) + 1
    extra = ", ".join(f"{v} {k}{'s' if v > 1 else ''}" for k, v in sorted(cnt.items()))
    print(f"{name}: {len(foot)} buildings{', ' + extra if extra else ''}, {len(roads)} roads/paths, "
          f"ground {gz.min():+.1f} to {gz.max():+.1f} m around {alt0:.1f} m a.s.l. -> {out}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name")
    ap.add_argument("lat", type=float)
    ap.add_argument("lon", type=float)
    ap.add_argument("--radius", type=float, default=None, help="metres (default 200, or the existing site's)")
    ap.add_argument("--level-m", type=float, default=3.5, help="height of one floor when OSM gives only levels")
    ap.add_argument("--new-origin", action="store_true", help="move an existing site's origin to LAT LON")
    ap.add_argument("--cached", action="store_true", help="reuse the downloaded OSM and aerial photo")
    a = ap.parse_args()
    build(a.name, a.lat, a.lon, a.radius, a.level_m, a.new_origin, a.cached)
    print(end="", flush=True)
    os._exit(0)   # don't wait for the slower Overpass servers still answering


if __name__ == "__main__":
    main()
