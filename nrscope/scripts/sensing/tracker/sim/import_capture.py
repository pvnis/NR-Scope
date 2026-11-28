#!/usr/bin/env python3
"""Turn a real NR-Scope capture and a ground truth of its targets into a run folder that the
GUI's Results tab (and analyze.py) open like a simulated run.

    import_capture.py CAPTURE_DIR --truth walk1.gpx --truth cars.csv --gnb-h 25 --rx-h 1.5 \
        [--origin LAT,LON] [--gnb LAT,LON --rx LAT,LON --boresight DEG] [--t-offset S] [--name NAME]

CAPTURE_DIR holds what NR-Scope wrote: map2d.csv, map2d.csv.time.csv (the UTC time of every
map) and the nrscope_config.yaml it ran with (or pass --config).

Geometry: the YAML's sensing.geometry (gnb, rx in metres, boresight_deg), or --gnb/--rx in
lat,lon with --origin (then everything is in metres east/north of the origin, and
--boresight is the array normal's bearing counterclockwise from east). Heights above the
ground are not in the YAML: --gnb-h and --rx-h.

Ground truth, one or more files, each either
  GPX    a track of one target (GPS logger app): time, lat, lon; its id is the file name
  CSV    columns: a time (t, time, unix_s or timestamp: UTC seconds or ISO 8601), id
         (optional, else the file name), x,y in metres or lat,lon, cls (optional:
         person, car, ...; default --cls)
Positions are resampled every 0.1 s and the velocity taken from them; a target is absent
where its log has a gap over 2 s. --t-offset (s) is added to the ground truth's times, to
correct a clock offset.

Known path + timed markers (no GPS): --marks marks.csv lists the marks on the ground
(columns mark, x,y or lat,lon), and a truth CSV with a "mark" column instead of positions
gives the UTC time at which the target passed each mark (a tap-logger on the walker's
phone). Between two taps less than MAX_TAP_GAP_S apart the target moves at constant speed
on the straight line between the marks; a longer gap ends the walk.

The ground truth is written as the simulator would: per target one echo at its bistatic
excess range, speed and angle from the geometry (no powers: every target counts with the
same weight). The run's settings carry map_clock = "utc", so analyze.py times the maps
with map2d.csv.time.csv instead of the cell's (frame, slot).
"""
import argparse
import csv
import datetime
import json
import math
import os
import shutil
import sys
import xml.etree.ElementTree as ET

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import track  # noqa: E402

R_EARTH = 6371008.8
DT = 0.1           # ground-truth sample period (s)
MAX_GAP_S = 2.0    # longer gaps in a log: the target is absent
Z = dict(person=1.0, car=0.75)   # echo height above the ground (m), as the simulator's targets


def enu(latlon, origin):
    """lat,lon (deg) to metres east/north of origin (equirectangular, fine over a few km)."""
    lat, lon = np.radians(np.asarray(latlon, float).T)
    lat0, lon0 = np.radians(origin)
    return np.stack([(lon - lon0) * math.cos(lat0) * R_EARTH, (lat - lat0) * R_EARTH], -1)


def parse_time(s):
    try:
        return float(s)
    except ValueError:
        d = datetime.datetime.fromisoformat(s.strip().replace("Z", "+00:00"))
        return (d if d.tzinfo else d.replace(tzinfo=datetime.timezone.utc)).timestamp()


MAX_TAP_GAP_S = 60.0


def read_marks(path, origin):
    marks = {}
    for r in csv.DictReader(open(path)):
        r = {k.strip().lower(): v for k, v in r.items()}
        if "x" in r:
            marks[r["mark"].strip()] = (float(r["x"]), float(r["y"]))
        else:
            if origin is None:
                sys.exit(f"{path}: lat/lon needs --origin")
            marks[r["mark"].strip()] = tuple(enu([(float(r["lat"]), float(r["lon"]))], origin)[0])
    return marks


def densify(a):
    """Taps (t, x, y) -> samples every DT along the straight segments, so resample sees no gap."""
    out = [a[:1]]
    for p, q in zip(a[:-1], a[1:]):
        if q[0] - p[0] <= MAX_TAP_GAP_S:
            f = np.arange(DT, q[0] - p[0], DT)[:, None] / (q[0] - p[0])
            out.append(p + f * (q - p))
        out.append(q[None])
    return np.concatenate(out)


def read_truth_file(path, origin, cls, marks=None):
    """[(id, cls, t (n,), xy (n, 2))]"""
    stem = os.path.splitext(os.path.basename(path))[0]
    if path.lower().endswith(".gpx"):
        pts = []
        for el in ET.parse(path).iter():
            if el.tag.endswith("trkpt") or el.tag.endswith("wpt"):
                tm = [c.text for c in el if c.tag.endswith("time")]
                if tm:
                    pts.append((parse_time(tm[0]), float(el.get("lat")), float(el.get("lon"))))
        if origin is None:
            sys.exit(f"{path}: lat/lon needs --origin")
        pts.sort()
        t = np.array([p[0] for p in pts])
        return [(stem, cls, t, enu([(p[1], p[2]) for p in pts], origin))]
    rows = list(csv.DictReader(open(path)))
    if not rows:
        return []
    keys = {k.strip().lower(): k for k in rows[0]}
    tk = next((keys[k] for k in ("t", "time", "unix_s", "timestamp") if k in keys), None)
    if tk is None:
        sys.exit(f"{path}: no time column (t, time, unix_s or timestamp)")
    out = {}
    for r in rows:
        i = r.get(keys.get("id"), stem) if "id" in keys else stem
        c = r[keys["cls"]] if "cls" in keys else cls
        if "mark" in keys:
            if not marks:
                sys.exit(f"{path}: a mark column needs --marks")
            p = marks[r[keys["mark"]].strip()]
        elif "x" in keys:
            p = (float(r[keys["x"]]), float(r[keys["y"]]))
        elif "lat" in keys:
            if origin is None:
                sys.exit(f"{path}: lat/lon needs --origin")
            p = tuple(enu([(float(r[keys["lat"]]), float(r[keys["lon"]]))], origin)[0])
        else:
            sys.exit(f"{path}: needs x,y or lat,lon columns")
        out.setdefault(i, (c, []))[1].append((parse_time(r[tk]), *p))
    res = []
    for i, (c, pts) in out.items():
        pts.sort()
        a = np.array(pts)
        if "mark" in keys:
            a = densify(a)
        res.append((f"{stem}:{i}" if "id" in keys else stem, c, a[:, 0], a[:, 1:3]))
    return res


def resample(t, xy, grid):
    """Positions and velocities on the grid; NaN outside the log and across gaps."""
    pos = np.full((len(grid), 2), np.nan)
    ok = (grid >= t[0]) & (grid <= t[-1])
    j = np.clip(np.searchsorted(t, grid), 1, len(t) - 1)
    ok &= (t[j] - t[j - 1]) <= MAX_GAP_S
    for k in range(2):
        pos[ok, k] = np.interp(grid[ok], t, xy[:, k])
    # velocity: slope over +-0.5 s (a GPS fix is noisy, a camera box jitters)
    h = int(round(0.5 / DT))
    vel = np.full_like(pos, np.nan)
    for i in range(len(grid)):
        a, b = max(i - h, 0), min(i + h, len(grid) - 1)
        if ok[a] and ok[b] and b > a and ok[i]:
            vel[i] = (pos[b] - pos[a]) / (grid[b] - grid[a])
    return pos, vel


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("capture")
    ap.add_argument("--config", help="NR-Scope YAML (default: the capture's *.yaml)")
    ap.add_argument("--truth", action="append", default=[], help="ground-truth file (GPX or CSV), repeatable")
    ap.add_argument("--gnb-h", type=float, required=True, help="gNB antenna height above the ground (m)")
    ap.add_argument("--rx-h", type=float, required=True, help="NR-Scope array height above the ground (m)")
    ap.add_argument("--origin", help="LAT,LON of the local frame's origin (needed for any lat/lon input)")
    ap.add_argument("--gnb", help="LAT,LON of the gNB (instead of the YAML's geometry)")
    ap.add_argument("--rx", help="LAT,LON of NR-Scope (instead of the YAML's geometry)")
    ap.add_argument("--boresight", type=float, help="array normal, deg counterclockwise from east (instead of the YAML's)")
    ap.add_argument("--marks", help="CSV of marks on the ground (mark, x,y or lat,lon) for tap-logged truth files")
    ap.add_argument("--cls", default="person", help="class of targets whose file has no cls column")
    ap.add_argument("--t-offset", type=float, default=0.0, help="seconds added to the ground truth's times")
    ap.add_argument("--scene", help="site the run is at, e.g. site:epfl_bc (map background, ghost prediction)")
    ap.add_argument("--name", default="capture")
    ap.add_argument("--out", help="output run folder (default logs/<date>_real_<name>)")
    a = ap.parse_args()

    import yaml
    cfgp = a.config or next((os.path.join(a.capture, f) for f in sorted(os.listdir(a.capture)) if f.endswith(".yaml")), None)
    cfg = yaml.safe_load(open(cfgp)) if cfgp else {}
    geo = (cfg.get("sensing") or {}).get("geometry") or {}
    origin = tuple(map(float, a.origin.split(","))) if a.origin else None
    pick = lambda s: tuple(enu([tuple(map(float, s.split(",")))], origin)[0]) if s else None
    if (a.gnb or a.rx) and origin is None:
        sys.exit("--gnb/--rx in lat,lon need --origin")
    gnb2, rx2 = pick(a.gnb) or geo.get("gnb"), pick(a.rx) or geo.get("rx")
    bore = a.boresight if a.boresight is not None else geo.get("boresight_deg")
    if gnb2 is None or rx2 is None or bore is None:
        sys.exit("no geometry: set sensing.geometry in the YAML or pass --gnb, --rx, --boresight")
    gnb = [float(gnb2[0]), float(gnb2[1]), a.gnb_h]
    rx = [float(rx2[0]), float(rx2[1]), a.rx_h]
    n_ant = int((cfg.get("usrp_setting_0") or {}).get("nof_antennas") or 4)

    mp, tp = os.path.join(a.capture, "map2d.csv"), os.path.join(a.capture, "map2d.csv.time.csv")
    if not os.path.exists(tp):
        sys.exit(f"{tp} missing: the capture needs NR-Scope's map timestamps")
    unix = np.array([float(r["unix_s"]) for r in csv.DictReader(open(tp))])
    t0 = float(unix[0]) - 1.0
    import analyze
    maps = analyze.read_maps(mp)
    if len(maps) != len(unix):
        sys.exit(f"{len(maps)} maps but {len(unix)} timestamps")
    lam = maps[0]["lam"]

    targets = []
    marks = read_marks(a.marks, origin) if a.marks else None
    for f in a.truth:
        targets += read_truth_file(f, origin, a.cls, marks)
    if not targets:
        print("no ground truth given: the run will have tracks but no scores")
    grid = np.arange(unix[0] - maps[0]["t_span"] - 1.0, unix[-1] + 1.0, DT)
    T, Rx = np.array(gnb), np.array(rx)
    L = float(np.linalg.norm(Rx - T))
    bo = math.radians(bore)
    ax = np.array([-math.sin(bo), math.cos(bo), 0.0])

    out = a.out or os.path.join(os.path.dirname(HERE), "logs",
                                datetime.datetime.fromtimestamp(unix[0]).strftime("%Y%m%d_%H%M%S") + "_real_" + a.name)
    os.makedirs(out, exist_ok=True)
    shutil.copy(mp, os.path.join(out, "map2d.csv"))
    shutil.copy(tp, os.path.join(out, "map2d.csv.time.csv"))
    if cfgp:
        shutil.copy(cfgp, os.path.join(out, "nrscope_config.yaml"))

    tracks = [(i, c, *resample(t + a.t_offset, xy, grid)) for i, c, t, xy in targets]
    with open(os.path.join(out, "ground_truth.csv"), "w", newline="") as fg, \
            open(os.path.join(out, "ground_truth_cars.csv"), "w", newline="") as fc:
        g, c = csv.writer(fg), csv.writer(fc)
        g.writerow("t_s,car_x,car_y,car_z,car_vx,car_vy,car_vz,path,via_car,delay_s,bistatic_range_m,gain_db,doppler_hz,aoa_deg".split(","))
        c.writerow("t_s,car,x,y,z,vx,vy,vz".split(","))
        for n, tg in enumerate(grid):
            ts = f"{tg - t0:.6f}"
            first = tracks[0][2][n] if tracks else (np.nan, np.nan)
            fv = tracks[0][3][n] if tracks else (np.nan, np.nan)
            head = [ts, *(f"{v:.4f}" for v in (first[0], first[1], 0.0, fv[0], fv[1], 0.0))]
            g.writerow(head + [0, 0, f"{L / 299792458.0:.6e}", f"{L:.3f}", "0.00", "0.000", "0.00"])   # the direct path
            for k, (i, cl, pos, vel) in enumerate(tracks, 1):
                p, v = pos[n], vel[n]
                c.writerow([ts, k, *(f"{x:.4f}" for x in (p[0], p[1], Z.get(cl, 1.0), v[0], v[1], 0.0))])
                if np.isnan(p[0]) or np.isnan(v[0]):
                    continue
                pz = np.array([p[0], p[1], Z.get(cl, 1.0)])
                e, spd, _ = track.bistatic_meas(np.r_[p, v], T, Rx, L, bo, pz[2])
                az = math.degrees(math.asin(np.clip(((pz - Rx) / np.linalg.norm(pz - Rx)) @ ax, -1, 1)))
                g.writerow(head + [k, k, f"{(e + L) / 299792458.0:.6e}", f"{e + L:.3f}", "0.00",
                                   f"{2 * spd / lam:.3f}", f"{az:.2f}"])
    json.dump(dict(cls=[cl for _, cl, _, _ in tracks], ids=[i for i, _, _, _ in tracks]),
              open(os.path.join(out, "agents.json"), "w"))
    settings = dict(gnb=gnb, rx=rx, boresight_deg=float(bore), n_ant=n_ant, map_clock="utc", t0_unix=t0,
                    duration=float(unix[-1] - unix[0]), scenario_name="real: " + a.name, real_capture=True,
                    scene=a.scene, truth_files=[os.path.basename(f) for f in a.truth], t_offset_s=a.t_offset,
                    snr_db=None, traffic="real cell")
    json.dump(settings, open(os.path.join(out, "settings.json"), "w"), indent=1)
    print(f"{out}: {len(maps)} maps over {unix[-1] - unix[0]:.0f} s, {len(tracks)} ground-truth target(s): "
          + ", ".join(f"{i} ({cl})" for i, cl, _, _ in tracks))


if __name__ == "__main__":
    main()
