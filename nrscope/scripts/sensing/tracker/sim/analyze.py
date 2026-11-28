#!/usr/bin/env python3
"""Score NR-Scope's sensing maps against the Sionna ground truth of one run.

    analyze.py LOGDIR            prints the summary; the GUI imports summary() and one_map()

Inputs, both written by a run into LOGDIR:
  map2d.csv          NR-Scope range-Doppler maps (format: NR-Scope/scripts/sensing/plot_range_doppler.py)
  ground_truth.csv   every propagation path at every ray trace (sim/sionna_channel.py)
  ground_truth_cars.csv  every car's position at every ray trace (runs with one car may lack it)

Timing. The ZMQ chain is sample-synchronous and the gNB starts at SFN 0 on the bridge's first
sample, so a map stamped frame.slot ends at frame*10 ms + slot*0.5 ms of signal (unwrapped
over the 10.24 s SFN period). The map covers the t_span_s before that, and is compared with
the ground truth at the middle of its window.

Quantities, in NR-Scope's own convention:
  excess range  path length of the echo minus that of the direct path (m)
  speed         Doppler * lambda / 2 (m/s); a target moving away has negative Doppler

The truth for a map is one echo per car: the power-weighted range and speed of all the paths
that touched that car. A car counts as detected in a map when an estimate falls within
`gate_bins` range bins and `gate_bins` speed cells of its echo. Two estimators are scored:
  detector  NR-Scope's own targets (markers with verdict 0); only with tdd_detect markers.
            Scored per car and map; a marker near no car is a false target.
  peak      the strongest map cell away from zero Doppler (|speed| > guard): a model-free
            check, a hit when it lies on any car in view (one peak cannot find several cars)

Tracking (track.py): the detector's markers linked from map to map, scored per car
(coverage, identity switches, time to confirm) and with GOSPA.

Position (runs with several antennas and the geometry in settings.json): every AoA cell
localised in 2D as NR-Scope does (track.localize), scored against the car's centre, then
tracked in x-y (track.localize_and_score).

Two receivers (settings rx2, boresight2_deg): the second NR-Scope's maps are map2d_rx2.csv and
its truth ground_truth_rx2.csv. summary(rx=2) scores it alone, summary(rx="fused") tracks with
both receivers' maps at once (track.localize_and_score's detection-level fusion).

Angle (runs with several antennas): NR-Scope's AoA cells (the map's aoa fields) matched to
a car like the detections; the truth is the power-weighted cone angle of the car's paths
(ground_truth.csv aoa_deg), asin(sin(azimuth) cos(elevation)), what the linear array measures.
"""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import track  # noqa: E402
import pd_map  # noqa: E402

C = 299792458.0
SFN_PERIOD_S = 10.24
META = 21  # fields before the power payload (see the map header)


def read_maps(path):
    maps, wrap, last = [], 0.0, None
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.rstrip().split(",")
            frame, slot = int(v[0]), int(v[1])
            nb, nf, nm, na = int(v[12]), int(v[13]), int(v[14]), int(v[15])
            t = frame * 0.01 + slot * 0.0005
            if last is not None and t + wrap < last - 1.0:  # SFN wrapped
                wrap += SFN_PERIOD_S
            t += wrap
            last = t
            p0 = META
            power = np.array(v[p0:p0 + nb * nf], dtype=float).reshape(nb, nf)
            mk = np.array(v[p0 + nb * nf:p0 + nb * nf + 4 * nm], dtype=float).reshape(-1, 4)
            a0 = p0 + nb * nf + 4 * nm
            aoa = np.array(v[a0:a0 + 5 * na], dtype=float).reshape(-1, 5)  # range_m, speed, angle, power, quality
            fc = float(v[11])
            maps.append(dict(frame=frame, slot=slot, aarx=int(v[2]), ports=int(v[3]), layer=int(v[4]),
                             t_end=t, t_span=float(v[8]), m_per_bin=float(v[9]), f_max=float(v[10]),
                             carrier=fc, lam=C / fc, n_bins=nb, n_freq=nf, bin_los=float(v[16]),
                             noise_ref=float(v[20]), k_step=int(v[5]), n_snap=int(v[7]),
                             pilots=(int(v[17]), int(v[18])), power=power, markers=mk, aoa=aoa))
    for m in maps:
        m["t"] = m["t_end"] - m["t_span"] / 2
        m["speeds"] = (-m["f_max"] + np.arange(m["n_freq"]) * 2 * m["f_max"] / (m["n_freq"] - 1)) * m["lam"] / 2
        m["excess"] = (np.arange(m["n_bins"]) - m["bin_los"]) * m["m_per_bin"]
        m["speed_step"] = m["speeds"][1] - m["speeds"][0]
    return maps


def read_truth(path):
    """{t: dict(car=(x, y) of car 1, cars=[(x, y) per car], los_m=direct path length,
    paths=[(excess_m, doppler_hz, gain_db, car)])}, car = 0 for no car, k for car k"""
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.setdefault(float(r["t_s"]), []).append(r)
    cars = {}
    cars_path = path.rsplit(".", 1)[0] + "_cars.csv"
    if os.path.exists(cars_path):
        with open(cars_path) as f:
            for r in csv.DictReader(f):
                cars.setdefault(float(r["t_s"]), []).append((float(r["x"]), float(r["y"])))
    truth = {}
    for t, rs in rows.items():
        los = min(float(r["bistatic_range_m"]) for r in rs)
        paths = []
        for r in rs:
            paths.append((float(r["bistatic_range_m"]) - los, float(r["doppler_hz"]), float(r["gain_db"]),
                          int(r["via_car"]), float(r.get("aoa_deg") or "nan")))
        car1 = (float(rs[0]["car_x"]), float(rs[0]["car_y"]))
        truth[t] = dict(car=car1, cars=cars.get(t, [car1]), los_m=los, paths=paths)
    return truth


def truth_at(truth, times, t, lam):
    """Ground truth at signal time t in map units: one echo per car (targets), and the
    strongest of them at the top level."""
    i = int(np.clip(np.searchsorted(times, t), 0, len(times) - 1))
    if i > 0 and abs(times[i - 1] - t) < abs(times[i] - t):
        i -= 1
    g = truth[times[i]]
    allcar = [p for p in g["paths"] if p[3]]
    out = dict(t_gt=times[i], car=g["car"], cars=g["cars"], los_m=g["los_m"],
               car_paths=[(p[0], p[1] * lam / 2, p[2], p[3]) for p in allcar], targets=[])
    for k in sorted({p[3] for p in allcar}):
        car = [p for p in allcar if p[3] == k]
        # The echo as a whole: its power-weighted range and speed. The strongest single path is
        # kept too, but it can jump between parts of the car from one trace to the next.
        w = np.array([10 ** (p[2] / 10) for p in car])
        best = max(car, key=lambda p: p[2])
        ang = np.deg2rad([p[4] for p in car])
        out["targets"].append(dict(car=k, excess=float(np.dot(w, [p[0] for p in car]) / w.sum()),
                                   aoa=float(np.rad2deg(np.angle(np.dot(w, np.exp(1j * ang))))),
                                   speed=float(np.dot(w, [p[1] for p in car]) / w.sum()) * lam / 2,
                                   gain_db=float(10 * np.log10(w.sum())),
                                   strongest_excess=best[0], strongest_speed=best[1] * lam / 2))
    if out["targets"]:
        out.update({k: v for k, v in max(out["targets"], key=lambda g: g["gain_db"]).items() if k != "car"})
    return out


def estimates(m, guard_ms):
    det = [(r[0] - m["bin_los"] * m["m_per_bin"], r[1], r[2]) for r in m["markers"] if int(r[3]) == 0]
    p = m["power"].copy()
    p[:, np.abs(m["speeds"]) <= guard_ms] = -np.inf
    b, k = np.unravel_index(np.argmax(p), p.shape)
    peak_db = 10 * np.log10(max(m["power"][b, k], 1e-30) / max(np.median(m["power"]), 1e-30))
    return det, (m["excess"][b], m["speeds"][k], peak_db)


def load(logdir, suffix=""):
    maps = read_maps(os.path.join(logdir, f"map2d{suffix}.csv"))
    st = os.path.join(logdir, "settings.json")
    st = json.load(open(st)) if os.path.exists(st) else {}
    if st.get("map_clock") == "utc":   # a real capture (import_capture.py): the maps' UTC times, not the cell's (frame, slot)
        with open(os.path.join(logdir, "map2d.csv.time.csv")) as f:
            unix = [float(r["unix_s"]) for r in csv.DictReader(f)]
        for m, u in zip(maps, unix):
            m["t_end"] = u - float(st["t0_unix"])
            m["t"] = m["t_end"] - m["t_span"] / 2
    truth = read_truth(os.path.join(logdir, f"ground_truth{suffix}.csv"))
    return maps, truth, np.array(sorted(truth))


def grants(logdir, maps=None, carrier_prb=273):
    """The downlink grants the maps were built from (prb_stats.py's main numbers), from the
    run's DCI CSV (dci_*.csv, GUI runs record it), and the pilots per map from the maps."""
    import csv, glob, statistics
    from collections import Counter
    out = {}
    files = sorted(glob.glob(os.path.join(logdir, "dci_*.csv")))
    rows = [r for f in files for r in csv.DictReader(open(f)) if r["direction"] == "DL" and r["nof_prb"]]
    if rows:
        prb = [int(r["nof_prb"]) for r in rows]
        top = lambda key: [[k, n / len(rows)] for k, n in Counter(r[key] for r in rows).most_common(4)]
        med = statistics.median(prb)
        out.update(n=len(rows), carrier_prb=carrier_prb, prb_mean=statistics.mean(prb), prb_median=med,
                   prb_min=min(prb), prb_max=max(prb), full_band=sum(p == carrier_prb for p in prb) / len(prb),
                   n_sizes=len(set(prb)), n_allocs=len(set(r["prbs"] for r in rows)), sizes=top("nof_prb"),
                   mcs_median=statistics.median(int(r["dci_mcs"]) for r in rows), layers=top("nof_layers"),
                   dmrs_add_pos=top("dmrs_add_pos"), mapping=top("mapping"), symbols=top("time_length"),
                   range_res_m=3e8 / (2 * med * 12 * 30e3))
    if maps:
        out.update(snap_median=float(np.median([m["n_snap"] for m in maps])),
                   pilots_median=float(np.median([m["pilots"][1] for m in maps])),
                   k_step=int(np.median([m["k_step"] for m in maps])))
    return out or None


MICRO_DOPPLER_MS = 2.5   # a walking person's limbs reach this far from the body's speed: those echoes are the person


def site_walls(scene):
    """The site map's building walls, for the tracker's multipath ghost prediction (heights ignored)."""
    geo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sites", str(scene)[5:], "geometry.json")
    if not str(scene).startswith("site:") or not os.path.exists(geo):
        return None
    import nonideal
    W = nonideal.walls(json.load(open(geo))["buildings"])
    return None if W is None else (W[0], W[1], np.full(len(W[0]), -1e9), np.full(len(W[0]), 1e9))


def site_footprints(scene):
    """The site map's building footprints: no track starts inside one."""
    geo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sites", str(scene)[5:], "geometry.json")
    if not str(scene).startswith("site:") or not os.path.exists(geo):
        return []
    return [b["footprint"] for b in json.load(open(geo))["buildings"]]


def summary(logdir, gate_bins=2.0, guard_ms=None, object_m=None, ekf=None, rx=1, _more=None):
    """object_m: the tracker's object size (track.py), default settings' track_object_m or track.OBJECT_M.
    ekf: the 2D tracker's settings (track.EKF_DEFAULTS), default the settings' "tracking".
    rx: 1, 2 (a run with a second receiver: that one alone) or "fused" (both: the first's
    detector and angle scores, the tracks from both receivers' maps)."""
    settings = {}
    try:
        settings = json.load(open(os.path.join(logdir, "settings.json")))
    except (OSError, ValueError):
        pass
    if rx == "fused":   # the second receiver's rows, without tracks, then the first's with both
        two = summary(logdir, gate_bins, guard_ms, object_m, ekf, rx=2, _more=False)
        if not two.get("rows"):
            return dict(summary(logdir, gate_bins, guard_ms, object_m, ekf), error="no second receiver's maps")
        geo2 = (settings["gnb"], settings["rx2"], float(settings["boresight2_deg"]))
        return dict(summary(logdir, gate_bins, guard_ms, object_m, ekf, rx=1, _more=[(two["rows"], geo2)]), rx="fused")
    suffix = "_rx2" if rx == 2 else ""
    rx_settings = settings
    if rx == 2:   # this receiver's place and array normal stand in for the first's
        rx_settings = dict(settings, rx=settings.get("rx2"), boresight_deg=settings.get("boresight2_deg"))
    have_maps = os.path.exists(os.path.join(logdir, f"map2d{suffix}.csv"))
    have_truth = os.path.exists(os.path.join(logdir, f"ground_truth{suffix}.csv"))
    out = dict(logdir=logdir, settings=settings, have_maps=have_maps, have_truth=have_truth, rx=rx)
    if not (have_maps and have_truth):
        return out
    maps, truth, times = load(logdir, suffix)
    out["grants"] = grants(logdir, maps)
    if not maps or not len(times):
        return out
    try:   # the people get a wider speed gate: their legs and arms spread the echo in Doppler
        people = {k + 1 for k, c in enumerate(json.load(open(os.path.join(logdir, "agents.json")))["cls"]) if c == "person"}
    except (OSError, ValueError, KeyError):
        people = set()
    rows = []
    for i, m in enumerate(maps):
        guard = guard_ms if guard_ms is not None else 2 * abs(m["speed_step"])
        gt = truth_at(truth, times, m["t"], m["lam"])
        det, peak = estimates(m, guard)
        r = dict(i=i, t=m["t"], frame=m["frame"], slot=m["slot"], aarx=m["aarx"], car=gt["car"], cars=gt["cars"],
                 gt_excess=gt.get("excess"), gt_speed=gt.get("speed"), gt_gain_db=gt.get("gain_db"),
                 peak_excess=float(peak[0]), peak_speed=float(peak[1]), peak_db=float(peak[2]),
                 det=[[float(a), float(b), float(c)] for a, b, c in det],
                 aoa=[[float(a[0]), float(a[1]), float(a[2])] for a in m["aoa"]])
        gr, gs = gate_bins * m["m_per_bin"], gate_bins * max(abs(m["speed_step"]) * 2, track.SPEED_FLOOR_MS)
        near = lambda e, v, g: abs(e - g["excess"]) <= gr and abs(v - g["speed"]) <= (max(gs, MICRO_DOPPLER_MS) if g["car"] in people else gs)
        used = set()  # markers that lie on some car
        r["targets"] = []
        for g in gt["targets"]:
            tg = dict(car=g["car"], excess=g["excess"], speed=g["speed"], gain_db=g["gain_db"], aoa=g["aoa"])
            # NR-Scope's angle for this car: the AoA cell nearest its echo, inside the gate
            cells = [(a[0], a[1], a[2]) for a in m["aoa"]]  # range_m is already excess range
            near_cells = [c for c in cells if near(c[0], c[1], g)]
            if near_cells and np.isfinite(g["aoa"]):
                # two cells at one range and speed (NR-Scope's second arrival): the one in its direction
                c = min(near_cells, key=lambda c: abs(c[0] - g["excess"]) / gr + abs(c[1] - g["speed"]) / gs + abs(c[2] - g["aoa"]) / 1e3)
                tg["aoa_est"], tg["aoa_err"] = c[2], c[2] - g["aoa"]
            tg["in_view"] = 0 <= g["excess"] <= m["excess"][-1] and guard < abs(g["speed"]) <= m["speeds"][-1]
            tg["peak_hit"] = near(peak[0], peak[1], g)
            hits = [j for j, d in enumerate(det) if near(d[0], d[1], g)]
            used.update(hits)
            tg["det_hit"] = bool(hits)
            if hits:
                h = det[min(hits, key=lambda j: abs(det[j][0] - g["excess"]) / gr + abs(det[j][1] - g["speed"]) / gs)]
                tg["det_err"] = [h[0] - g["excess"], h[1] - g["speed"]]
            if tg["peak_hit"]:
                tg["peak_err"] = [peak[0] - g["excess"], peak[1] - g["speed"]]
            r["targets"].append(tg)
        view = [g for g in r["targets"] if g["in_view"]]
        r["in_view"] = bool(view)
        r["peak_hit"] = any(g["peak_hit"] for g in r["targets"])
        r["det_hit"] = any(g["det_hit"] for g in view)
        r["det_false"] = len(det) - len(used)
        rows.append(r)

    def score(key, car=None):
        """Detector: over every (map, car in view); peak: over the maps with a car in view."""
        if key == "peak":
            view = [r for r in rows if r["in_view"]]
            hits = [next(g for g in r["targets"] if g["peak_hit"]) for r in view if r["peak_hit"]]
        else:
            view = [g for r in rows for g in r["targets"] if g["in_view"] and car in (None, g["car"])]
            hits = [g for g in view if g["det_hit"]]
        errs = np.array([g[key + "_err"] for g in hits]) if hits else np.zeros((0, 2))
        return dict(n_in_view=len(view), n_hit=len(hits),
                    p_detect=len(hits) / len(view) if view else None,
                    rmse_range_m=float(np.sqrt((errs[:, 0] ** 2).mean())) if len(errs) else None,
                    rmse_speed_ms=float(np.sqrt((errs[:, 1] ** 2).mean())) if len(errs) else None,
                    bias_range_m=float(errs[:, 0].mean()) if len(errs) else None,
                    bias_speed_ms=float(errs[:, 1].mean()) if len(errs) else None)

    if object_m is None:
        object_m = float(settings.get("track_object_m") or track.OBJECT_M)
    track.group_measurements(rows, maps, object_m)
    angle = None
    if any(len(m["aoa"]) for m in maps):
        per = []
        for k in sorted({g["car"] for r in rows for g in r["targets"]}):
            view = [g for r in rows for g in r["targets"] if g["car"] == k and g["in_view"]]
            e = np.array([g["aoa_err"] for g in view if "aoa_err" in g])
            per.append(dict(car=k, n_in_view=len(view), n_angle=len(e),
                            rmse_deg=float(np.sqrt((e ** 2).mean())) if len(e) else None,
                            bias_deg=float(e.mean()) if len(e) else None,
                            median_abs_deg=float(np.median(np.abs(e))) if len(e) else None))
        e = np.array([g["aoa_err"] for r in rows for g in r["targets"] if g["in_view"] and "aoa_err" in g])
        angle = dict(cars=per, n_angle=len(e), rmse_deg=float(np.sqrt((e ** 2).mean())) if len(e) else None,
                     bias_deg=float(e.mean()) if len(e) else None,
                     median_abs_deg=float(np.median(np.abs(e))) if len(e) else None)
    classes = None
    try:  # scenario runs: the class of every agent, numbered as in the ground truth
        classes = json.load(open(os.path.join(logdir, "agents.json")))["cls"]
    except (OSError, ValueError, KeyError):
        pass
    half_len = None
    if classes:
        import scenario
        half_len = {k + 1: scenario.CLASSES[c]["size"][0] / 2 for k, c in enumerate(classes) if c in scenario.CLASSES}
    position = None
    if angle and settings.get("gnb") and settings.get("boresight_deg") is not None and _more is not False:
        cfg = track.ekf_cfg(ekf if ekf is not None else settings.get("tracking"))
        geo = (settings["gnb"], rx_settings["rx"], float(rx_settings["boresight_deg"]))
        scene = settings.get("scene") or (settings.get("scenario") or {}).get("scene")
        buildings = site_footprints(scene)
        pd = pd_map.site_grids(logdir, scene, [geo] + [g for _, g in _more or ()])   # each receiver's line of sight
        learned = None
        if cfg["ghost_map"] == 2:   # first pass without walls, learn them from its tracks, then track again
            import wall_learn
            track.localize_and_score(rows, geo, object_m, cfg=dict(cfg, ghost_map=0), half_len=half_len, more=_more or (), pd=pd, buildings=buildings)
            learned = wall_learn.learn([(r["t"], r["xy_tracks"]) for r in rows], settings["gnb"], rx_settings["rx"])
            walls = wall_learn.as_walls(learned)
        else:
            walls = site_walls(scene)
        position = track.localize_and_score(rows, geo, object_m, cfg=cfg, half_len=half_len, walls=walls, more=_more or (), pd=pd, buildings=buildings)
        if learned is not None:
            position["learned_walls"] = learned
    per_car = [dict(car=k, **score("det", k)) for k in range(1, len(truth[times[0]]["cars"]) + 1)]
    if classes:
        per_class = dict(
            detector=by_class(classes, [p for p in per_car if p["n_in_view"]], "n_in_view",
                              dict(p_detect=("p_detect", "rate"), n_hit=("n_hit", "sum"))),
            position=by_class(classes, position["cars"] if position else [], "maps_in_view",
                              dict(coverage=("coverage", "rate"), switches=("switches", "sum"),
                                   fair_coverage=("fair_coverage", "rate"), fair_switches=("fair_switches", "sum"),
                                   alone_share=("alone_share", "rate"),
                                   track_median_m=("track_median_m", "median"), pos_median_m=("pos_median_m", "median"))),
            angle=by_class(classes, angle["cars"] if angle else [], "n_in_view",
                           dict(n_angle=("n_angle", "sum"), median_abs_deg=("median_abs_deg", "median"))))
        out.update(classes=classes, per_class=per_class)
    m0 = maps[0]
    out.update(
        n_maps=len(maps), t_first=maps[0]["t"], t_last=maps[-1]["t"], sim_time_s=float(times[-1]),
        resolution=dict(range_m=m0["m_per_bin"], speed_ms=abs(m0["speed_step"]) * 2,
                        max_excess_m=float(m0["excess"][-1]), max_speed_ms=float(m0["speeds"][-1]),
                        window_s=float(np.mean([m["t_span"] for m in maps]))),
        gate=dict(bins=gate_bins), has_markers=any(len(m["markers"]) for m in maps),
        detector=score("det"), peak=score("peak"), angle=angle, position=position,
        n_cars=max([len(truth[t]["cars"]) for t in times[:1]] + [0]),
        per_car=per_car,
        false_per_map=float(np.mean([r["det_false"] for r in rows])),
        track=[dict(t=float(t), car=truth[t]["car"], cars=truth[t]["cars"]) for t in times[::max(1, len(times) // 400)]],
        rows=rows)
    return out


def by_class(classes, items, n_key, extra):
    """Scores per target class (scenario runs): items are per-target dicts with "car"; n_key
    counts the maps; extra = {name: (key, how)}, how "rate" (a fraction of n_key, weighted),
    "sum", or "median" (the n_key-weighted mean of the targets' values)."""
    out = {}
    for it in items or []:
        c = classes[it["car"] - 1] if 0 < it["car"] <= len(classes) else "?"
        out.setdefault(c, []).append(it)
    res = []
    for c, its in sorted(out.items()):
        n = sum(it.get(n_key) or 0 for it in its)
        r = dict(cls=c, targets=len(its), maps=n)
        for name, (key, how) in extra.items():
            v = [(it.get(key), it.get(n_key) or 0) for it in its if it.get(key) is not None]
            if how == "sum":
                r[name] = sum(x for x, _ in v)
            elif not v or sum(w for _, w in v) == 0:
                r[name] = None
            else:
                r[name] = sum(x * w for x, w in v) / sum(w for _, w in v)
        res.append(r)
    return res


def one_map(logdir, i, floor_db=40.0):
    return one_map_loaded(load(logdir), i, floor_db)


def one_map_loaded(loaded, i, floor_db=40.0):
    maps, truth, times = loaded
    m = maps[int(i)]
    gt = truth_at(truth, times, m["t"], m["lam"])
    # dB re the median non-zero cell, as plot_range_doppler.py --norm median; the page can
    # switch to --norm peak by subtracting peak_db
    pw = m["power"]
    med = np.median(pw[pw > 0]) if np.any(pw > 0) else 1.0
    db = 10 * np.log10(np.maximum(pw, 1e-30) / med)
    return dict(i=int(i), n=len(maps), t=m["t"], frame=m["frame"], slot=m["slot"],
                excess=m["excess"].round(3).tolist(), speeds=m["speeds"].round(4).tolist(),
                db=np.maximum(db, -100).round(2).tolist(), peak_db=float(db.max()),
                markers=[[float(r[0] - m["bin_los"] * m["m_per_bin"]), float(r[1]), float(r[2]), int(r[3])]
                         for r in m["markers"]],
                gt=gt)


def main():
    s = summary(sys.argv[1])
    rows = s.pop("rows", [])
    s.pop("track", None)
    print(json.dumps(s, indent=1, default=float))
    for r in rows[:40]:
        print(f"t {r['t']:6.3f}  car ({r['car'][0]:6.1f},{r['car'][1]:5.1f})  truth "
              + ("  ".join(f"car{g['car']} {g['excess']:6.1f} m {g['speed']:6.2f} m/s" for g in r["targets"]) or "   no car path   ")
              + f"  peak {r['peak_excess']:6.1f} m {r['peak_speed']:6.2f} m/s {r['peak_db']:5.1f} dB"
              + f"  {'HIT' if r.get('peak_hit') else '   '}  view {int(r['in_view'])}  markers {len(r['det'])}")


if __name__ == "__main__":
    main()
