"""Real-tab replay of a scored run: per map (both streams, in time order) the accepted detections (range, speed, SNR,
position when they have an angle), the confirmed tracks (position, velocity) and the truth movers' bistatic ranges;
per track its life, speed, share of life on a mover and whether it is false (the paper's rule: on a mover less than
half its life); per mover the paper KPIs of the scored setup. Written as tracker_<name>.json in the lidar recording
folder, where the Real tab's Tracker selector lists it.

  <paper python> replay.py RUN_DIR LIDAR_DIR NAME"""
import collections as Cl, csv, json, os, pickle, sys
import numpy as np


def main(run_dir, lidar_dir, name):
    s = json.load(open(os.path.join(run_dir, "settings.json"))); geo = json.load(open(os.path.join(run_dir, "geometry_origin.json")))
    rep = pickle.load(open(os.path.join(run_dir, "replay.pkl"), "rb")); sc = json.load(open(os.path.join(run_dir, "scores.json")))
    o = np.radians(geo["origin"]); RE = 6371008.8   # the importer's enu(): east = dlon cos(lat0) R, north = dlat R
    ll = lambda x, y: [round(float(np.degrees(o[0] + y / RE)), 8), round(float(np.degrees(o[1] + x / (RE * np.cos(o[0])))), 8)]
    t0 = float(s["t0_unix"]); T, R = np.asarray(s["gnb"], float), np.asarray(s["rx"], float); LTR = float(np.linalg.norm(T - R))
    exc = lambda x, y: float(np.linalg.norm([x - T[0], y - T[1], 1.0 - T[2]]) + np.linalg.norm([x - R[0], y - R[1], 1.0 - R[2]]) - LTR)
    tru = Cl.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(run_dir, "ground_truth_cars.csv"))):
        if r["x"] != "nan":
            tru[float(r["t_s"])].append((int(r["car"]), float(r["x"]), float(r["y"])))
    tt = np.array(sorted(tru))
    def truth_at(t):
        i = int(np.clip(np.searchsorted(tt, t), 0, len(tt) - 1))
        return tru[tt[i]] if len(tt) and abs(tt[i] - t) < 0.3 else []
    maps, life, spd, on, mov, first, last = [], Cl.Counter(), Cl.defaultdict(list), Cl.Counter(), Cl.defaultdict(Cl.Counter), {}, {}
    for fr, trks in zip(rep["frames"], rep["out"]):
        t = float(fr[0]); rx = fr[5] if len(fr) > 5 else 1; tr = truth_at(t)
        dets = [dict(e=round(float(m["e"]), 2), v=round(float(m["v"]), 2), snr=None if m.get("p_db") is None else round(float(m["p_db"]), 1),
                     ll=None if m.get("xy") is None else ll(*m["xy"]), stream="trs" if rx == 1 else "dmrs") for m in fr[1]]
        tk = []
        for q in trks:
            i, x, y, vx, vy = int(q[0]), *map(float, q[1:5])
            tk.append(dict(id=i, ll=ll(x, y), v=[round(vx, 2), round(vy, 2)]))
            life[i] += 1; spd[i].append(np.hypot(vx, vy)); first.setdefault(i, t + t0); last[i] = t + t0
            d = [(np.hypot(x - a, y - b), k) for k, a, b in tr]
            if d and min(d)[0] <= 3.0:
                on[i] += 1; mov[i][min(d)[1]] += 1
        maps.append(dict(t=round(t + t0, 3), stream="trs" if rx == 1 else "dmrs", dets=dets, trk=tk, truth=[dict(id=k, e=round(exc(a, b), 2)) for k, a, b in tr]))
    tracks = [dict(id=i, maps=life[i], t0=round(first[i], 2), t1=round(last[i], 2), speed=round(float(np.median(spd[i])), 2),
                   near_pct=round(100.0 * on[i] / life[i]), mover=(mov[i].most_common(1)[0][0] if mov[i] else None), in_fov=True, false=on[i] < life[i] / 2)
              for i in sorted(life)]
    P = sc[rep["setup"]]["position"]; minutes = max((maps[-1]["t"] - maps[0]["t"]) / 60.0, 1e-9) if maps else 1.0
    kpi = dict(movers=[dict(id=c["car"], maps=c["maps_in_view"], fov_maps=c["maps_in_view"], coverage_strict=c["coverage"], coverage_fair=c["fair_coverage"],
                            switches=c["switches"], error_m=None if c.get("track_median_m") is None else round(c["track_median_m"], 1), miss_rate=None) for c in P["cars"]],
               false_tracks_fov=P["false_tracks"], false_tracks_per_min=round(P["false_tracks"] / minutes, 2), false_targets_per_map=round(P["false_per_map"], 2),
               n_tracks=P["n_tracks"], minutes=round(minutes, 2))
    c1 = [c for c in P["cars"] if c["car"] == 1]
    score = dict(maps=len(maps), dets_per_map=round(float(np.mean([len(m["dets"]) for m in maps])), 2) if maps else 0,
                 with_angle=round(float(np.mean([sum(d["ll"] is not None for d in m["dets"]) for m in maps])), 2) if maps else 0,
                 tracks=P["n_tracks"], false_tracks=P["false_tracks"], coverage_pct=round(100 * c1[0]["fair_coverage"], 1) if c1 else None,
                 gospa_m=round(P["gospa_m"], 2) if P.get("gospa_m") is not None else None, angle_err_deg_median=None, angle_err_deg_mad=None, setup=rep["setup"])
    out = os.path.join(lidar_dir, "tracker_%s.json" % name)
    json.dump(dict(label="%s (%s)" % (name, rep["setup"]), gnb=[*geo["geometry"]["gnb"], geo["geometry"]["gnb_h"]], rx=[*geo["geometry"]["rx"], geo["geometry"]["rx_h"]],
                   bore_deg=geo["geometry"]["boresight"], score=score, kpi=kpi, tracks=tracks, maps=maps, params=rep["params"]), open(out, "w"))
    print(out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
