"""Tracking and scoring of an imported run with the paper tracker (harness OURS_full), run under the paper's python
from its folder. Each setup in tracker.setups; optional self-calibration (per-stream p_d and clutter from a first fused
pass, no truth). Two adaptations for real runs, nothing in the paper code edited: each receiver's maps are timed by
its own UTC time file, and the run's line-of-sight map (pd_grids.npz) serves every receiver (both streams are one
receiver). The replay setup's frames and confirmed tracks are saved (replay.pkl) for replay.py.

  <paper python> score.py RUN_DIR ARGS_JSON        (ARGS_JSON: cfg["tracker"] + "replay_setup"; tracker.t_window_unix cuts maps and truth)"""
import collections as Cl, csv, json, os, pickle, sys
import numpy as np


def main(run_dir, A):
    os.chdir(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tracker")); sys.path[:0] = ["sim", "tracking_baselines"]   # the paper tracker (sensing/tracker)
    import analyze, pd_map, track as TK
    if A.get("hit_m"):   # diagnostic only: the radius within which a track counts as on the target (paper: 3 m)
        TK.XY_HIT_M = float(A["hit_m"])
    analyze.MICRO_DOPPLER_MS = A["micro_doppler_ms"]
    _load = analyze.load
    def load(logdir, suffix=""):
        maps, truth, times = _load(logdir, suffix); st = json.load(open(os.path.join(logdir, "settings.json")))
        if suffix and st.get("map_clock") == "utc":
            unix = [float(r["unix_s"]) for r in csv.DictReader(open(os.path.join(logdir, "map2d%s.csv.time.csv" % suffix)))]
            for m, u in zip(maps, unix):
                m["t_end"] = u - float(st["t0_unix"]); m["t"] = m["t_end"] - m["t_span"] / 2
        W = A.get("t_window_unix")
        if W:   # track and score only inside the window (cfg tracker.t_window_unix): maps and truth cut alike
            lo, hi = W[0] - float(st["t0_unix"]), W[1] - float(st["t0_unix"])
            maps = [m for m in maps if lo <= m["t"] <= hi]; truth = {k: v for k, v in truth.items() if lo <= k <= hi}
            sh = float(A.get("truth_shift_s") or 0)
            if sh:   # chance control: the truth moved in time (cyclic inside the window), detections and tracks unchanged
                ks = sorted(truth); new = {}
                for k in ks:
                    tgt = lo + ((k - lo + sh) % (hi - lo)); src = ks[int(np.argmin([abs(q - tgt) for q in ks]))]; new[k] = truth[src]
                truth = new
            times = np.array(sorted(truth))
        return maps, truth, times
    analyze.load = load
    z = np.load(os.path.join(run_dir, "pd_grids.npz"))
    zs = [z] + ([np.load(os.path.join(run_dir, "pd_grids_rx2.npz"))] if os.path.exists(os.path.join(run_dir, "pd_grids_rx2.npz")) else [])
    def site_grids(logdir, scene, geos):   # one line-of-sight map per receiver geometry (pd_grids_rx2.npz: a second receiver elsewhere)
        out = []
        for g in geos:
            k = [float(v) for v in (*g[0], *g[1], g[2])]; m = [q for q in zs if np.allclose(k, q["key"])]
            assert m, "line-of-sight map made for another geometry"
            out.append((*m[0]["o0"], pd_map.STEP_M, m[0]["g0"]))
        return out
    pd_map.site_grids = site_grids
    CAP = {}; _run = TK.run_ekf_tracker
    def run_(frames, geo, *a, **k):
        out = _run(frames, geo, *a, **k); CAP["frames"], CAP["out"] = frames, out; return out
    TK.run_ekf_tracker = run_
    if A.get("external_tracks"):   # tracks made elsewhere (tbd.py paths): per frame [id, x, y, vx, vy], scored as the tracker's
        X = json.load(open(A["external_tracks"])); t0 = float(json.load(open(os.path.join(run_dir, "settings.json")))["t0_unix"])
        tm = np.array([u - t0 - sp / 2 for u, sp in zip(X["unix"], X["span"])])
        def ext(frames, geo, *a, **k):
            out, prev = [], {}
            for f in frames:
                j = int(np.argmin(np.abs(tm - f[0]))); T_, R_, _, bore = f[4]; cur = []
                for i, e, ang, v in (X["tracks"][j] if abs(tm[j] - f[0]) < 0.05 else []):
                    p = TK.localize(e, ang, T_, R_, np.degrees(bore))
                    if p is not None:
                        q = prev.get(i); dt = f[0] - q[0] if q else 0
                        cur.append([i, p[0], p[1], (p[0] - q[1]) / dt if q and dt > 0 else 0.0, (p[1] - q[2]) / dt if q and dt > 0 else 0.0]); prev[i] = (f[0], p[0], p[1])
                out.append(cur)
            CAP["frames"], CAP["out"] = frames, out; return out
        TK.run_ekf_tracker = ext
    import harness
    name = os.path.basename(run_dir.rstrip("/")); res = {}
    def score(setup, params):
        r = harness.score(name, "OURS_full", setup if setup == "fused" else int(setup), params)
        return dict(position={k: v for k, v in r["position"].items() if k != "ekf"}, per_class=r.get("per_class"))
    params = dict(A["params"])
    for s in A["setups"]:
        res[s] = score(s, params)
        if s == A["replay_setup"]:
            rep = dict(frames=CAP["frames"], out=CAP["out"], setup=s, params=params)
    if A.get("selfcal"):
        score("fused", params); hit, tot, clut = Cl.Counter(), Cl.Counter(), Cl.defaultdict(list)
        for fr, trks in zip(CAP["frames"], CAP["out"]):
            rx = fr[5] if len(fr) > 5 else 1; T, R, L, bore = fr[4]; rb, sc = fr[2], fr[3]; used = set()
            for tq in trks:
                x = np.array(tq[1:5]); a = np.arctan2(x[1] - R[1], x[0] - R[0]) - bore
                if abs((a + np.pi) % (2 * np.pi) - np.pi) > np.deg2rad(pd_map.FOV_DEG):
                    continue
                h = TK.bistatic_h(x, T, R, L, bore)[0]; tot[rx] += 1
                near = [j for j, m in enumerate(fr[1]) if abs(m["e"] - h[0]) <= 2 * rb and abs(m["v"] - h[1]) <= 2 * max(sc, TK.SPEED_FLOOR_MS)]
                if near:
                    hit[rx] += 1; used.update(near)
            clut[rx].append(len(fr[1]) - len(used))
        cl = lambda v: float(min(0.95, max(0.05, v)))
        est = dict(p_d=cl(hit[1] / max(tot[1], 1)), p_d_rx2=cl(hit[2] / max(tot[2], 1)), clutter=max(0.1, float(np.mean(clut[1]))), clutter_rx2=max(0.1, float(np.mean(clut[2]))))
        res["selfcal_estimates"] = est; res["fused_selfcal"] = score("fused", dict(params, **est))
        if A["replay_setup"] == "fused_selfcal":
            rep = dict(frames=CAP["frames"], out=CAP["out"], setup="fused_selfcal", params=dict(params, **est))
    if "rep" in locals():
        pickle.dump(rep, open(os.path.join(run_dir, "replay.pkl"), "wb"))
    json.dump(res, open(os.path.join(run_dir, "scores.json"), "w"), indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), json.loads(sys.argv[2]))
