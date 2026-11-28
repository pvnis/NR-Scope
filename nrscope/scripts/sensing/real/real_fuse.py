#!/usr/bin/env python3
"""Real-data sensing + tracking of one recording: the TRS and the DM-RS as two fused sensors, the paper tracker.

  python3 real_fuse.py --config configs/<recording>.json [--set key=value ...] [--name NAME] [--maps-only]
  --maps-only: calibration, estimates, sensing and the maps (cap_<NAME>_<stream>) only, e.g. on a second receiver's host
  (session = a copy whose "lidar" is that receiver's pose); its truth and scoring are made on the lidar host.

Steps (each a module of this folder, also runnable alone): truth (capture.py, from the lidar movers unless given) ->
angle calibration (calib.py) -> per-stream estimates (estimates.py, kept) -> sensing (sense.py) -> captures (capture.py)
-> import as two receivers + line-of-sight map (importer.py, los.py) -> tracking and scoring (score.py) ->
Real-tab replay (replay.py). Every parameter: cfg.py (python3 cfg.py --config ... prints them all).
Results: <data>/<recording>/results.jsonl, the run in paper/logs/real_<NAME>, the replay tracker_<NAME>.json next to
the lidar recording (Real tab, Tracker selector)."""
import argparse, json, os, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cfg as C, calib, estimates, sense, capture, importer


def run(c, name, maps_only=False):
    rec = c["recording"]; D = C.recording_dir(c); os.makedirs(D, exist_ok=True); geo = C.geometry(c); res = dict(name=name, time=time.strftime("%Y-%m-%d %H:%M:%S"), config=c, geometry=geo)
    truth = rec["truth"] or os.path.join(D, "truth_lidar.csv")
    if not rec["truth"] and not maps_only:
        res["truth_rows"] = capture.truth(c, truth)
    res["aoa_cal_deg"], res["aoa_cal_source"] = calib.get(c)
    caps = {}
    for s in ("trs", "dmrs"):
        e, ok, n = estimates.build(c, s); res["estimates_" + s] = "%d/%d" % (ok, n)
        sd = os.path.join(D, "sense_%s_%s" % (name, s)); shutil.rmtree(sd, ignore_errors=True)
        y, ok, n = sense.run(c, s, res["aoa_cal_deg"], sd); res["sensed_" + s] = "%d/%d" % (ok, n)
        caps[s] = os.path.join(D, "cap_%s_%s" % (name, s)); res["maps_" + s] = capture.maps(c, sd, caps[s]); shutil.rmtree(sd)
    if maps_only:
        res["caps"] = caps; return res
    run_dir, seen = importer.run(c, geo, caps["trs"], caps["dmrs"], truth, "real_" + name); res["run"] = run_dir; res["los_seen_share"] = seen
    args = dict(c["tracker"], replay_setup=c["replay"]["setup"]); P = c["paths"]
    r = subprocess.run([P["python"], os.path.join(HERE, "score.py"), run_dir, json.dumps(args)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("scoring failed:\n" + r.stderr[-1500:])
    res["scores"] = json.loads(r.stdout.strip().splitlines()[-1])
    r = subprocess.run([P["python"], os.path.join(HERE, "replay.py"), run_dir, rec["lidar"], name], capture_output=True, text=True)
    res["replay"] = r.stdout.strip() if r.returncode == 0 else "FAILED: " + r.stderr[-500:]
    return res


def summary(res):
    lines = ["%s: aoa_cal %s; maps TRS %s, DM-RS %s; line of sight %.0f %% of the grid" % (res["name"], res.get("aoa_cal_deg"), res.get("maps_trs"), res.get("maps_dmrs"), 100 * res.get("los_seen_share", 0))]
    for k, v in res.get("scores", {}).items():
        if not isinstance(v, dict) or "position" not in v:
            lines.append("  %-14s %s" % (k, json.dumps(v))); continue
        p = v["position"]
        movers = ", ".join("#%d %.0f%% (strict %.0f%%)" % (c["car"], 100 * c["fair_coverage"], 100 * c["coverage"]) for c in p["cars"] if c["maps_in_view"])
        lines.append("  %-14s tracks %3d, false %3d, false/map %.2f | coverage %s" % (k, p["n_tracks"], p["false_tracks"], p["false_per_map"], movers))
    lines.append("  replay: %s" % res.get("replay"))
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[]); ap.add_argument("--name", default="fused"); ap.add_argument("--maps-only", action="store_true")
    a = ap.parse_args(); c = C.load(a.config, a.set)
    try:
        res = run(c, a.name, a.maps_only)
    except Exception as e:
        res = dict(name=a.name, config=c, FAILED=str(e)[-1500:])
    with open(os.path.join(C.recording_dir(c), "results.jsonl"), "a") as f:
        f.write(json.dumps(res) + "\n")
    print(summary(res) if "FAILED" not in res else "FAILED: " + res["FAILED"])
