"""Runs trackers on one run's NR-Scope maps and scores them exactly as our tracker is scored.

analyze.summary does everything as usual (detections, grouping, localisation, scoring); only the
tracker it calls (track.run_ekf_tracker) is swapped for the one under test. So every tracker gets
the same measurements and the same scoring (GOSPA, false tracks, coverage, switches...).

python tracking_baselines/harness.py <run> [tracker ...] [--setups 1,2,fused] [--set key=value ...]
  trackers: names from trackers.REGISTRY (default: all); --set: the tracker's parameters
Writes tracking_baselines/results/<run>/<tracker>[__<tag>]__<setup>.json, then summary.py tables them.
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "sim"))
sys.path.insert(0, HERE)
import analyze  # noqa: E402
import track  # noqa: E402
from trackers import REGISTRY  # noqa: E402

ORIGINAL = track.run_ekf_tracker


def score(run, name, setup, params):
    fn = REGISTRY[name]
    if fn is None:   # ours: the original tracker, params = its settings (track.EKF_DEFAULTS keys)
        track.run_ekf_tracker = ORIGINAL
    else:
        track.run_ekf_tracker = lambda frames, geo=None, cfg=None, walls=None, pd=None, buildings=(): fn(frames, **params)
    try:
        s = analyze.summary(os.path.join(ROOT, "logs", run), rx=setup, ekf=params if fn is None else {})
    finally:
        track.run_ekf_tracker = ORIGINAL
    pos = s.get("position") or {}
    return dict(run=run, tracker=name, setup=str(setup), params=params,
                classes=s.get("classes"), per_class=(s.get("per_class") or {}).get("position"),
                position={k: v for k, v in pos.items() if k != "ekf"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("trackers", nargs="*")
    ap.add_argument("--setups", default="1,2,fused")
    ap.add_argument("--set", nargs="*", default=[])
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    params = {}
    for kv in a.set:
        k, v = kv.split("=", 1)
        params[k] = json.loads(v)
    out = os.path.join(HERE, "results", a.run)
    os.makedirs(out, exist_ok=True)
    for name in a.trackers or list(REGISTRY):
        for setup in a.setups.split(","):
            setup = setup if setup == "fused" else int(setup)
            t0 = time.time()
            r = score(a.run, name, setup, params)
            r["seconds"] = round(time.time() - t0, 1)
            p = r["position"]
            print(f"{name:42s} {str(setup):6s} GOSPA {p.get('gospa_m') or float('nan'):.2f} m  false/map "
                  f"{p.get('false_per_map') or 0:.2f}  missed/map {p.get('missed_per_map') or 0:.2f}  "
                  f"tracks {p.get('n_tracks')} false tracks {p.get('false_tracks')}  ({r['seconds']} s)", flush=True)
            fn = f"{name}{'__' + a.tag if a.tag else ''}__{setup}.json"
            json.dump(r, open(os.path.join(out, fn), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
