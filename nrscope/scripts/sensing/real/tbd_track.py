#!/usr/bin/env python3
"""Track-before-detect run on a recording's existing captures (real_fuse.py's cap_<CAPS>_trs / cap_<CAPS>_dmrs, whose
map windows it reuses): cubes from the estimates -> TBD recursion -> TBD capture -> import (TBD as receiver 1) ->
paper tracker, setup 1 -> scores + Real-tab replay. Results appended to <data>/<recording>/results.jsonl.

  python3 tbd_track.py --config configs/<recording>.json --caps CAPS --name NAME [--as-rx2] [--set key=value ...]"""
import argparse, json, os, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cfg as C, calib, importer, tbd, real_fuse


def run(c, caps_name, name, as_rx2=False):
    rec = c["recording"]; D = C.recording_dir(c); geo = C.geometry(c)
    res = dict(name=name, kind="tbd", caps=caps_name, time=time.strftime("%Y-%m-%d %H:%M:%S"), config=c, geometry=geo)
    truth = rec["truth"] or os.path.join(D, "truth_lidar.csv")
    res["aoa_cal_deg"], res["aoa_cal_source"] = calib.get(c)
    caps = {s: os.path.join(D, "cap_%s_%s" % (caps_name, s)) for s in ("trs", "dmrs")}
    out = os.path.join(D, "cap_%s_tbd" % name); shutil.rmtree(out, ignore_errors=True)
    res["tbd"] = tbd.run(c, caps, res["aoa_cal_deg"], out); res["maps_trs"] = res["tbd"]["maps"]; res["maps_dmrs"] = 0
    if as_rx2:   # NR-Scope's TRS stream as receiver 1, this stream as receiver 2, fused as in real_fuse.py
        run_dir, seen = importer.run(c, geo, caps["trs"], out, truth, "real_" + name); args = dict(c["tracker"], replay_setup="fused")
    else:
        run_dir, seen = importer.run(c, geo, out, caps["dmrs"], truth, "real_" + name); args = dict(c["tracker"], setups=["1"], replay_setup="1")
    res["run"] = run_dir; res["los_seen_share"] = seen; P = c["paths"]
    if c["tbd"]["output"] == "paths":
        args["external_tracks"] = os.path.join(out, "tbd_paths.json")
    r = subprocess.run([P["python"], os.path.join(HERE, "score.py"), run_dir, json.dumps(args)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("scoring failed:\n" + r.stderr[-1500:])
    res["scores"] = json.loads(r.stdout.strip().splitlines()[-1])
    r = subprocess.run([P["python"], os.path.join(HERE, "replay.py"), run_dir, rec["lidar"], name], capture_output=True, text=True)
    res["replay"] = r.stdout.strip() if r.returncode == 0 else "FAILED: " + r.stderr[-500:]
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[]); ap.add_argument("--caps", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--as-rx2", action="store_true", help="fuse this stream (as receiver 2) with NR-Scope's TRS capture")
    a = ap.parse_args(); c = C.load(a.config, a.set)
    try:
        res = run(c, a.caps, a.name, a.as_rx2)
    except Exception as e:
        import traceback; res = dict(name=a.name, kind="tbd", config=c, FAILED=traceback.format_exc()[-1500:])
    with open(os.path.join(C.recording_dir(c), "results.jsonl"), "a") as f:
        f.write(json.dumps(res) + "\n")
    print(real_fuse.summary(res) + "\n  tbd %s" % json.dumps(res.get("tbd")) if "FAILED" not in res else "FAILED: " + res["FAILED"])
