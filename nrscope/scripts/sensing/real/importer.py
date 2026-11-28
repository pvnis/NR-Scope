"""Import the two streams' captures as one run of the paper chain: the TRS capture as receiver 1 and the DM-RS capture
as receiver 2 at the same place (sim/import_capture.py, then map2d_rx2 / ground_truth_rx2 / rx2 settings), the origin
saved (geometry_origin.json), and the line-of-sight map (los.py) written into the run as the tracker's cache.

  python3 importer.py [--config ...] TRS_CAPTURE DMRS_CAPTURE TRUTH.csv RUN_NAME"""
import argparse, json, os, shutil, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C, los


def run(cfg, geo, cap_trs, cap_dmrs, truth, name):
    P = cfg["paths"]; out = os.path.join(P["paper"], "logs", name); tmp = out + "_rx2tmp"
    G = ["--origin", "%.8f,%.8f" % geo["origin"], "--gnb", "%.8f,%.8f" % geo["gnb"], "--rx", "%.8f,%.8f" % geo["rx"],
         "--boresight", "%.3f" % geo["boresight"], "--gnb-h", "%g" % geo["gnb_h"], "--rx-h", "%g" % geo["rx_h"], "--scene", "site:" + cfg["recording"]["site"]]
    for d, cap in ((out, cap_trs), (tmp, cap_dmrs)):
        shutil.rmtree(d, ignore_errors=True)
        r = subprocess.run([P["python"], "sim/import_capture.py", cap, "--truth", truth, *G, "--name", name, "--out", d], cwd=P["paper"], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError("import_capture failed for %s:\n%s" % (cap, (r.stdout + r.stderr)[-1500:]))
    for a, b in (("map2d.csv", "map2d_rx2.csv"), ("map2d.csv.time.csv", "map2d_rx2.csv.time.csv"),
                 ("ground_truth.csv", "ground_truth_rx2.csv"), ("ground_truth_cars.csv", "ground_truth_rx2_cars.csv")):
        shutil.copy(os.path.join(tmp, a), os.path.join(out, b))
    shutil.rmtree(tmp)
    s = json.load(open(os.path.join(out, "settings.json"))); s["rx2"] = s["rx"]; s["boresight2_deg"] = s["boresight_deg"]
    json.dump(s, open(os.path.join(out, "settings.json"), "w"), indent=1)
    json.dump(dict(origin=list(geo["origin"]), geometry=geo), open(os.path.join(out, "geometry_origin.json"), "w"), indent=1)
    g = los.write(cfg, out)
    return out, float(g.mean())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cap_trs"); ap.add_argument("cap_dmrs"); ap.add_argument("truth"); ap.add_argument("name")
    ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[]); a = ap.parse_args(); c = C.load(a.config, a.set)
    out, seen = run(c, C.geometry(c), a.cap_trs, a.cap_dmrs, a.truth, a.name); print("%s (line of sight %.1f %% of the grid)" % (out, 100 * seen))
