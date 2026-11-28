"""Captures and ground truth for the import:
- capture: every chunk's maps of a sensing folder in time order (map2d.csv) with their UTC end times
  (map2d.csv.time.csv): t_rec + iq_offset_s + k * 2.5 s + ((slot - a0_k) mod 20480) * 0.5 ms, a0_k from the
  recording's occupancy grid (grid/b<k>.json);
- truth: the lidar movers (movers.json) as a CSV unix_s,id,lat,lon,cls (lidar frame -> east/north by the heading); labels.json applied; only the lidar front half (azimuth within +-90 deg).

  python3 capture.py [--config ...] maps SENSE_DIR OUT_DIR      python3 capture.py [--config ...] truth OUT.csv"""
import argparse, datetime, glob, json, os, re, shutil, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C


def t_rec(cfg):
    return datetime.datetime.strptime(cfg["recording"]["t_rec"], "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=datetime.timezone.utc).timestamp()


def maps(cfg, sense_dir, out):
    rec = cfg["recording"]; t0 = t_rec(cfg) + rec["iq_offset_s"]; rows, head = [], None
    for f in glob.glob(os.path.join(sense_dir, "b*.map2d.csv")):
        k = int(re.match(r"b(\d+)", os.path.basename(f)).group(1)); a0 = json.load(open(os.path.join(rec["lidar"], "grid", "b%d.json" % k)))["a0"]
        for line in open(f):
            if line.startswith("#"):
                head = head or line; continue
            c = line.split(",", 2); a = int(c[0]) * 20 + int(c[1])
            rows.append((t0 + k * 2.5 + ((a - a0) % 20480) * 5e-4, line))
    rows.sort(); shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
    with open(os.path.join(out, "map2d.csv"), "w") as f:
        f.write(head or "")
        for _, l in rows:
            f.write(l)
    with open(os.path.join(out, "map2d.csv.time.csv"), "w") as f:
        f.write("unix_s\n")
        for t, _ in rows:
            f.write("%.4f\n" % t)
    y = os.path.join(sense_dir, "sensing.yaml")
    if os.path.exists(y):
        shutil.copy(y, os.path.join(out, "nrscope_config.yaml"))
    return len(rows)


def truth(cfg, out):
    rec = cfg["recording"]; L = json.load(open(rec["session"]))["lidar"]
    KY = 6378137 * np.pi / 180; KX = KY * np.cos(np.radians(L["lat"])); h = np.radians(L["heading"]); n_ = 0
    with open(out, "w") as f:
        f.write("unix_s,id,lat,lon,cls\n")
        # the user's labels (Real tab, labels.json): tracks labelled ignore are left out, a label's class replaces the
        # movers' one, and tracks with the same name (a trailing number or "(n)" stripped: "GIO", "GIO 2") are one target
        lab = os.path.join(rec["lidar"], "labels.json"); lab = json.load(open(lab)).get("tracks", {}) if os.path.exists(lab) else {}
        base = lambda nm: re.sub(r"[\s(]*\d+\)?\s*$", "", nm.strip()).lower()
        named = {}
        for k, v in lab.items():
            if v.get("cls") != "ignore" and v.get("name", "").strip():
                named.setdefault(base(v["name"]), []).append(int(k))
        merge = {i: min(ids) for ids in named.values() for i in ids}
        for tr in json.load(open(os.path.join(rec["lidar"], "movers.json")))["tracks"]:
            lv = lab.get(str(tr["id"]), {})
            if lv.get("cls") == "ignore":
                continue
            cls = lv.get("cls") or tr.get("cls") or "person"; tid = merge.get(tr["id"], tr["id"])
            for t, x, y in tr["path"]:
                if y < 0:   # only the lidar's front half (azimuth 0 = +y, within +-90 deg): points behind it are not truth
                    continue
                e, n = x * np.cos(h) + y * np.sin(h), -x * np.sin(h) + y * np.cos(h)
                f.write("%.3f,%d,%.8f,%.8f,%s\n" % (t, tid, L["lat"] + n / KY, L["lon"] + e / KX, cls)); n_ += 1
    return n_


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["maps", "truth"]); ap.add_argument("paths", nargs="+")
    ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[]); a = ap.parse_args(); c = C.load(a.config, a.set)
    if a.what == "maps":
        print("%d maps -> %s" % (maps(c, a.paths[0], a.paths[1]), a.paths[1]))
    else:
        print("%d truth rows -> %s" % (truth(c, a.paths[0]), a.paths[0]))
