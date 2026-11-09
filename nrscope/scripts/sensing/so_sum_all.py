"""Aggregate run_eval_all.sh: per method, target speed class and level, targets found over all bursts."""
import sys, os, csv, glob, numpy as np, collections as C
sys.argv = ["x"]; exec(open("/home/sens/field/so_score.py").read().split("import os")[0])
O = "/home/sens/field/so/all"; R = C.defaultdict(lambda: [0, 0, 0, 0, [], 0, []])   # maps, found, within2bins, n, err, alias, other
for d in glob.glob(O + "/m_b*"):
    _, b, v, tg, lv = os.path.basename(d).split("_"); p = d + "/map2d.csv"; k = (v, tg, lv); r = R[k]
    if not os.path.exists(p): continue
    T = {}
    for row in csv.DictReader(open(p + ".truth.csv")): T.setdefault((int(row["frame"]), int(row["slot"])), []).append((float(row["range_m"]), float(row["speed_ms"])))
    for m in maps(p):
        r[0] += 1; acc = [x for x in m["markers"] if x[3] == 0]; u = 0
        for rr, ss in T.get((m["frame"], m["slot"]), []):
            r[3] += 1; ok = [x for x in acc if abs(x[1] - ss) <= 2 * m["res"] and abs(x[0] - rr) <= 10]
            if ok:
                e = min(ok, key=lambda x: abs(x[0] - rr))[0] - rr; r[1] += 1; u += 1; r[2] += abs(e) <= 4.9; r[4].append(e)
            elif any(abs(x[0] - rr) <= 5 for x in acc): r[5] += 1
        r[6].append(len(acc) - u)
lab = {"raw": "tablet, raw", "proc": "tablet, stitching", "trs": "tablet, TRS-referenced", "light": "no tablet: TRS + DM-RS", "trsonly": "no tablet: TRS only"}
print("%-24s %-5s %-4s %6s %16s %14s %10s %12s" % ("method", "speed", "dB", "maps", "found", "within 2 bins", "alias", "other/map"))
for v in ("raw", "proc", "trs", "light", "trsonly"):
    for tg in ("slow", "fast"):
        for lv in ("-10", "-20"):
            r = R.get((v, tg, lv))
            if not r or not r[3]: continue
            print("%-24s %-5s %-4s %6d %7d/%-4d %3.0f%% %8d %13d %9.1f" % (lab[v], tg, lv, r[0], r[1], r[3], 100 * r[1] / r[3], r[2], r[5], np.mean(r[6])))
