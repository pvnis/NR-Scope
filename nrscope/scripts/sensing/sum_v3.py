"""Main table: per variant / speed class / level, share of injected targets found (speed within 2 res, range within
10 m) and, of those, share within one range bin (2.44 m). usage: sum_v3.py <dir> [variants...]"""
import sys, os, csv, glob, numpy as np, collections as C
sys.argv, ARGS = ["x"], sys.argv[1:]; exec(open("/home/sens/field/so_score.py").read().split("import os")[0])
O = ARGS[0]; VS = ARGS[1:]; R = C.defaultdict(lambda: [0, 0, 0])
for d in glob.glob(O + "/m_b*"):
    parts = os.path.basename(d).split("_")
    if len(parts) != 5: continue
    _, b, v, tg, lv = parts; p = d + "/map2d.csv"
    if VS and v not in VS or not os.path.exists(p + ".truth.csv"): continue
    T = {}
    for row in csv.DictReader(open(p + ".truth.csv")): T.setdefault((int(row["frame"]), int(row["slot"])), []).append((float(row["range_m"]), float(row["speed_ms"])))
    r = R[(v, tg, lv)]
    for m in maps(p):
        acc = [x for x in m["markers"] if x[3] == 0]
        for rr, ss in T.get((m["frame"], m["slot"]), []):
            r[0] += 1; ok = [x for x in acc if abs(x[1] - ss) <= 2 * m["res"] and abs(x[0] - rr) <= 10]
            if ok: e = min(ok, key=lambda x: abs(x[0] - rr))[0] - rr; r[1] += 1; r[2] += abs(e) <= 2.44
print("%-12s" % "variant" + "".join("%15s" % f"{tg} {lv}" for tg in ("slow", "fast", "faster") for lv in ("-10", "-20")))
for v in sorted(set(k[0] for k in R), key=lambda x: (VS.index(x) if x in VS else 99, x)):
    s = "%-12s" % v
    for tg in ("slow", "fast", "faster"):
        for lv in ("-10", "-20"):
            r = R.get((v, tg, lv)); s += "%15s" % (f"{100*r[1]/r[0]:.0f}% ({100*r[2]/max(r[1],1):.0f}%)" if r and r[0] else "-")
    print(s)
