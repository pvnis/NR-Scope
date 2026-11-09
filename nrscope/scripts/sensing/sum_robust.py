import sys, os, csv, glob, numpy as np, collections as C
sys.argv = ["x"]; exec(open("/home/sens/field/so_score.py").read().split("import os")[0])
O = sys.argv[1] if len(sys.argv) > 1 else "/home/sens/field/so/robust"
O = "/home/sens/field/so/robust"
lab = {"raw": "tablet, raw", "trs": "tablet, TRS-referenced", "ltab": "tablet + TRS + CSI-RS + DM-RS", "lnotab": "NO tablet: TRS + CSI-RS + DM-RS"}
tgt = {0: "25 m, +2.5 m/s, -20 dB", 1: "60 m, -5 m/s, -22 dB", 2: "100 m, +8 m/s, -25 dB"}
R = C.defaultdict(lambda: [0, 0, []]); maps_n = C.Counter(); oth = C.defaultdict(list); nb = C.Counter()
for d in glob.glob(O + "/m_b*"):
    _, b, v = os.path.basename(d).split("_"); p = d + "/map2d.csv"
    if not os.path.exists(p) or not os.path.exists(p + ".truth.csv"): continue
    nb[v] += 1; T = {}
    for r in csv.DictReader(open(p + ".truth.csv")): T.setdefault((int(r["frame"]), int(r["slot"])), []).append((int(r["target"]), float(r["range_m"]), float(r["speed_ms"])))
    for m in maps(p):
        maps_n[v] += 1; acc = [x for x in m["markers"] if x[3] == 0]; u = 0
        for tid, rr, ss in T.get((m["frame"], m["slot"]), []):
            k = (v, tid); R[k][1] += 1; ok = [x for x in acc if abs(x[1] - ss) <= 2 * m["res"] and abs(x[0] - rr) <= 10]
            if ok: R[k][0] += 1; u += 1; R[k][2].append(min(ok, key=lambda x: abs(x[0] - rr))[0] - rr)
        oth[v].append(len(acc) - u)
for v in ("raw", "trs", "ltab", "lnotab"):
    if not nb[v]: continue
    print("%s  (%d bursts, %d maps, other detections/map %.1f)" % (lab[v], nb[v], maps_n[v], np.mean(oth[v]) if oth[v] else 0))
    for tid in (0, 1, 2):
        f, n, e = R[(v, tid)]
        if n: e = np.array(e); print("    %-24s found %4d/%-4d %3.0f%%  | within 1 bin %3.0f%% of found | range error median %+.1f m, RMS %.1f m" % (tgt[tid], f, n, 100 * f / n, 100 * np.mean(np.abs(e) <= 2.44) if len(e) else 0, np.median(e) if len(e) else np.nan, np.sqrt(np.mean(e ** 2)) if len(e) else np.nan))
