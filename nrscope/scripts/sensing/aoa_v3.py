"""Three raw-injected targets: per variant and target, share of maps where it is detected (range +-1 bin, speed
+-1.5 res, strongest match) and, of those, angle within 6 deg; median angle and median |error|. usage: aoa_v3.py <glob>"""
import sys, glob, os, csv, numpy as np, collections as C
G = sys.argv[1]; sys.argv = ["x", "/nonexistent"]; sys.path.insert(0, "/home/sens/field"); from track_real import read_maps
S = C.defaultdict(lambda: [0, []]); TA = {}
for d in sorted(glob.glob(G)):
    p = d + "/map2d.csv"
    if not os.path.exists(p + ".truth.csv"): continue
    v = os.path.basename(d).split("_", 2)[2]; TR = C.defaultdict(list)
    for r in csv.DictReader(open(p + ".truth.csv")):
        TR[(int(r["frame"]), int(r["slot"]))].append((int(r["target"]), float(r["range_m"]), float(r["speed_ms"]), float(r["angle_deg"])))
    for m in read_maps(p):
        for (k, rr, vv, aa) in TR.get((m["fr"], m["sl"]), []):
            TA[k] = aa; s = S[(v, k)]; s[0] += 1
            c = [x for x in m["aoa"] if abs(x[0] - rr) <= m["mpb"] and abs(x[1] - vv) <= 1.5 * m["res"]]
            if c: s[1].append(max(c, key=lambda x: x[3])[2])
print("%-11s %-26s %-26s %-26s" % ("variant", "25 m @+15", "60 m @-25", "100 m @+30"))
for v in sorted(set(k[0] for k in S)):
    out = "%-11s" % v
    for k in (0, 1, 2):
        n, L = S.get((v, k), [0, []]); L = np.array(L)
        out += " %-26s" % (f"det {100*len(L)/max(n,1):3.0f}% ok {100*np.mean(abs(L-TA[k])<6) if len(L) else 0:3.0f}% med {np.median(L) if len(L) else float('nan'):+5.1f}")
    print(out)
