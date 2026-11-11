"""KPI gaps on target-free data (/mnt/mydrive/field/kpi): per variant
 - detections per map and per range-Doppler cell (accepted markers, verdict 0), clutter removal on (Pfa upper bound:
   real movers count as false alarms here), also split |v| < 2 res (residual clutter) vs moving
 - clutter-to-noise ratio (CNR): zero-Doppler cells (+-1 bin) at the direct path (strongest) and over all range bins,
   against the median of the cells with |f| > f_max/2, with and without clutter removal; suppression = CNR drop (dB).
usage: kpi_pfa.py [dir]"""
import sys, glob, os, numpy as np, collections as C
O = sys.argv[1] if len(sys.argv) > 1 else "/mnt/mydrive/field/kpi"
def maps(p):
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split(","); nb, nf, nm = int(f[12]), int(f[13]), int(f[14])
        P = np.array(f[21:21 + nb * nf], float).reshape(nb, nf); p0 = 21 + nb * nf
        mk = [tuple(float(f[p0 + 4 * i + j]) for j in range(4)) for i in range(nm)]
        tsp, fc = float(f[8]), float(f[11]); res = 3e8 / fc / 2 / max(tsp, 1e-3)
        yield dict(fr=int(f[0]), sl=int(f[1]), P=P, mk=mk, res=res, blos=float(f[16]))
S = C.defaultdict(lambda: C.defaultdict(list))
for d in sorted(glob.glob(O + "/*_b*_*")):
    cl, b, v = os.path.basename(d).split("_", 2); p = d + "/map2d.csv"
    if not os.path.exists(p): continue
    for m in maps(p):
        P = m["P"]; nb, nf = P.shape; j0 = nf // 2; noise = np.median(P[:, np.r_[0:nf // 4, 3 * nf // 4:nf]])
        z = P[:, j0 - 1:j0 + 2].max(1); bl = int(round(m["blos"]))
        s = S[(v, cl)]; s["cnr_los"].append(10 * np.log10(z[min(bl, nb - 1)] / noise)); s["cnr_all"].append(10 * np.log10(z.sum() / (noise * nb)))
        acc = [x for x in m["mk"] if x[3] == 0]; s["det"].append(len(acc)); s["cells"].append(nb * nf)
        s["det_slow"].append(sum(abs(x[1]) < 2 * m["res"] for x in acc)); s["maps_b"].append(b)
print("%-9s %6s %9s %11s %10s %9s | %s" % ("variant", "maps", "dets/map", "Pfa/cell", "slow/map", "moving", "CNR at LoS: no removal -> removal (suppression); all bins"))
for v in sorted({k[0] for k in S}):
    a, n = S[(v, "clut")], S[(v, "noclut")]
    if not a["det"]: continue
    dm = np.mean(a["det"]); pc = np.sum(a["det"]) / np.sum(a["cells"])
    line = "%-9s %6d %9.2f %11.1e %10.2f %9.2f" % (v, len(a["det"]), dm, pc, np.mean(a["det_slow"]), dm - np.mean(a["det_slow"]))
    if n["cnr_los"]:
        line += " | %5.1f -> %5.1f dB (%4.1f dB); %5.1f -> %5.1f dB (%4.1f dB)" % (np.median(n["cnr_los"]), np.median(a["cnr_los"]), np.median(n["cnr_los"]) - np.median(a["cnr_los"]),
                                                                              np.median(n["cnr_all"]), np.median(a["cnr_all"]), np.median(n["cnr_all"]) - np.median(a["cnr_all"]))
    print(line)
