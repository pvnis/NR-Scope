import sys, csv, glob, os, collections as C, numpy as np
# per target: power at the expected cell (max over +-1 bin, +-2 Doppler) over the map median, in dB; maps of a glob
def maps(p):
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split(","); nb, nf = int(f[12]), int(f[13])
        P = np.array(f[21:21 + nb * nf], float).reshape(nb, nf)
        yield int(f[0]), int(f[1]), float(f[9]), float(f[10]), float(f[11]), float(f[16]), P
for v in sys.argv[2:]:
    R = C.defaultdict(list)
    for d in sorted(glob.glob(sys.argv[1] % v)):
        p = d + "/map2d.csv"
        if not os.path.exists(p + ".truth.csv"): continue
        TR = C.defaultdict(list)
        for r in csv.DictReader(open(p + ".truth.csv")): TR[(int(r["frame"]), int(r["slot"]))].append((int(r["target"]), float(r["range_m"]), float(r["speed_ms"])))
        for fr, sl, mpb, fmax, fc, blos, P in maps(p):
            nb, nf = P.shape; med = np.median(P)
            for k, rr, vv in TR.get((fr, sl), []):
                b = int(round(blos + rr / mpb)); fq = 2 * vv / (3e8 / fc); j = int(round((fq + fmax) / (2 * fmax) * (nf - 1)))
                if not (1 <= b < nb - 1 and 2 <= j < nf - 2): continue
                R[k].append(10 * np.log10(P[b - 1:b + 2, j - 2:j + 3].max() / med))
    print(v, "  ".join("t%d med %.1f dB (n %d)" % (k, np.median(R[k]), len(R[k])) for k in sorted(R)))
