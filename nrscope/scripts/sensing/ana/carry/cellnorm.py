import csv, json, sys, numpy as np
run, W0, W1, vmin = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
st = json.load(open(run + "/settings.json")); t0 = float(st["t0_unix"]); W = (W0 - t0, W1 - t0); lam = 299792458 / 3.75e9
def truth(fn):
    r = list(csv.DictReader(open(fn))); los = float([x for x in r if x["path"] == "0"][0]["bistatic_range_m"])
    w = [x for x in r if x["via_car"] == "1"]; t = np.array([float(x["t_s"]) for x in w])
    return t, np.array([float(x["bistatic_range_m"]) - los for x in w]), np.array([float(x["doppler_hz"]) for x in w])
def maps(fn):
    unix = [float(r["unix_s"]) for r in csv.DictReader(open(fn + ".time.csv"))]; T, P = [], []
    for l, u in zip([l for l in open(fn) if not l.startswith("#")], unix):
        v = l.strip().split(","); nb, nf = int(v[12]), int(v[13]); mpb, fmax = float(v[9]), float(v[10])
        T.append(u - t0 - float(v[8]) / 2); P.append((np.array(v[21:21 + nb * nf], float).reshape(nb, nf), fmax))
    from collections import Counter
    shp = Counter(p.shape for p, _ in P).most_common(1)[0][0]; keep = [i for i, (p, _) in enumerate(P) if p.shape == shp]
    print("  %s: %d maps, %d of the main shape %s" % (fn.split("/")[-1], len(P), len(keep), shp))
    return np.array(T)[keep], np.array([P[i][0] for i in keep]), mpb, P[keep[0]][1]
for name, mf, tf in (("TRS", "map2d.csv", "ground_truth.csv"), ("DMRS", "map2d_rx2.csv", "ground_truth_rx2.csv")):
    tt, ex, fd = truth(run + "/" + tf); T, P, mpb, fmax = maps(run + "/" + mf); nb, nf = P.shape[1:]
    N = np.median(P, 0) + 1e-12; Z = P / N
    def val(i, e, v):
        ri = int(round(e / mpb)); fi = int(round(((2 * v / lam) + fmax) / (2 * fmax) * (nf - 1)))
        if not (1 <= ri < nb - 1 and 1 <= fi < nf - 1): return np.nan
        return Z[i, ri - 1:ri + 2, fi - 1:fi + 2].max()
    sel = [i for i in range(len(T)) if W[0] <= T[i] <= W[1]]
    res = {"true": [], "reversed": [], "+15 m": [], "shift": []}
    for i in sel:
        k = np.argmin(abs(tt - T[i])); e, v = ex[k], fd[k] * lam / 2
        if abs(v) < vmin: continue
        res["true"].append(val(i, e, v)); res["reversed"].append(val(i, e, -v)); res["+15 m"].append(val(i, e + 15, v))
        for s in (8, 16, 24, 32, 40):
            ts = W[0] + ((T[i] + s - W[0]) % (W[1] - W[0])); k2 = np.argmin(abs(tt - ts)); res["shift"].append(val(i, ex[k2], fd[k2] * lam / 2))
    out = " | ".join("%s %.2f dB" % (k, 10 * np.log10(np.nanmean(v))) for k, v in res.items())
    print("%-4s maps %3d (|v|>=%.1f): mean cell-normalised power: %s" % (name, len(res["true"]), vmin, out))
