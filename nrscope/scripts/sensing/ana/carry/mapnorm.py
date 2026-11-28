import csv, json, sys, numpy as np
W0, W1 = float(sys.argv[1]), float(sys.argv[2]); runs = sys.argv[3:]; lam = 299792458 / 3.75e9
def truth(fn):
    r = list(csv.DictReader(open(fn))); los = float([x for x in r if x["path"] == "0"][0]["bistatic_range_m"])
    w = [x for x in r if x["via_car"] == "1"]; t = np.array([float(x["t_s"]) for x in w])
    return t, np.array([float(x["bistatic_range_m"]) - los for x in w]), np.array([float(x["doppler_hz"]) for x in w])
for run in runs:
    st = json.load(open(run + "/settings.json")); t0 = float(st["t0_unix"]); W = (W0 - t0, W1 - t0)
    for name, mf, tf in (("TRS", "map2d.csv", "ground_truth.csv"), ("DMRS", "map2d_rx2.csv", "ground_truth_rx2.csv")):
        if name == "TRS" and run != runs[0]: continue
        tt, ex, fd = truth(run + "/" + tf); unix = [float(r["unix_s"]) for r in csv.DictReader(open(run + "/" + mf + ".time.csv"))]
        res = {k: [] for k in ("true", "reversed", "+15 m", "+30 m", "shift")}; resm = {k: [] for k in res}
        for l, u in zip([l for l in open(run + "/" + mf) if not l.startswith("#")], unix):
            v = l.strip().split(","); nb, nf = int(v[12]), int(v[13]); mpb, fmax = float(v[9]), float(v[10]); t = u - t0 - float(v[8]) / 2
            if not (W[0] <= t <= W[1]): continue
            P = np.array(v[21:21 + nb * nf], float).reshape(nb, nf); r4 = int(np.ceil(4 / mpb))
            Z = P / (np.median(P[r4:]) + 1e-12)
            def val(e, sp):
                ri = int(round(e / mpb)); fi = int(round(((2 * sp / lam) + fmax) / (2 * fmax) * (nf - 1)))
                df = max(1, int(round(0.25 / (lam / 2) / (2 * fmax / (nf - 1)))))   # +-0.25 m/s in bins
                if not (r4 <= ri < nb - 1 and df <= fi < nf - df): return np.nan
                return Z[ri - 1:ri + 2, fi - df:fi + df + 1].max()
            k = np.argmin(abs(tt - t)); e, sp = ex[k], fd[k] * lam / 2
            d = res if True else None
            for D, ok in ((res, True), (resm, abs(sp) >= 0.6)):
                if not ok: continue
                D["true"].append(val(e, sp)); D["reversed"].append(val(e, -sp)); D["+15 m"].append(val(e + 15, sp)); D["+30 m"].append(val(e + 30, sp))
                for s in (8, 16, 24, 32, 40):
                    ts = W[0] + ((t + s - W[0]) % (W[1] - W[0])); k2 = np.argmin(abs(tt - ts)); D["shift"].append(val(ex[k2], fd[k2] * lam / 2))
        for lab, D in (("all", res), ("|v|>=0.6", resm)):
            print("%-22s %-4s %-8s n %3d | " % (run.split("/")[-1], name, lab, len(D["true"])) + " | ".join("%s %.1f dB" % (k, 10 * np.log10(np.nanmedian(x))) for k, x in D.items()))
