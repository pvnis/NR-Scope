import csv, json, sys, numpy as np
run, W0, W1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]); lam = 299792458 / 3.75e9
st = json.load(open(run + "/settings.json")); t0 = float(st["t0_unix"]); W = (W0 - t0, W1 - t0)
for name, mf, tf in (("TRS", "map2d.csv", "ground_truth.csv"), ("DMRS", "map2d_rx2.csv", "ground_truth_rx2.csv")):
    r = list(csv.DictReader(open(run + "/" + tf))); los = float([x for x in r if x["path"] == "0"][0]["bistatic_range_m"])
    w = [x for x in r if x["via_car"] == "1"]; tt = np.array([float(x["t_s"]) for x in w]); ex = np.array([float(x["bistatic_range_m"]) - los for x in w]); fd = np.array([float(x["doppler_hz"]) for x in w])
    unix = [float(q["unix_s"]) for q in csv.DictReader(open(run + "/" + mf + ".time.csv"))]; L = []; B = []; C = []
    Wt = (max(W[0], tt.min()), min(W[1], tt.max()))
    for l, u in zip([l for l in open(run + "/" + mf) if not l.startswith("#")], unix):
        v = l.strip().split(","); nb, nf = int(v[12]), int(v[13]); mpb, fmax = float(v[9]), float(v[10]); t = u - t0 - float(v[8]) / 2
        if not W[0] <= t <= W[1]: continue
        P = np.array(v[21:21 + nb * nf], float).reshape(nb, nf); k = np.argmin(abs(tt - t)); e, sp = ex[k], fd[k] * lam / 2
        ri = int(round(e / mpb)); fi = int(round(((2 * sp / lam) + fmax) / (2 * fmax) * (nf - 1))); df = max(1, int(round(0.25 / (lam / 2) / (2 * fmax / (nf - 1)))))
        if abs(sp) < 0.6 or not (2 <= ri < nb - 1 and df <= fi < nf - df): continue
        direct = P[:2, :].max()                                  # direct-path peak (range 0-1 bins, any Doppler)
        cell = P[ri - 1:ri + 2, fi - df:fi + df + 1].max(); mir = P[ri - 1:ri + 2, nf - 1 - fi - df:nf - fi + df].max()
        L.append(10 * np.log10(cell / direct)); B.append(10 * np.log10(np.median(P[ri - 1:ri + 2, :]) / direct))
        for sft in (7, 14, 21, 28):   # control: the same neighbourhood statistic at the truth of another time
            ts = Wt[0] + ((t + sft - Wt[0]) % (Wt[1] - Wt[0])); k2 = np.argmin(abs(tt - ts)); e2, s2 = ex[k2], fd[k2] * lam / 2
            r2 = int(round(e2 / mpb)); f2 = int(round(((2 * s2 / lam) + fmax) / (2 * fmax) * (nf - 1)))
            if 2 <= r2 < nb - 1 and df <= f2 < nf - df: C.append(10 * np.log10(P[r2 - 1:r2 + 2, f2 - df:f2 + df + 1].max() / np.median(P[r2 - 1:r2 + 2, :])))
    L = np.array(L); B = np.array(B)
    print("%-4s moving maps %d: walker above its range background %.1f dB | control (time-shifted truth) %.1f dB | excess %.1f dB" % (name, len(L), np.median(L - B), np.median(C), np.median(L - B) - np.median(C)))
