#!/usr/bin/env python3
"""Detection statistics per range x speed cell, at one false-alarm rate: chain 0's power (NR-Scope's detector),
the best whitened 4-chain beam over angle, the whitened energy D^H R^-1 D (no steering), the plain 4-chain beam.
Each normalised per range bin (median over moving cells); threshold = the pfa quantile of all moving cells;
Pd = share of in-view maps where a cell within +-1 range cell and +-speed_tol of the walker's truth crosses it.

  python3 detstat.py --config ... --caps CAPS --run RUN_DIR [--stream trs] [--pfa 1e-3] [--set ...]"""
import argparse, csv, json, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C, calib, cube, est_read


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--config"); ap.add_argument("--caps", required=True); ap.add_argument("--run", required=True)
    ap.add_argument("--stream", default="trs"); ap.add_argument("--pfa", type=float, default=1e-3); ap.add_argument("--speed-tol", type=float, default=0.6)
    ap.add_argument("--set", nargs="+", action="extend", default=[]); a = ap.parse_args(); c = C.load(a.config, a.set); T_ = c["tbd"]
    cal = calib.get(c)[0]; r, v, s = cube.grids(T_); cap = os.path.join(C.recording_dir(c), "cap_%s_%s" % (a.caps, a.stream)); L, U = cube.rows(cap); rec = c["recording"]
    a0s = {int(re.search(r"b(\d+)", f).group(1)): json.load(open(os.path.join(rec["lidar"], "grid", f)))["a0"] for f in os.listdir(os.path.join(rec["lidar"], "grid")) if re.match(r"b\d+\.json$", f)}
    by = {}
    for i, (l, u) in enumerate(zip(L, U)):
        k, an, span = cube.chunk_of(c, l, u, a0s)
        if k is not None:
            by.setdefault(k, []).append((i, an, span))
    mv = np.abs(v) >= c["streams"][a.stream]["target_dv"]; S = np.exp(-1j * np.pi * np.arange(4)[:, None] * s[None, :]); stats = {}
    G, LOS = {}, {}
    for q in csv.DictReader(open(os.path.join(a.run, "ground_truth.csv"))):
        t = float(q["t_s"])
        if q["via_car"] == "0": LOS[t] = float(q["bistatic_range_m"])
        elif q["via_car"] == "1": G[t] = (float(q["bistatic_range_m"]), float(q["aoa_deg"]), float(q["doppler_hz"]))
    tg = np.array(sorted(G)); t0 = float(json.load(open(os.path.join(a.run, "settings.json")))["t0_unix"])
    for k, items in sorted(by.items()):
        T, X = est_read.snapshots(os.path.join(C.est_dir(c, a.stream), "b%d.bin" % k))
        for i, an, span in items:
            D = cube.transform(X, T, an, span, r, v, cal)
            if D is None:
                continue
            R = cube.noise_cov(D, mv); Rn = R / (np.trace(R).real / 4) + 0.01 * np.eye(4); Ri = np.linalg.inv(Rn)
            W = Ri @ S.conj(); W = W / np.sqrt(np.real(np.sum(S * W, 0)))[None, :]
            st = {"chain0": np.abs(D[0]) ** 2, "beam4": (np.abs(np.einsum("arv,ag->rvg", D, S)) ** 2).max(2),
                  "white_beam": (np.abs(np.einsum("arv,ag->rvg", D, W.conj())) ** 2).max(2), "white_energy": np.real(np.einsum("arv,ab,brv->rv", D.conj(), Ri, D))}
            tt = float(U[i]) - t0 - span / 2; j = np.argmin(abs(tg - tt))
            truth = None
            if abs(tg[j] - tt) < 0.3 and tg[j] in LOS:
                b_, ang, fd = G[tg[j]]; ex = b_ - LOS[tg[j]]
                if abs(ang) <= 60 and r[0] <= ex <= r[-1]:
                    truth = (ex, fd * (299792458.0 / 3.75e9) / 2)
            for name, P in st.items():
                n = np.median(P[:, mv], 1); Z = P / np.maximum(n, 1e-30)[:, None]; Z[:, ~mv] = 0
                d = stats.setdefault(name, dict(noise=[], walker=[], walker_mirror=[]))
                d["noise"].append(Z[:, mv].ravel()[::7])
                if truth:
                    ir = np.argmin(abs(r - truth[0])); w = (slice(max(0, ir - 1), ir + 2))
                    for key, vv in (("walker", truth[1]), ("walker_mirror", -truth[1])):   # control: the opposite speed
                        sv = mv & (np.abs(v - vv) <= a.speed_tol)
                        d[key].append(Z[w][:, sv].max() if sv.any() else 0.0)
    out = {}
    for name, d in stats.items():
        th = np.quantile(np.concatenate(d["noise"]), 1 - a.pfa); w = np.array(d["walker"]); wm = np.array(d["walker_mirror"])
        out[name] = dict(maps=len(w), pd=round(float(np.mean(w > th)), 3), pd_control=round(float(np.mean(wm > th)), 3), walker_z_db=round(float(10 * np.log10(np.median(w))), 1))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
