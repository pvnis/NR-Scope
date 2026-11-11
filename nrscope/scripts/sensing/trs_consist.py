"""Is the TRS a fixed reference? Compare snapshot pairs of the TRS and of the tablet's DM-RS at the same time lag, after
the same corrections as the pipeline (window moves undone, drift removed, range zero on the direct path).
Per 4-PRB group b of a pair (X1, X2): z_b = <X1_b, X2_b> over 4 antennas x 12 pilots.
 - shape  rho_b = |z_b| / (|X1_b| |X2_b|): same channel shape on the group (1 = identical up to one number)
 - phase  delta_b = angle(z_b) - common phase of the pair: a precoder change shows up as a different phase per group
 - coherence |sum z_b| / sum |z_b|: 1 when every group changed by the same phase (only the channel moved)
usage: trs_consist.py <burst dir> <out.json>"""
import sys, os, json, numpy as np, collections as C
here = "/home/sens/field/"; BD, OUTJ = sys.argv[1], sys.argv[2]
src = open(here + "so_light.py").read().split("# CSI-RS as extra TRS-beam samples")[0]
sys.argv = ["so_light.py", BD, "/dev/null"]; os.environ.pop("INJ", None); __file__ = here + "so_light.py"; exec(src)
def metr(X1, X2, pm):
    out = []
    for b in np.flatnonzero(pm):
        sl = slice(12 * b, 12 * b + 12); z = np.vdot(X1[:, sl], X2[:, sl]); n = np.linalg.norm(X1[:, sl]) * np.linalg.norm(X2[:, sl])
        if n > 0: out.append((z, abs(z) / n))
    if len(out) < 4: return None
    z = np.array([o[0] for o in out]); rho = np.array([o[1] for o in out]); ph = np.angle(np.sum(z))
    d = np.degrees(np.abs(np.angle(z * np.exp(-1j * ph))))
    return float(np.median(rho)), float(np.median(d)), float(abs(z.sum()) / np.abs(z).sum())
res = C.defaultdict(list)
T = {(a, l): R for t, a, l, R in TRS}
for (a, l), R in T.items():
    for lag, name in ((40, "20ms"), (1, "0.5ms"), (80, "40ms")):
        if (a + lag, l) in T:
            m = metr(R, T[(a + lag, l)], np.ones(NG, bool))
            if m: res["trs_" + name].append(m + (int(shifts(a) != shifts(a + lag)),))
Dd = [load(f) for f in files]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
cnt_ = C.Counter(k[1] for k in keys if k[1] != 0); tab = cnt_.most_common(1)[0][0] if cnt_ else None
if tab is not None:
    S = {}
    for k in keys:
        a, u, l = k
        if u != tab: continue
        pm = np.logical_and.reduce([d[k][0] for d in Dd]); t = a * 0.5e-3 + l / 28e3
        Y = np.array([d[k][1] for d in Dd]).astype(complex); S[(a, l)] = (pm, dedrift(t, undo_wsh(Y, shifts(a), -1.0)) * ramp)
    del Dd
    for (a, l), (pm, Y) in S.items():
        for lag, name in ((40, "20ms"), (1, "0.5ms"), (80, "40ms")):
            if (a + lag, l) in S:
                pm2, Y2 = S[(a + lag, l)]; m = metr(Y, Y2, pm & pm2)
                if m: res["dmrs_" + name].append(m + (int(shifts(a) != shifts(a + lag)),))
summ = {k: {"n": len(v), "rho": float(np.median([x[0] for x in v])), "phase_deg": float(np.median([x[1] for x in v])),
            "coh": float(np.median([x[2] for x in v])), "coh_p10": float(np.percentile([x[2] for x in v], 10))} for k, v in res.items()}
json.dump(summ, open(OUTJ, "w"), indent=1); print(json.dumps(summ))
json.dump({k: v for k, v in res.items()}, open(OUTJ + ".raw", "w"))
for k in ("trs_20ms", "dmrs_20ms", "trs_0.5ms", "dmrs_0.5ms"):
    v = np.array(res[k]); mv = v[:, 3] == 1
    print("%-10s n %4d | coh<0.9 %.0f%% | phase>20deg %.0f%% | with a window move between: n %d coh med %.3f, without: coh med %.3f, coh<0.9 %.0f%%" % (
        k, len(v), 100 * np.mean(v[:, 2] < 0.9), 100 * np.mean(v[:, 1] > 20), mv.sum(), np.median(v[mv, 2]) if mv.any() else np.nan,
        np.median(v[~mv, 2]), 100 * np.mean(v[~mv, 2] < 0.9)))
