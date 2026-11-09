"""Gluing precoder epochs in time (one unknown complex number per epoch), on the tablet.
Series: per (antenna, PRG) the mean of the epoch-stitched DM-RS over the PRG's 12 samples, at each snapshot time.
Controlled: a long epoch (>= 40 ms) cut at its middle, second half times a random number, optional time gap g removed
around the cut. Methods for the number:
  anchor   nearest snapshot pair across the cut, one number over shared PRGs and antennas (what gs_all does for UEs)
  mean     ratio of the time-averaged series (static channel assumption), shared PRGs
  doppler  global phase search, simultaneous OMP over a Doppler grid (-300..300 Hz, 1 Hz) with P tones shared by all
           series, smallest residual wins, then LS amplitude
Real: consecutive real epochs (gap <= 20 ms): the three methods, their agreement, and the fit residual of the static model.
usage: gs_omp_t.py <rec dir> <out.json> [bursts]   env CH (0,1,2,3)"""
import sys, os, json, numpy as np, collections as C
REC, OUTF = sys.argv[1], sys.argv[2]; ONLY = [int(x) for x in sys.argv[3:]]; saved = list(sys.argv)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", REC, "/dev/null"]; exec(src); sys.argv = saved
rng = np.random.default_rng(5); FD = np.arange(-300, 301, 1.0)

def series(XA, idx):
    """-> times (n), prg mask (NG), V: n x (NA*NG) complex (nan where absent)"""
    t = np.array([XA[i][0] for i in idx]); V = np.full((len(idx), NA, NG), np.nan + 0j)
    for j, i in enumerate(idx):
        _, _, pm, X = XA[i]; v = X.reshape(NA, NG, 12).mean(2); V[j][:, pm] = v[:, pm]
    return t, V

def fit_c(t1, V1, t2, V2, method, P=4):
    ok = np.isfinite(V1).mean(0) >= 0.6
    ok &= np.isfinite(V2).mean(0) >= 0.6
    if ok.sum() < 8: return np.nan, np.nan
    if method == "anchor":
        a = V1[-1]; b = V2[0]; m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 8: return np.nan, np.nan
        return np.vdot(a[m], b[m]) / np.vdot(a[m], a[m]), np.nan
    if method == "mean":
        a = np.nanmean(V1, 0)[ok]; b = np.nanmean(V2, 0)[ok]
        c = np.vdot(a, b) / np.vdot(a, a)
        res = (np.nansum(np.abs(V1[:, ok] - a) ** 2) + np.nansum(np.abs(V2[:, ok] / c - a) ** 2)) / (np.nansum(np.abs(V1[:, ok]) ** 2) + np.nansum(np.abs(V2[:, ok] / c) ** 2))
        return c, float(res)
    # doppler: columns (ant, PRG) present at every snapshot of both halves; P tones shared by all columns
    ok = (np.isfinite(V1).mean(0) >= 0.8) & (np.isfinite(V2).mean(0) >= 0.8)
    if ok.sum() < 8: return np.nan, np.nan
    def fillc(t_, V_):
        V_ = V_[:, ok].copy()
        for j in range(V_.shape[1]):
            m_ = np.isfinite(V_[:, j])
            if not m_.all(): V_[~m_, j] = np.interp(t_[~m_], t_[m_], V_[m_, j].real) + 1j * np.interp(t_[~m_], t_[m_], V_[m_, j].imag)
        return V_
    W1 = fillc(t1, V1); W2 = fillc(t2, V2); amp0 = np.sqrt(np.mean(np.abs(W2) ** 2) / np.mean(np.abs(W1) ** 2)); W2 = W2 / amp0
    t = np.r_[t1, t2] - t1[0]; best = None; Aall = np.exp(2j * np.pi * np.outer(t, FD))
    for ph in np.linspace(-np.pi, np.pi, 72, endpoint=False):
        Vz = np.r_[W1, W2 * np.exp(-1j * ph)]; sel = []; R = Vz.copy()
        for _ in range(P):
            k = int(np.argmax(np.sum(np.abs(Aall.conj().T @ R) ** 2, 1))); sel.append(k)
            A = Aall[:, sel]; x, *_ = np.linalg.lstsq(A, Vz, rcond=None); Mod = A @ x; R = Vz - Mod
        r = np.sum(np.abs(R) ** 2) / np.sum(np.abs(Vz) ** 2)
        if best is None or r < best[0]: best = (r, ph, Mod)
    r, ph, Mod = best
    return amp0 * np.exp(1j * ph), float(r)

def err(est, true):
    if not np.isfinite(est): return None
    r = est / true; return abs(float(np.degrees(np.angle(r)))), abs(float(20 * np.log10(abs(r))))

BURSTS = sorted(glob.glob(REC + "/b*/"), key=lambda s: int(s.rstrip("/").split("/b")[-1]))
Rz = C.defaultdict(list)
for bi, bd in enumerate(BURSTS):
    if ONLY and bi not in ONLY: continue
    fs = [bd + "replay_ch%d/ue_dmrs.bin" % c for c in CH]
    if not all(os.path.exists(f) and os.path.getsize(f) > 20e6 for f in fs): continue
    Dd = [load(f) for f in fs]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
    S = C.defaultdict(list)
    for k in sorted(keys):
        pm = np.logical_and.reduce([Dd[c][k][0] for c in range(NA)]); Y = np.array([Dd[c][k][1] for c in range(NA)]).astype(complex)
        if pm.sum() and np.all(np.isfinite(Y[:, np.repeat(pm, 12)])): S[k[1]].append((k[0] * 0.5e-3 + k[2] / 28e3, k[0], pm, Y))
    del Dd
    UA = max(S, key=lambda u: len(S[u])); A = sorted(S[UA], key=lambda x: x[0])
    # per-subframe phase step: consecutive snapshots crossing exactly one subframe boundary vs none
    st1 = []
    for (t1_, a1, p1, Y1), (t2_, a2, p2, Y2) in zip(A[:-1], A[1:]):
        m = np.repeat(p1 & p2, 12)
        if m.sum() < 120 or a2 - a1 > 2: continue
        nb = a2 // 2 - a1 // 2
        if nb == 1: st1.append(np.vdot(Y1[:, m], Y2[:, m]) / abs(np.vdot(Y1[:, m], Y2[:, m])))
    phs = float(np.angle(np.mean(st1))) if st1 else 0.0
    if os.environ.get("NOSF") != "1":
        A = [(t, a, pm, Y * np.exp(-1j * phs * (a // 2))) for t, a, pm, Y in A]
    Rz[("sf step deg",)].append(float(np.degrees(phs)))
    EP, TPL = epochs(A); EG = {k: joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0) for k, (S2, C2) in TPL.items()}
    XA = []
    for i, (t, a, pm, Y) in enumerate(A):
        X = Y / np.repeat(EG[EP[i]], 12)[None, :]; X[:, ~np.repeat(pm, 12)] = 0; XA.append((t, a, pm, X))
    byep = C.defaultdict(list)
    for i, e in enumerate(EP): byep[e].append(i)
    ne = 0
    # ---- controlled split of long epochs
    for e, idx in byep.items():
        T = np.array([XA[i][0] for i in idx])
        if T[-1] - T[0] < 0.04 or len(idx) < 30: continue
        ne += 1
        for g in (0.0, 0.005, 0.02):
            tm = (T[0] + T[-1]) / 2; i1 = [i for i in idx if XA[i][0] < tm - g / 2]; i2 = [i for i in idx if XA[i][0] > tm + g / 2]
            if len(i1) < 8 or len(i2) < 8: continue
            cc = rng.uniform(0.5, 2) * np.exp(1j * rng.uniform(-np.pi, np.pi))
            t1, V1 = series(XA, i1); t2, V2 = series(XA, i2); V2 = V2.reshape(len(i2), -1) * cc; V1 = V1.reshape(len(i1), -1)
            for meth in ("anchor", "mean", "doppler"):
                c, res = fit_c(t1, V1, t2, V2, meth)
                ee = err(c, cc)
                if ee: Rz[("ctl", g, meth, "deg")].append(ee[0]); Rz[("ctl", g, meth, "db")].append(ee[1])
                if np.isfinite(res): Rz[("ctl", g, meth, "res")].append(res)
    # ---- real consecutive epochs
    eps = sorted(byep)
    for e1, e2 in zip(eps[:-1], eps[1:]):
        i1, i2 = byep[e1], byep[e2]
        if len(i1) < 6 or len(i2) < 6: continue
        gap = XA[i2[0]][0] - XA[i1[-1]][0]
        if gap > 0.02: Rz[("real", "gap > 20 ms")].append(1); continue
        t1, V1 = series(XA, i1); t2, V2 = series(XA, i2); V1 = V1.reshape(len(i1), -1); V2 = V2.reshape(len(i2), -1)
        cs = {}
        for meth in ("anchor", "mean", "doppler"):
            c, res = fit_c(t1, V1, t2, V2, meth); cs[meth] = c
            if np.isfinite(res): Rz[("real", meth, "res")].append(res)
        if all(np.isfinite(list(cs.values()))):
            Rz[("real", "anchor vs mean", "deg")].append(err(cs["anchor"], cs["mean"])[0])
            Rz[("real", "doppler vs mean", "deg")].append(err(cs["doppler"], cs["mean"])[0])
            Rz[("real", "doppler vs anchor", "deg")].append(err(cs["doppler"], cs["anchor"])[0])
        Rz[("real", "gap ms")].append(gap * 1e3)
    print("burst %d: %d long epochs, %d epochs" % (bi, ne, len(eps)), flush=True)
json.dump({"|".join(map(str, k)): v for k, v in Rz.items()}, open(OUTF, "w"))
def q(k):
    v = np.array(Rz.get(k, []), float); return "%5.1f / %5.1f (n %d)" % (np.median(v), np.percentile(v, 90), len(v)) if len(v) else "-"
for g in (0.0, 0.005, 0.02):
    print("ctl gap %g ms:" % (g * 1e3), " | ".join("%s %s deg, %s dB" % (m, q(("ctl", g, m, "deg")), q(("ctl", g, m, "db"))) for m in ("anchor", "mean", "doppler")))
    print("   residual: mean", q(("ctl", g, "mean", "res")), "doppler", q(("ctl", g, "doppler", "res")))
print("subframe step per burst (deg):", [round(x, 1) for x in Rz[("sf step deg",)]])
print("real: gap ms", q(("real", "gap ms")), "skipped >20 ms", len(Rz[("real", "gap > 20 ms")]))
print("real residual: mean", q(("real", "mean", "res")), "doppler", q(("real", "doppler", "res")))
print("real agreement: anchor vs mean", q(("real", "anchor vs mean", "deg")), "doppler vs mean", q(("real", "doppler vs mean", "deg")), "doppler vs anchor", q(("real", "doppler vs anchor", "deg")))
