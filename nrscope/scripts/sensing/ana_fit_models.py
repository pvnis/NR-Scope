"""Real 18 Nov chunks: how well does DM-RS = TRS x (model) fit, per PRG of the tablet? Models: one complex number shared by
the 4 antennas; one per antenna; per antenna + linear slope across the PRG; per antenna + quadratic. Relative residual
||Y - fit|| / ||Y|| per PRG (median). Noise floor from the two TRS symbols of a slot (4 and 8). Also by TRS-DM-RS time gap.
usage: fit_models.py fig_b*.npz"""
import sys, numpy as np, collections
xs = (12 * (np.arange(12) // 3) + 4 * (np.arange(12) % 3)).astype(float); xs = (xs - xs.mean()) / (xs.max() - xs.mean())
V = {k: np.vander(xs, k + 1, increasing=True) for k in (0, 1, 2)}          # [12, k+1]
def fit_res(Y, R, k, shared=False):
    if shared:   # one number for all antennas
        c = np.vdot(R.ravel(), Y.ravel()) / np.vdot(R.ravel(), R.ravel()).real; E = Y - c * R
    else:
        E = np.empty_like(Y)
        for q in range(Y.shape[0]):
            A = R[q][:, None] * V[k]; sol, *_ = np.linalg.lstsq(A, Y[q], rcond=None); E[q] = Y[q] - A @ sol
    return np.linalg.norm(E) / np.linalg.norm(Y)
ALL = collections.defaultdict(list)
for f in sys.argv[1:]:
    Z = np.load(f); tt, tl, th, ta = Z["trs_t"], Z["trs_l"], Z["trs_h"], Z["trs_a"]
    # noise floor: TRS symbols 4 and 8 of the same slot, per PRG and antenna
    nz = []
    for a in np.unique(ta):
        i4 = np.flatnonzero((ta == a) & (tl == 4)); i8 = np.flatnonzero((ta == a) & (tl == 8))
        if len(i4) and len(i8):
            A, B = th[i4[0]], th[i8[0]]
            for b in range(0, A.shape[1] // 12, 4):
                sl = slice(12 * b, 12 * b + 12); nz.append(np.linalg.norm(A[:, sl] - B[:, sl]) / np.sqrt(2) / np.linalg.norm(A[:, sl]))
    u = Z["d_u"]; tab = collections.Counter(u[u != 0].tolist()).most_common(1)[0][0]
    res = collections.defaultdict(list); n = 0
    for t, uu, pm, y in zip(Z["d_t"], u, Z["d_pm"], Z["d_y"]):
        if uu != tab: continue
        j = int(np.argmin(np.abs(tt - t))); dt = abs(tt[j] - t)
        if dt > 0.01: continue
        R = th[j]; gb = "<=2ms" if dt <= 0.002 else "2-10ms"
        for b in np.flatnonzero(pm)[::2]:
            sl = slice(12 * b, 12 * b + 12); Yb, Rb = y[:, sl], R[:, sl]
            if not np.all(np.abs(Rb) > 0) or np.linalg.norm(Yb) == 0: continue
            for name, k, sh in (("shared", 0, True), ("per antenna", 0, False), ("per antenna + slope", 1, False), ("per antenna + quadratic", 2, False)):
                r = fit_res(Yb, Rb, k, sh); res[name].append(r); res[name + " " + gb].append(r)
            n += 1
    nfl = np.median(nz)
    print("%s: tablet 0x%x, %d PRGs; noise floor (TRS 4 vs 8) %.1f %%" % (f.split("/")[-1], tab, n, 100 * nfl))
    for name in ("shared", "per antenna", "per antenna + slope", "per antenna + quadratic"):
        k = {"shared": 0, "per antenna": 0, "per antenna + slope": 1, "per antenna + quadratic": 2}[name]
        print("   %-24s residual %5.1f %%   (noise alone would leave ~%.1f %%)   <=2 ms: %5.1f %%   2-10 ms: %5.1f %%" % (
            name, 100 * np.median(res[name]), 100 * nfl * np.sqrt((12 - k - 1) / 12) * np.sqrt(2),
            100 * np.median(res[name + " <=2ms"]) if res[name + " <=2ms"] else np.nan, 100 * np.median(res[name + " 2-10ms"]) if res[name + " 2-10ms"] else np.nan))
        ALL[name] += res[name]
    ALL["noise"] += nz
print("ALL chunks: noise floor %.1f %%; " % (100 * np.median(ALL["noise"])) + "; ".join("%s %.1f %%" % (k, 100 * np.median(v)) for k, v in ALL.items() if k != "noise"))
