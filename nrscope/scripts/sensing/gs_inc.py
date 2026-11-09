"""Coherent merge vs incoherent sum of range-Doppler maps, per 150 ms window, on the tablet bursts.
Truth T: the tablet's own full-band snapshots (>= 60 of 68 PRGs), subframe step removed, epochs stitched in frequency
and glued in time by the nearest-snapshot anchor -> one coherent map over the full band.
Part 1, controlled (truth known): the same snapshots cut into pieces that act as different "UEs":
  split2   two touching halves       split4  four touching quarters     gap  30 PRG | 8 PRG hole | 30 PRG
  every piece but the first gets an unknown complex number (constant in the window) and a beam mismatch:
  a per-PRG complex perturbation of power sigma2 (0 = exact one-number model, 0.15 ~ measured between UEs, 0.5 ~ per-PRG UEs)
  methods: oracle (true number, mismatch left), edge (number per snapshot from the boundary, touching cases only),
           raw (no correction), inc (one map per piece, powers added), single (first piece only)
Part 2, real UEs: each other UE's snapshots in the window, self-stitched, split by its own epochs, one map per
  (UE, epoch) segment, powers added (inc_real). Baseline inc_tab: the tablet's data on exactly the same (slot, PRB)
  support and segments (same holes, same short times, but the tablet's beam, from its nearest snapshot within 1 ms). Compared to T.
Metrics against T (power summed over 4 antennas, +-64 range bins around T's peak, Doppler -300..300 Hz):
  corr      correlation of dB maps clipped 30 dB under the peak
  dyn       same on |nu| >= 15 Hz only (moving things)
  recall    of T's 10 strongest moving peaks (within 20 dB), a test peak within +-2 range bins and +-4 Hz
  false     share of the test map's 10 strongest moving peaks that are not near a T peak
  width     range width of the main peak (m, -3 dB)
usage: gs_inc.py <rec dir> <out.json> [bursts]   env CH (0,3 for other UEs on 5 Nov), PART (12, 1 or 2)"""
import sys, os, json, numpy as np, collections as C
REC, OUTF = sys.argv[1], sys.argv[2]; ONLY = [int(x) for x in sys.argv[3:]]; saved = list(sys.argv)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", REC, "/dev/null"]; exec(src); sys.argv = saved
rng = np.random.default_rng(11)
FD = np.arange(-300, 301, 2.0); DYN = np.abs(FD) >= 15; HW = 64; WIN = 0.15

def sfcorr(A):
    st = []
    for (t1, a1, p1, Y1), (t2, a2, p2, Y2) in zip(A[:-1], A[1:]):
        m = np.repeat(p1 & p2, 12)
        if m.sum() < 120 or a2 - a1 > 2: continue
        if a2 // 2 - a1 // 2 == 1:
            v = np.vdot(Y1[:, m], Y2[:, m]); st.append(v / abs(v))
    return float(np.angle(np.mean(st))) if st else 0.0

def stitch(A):
    EP, TPL = epochs(A); EG = {k: joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0) for k, (S2, C2) in TPL.items()}
    X0 = []
    for i, (t, a, pm, Y) in enumerate(A):
        Z = Y / np.repeat(EG[EP[i]], 12)[None, :]; Z[:, ~np.repeat(pm, 12)] = 0; X0.append(Z)
    # glue: one factor per epoch, chained through the last snapshot of the previous (already glued) epoch
    X = [X0[0]]; g = 1 + 0j
    for i in range(1, len(A)):
        if EP[i] != EP[i - 1]:
            f, n = frame(X[i - 1], A[i - 1][2], X0[i], A[i][2])
            g = f if (np.isfinite(f) and n >= 4 and f != 0) else 1 + 0j
        X.append(X0[i] / g)
    return X, EP

def rd(ts, Xs, masks, t0):
    """ts: times, Xs: NA x L arrays (zero outside), masks: L bool support for the range window -> power map (FD x NF)"""
    if len(ts) < 3: return None
    R = []
    for X, m in zip(Xs, masks):
        idx = np.flatnonzero(m); w = np.zeros(L)
        if len(idx) < 24: R.append(np.zeros((NA, NF), complex)); continue
        w[idx[0]:idx[-1] + 1] = np.hanning(idx[-1] - idx[0] + 1)
        R.append(np.fft.ifft(X * (w * m)[None, :], NF, axis=1))
    R = np.array(R)                                       # n x NA x NF
    tt = np.array(ts) - t0; wt = np.hanning(len(tt) + 2)[1:-1] if len(tt) > 2 else np.ones(len(tt))
    E = np.exp(-2j * np.pi * np.outer(FD, tt)) * wt[None, :]
    M = np.einsum("fn,nak->fak", E, R)
    return np.sum(np.abs(M) ** 2, 1) / max(len(tt), 1)

def peaks(P, k0, mask_f, n=10, rel=20):
    sub = P[:, k0 - HW:k0 + HW]; out = []
    if sub.max() <= 0: return out
    pm = np.zeros_like(sub, bool)
    pm[1:-1, 1:-1] = (sub[1:-1, 1:-1] >= np.maximum.reduce([sub[:-2, 1:-1], sub[2:, 1:-1], sub[1:-1, :-2], sub[1:-1, 2:]]))
    pm &= mask_f[:, None]
    ii = np.argwhere(pm); v = sub[pm]
    if not len(v): return out
    o = np.argsort(-v)[:n]; top = v[o[0]]
    return [tuple(ii[j]) for j in o if v[j] >= top * 10 ** (-rel / 10)]

def metrics(P, T, k0):
    if P is None: return None
    a = P[:, k0 - HW:k0 + HW]; b = T[:, k0 - HW:k0 + HW]
    def dbc(x, sel=None):
        x = x if sel is None else x[sel]
        return 10 * np.log10(np.maximum(x / max(x.max(), 1e-30), 1e-3))
    c = float(np.corrcoef(dbc(a).ravel(), dbc(b).ravel())[0, 1])
    cd = float(np.corrcoef(dbc(a, DYN).ravel(), dbc(b, DYN).ravel())[0, 1])
    pt = peaks(T, k0, DYN); pp = peaks(P, k0, DYN)
    near = lambda p, lst: any(abs(p[0] - q[0]) <= 2 and abs(p[1] - q[1]) <= 2 for q in lst)
    rec = float(np.mean([near(p, pp) for p in pt])) if pt else np.nan
    fal = float(np.mean([not near(p, pt) for p in pp])) if pp else np.nan
    prof = a.sum(0); k = int(np.argmax(prof)); h = prof >= prof[k] / 2
    lo = k; hi = k
    while lo > 0 and h[lo - 1]: lo -= 1
    while hi < len(h) - 1 and h[hi + 1]: hi += 1
    return dict(corr=c, dyn=cd, recall=rec, false=fal, width=(hi - lo + 1) * BIN)

CASES = {"split2": [(0, 34), (34, 68)], "split4": [(0, 17), (17, 34), (34, 51), (51, 68)], "gap": [(0, 30), (38, 68)]}
SIG = (0.0, 0.15, 0.5)
Rz = C.defaultdict(list)
BURSTS = sorted(glob.glob(REC + "/b*/"), key=lambda s: int(s.rstrip("/").split("/b")[-1]))
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
    UA = max(S, key=lambda u: len(S[u])); phs = sfcorr(sorted(S[UA], key=lambda x: x[0]))
    for u in S: S[u] = sorted([(t, a, pm, Y * np.exp(-1j * phs * (a // 2))) for t, a, pm, Y in S[u]], key=lambda x: x[0])
    A = S[UA]; XA, EPA = stitch(A)
    TA = np.array([x[0] for x in A]); slotA = C.defaultdict(list)
    for i, x in enumerate(A): slotA[x[1]].append(i)
    nw = 0
    for w0 in np.arange(TA[0], TA[-1] - WIN, WIN):
        full = [i for i in np.flatnonzero((TA >= w0) & (TA < w0 + WIN)) if A[i][2].sum() >= 60]
        if len(full) < 100: continue
        nw += 1
        ts = [A[i][0] for i in full]; Xf = [XA[i] for i in full]; Mf = [np.repeat(A[i][2], 12) for i in full]
        T = rd(ts, Xf, Mf, w0); k0 = int(np.argmax(T.sum(0))); k0 = min(max(k0, HW), NF - HW)
        # ---- part 1
        for cname, pieces in (CASES.items() if os.environ.get("PART", "12").find("1") >= 0 else []):
            pmask = [np.zeros(L, bool) for _ in pieces]
            for j, (s, e) in enumerate(pieces): pmask[j][12 * s:12 * e] = True
            for sg in SIG:
                cs = [1 + 0j] + [rng.uniform(0.5, 2) * np.exp(1j * rng.uniform(-np.pi, np.pi)) for _ in pieces[1:]]
                pert = [np.ones(L, complex)] + [np.repeat(1 + np.sqrt(sg / 2) * (rng.standard_normal(NG) + 1j * rng.standard_normal(NG)), 12) for _ in pieces[1:]]
                obs = []                                   # per snapshot: list of piece arrays as observed
                for X, M in zip(Xf, Mf):
                    obs.append([X * (pmask[j] & M)[None, :] * (cs[j] * pert[j])[None, :] for j in range(len(pieces))])
                maps = {}
                sup = [M & np.logical_or.reduce(pmask) for M in Mf]
                maps["oracle"] = rd(ts, [sum(o[j] / cs[j] for j in range(len(pieces))) for o in obs], sup, w0)
                maps["raw"] = rd(ts, [sum(o) for o in obs], sup, w0)
                if cname != "gap":
                    Xe = []
                    for o, M in zip(obs, Mf):
                        Z = o[0].copy(); mL = pmask[0] & M; ok = True
                        for j in range(1, len(pieces)):
                            mR = pmask[j] & M
                            if mL.sum() < 24 or mR.sum() < 24: ok = False; break
                            f = edge_number(Z, mL, o[j], mR, 24)
                            if not np.isfinite(f) or f == 0: ok = False; break
                            Z = Z + o[j] / f; mL = mL | mR
                        Xe.append(Z if ok else o[0])
                    maps["edge"] = rd(ts, Xe, sup, w0)
                inc = None
                for j in range(len(pieces)):
                    Pj = rd(ts, [o[j] for o in obs], [pmask[j] & M for M in Mf], w0)
                    if Pj is not None: inc = Pj if inc is None else inc + Pj
                maps["inc"] = inc
                maps["single"] = rd(ts, [o[0] for o in obs], [pmask[0] & M for M in Mf], w0)
                for mname, P in maps.items():
                    r = metrics(P, T, k0)
                    if r:
                        for kk, v in r.items():
                            if np.isfinite(v): Rz[("ctl", cname, sg, mname, kk)].append(v)
        rT = metrics(T, T, k0); Rz[("truth", "width")].append(rT["width"])
        # ---- part 2: real other UEs on their own support, incoherent per (UE, epoch)
        inc_r = None; inc_t = None; cover = np.zeros((len(range(int(WIN / 0.5e-3))), NG), bool); nseg = 0
        for u, lst in S.items():
            if u == UA: continue
            B = [x for x in lst if w0 <= x[0] < w0 + WIN]
            if len(B) < 3: continue
            EPb, _ = epochs(B); segs = C.defaultdict(list)
            for i, e in enumerate(EPb): segs[e].append(i)
            for e, idx in segs.items():
                tr, Xr, Mr, Xt, Mt = [], [], [], [], []
                for i in idx:
                    t, a, pb, Yb = B[i]
                    # same PRBs are never shared in one slot: nearest tablet snapshot within 1 ms covering B's PRBs
                    near = [k for k in np.flatnonzero(np.abs(TA - t) <= 0.001) if A[k][2][pb].all()]
                    if not near: continue
                    ia = [min(near, key=lambda k: abs(TA[k] - t))]
                    g = joint_gains(Yb, pb); Z = Yb / np.repeat(g, 12)[None, :]; m = np.repeat(pb, 12); Z[:, ~m] = 0
                    tr.append(t); Xr.append(Z); Mr.append(m); Xt.append(XA[ia[0]] * m[None, :]); Mt.append(m)
                if len(tr) < 3: continue
                nseg += 1
                P1 = rd(tr, Xr, Mr, w0); P2 = rd(tr, Xt, Mt, w0)
                # scale every segment to the power the tablet has there (heights are not comparable across beams anyway)
                s1 = P1.sum(); s2 = P2.sum()
                if s1 > 0 and s2 > 0:
                    P1 = P1 / s1 * s2
                    inc_r = P1 if inc_r is None else inc_r + P1; inc_t = P2 if inc_t is None else inc_t + P2
        if inc_r is not None:
            Rz[("real", "segments")].append(nseg)
            for mname, P in (("inc_real", inc_r), ("inc_tab", inc_t)):
                r = metrics(P, T, k0)
                for kk, v in r.items():
                    if np.isfinite(v): Rz[("real", mname, kk)].append(v)
    print("burst %d: %d windows" % (bi, nw), flush=True)
json.dump({"|".join(map(str, k)): v for k, v in Rz.items()}, open(OUTF, "w"))
def q(k):
    v = np.array(Rz.get(k, []), float); return "%5.2f" % np.median(v) if len(v) else "  -  "
print("truth range width m", q(("truth", "width")))
print("%-7s %-5s %-7s %6s %6s %6s %6s %6s" % ("case", "sig2", "method", "corr", "dyn", "recall", "false", "width"))
for cname in CASES:
    for sg in SIG:
        for m in ("oracle", "edge", "raw", "inc", "single"):
            if ("ctl", cname, sg, m, "corr") not in Rz: continue
            print("%-7s %-5s %-7s %s %s %s %s %s" % (cname, sg, m, *[q(("ctl", cname, sg, m, kk)) for kk in ("corr", "dyn", "recall", "false", "width")]))
print("real: windows", len(Rz.get(("real", "segments"), [])), "median segments", q(("real", "segments")))
for m in ("inc_real", "inc_tab"):
    print("%-9s %s %s %s %s %s" % (m, *[q(("real", m, kk)) for kk in ("corr", "dyn", "recall", "false", "width")]))
