"""Sparse multi-path model vs edge matching, controlled on full-band tablet snapshots (truth known).
Model: piece k (a run of PRGs) of antenna a:  X_a[f] = c_k * sum_p alpha_{p,a} exp(-j 2 pi f tau_p)
  common delays tau_p for all antennas and pieces, own amplitude per antenna, one unknown complex number c_k per piece.
Solver: alternate (1) simultaneous OMP over antennas on the data divided by c (zero-filled, 4x oversampled FFT grid),
LS refit of all amplitudes, (2) LS update of each c_k from the model; c_1 = 1 fixes the gauge.
Truth: epoch-stitched full-band snapshot (>= 64 of 68 PRGs). Cases:
  A  one hole of g PRGs in the middle, right side times a random number
  B  two 'UEs' of w PRGs each, s PRGs apart, the second times a random number (no shared PRGs: no anchor possible)
  C  raw snapshot (precoder jump every PRG), one unknown number per PRG (expected to fail: gauge freedom)
Metrics: number error (phase deg, amp dB) for edge vs OMP; OMP fill of the hole vs truth (nmse).
usage: gs_omp.py <rec dir> <out.json> [bursts]   env CH (0,1,2,3), EVERY (15)"""
import sys, os, json, numpy as np, collections as C
REC, OUTF = sys.argv[1], sys.argv[2]; ONLY = [int(x) for x in sys.argv[3:]]; saved = list(sys.argv)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", REC, "/dev/null"]; exec(src); sys.argv = saved
EVERY = int(os.environ.get("EVERY", "15")); OS = 4; G = OS * L; PMAX = 12
rng = np.random.default_rng(11)
fidx = np.arange(L)

def somp(Z, m, pmax=PMAX, tol=0.01):
    """Z: NA x L (zero where m False). returns delays (grid idx) and amplitudes NA x P, model NA x L"""
    R = Z * m[None, :]; sel = []; E0 = np.sum(np.abs(R) ** 2); prevE = E0
    A = None; amp = None
    for _ in range(pmax):
        cor = np.sum(np.abs(np.fft.fft(R, G, axis=1)) ** 2, 0)      # sum_a |sum_f R e^{-j2pi f k/G}|^2
        k = int(np.argmax(cor))
        if k in sel: break
        sel.append(k)
        A = np.exp(2j * np.pi * np.outer(fidx[m], sel) / G)           # Mm x P
        amp, *_ = np.linalg.lstsq(A, Z[:, m].T, rcond=None)            # P x NA
        R = np.zeros_like(Z); R[:, m] = Z[:, m] - (A @ amp).T
        E = np.sum(np.abs(R) ** 2)
        if (prevE - E) < tol * E0: break
        prevE = E
    Af = np.exp(2j * np.pi * np.outer(fidx, sel) / G)
    return sel, (Af @ amp).T

def omp_joint(Y, pieces, iters=8, c0=None):
    """Y: NA x L, pieces: list of bool masks (len L). returns c (len K, c[0]=1), model (unit-gain) NA x L"""
    K = len(pieces); c = np.ones(K, complex) if c0 is None else c0.copy(); m = np.logical_or.reduce(pieces)
    for _ in range(iters):
        Z = np.zeros_like(Y)
        for k, pk in enumerate(pieces): Z[:, pk] = Y[:, pk] / c[k]
        sel, M = somp(Z, m)
        for k in range(1, K):
            pk = pieces[k]; d = np.vdot(M[:, pk], M[:, pk])
            if abs(d) > 0: c[k] = np.vdot(M[:, pk], Y[:, pk]) / d
    return c, M

def omp_search(Y, pieces, P, nph=72):
    """global search of piece 2's phase: fixed P paths, keep the phase with the smallest residual, then LS amplitude"""
    best = None; m = np.logical_or.reduce(pieces)
    for ph in np.linspace(-np.pi, np.pi, nph, endpoint=False):
        Z = Y.copy(); Z[:, pieces[1]] *= np.exp(-1j * ph)
        sel, M = somp(Z, m, pmax=P, tol=0)
        r = np.sum(np.abs((Z - M)[:, m]) ** 2)
        if best is None or r < best[0]: best = (r, ph, M)
    r, ph, M = best; pk = pieces[1]
    c = np.vdot(M[:, pk], Y[:, pk] * np.exp(-1j * ph)) / np.vdot(M[:, pk], M[:, pk]) * np.exp(1j * ph)
    return c * np.exp(1j*0), M

def err(est, true):
    r = est / true; return abs(float(np.degrees(np.angle(r)))), abs(float(20 * np.log10(abs(r))))

def rnd(): return rng.uniform(0.5, 2) * np.exp(1j * rng.uniform(-np.pi, np.pi))

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
    EP, TPL = epochs(A); EG = {k: joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0) for k, (S2, C2) in TPL.items()}
    n = 0
    for i in range(0, len(A), EVERY):
        t, a, pm, Yraw = A[i]
        if pm.sum() < 64: continue
        X = Yraw / np.repeat(EG[EP[i]], 12)[None, :]
        r = runs(pm); s0, e0 = max(r, key=lambda x: x[1] - x[0])
        if e0 - s0 < 56: continue
        n += 1
        P = lambda lo, hi: np.r_[np.zeros(12 * lo, bool), np.ones(12 * (hi - lo), bool), np.zeros(L - 12 * hi, bool)]
        # ---- A: hole of g PRGs in the middle
        for g in (0, 1, 2, 4, 8, 16):
            mid = (s0 + e0 - g) // 2; pL = P(s0, mid); pR = P(mid + g, e0); cc = rnd()
            Y = X.copy(); Y[:, pR] *= cc; Y[:, ~(pL | pR)] = 0
            fe = edge_number(Y, pL, Y, pR, 24)
            if np.isfinite(fe): e = err(fe, cc); Rz[("A", g, "edge", "deg")].append(e[0]); Rz[("A", g, "edge", "db")].append(e[1])
            c, M = omp_joint(Y, [pL, pR])
            e = err(c[1], cc); Rz[("A", g, "omp", "deg")].append(e[0]); Rz[("A", g, "omp", "db")].append(e[1])
            c2, M2 = omp_joint(Y, [pL, pR], c0=np.array([1, fe if np.isfinite(fe) else 1]))
            e = err(c2[1], cc); Rz[("A", g, "omp+edge init", "deg")].append(e[0]); Rz[("A", g, "omp+edge init", "db")].append(e[1])
            for Pn in (3, 6, 12):
                cs_, _ = omp_search(Y, [pL, pR], Pn)
                e = err(cs_, cc); Rz[("A", g, "search P%d" % Pn, "deg")].append(e[0])
            if g:
                Yo = X.copy(); Yo[:, ~(pL | pR)] = 0; h = P(mid, mid + g)
                for Pn in (6, 12, 24, 48):
                    sel, Mo = somp(Yo, pL | pR, pmax=Pn, tol=0); Rz[("A", g, "fill true number P%d" % Pn)].append(nmse(Mo[:, h], X[:, h]))
                h = P(mid, mid + g); sc = np.vdot(M2[:, pL], X[:, pL]) / np.vdot(M2[:, pL], M2[:, pL])
                Rz[("A", g, "fill nmse")].append(nmse(M2[:, h] * sc, X[:, h]))
        # ---- B: two small UEs, no shared PRGs
        for w in (4, 8):
            for sgap in (0, 2, 8, 16, 32):
                if 2 * w + sgap > e0 - s0: continue
                lo = s0 + rng.integers(0, e0 - s0 - 2 * w - sgap + 1); p1 = P(lo, lo + w); p2 = P(lo + w + sgap, lo + 2 * w + sgap); cc = rnd()
                Y = np.zeros_like(X); Y[:, p1] = X[:, p1]; Y[:, p2] = X[:, p2] * cc
                fe = edge_number(Y, p1, Y, p2, min(24, 6 * w))
                if np.isfinite(fe): e = err(fe, cc); Rz[("B", w, sgap, "edge", "deg")].append(e[0]); Rz[("B", w, sgap, "edge", "db")].append(e[1])
                c, M = omp_joint(Y, [p1, p2], c0=np.array([1, fe if np.isfinite(fe) else 1]))
                e = err(c[1], cc); Rz[("B", w, sgap, "omp", "deg")].append(e[0]); Rz[("B", w, sgap, "omp", "db")].append(e[1])
                for Pn in (3, 6):
                    cs_, _ = omp_search(Y, [p1, p2], Pn); e = err(cs_, cc); Rz[("B", w, sgap, "search P%d" % Pn, "deg")].append(e[0])
        # ---- C: raw snapshot, one number per PRG (on 1 in 3 snapshots, slow)
        if n % 3 == 1:
            pcs = [P(b, b + 1) for b in range(s0, e0)]
            gtrue = EG[EP[i]][s0:e0] / EG[EP[i]][s0]          # the numbers stitching found = reference
            Y = Yraw.copy(); Y[:, ~P(s0, e0)] = 0
            c, M = omp_joint(Y, pcs, iters=6)
            ph = np.degrees(np.angle(c[1:] / gtrue[1:])); Rz[("C", "omp per PRG", "deg")].append(float(np.median(np.abs(ph))))
            far = np.degrees(np.angle(c[-1] / gtrue[-1])); Rz[("C", "omp per PRG far end", "deg")].append(abs(float(far)))
    print("burst %d: %d snapshots tested" % (bi, n), flush=True)
out = {"|".join(map(str, k)): v for k, v in Rz.items()}
json.dump(out, open(OUTF, "w"))
def q(k):
    v = np.array(Rz.get(k, []), float)
    return "%5.1f / %5.1f (n %d)" % (np.median(v), np.percentile(v, 90), len(v)) if len(v) else "-"
print("A hole: g PRGs | edge deg | omp deg | omp+edge-init deg | omp+edge dB | fill nmse")
for g in (0, 1, 2, 4, 8, 16):
    print(g, "|", q(("A", g, "edge", "deg")), "|", q(("A", g, "omp", "deg")), "|", q(("A", g, "omp+edge init", "deg")), "|", q(("A", g, "omp+edge init", "db")), "|", q(("A", g, "fill nmse")))
print("A search: g | P3 | P6 | P12 | fill(true number) P6 | P12 | P24 | P48")
for g in (0, 1, 2, 4, 8, 16):
    print(g, "|", q(("A", g, "search P3", "deg")), "|", q(("A", g, "search P6", "deg")), "|", q(("A", g, "search P12", "deg")), "|", " | ".join(q(("A", g, "fill true number P%d" % n)) for n in (6, 12, 24, 48)))
print("B search: w s | P3 | P6")
for w in (4, 8):
    for s_ in (0, 2, 8, 16, 32): print(w, s_, "|", q(("B", w, s_, "search P3", "deg")), "|", q(("B", w, s_, "search P6", "deg")))
print("B two UEs: w, sep | edge deg | omp deg | omp dB")
for w in (4, 8):
    for s in (0, 2, 8, 16, 32): print(w, s, "|", q(("B", w, s, "edge", "deg")), "|", q(("B", w, s, "omp", "deg")), "|", q(("B", w, s, "omp", "db")))
print("C per-PRG numbers on raw:", q(("C", "omp per PRG", "deg")), "far end", q(("C", "omp per PRG far end", "deg")))
