"""6f: TRS and CSI-RS as references for DM-RS sensing (5 Nov recording, replays with golden_sniffer.dump_blind).
 A  TRS stability: normalised correlation between TRS snapshots vs time lag (same symbol), 1 - |rho|^2 = what one
    complex number cannot explain. Lags: 4 symbols (same slot), 1 slot, 20 ms, 40 ms, 100 ms, 500 ms.
 B  TRS as anchor: each DM-RS snapshot (self-stitched across its PRGs) against the nearest TRS snapshot within BR,
    one complex number on its PRGs, residual eps (the gate quantity, gate 0.3). By UE kind: tablet, other placed UE,
    blind (6e), and by time distance. Reference: the same with the nearest tablet snapshot (today's anchor).
 C  Light network (tablet removed): share of DM-RS PRG-snapshots that join coherently through the TRS (anchor within
    BR and eps < GATE), and median band per 150 ms window: DM-RS joined via TRS, plus the TRS symbols themselves.
 D  CSI-RS (8 ports): stability vs lag per port in A; and (4 antennas only) whether a tablet DM-RS snapshot lies in the span of the 8 ports
    (per PRG, 8 unknown weights, LS residual), against a null with the CSI-RS taken 500 ms away.
    SLOPE=1 (default): the anchor fits also a relative delay slope (NR-Scope timing-window jumps between snapshots,
    up to ~16 samples between TRS occasions); SLOPE=0: one complex number only.
usage: gs_trs.py <rec dir> <out.json> [bursts]   env CH (0,3), DIR (replay_blind_ch), BR (0.01), GATE (0.3), SLOPE (1)"""
import sys, os, json, numpy as np, collections as C
REC, OUTF = sys.argv[1], sys.argv[2]; ONLY = [int(x) for x in sys.argv[3:]]; saved = list(sys.argv)
os.environ.setdefault("CH", "0,3")
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", REC, "/dev/null"]; exec(src); sys.argv = saved
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from gs_trs_lib import *
DIR = os.environ.get("DIR", "replay_blind_ch"); BR = float(os.environ.get("BR", "0.01")); GATE = float(os.environ.get("GATE", "0.3")); SLOPE = os.environ.get("SLOPE", "1") == "1"
TRS_C = 2; WIN = 300

def sfstep(A):
    st = []
    for (t1, a1, p1, Y1), (t2, a2, p2, Y2) in zip(A[:-1], A[1:]):
        m = np.repeat(p1 & p2, 12)
        if m.sum() < 120 or a2 - a1 > 2: continue
        if a2 // 2 - a1 // 2 == 1:
            v = np.vdot(Y1[:, m], Y2[:, m]); st.append(v / abs(v))
    return float(np.angle(np.mean(st))) if st else 0.0

def rho2(x, y):
    m = (x != 0) & (y != 0)
    if m.sum() < 24: return np.nan
    a, b = x[m], y[m]; return float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real))

def fit_eps(ref, x, pm):
    """one number mapping ref onto x on PRGs pm (NA x 816 each): eps = |x - c ref|^2 / |x|^2"""
    m = np.repeat(pm, 12)
    if m.sum() < 48: return np.nan
    r = ref[:, m]; v = x[:, m]; den = np.vdot(r, r).real
    if den <= 0: return np.nan
    c = np.vdot(r, v) / den; return float(np.sum(np.abs(v - c * r) ** 2) / np.sum(np.abs(v) ** 2))

NFS = 16384
def fit_eps_slope(ref, x, pm):
    """one number and a relative delay slope (timing-window jumps between the two snapshots): x ~ c e^{-j s q} ref.
    s from the peak of the zero-padded FFT of sum_a conj(ref) x, parabolic refinement; -> eps"""
    m = np.repeat(pm, 12)
    if m.sum() < 48: return np.nan
    z = np.sum(np.conj(ref) * x * m[None, :], 0)
    F = np.fft.fft(z, NFS); k = int(np.argmax(np.abs(F)))
    a_, b_, c_ = np.abs(F[k - 1]), np.abs(F[k]), np.abs(F[(k + 1) % NFS]); d = 0.5 * (a_ - c_) / (a_ - 2 * b_ + c_ + 1e-30)
    s_ = -2 * np.pi * (k + d) / NFS
    r = ref[:, m] * np.exp(-1j * s_ * q[m])[None, :]; v = x[:, m]; den = np.vdot(r, r).real
    if den <= 0: return np.nan
    cc = np.vdot(r, v) / den; return float(np.sum(np.abs(v - cc * r) ** 2) / np.sum(np.abs(v) ** 2))

LAGS = [("4 symbols", 0), ("1 slot", 1), ("20 ms", 40), ("40 ms", 80), ("100 ms", 200), ("500 ms", 1000)]   # in slots
Rz = C.defaultdict(list)
BURSTS = sorted(glob.glob(REC + "/b*/"), key=lambda s: int(s.rstrip("/").split("/b")[-1]))
for bi, bd in enumerate(BURSTS):
    if ONLY and bi not in ONLY: continue
    fs = [bd + DIR + "%d/ue_dmrs.bin" % c for c in CH]
    if not all(os.path.exists(f) and os.path.getsize(f) > 20e6 for f in fs): continue
    Dd = [load(f) for f in fs]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
    S = C.defaultdict(list)
    for k in sorted(keys):
        pm = np.logical_and.reduce([Dd[c][k][0] for c in range(NA)]); Y = np.array([Dd[c][k][1] for c in range(NA)]).astype(complex)
        if pm.sum() and np.all(np.isfinite(Y[:, np.repeat(pm, 12)])): S[k[1]].append((k[0] * 0.5e-3 + k[2] / 28e3, k[0], pm, Y))
    del Dd
    real = sorted((u for u in S if u != 0), key=lambda u: -len(S[u]))
    if not real: continue
    TAB = real[0]
    if 0 in S:   # blind minus placed PRGs in the slot
        cov0 = C.defaultdict(lambda: np.zeros(NG, bool))
        for u in real:
            for x in S[u]: cov0[x[1]] |= x[2]
        S[0] = [(t, a, pm & ~cov0[a], Y) for t, a, pm, Y in S[0] if (pm & ~cov0[a]).sum() >= 1]
    phs = sfstep(sorted(S[TAB], key=lambda x: x[0]))
    for u in S: S[u] = sorted([(t, a, pm, Y * np.exp(-1j * phs * (a // 2))) for t, a, pm, Y in S[u]], key=lambda x: x[0])
    # TRS and CSI-RS snapshots (NA x 816, on the DM-RS pair grid: 3 values per PRB, first 272 PRBs)
    Gs = [load_gsc(f) for f in fs]; gk = set(Gs[0]).intersection(*[set(g) for g in Gs[1:]])
    TR = []; CS = []
    for (a, l) in sorted(gk):
        t = a * 0.5e-3 + l / 28e3; rot = np.exp(-1j * phs * (a // 2))
        if l in (4, 8):
            H = np.array([trs_despread(*Gs[c][(a, l)], a % 20, l, TRS_C)[:N].ravel() for c in range(NA)]) * rot
            pmt = (np.abs(H).reshape(NA, N, 3).mean(2).min(0) > 0); pmt = pmt.reshape(NG, 4).all(1)
            TR.append((t, a, l, pmt, H))
        elif l == 13:
            P = np.array([csi8_despread(*Gs[c][(a, l)], a % 20)[:N] for c in range(NA)]) * rot   # NA x N x 8
            CS.append((t, a, P))
    del Gs
    TT = np.array([x[0] for x in TR]); TC = np.array([x[0] for x in CS])
    # ---- A: stability vs lag, exact slot pairs (same symbol; "4 symbols" = symbol 4 -> 8 of one slot)
    TRk = {(x[1], x[2]): x[4] for x in TR}; CSk = {x[1]: x[2] for x in CS}
    for name, ls in LAGS:
        for (a, l), H in TRk.items():
            k2 = (a, 8) if ls == 0 else (a + ls, l)
            if (ls == 0 and l != 4) or k2 not in TRk: continue
            r = rho2(H.ravel(), TRk[k2].ravel())
            if np.isfinite(r): Rz[("A", "trs", name)].append(1 - r)
        if ls == 0: continue
        for a, P in CSk.items():
            if a + ls not in CSk: continue
            for p_ in range(8):
                r = rho2(P[:, :, p_].ravel(), CSk[a + ls][:, :, p_].ravel())
                if np.isfinite(r): Rz[("A", "csi port %d" % p_, name)].append(1 - r)
    # ---- B/C: anchors
    tabX = []
    for t, a, pm, Y in S[TAB]:
        g = joint_gains(Y, pm); X = Y / np.repeat(g, 12)[None, :]; X[:, ~np.repeat(pm, 12)] = 0; tabX.append((t, pm, X))
    TA = np.array([x[0] for x in tabX])
    kinds = {TAB: "tablet"}
    for u in S: kinds.setdefault(u, "blind" if u == 0 else "other")
    win_trs = C.Counter(); win_dm = C.Counter(); win_dl = C.Counter()
    for t, a, l, pmt, H in TR: win_trs[a // WIN] += 4 * int(pmt.sum()) / 2      # PRB-slots (2 TRS symbols per slot)
    surv = set()
    for line in list(open(bd + DIR + "0/survey.csv"))[1:]:
        r = line.split(",")
        if len(r) > 2: surv.add(int(r[0]) * 20 + int(r[1]))
    a_min = min(x[1] for x in S[TAB])
    for u, lst in S.items():
        kd = kinds[u]
        for t, a, pm, Y in lst:
            if pm.sum() < 1: continue
            g = joint_gains(Y, pm); X = Y / np.repeat(g, 12)[None, :]; X[:, ~np.repeat(pm, 12)] = 0
            n = int(pm.sum())
            j = int(np.argmin(np.abs(TT - t))) if len(TT) else -1
            dt = abs(TT[j] - t) if j >= 0 else 9
            e = (fit_eps_slope if SLOPE else fit_eps)(TR[j][4], X, pm & TR[j][3]) if dt <= BR else np.nan
            if np.isfinite(e):
                b = "0-2 ms" if dt <= 0.002 else "2-5 ms" if dt <= 0.005 else "5-10 ms"
                Rz[("B", "trs", kd, b)].append(e); Rz[("B", "trs", kd, "all")].append(e)
            if kd != "tablet":
                jt = int(np.argmin(np.abs(TA - t)))
                if abs(TA[jt] - t) <= BR:
                    e2 = (fit_eps_slope if SLOPE else fit_eps)(tabX[jt][2], X, pm & tabX[jt][1])
                    if np.isfinite(e2): Rz[("B", "tablet", kd, "all")].append(e2)
            if kd != "tablet":   # C: light network, tablet removed
                ok = np.isfinite(e) and e < GATE
                Rz[("C", kd, "PRG-snapshots")].append(n); Rz[("C", kd, "joined via TRS")].append(n if ok else 0)
                if ok: win_dm[a // WIN] += 4 * n / 2
    for a in surv:
        if a >= a_min: win_dl[a // WIN] += 1
    for w, ndl in win_dl.items():
        if ndl >= 50:
            Rz[("C", "band", "DM-RS via TRS")].append(win_dm[w] / (273 * ndl))
            Rz[("C", "band", "TRS itself")].append(win_trs[w] / (273 * ndl))
            Rz[("C", "band", "both")].append((win_dm[w] + win_trs[w]) / (273 * ndl))
    # ---- D: tablet DM-RS in the span of the 8 CSI-RS ports, per PRG
    for t, pm, X in (tabX[::3] if NA >= 4 else []):   # needs 4 antennas: with 2, 8 weights fit 8 independent values
        j = int(np.argmin(np.abs(TC - t))) if len(TC) else -1
        if j < 0 or abs(TC[j] - t) > 0.001: continue
        jn = int(np.argmin(np.abs(TC - t - 0.5)))
        for tag, P in (("near", CS[j][2]), ("500 ms away", CS[jn][2] if abs(TC[jn] - t - 0.5) < 0.05 else None)):
            if P is None: continue
            res = []
            for b in np.flatnonzero(pm):
                y = X[:, 12 * b:12 * b + 12].reshape(NA, 4, 3)           # antenna, PRB, pair
                A_ = np.repeat(P[:, 4 * b:4 * b + 4, :][:, :, None, :], 3, 2).reshape(-1, 8); yv = y.reshape(-1)
                w, *_ = np.linalg.lstsq(A_, yv, rcond=None); res.append(np.sum(np.abs(yv - A_ @ w) ** 2) / np.sum(np.abs(yv) ** 2))
            if res: Rz[("D", tag)].append(float(np.median(res)))
    print("burst %d: %d TRS, %d CSI-RS snapshots, subframe step %.1f deg, UEs %d" % (bi, len(TR), len(CS), np.degrees(phs), len(S)), flush=True)

json.dump({"|".join(map(str, k)): v for k, v in Rz.items()}, open(OUTF, "w"))
def q(k, f="%.2f"):
    v = np.array(Rz.get(k, []), float); v = v[np.isfinite(v)]
    return (f % np.median(v) + " (n %d)" % len(v)) if len(v) else "-"
print("A  1-|rho|^2 vs lag      TRS                 CSI-RS port 0          port 3               port 7")
for name, _ in LAGS: print("   %-10s %20s %20s %20s %20s" % (name, q(("A", "trs", name), "%.3f"), q(("A", "csi port 0", name), "%.3f"), q(("A", "csi port 3", name), "%.3f"), q(("A", "csi port 7", name), "%.3f")))
print("B  eps against the TRS (median, gate %.1f), and pass rate" % GATE)
for kd in ("tablet", "other", "blind"):
    v = np.array(Rz.get(("B", "trs", kd, "all"), [])); vt = np.array(Rz.get(("B", "tablet", kd, "all"), []))
    print("   %-7s TRS %s pass %s | " % (kd, q(("B", "trs", kd, "all")), "%.0f%%" % (100 * np.mean(v < GATE)) if len(v) else "-") +
          " ".join("%s %s" % (b, q(("B", "trs", kd, b))) for b in ("0-2 ms", "2-5 ms", "5-10 ms")) +
          (" | vs tablet anchor %s pass %.0f%%" % (q(("B", "tablet", kd, "all")), 100 * np.mean(vt < GATE)) if len(vt) else ""))
for kd in ("other", "blind"):
    a_ = sum(Rz.get(("C", kd, "PRG-snapshots"), [])); b_ = sum(Rz.get(("C", kd, "joined via TRS"), []))
    if a_: print("C  %-6s PRG-snapshots joined via TRS: %.0f%% of %d" % (kd, 100 * b_ / a_, a_))
for k in ("DM-RS via TRS", "TRS itself", "both"):
    v = np.array(Rz.get(("C", "band", k), []))
    if len(v): print("C  band per 150 ms window, %-14s median %.1f%%  p90 %.1f%%  zero %.0f%%" % (k, 100 * np.median(v), 100 * np.percentile(v, 90), 100 * np.mean(v < 1e-3)))
print("D  tablet DM-RS per PRG in the span of the 8 CSI-RS ports: residual near %s, CSI-RS 500 ms away %s; random-subspace level %.2f (4 antennas only)" %
      (q(("D", "near")), q(("D", "500 ms away")), 1 - 8 / (NA * 4)))
