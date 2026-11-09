"""Simulation with known targets on the REAL occupancy of the 5 Nov recording: does each way of combining signals help
a 150 ms range-Doppler map?
Scene (bistatic): gNB with an 8-port ULA (the CSI-RS has 8 ports), receiver with 4 antennas (ULA). Paths: direct path,
20 static scatterers (CLO..CHI dB = -45..-25, 10..CR = 300 m excess path; calibrated so that PRG stitching errs ~2 deg
per edge, measured 1.6 deg on the cell), 2 moving targets (-25 dB, 20-300 m, Doppler 15-250 Hz,
random sign). y[rx, k] = sum_p alpha_p a_rx(theta_p) (a_tx(phi_p)^H w) e^{-j 2 pi f_k tau_p} e^{j 2 pi nu_p t} + noise.
Signals, exactly where the recording had them (usage-map JSON of each burst: slots, PRBs, UE, TRS slots, CSI-RS slots):
  DM-RS  per UE: w = beam to the UE's direction plus a second beam mixed with a random complex weight, new weight every
         precoder epoch (~19 ms) and, for per-PRG UEs (tablet, per-PRG UEs, blind UEs), per 4-PRB group as well;
         blind runs are given to 6 hidden UEs by frequency position; symbols 2 and 11
  TRS    fixed wide beam (port 0), full band, symbols 4 and 8 of the TRS slots
  CSI-RS 8 ports unprecoded (one per port), one value per PRB, symbol 13 of the CSI-RS slots
Variants (each gives one map; power summed over antennas, and ports for the CSI-RS):
  ideal          full band every DL slot, fixed beam (upper bound)
  TRS            TRS only
  TRS+DMRS@TRS   TRS + DM-RS pieces self-stitched and put on the nearest TRS (<= 10 ms) with one number, gate 0.3
  TRS+incoh      TRS map + incoherent sum of DM-RS groups (consecutive coherent pieces, epoch detector)
  CSI-RS         CSI-RS only (8 port maps, powers added)
  DMRS today     (with tablet only) tablet epochs stitched + glued, other UEs merged with the tablet anchor, gate 0.3
  DMRS today+TRS (with tablet only) the same plus the TRS on the TRS' own reference (one number fit to the tablet)
Metrics per window: detections (moving peaks >= 12 dB over the floor and within 15 dB of the strongest), targets
found (detection within 2 range bins and 4 Hz), false detections per window (Doppler aliases count as false), and the
target peak over the median moving floor (dB).
usage: sim_sense.py <usemap json dir> <out.json> <bursts...>   env SC (notab | tab), SNR (25 dB, per RE on the direct path), MAXW, SEED,
  MIX (0.2): weight of the second beam (0 = one beam per UE: precoder changes are pure complex numbers per PRG and epoch)"""
import sys, os, json, numpy as np, collections as C
UM, OUTF = sys.argv[1], sys.argv[2]; BUR = [int(x) for x in sys.argv[3:]]; saved = list(sys.argv)
os.environ["CH"] = "0,1,2,3"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
SC = os.environ.get("SC", "notab"); MIX = float(os.environ.get("MIX", "0.2")); CLO, CHI, CR = float(os.environ.get("CLO", "-45")), float(os.environ.get("CHI", "-25")), float(os.environ.get("CR", "300")); SNR = float(os.environ.get("SNR", "25"))
rng = np.random.default_rng(int(os.environ.get("SEED", "3")))
NT = 8; FD = np.arange(-300, 301, 2.0); MOV = np.abs(FD) >= 10; DF_SC = 30e3
FK = (12 * np.repeat(np.arange(N), 3) + np.tile([1, 5, 9], N)) * DF_SC          # DM-RS pair grid (816)
FT = (12 * np.repeat(np.arange(N), 3) + np.tile([2, 6, 10], N)) * DF_SC         # TRS grid (816)
FC = (12 * np.arange(N) + 3.5) * DF_SC                                          # CSI-RS, one per PRB
c0 = 3e8
def atx(phi): return np.exp(1j * np.pi * np.arange(NT) * np.sin(phi)) / np.sqrt(NT)
def arx(th): return np.exp(1j * np.pi * np.arange(NA) * np.sin(th))

def scene():
    P = []
    P.append((1.0, 0.0, 0.0, rng.uniform(-0.3, 0.3), rng.uniform(-0.5, 0.5)))           # direct: amp, tau, nu, aod, aoa
    for _ in range(20):   # calibrated: stitching error per PRG edge ~2 deg (measured on the cell: 1.6 deg)
        P.append((10 ** (rng.uniform(CLO, CHI) / 20) * np.exp(1j * rng.uniform(0, 2 * np.pi)), rng.uniform(10, CR) / c0, 0.0,
                  rng.uniform(-1.2, 1.2), rng.uniform(-1.2, 1.2)))
    tg = []
    for _ in range(2):
        r = rng.uniform(20, 300); nu = rng.uniform(15, 250) * rng.choice([-1, 1])
        P.append((10 ** (-25 / 20) * np.exp(1j * rng.uniform(0, 2 * np.pi)), r / c0, nu, rng.uniform(-1.2, 1.2), rng.uniform(-1.2, 1.2)))
        tg.append((r, nu))
    return P, tg

class Synth:
    def __init__(s, P):
        s.amp = np.array([p[0] for p in P]); s.tau = np.array([p[1] for p in P]); s.nu = np.array([p[2] for p in P])
        s.At = np.array([atx(p[3]) for p in P]); s.Ar = np.array([arx(p[4]) for p in P]).T      # P x NT, NA x P
        s.E = {nm: np.exp(-2j * np.pi * np.outer(s.tau, f)) for nm, f in (("dm", FK), ("trs", FT), ("csi", FC))}
        s.sig = 10 ** (-SNR / 20) / np.sqrt(2)
    def obs(s, grid, t, Wk):
        """Wk: NT x K precoder per sample (K = len of grid) -> NA x K"""
        g = s.amp * np.exp(2j * np.pi * s.nu * t)                                  # P
        B = (np.conj(s.At) @ Wk)                                                    # P x K : a_tx^H w
        M = s.E[grid] * B * g[:, None]
        y = s.Ar @ M
        return y + s.sig * (rng.standard_normal(y.shape) + 1j * rng.standard_normal(y.shape))

def rd(ts, Xs, masks, t0, mti=True):
    """static clutter removal (mean over the snapshots where each sample is present), range IFFT per snapshot (Hann on
    its own span), NUDFT over time (Hann), power summed over antennas"""
    if len(ts) < 2: return None
    if mti:
        Ms = np.array(masks, float); Xa = np.array(Xs)
        mean = (Xa * Ms[:, None, :]).sum(0) / np.maximum(Ms.sum(0), 1)[None, :]
        Xs = [(X - mean) * m[None, :] for X, m in zip(Xs, masks)]
    R = []
    for X, m in zip(Xs, masks):
        idx = np.flatnonzero(m); w = np.zeros(L)
        if len(idx) < 24: R.append(np.zeros((X.shape[0], NF), complex)); continue
        w[idx[0]:idx[-1] + 1] = np.hanning(idx[-1] - idx[0] + 1)
        R.append(np.fft.ifft(X * (w * m)[None, :], NF, axis=1))
    R = np.array(R); tt = np.array(ts) - t0; wt = np.hanning(len(tt) + 2)[1:-1]
    E = np.exp(-2j * np.pi * np.outer(FD, tt)) * wt[None, :]
    M = np.einsum("fn,nak->fak", E, R); return np.sum(np.abs(M) ** 2, 1)

def score(Pm, tg):
    """detections = local maxima on |nu| >= 10 Hz, 0-390 m, at least 12 dB over the median moving floor and within
    15 dB of the strongest moving peak. found = targets with a detection within 2 range bins and 4 Hz;
    false = detections not near a target (Doppler aliases of a target count as false)"""
    if Pm is None: return None
    sub = Pm[:, :160]
    pk = np.zeros_like(sub, bool)
    pk[1:-1, 1:-1] = sub[1:-1, 1:-1] >= np.maximum.reduce([sub[:-2, 1:-1], sub[2:, 1:-1], sub[1:-1, :-2], sub[1:-1, 2:], sub[:-2, :-2], sub[2:, 2:]])
    pk &= MOV[:, None]; fl = np.median(sub[MOV]) + 1e-30
    v = sub[pk]; ii = np.argwhere(pk)
    if not len(v): return dict(hit=0.0, false=0, snr=0.0)
    keep = (v >= fl * 10 ** 1.2) & (v >= v.max() / 10 ** 1.5); det = ii[keep]
    hit = 0; used = np.zeros(len(det), bool); snr = []
    for r, nu in tg:
        kb = r / (c0 / (NF * 120e3)); fb = (nu + 300) / 2
        ok = [j for j, (f, k) in enumerate(det) if abs(k - kb) <= 2 and abs(f - fb) <= 2]
        if ok: hit += 1; used[ok] = True
        f0 = int(round(fb)); k0 = int(round(kb)); snr.append(10 * np.log10(sub[max(f0 - 2, 0):f0 + 3, max(k0 - 2, 0):k0 + 3].max() / fl))
    return dict(hit=hit / len(tg), false=int((~used).sum()), snr=float(np.mean(snr)))

def fit1(ref, x, pm):
    m = np.repeat(pm, 12); r = ref[:, m]; v = x[:, m]; den = np.vdot(r, r).real
    if den <= 0 or m.sum() < 48: return np.nan, np.nan
    c = np.vdot(r, v) / den; return c, float(np.sum(np.abs(v - c * r) ** 2) / np.sum(np.abs(v) ** 2))

Rz = C.defaultdict(list)
for bi in BUR:
    J = json.load(open(os.path.join(UM, "um_b%d.json" % bi))); seg = np.array(J["seg"][SC]).reshape(-1, 5)
    rn = J["rnti"]; slots = {s[0]: s for s in J["slots"]}
    kind = {i: r[1] for i, r in enumerate(rn)}
    nw = 0
    for w in J["win"][SC]:
        w0 = w[0]; wsl = [a for a in slots if w0 <= a < w0 + 300]
        if len(wsl) < 150: continue
        nw += 1
        if nw > int(os.environ.get("MAXW", "999")): break
        P, tg = scene(); syn = Synth(P)
        # precoders: per UE two beams; epochs ~ exponential(19 ms); per-PRG UEs a weight per 4-PRB group
        ues = {}
        def W_of(u, a, typ):
            if u not in ues:
                phi = rng.uniform(-1.0, 1.0); ues[u] = dict(b1=atx(phi), b2=atx(phi + rng.uniform(0.15, 0.4)), next=-1, g=None, typ=typ)
            U = ues[u]
            if a >= U["next"]:
                U["next"] = a + int(rng.exponential(38)) + 1
                U["g"] = (rng.standard_normal(NG) + 1j * rng.standard_normal(NG)) * (MIX if U["typ"] != "wideband" else 0)
                U["g"] += MIX * (rng.standard_normal() + 1j * rng.standard_normal())
                U["sc"] = np.exp(1j * rng.uniform(0, 2 * np.pi)) * rng.uniform(0.5, 2)   # epoch: one complex number
                U["ph"] = np.exp(1j * rng.choice([np.pi / 2, -np.pi / 2], NG)) if U["typ"] != "wideband" else np.ones(NG)
            Wb = (U["b1"][:, None] + U["b2"][:, None] * U["g"][None, :]) * U["ph"][None, :]
            Wb = Wb / np.linalg.norm(Wb, axis=0, keepdims=True) * U["sc"]
            return np.repeat(Wb, 12, 1)                                         # NT x 816
        # DM-RS pieces
        pieces = []; ORA = []; ORU = []   # (t, a, ue, pm, Y); oracle copies
        m = (seg[:, 0] >= w0) & (seg[:, 0] < w0 + 300) & np.isin(seg[:, 3], [0, 1, 3, 4, 5, 6, 9])
        byslot = C.defaultdict(list)
        for a, lo, hi, cat, u in seg[m]: byslot[int(a)].append((int(lo), int(hi), int(u)))
        for a, lst in byslot.items():
            byue = C.defaultdict(lambda: np.zeros(NG, bool))
            for lo, hi, u in lst:
                if kind.get(u) == "blind":
                    u = 1000 + min(lo // 46, 5)                                 # 6 hidden UEs by frequency position
                byue[u][lo // 4:hi // 4 + 1] = True
            for u, pm in byue.items():
                typ = "per-PRG" if u >= 1000 else ("per-PRG" if kind.get(u) in ("tablet", "per-PRG") else "wideband")
                Wk = W_of(u, a, typ); mk = np.repeat(pm, 12)
                for l in (2, 11):
                    t = a * 0.5e-3 + l / 28e3; Y = syn.obs("dm", t, Wk); Y[:, ~mk] = 0
                    Yo = syn.obs("dm", t, np.repeat(np.eye(NT)[:, :1], L, 1)); Yo[:, ~mk] = 0   # oracle: same piece, one fixed beam for all
                    Yu = syn.obs("dm", t, np.repeat(ues[u]["b1"][:, None], L, 1)); Yu[:, ~mk] = 0  # the UE's own beam, no precoder changes
                    pieces.append((t, a, u, pm, Y)); ORA.append(Yo); ORU.append(Yu)
        o_ = np.argsort([x[0] for x in pieces], kind='stable'); pieces = [pieces[i] for i in o_]; ORA = [ORA[i] for i in o_]; ORU = [ORU[i] for i in o_]
        wtrs = np.repeat(np.eye(NT)[:, :1], L, 1)
        TRS = []
        for a in wsl:
            if slots[a][1]:
                for l in (4, 8):
                    t = a * 0.5e-3 + l / 28e3; TRS.append((t, syn.obs("trs", t, wtrs)))
        CSI = []
        for a in wsl:
            if slots[a][2]:
                t = a * 0.5e-3 + 13 / 28e3; CSI.append((t, [syn.obs("csi", t, np.repeat(np.eye(NT)[:, p:p + 1], N, 1)) for p in range(NT)]))
        t0 = w0 * 0.5e-3; full = np.ones(L, bool); maps = {}
        # ideal: every DL slot, full band, fixed beam
        ts = [a * 0.5e-3 + l / 28e3 for a in sorted(wsl) for l in (2, 11)]
        maps["ideal"] = rd(ts, [syn.obs("dm", t, wtrs) for t in ts], [full] * len(ts), t0)
        maps["TRS"] = rd([x[0] for x in TRS], [x[1] for x in TRS], [full] * len(TRS), t0)
        if pieces:
            maps["DMRS oracle"] = rd([x[0] for x in pieces], ORA, [np.repeat(x[3], 12) for x in pieces], t0)
        # DM-RS pieces self-stitched
        XP = []
        for t, a, u, pm, Y in pieces:
            g = joint_gains(Y, pm); X = Y / np.repeat(g, 12)[None, :]; X[:, ~np.repeat(pm, 12)] = 0; XP.append(X)
        TT = np.array([x[0] for x in TRS]); add_t, add_X, add_m = [], [], []
        for (t, a, u, pm, Y), X in zip(pieces, XP):
            if not len(TT): break
            j = int(np.argmin(np.abs(TT - t)))
            if abs(TT[j] - t) > 0.01: continue
            c, e = fit1(TRS[j][1], X, pm)
            if np.isfinite(e) and e < 0.3: add_t.append(t); add_X.append(X / c); add_m.append(np.repeat(pm, 12))
        Rz[(SC, "joined via TRS share")].append(len(add_t) / max(len(pieces), 1))
        maps["TRS+DMRS@TRS"] = rd([x[0] for x in TRS] + add_t, [x[1] for x in TRS] + add_X, [full] * len(TRS) + add_m, t0)
        # incoherent: groups of consecutive coherent pieces (epoch detector over all pieces in time order)
        inc = maps["TRS"].copy() if maps["TRS"] is not None else None
        if pieces:
            EP, _ = epochs([(t, a, pm, X) for (t, a, u, pm, Y), X in zip(pieces, XP)]); grp = C.defaultdict(list)
            for i, e in enumerate(EP): grp[e].append(i)
            for e, idx in grp.items():
                if len(idx) < 2: continue
                Pg = rd([pieces[i][0] for i in idx], [XP[i] for i in idx], [np.repeat(pieces[i][3], 12) for i in idx], t0)
                if Pg is not None and inc is not None: inc = inc + Pg
        maps["TRS+incoh"] = inc
        if CSI:
            cm = None
            for p in range(NT):
                Pp = rd([x[0] for x in CSI], [np.repeat(x[1][p], 3, 1)[:, :L] for x in CSI], [full] * len(CSI), t0)
                cm = Pp if cm is None else cm + Pp
            maps["CSI-RS"] = cm
        if SC == "tab":
            # today's pipeline: reference = UE with most PRG-snapshots, epochs stitched + glued; others anchor-merged
            cnt = C.Counter()
            for t, a, u, pm, Y in pieces: cnt[u] += int(pm.sum())
            if cnt:
                R = cnt.most_common(1)[0][0]; A = [(t, a, pm, Y) for t, a, u, pm, Y in pieces if u == R]
                EP, TPL = epochs(A); EG = {k: joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0) for k, (S2, C2) in TPL.items()}
                X0 = []
                for i, (t, a, pm, Y) in enumerate(A):
                    Z = Y / np.repeat(EG[EP[i]], 12)[None, :]; Z[:, ~np.repeat(pm, 12)] = 0; X0.append(Z)
                XA = [X0[0]]; gg = 1 + 0j
                for i in range(1, len(A)):
                    if EP[i] != EP[i - 1]:
                        f, n = frame(XA[i - 1], A[i - 1][2], X0[i], A[i][2]); gg = f if (np.isfinite(f) and n >= 4 and f != 0) else 1 + 0j
                    XA.append(X0[i] / gg)
                TA = np.array([x[0] for x in A]); ct, cX, cm_ = [x[0] for x in A], list(XA), [np.repeat(x[2], 12) for x in A]
                for (t, a, u, pm, Y), X in zip(pieces, XP):
                    if u == R: continue
                    j = int(np.argmin(np.abs(TA - t)))
                    if abs(TA[j] - t) > 0.01: continue
                    c, e = fit1(XA[j], X, pm & A[j][2])
                    if np.isfinite(e) and e < 0.3: ct.append(t); cX.append(X / c); cm_.append(np.repeat(pm, 12))
                maps["DMRS today"] = rd(ct, cX, cm_, t0)
                Ro = [ORU[i] for i, x in enumerate(pieces) if x[2] == R]
                num = sum(np.vdot(o, x) for o, x in zip(Ro, XA)); den = sum(np.vdot(o, o).real for o in Ro)
                cg = num / den; Rz[(SC, "diag", "reference vs oracle, one number over the window")].append(
                    float(sum(np.sum(np.abs(x - cg * o) ** 2) for o, x in zip(Ro, XA)) / sum(np.sum(np.abs(x) ** 2) for x in XA)))
                ers = []
                for o, x, (t, a, pm, Y) in zip(Ro, XA, A):
                    m_ = np.repeat(pm, 12); c1 = np.vdot(o[:, m_], x[:, m_]) / np.vdot(o[:, m_], o[:, m_]); ers.append(float(np.sum(np.abs(x[:, m_] - c1 * o[:, m_]) ** 2) / np.sum(np.abs(x[:, m_]) ** 2)))
                Rz[(SC, "diag", "reference vs oracle, one number per snapshot (stitching only)")].append(float(np.median(ers)))
                # + TRS on the reference's frame: one number from the nearest reference snapshot
                tt2, tX2, tm2 = list(ct), list(cX), list(cm_)
                for t, H in TRS:
                    j = int(np.argmin(np.abs(TA - t)))
                    if abs(TA[j] - t) > 0.01: continue
                    c, e = fit1(H, XA[j], A[j][2])
                    if np.isfinite(e) and e < 0.3: tt2.append(t); tX2.append(H * c); tm2.append(full)
                maps["DMRS today+TRS"] = rd(tt2, tX2, tm2, t0)
        for nm, Pm in maps.items():
            r = score(Pm, tg)
            if r:
                for k, v in r.items(): Rz[(SC, nm, k)].append(v)
    print("burst %d done (%s)" % (bi, SC), flush=True)
json.dump({"|".join(k): v for k, v in Rz.items()}, open(OUTF, "w"))
print("%-16s %8s %14s %10s %6s" % ("variant (" + SC + ")", "found", "false/window", "SNR dB", "n"))
for k, v in Rz.items():
    if k[1] == "diag": print("diag:", k[2], "median %.3f" % np.median(v))
for nm in ("ideal", "DMRS oracle", "TRS", "TRS+DMRS@TRS", "TRS+incoh", "CSI-RS", "DMRS today", "DMRS today+TRS"):
    v = Rz.get((SC, nm, "hit"))
    if v: print("%-16s %7.0f%% %14.1f %10.1f %6d" % (nm, 100 * np.mean(v), np.mean(Rz[(SC, nm, "false")]), np.median(Rz[(SC, nm, "snr")]), len(v)))
print("DM-RS pieces joined via TRS: %.0f%%" % (100 * np.mean(Rz[(SC, "joined via TRS share")])))
