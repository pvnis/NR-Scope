"""Walker echo in the tablet DM-RS with the co-phasing model. Per chunk (BA wide, tablet only, in-process):
 1. phi_b per 4-PRB group: each transmission's G_b = <R_b, Y_b>/|R_b|^2 fitted as B1(b) + phi_b B2(b) (smooth B, phi in
    {1, j, -1, -j}); transmissions aligned to a reference by the best global rotation; majority vote per group.
 2. DM-RS dynamic part d = (Y e^{-j theta} - S_local) / scale (one real scale per beam: NO per-group division).
 3. Per 0.2 s window and antenna, the walker echo fitted at the lidar range/Doppler as (A1 + phi_b A2) m (2 unknowns),
    and, for comparison, as A m (1 unknown, the old model); the same for ranges +12..+40 m (nulls).
 Output per model: echo-to-null ratio (projection energy at the truth over the mean at the null ranges), TRS alongside.
usage: cophase_mf.py <chunk list> <offset s>"""
import sys, os, io, json, contextlib, collections, numpy as np
KS = [int(x) for x in sys.argv[1].split(",")]; OFF = float(sys.argv[2]); HP = 0.25; WIN = 0.2; DEG = 2; GATE = 10 ** (-12 / 10)
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
GD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth/grid"
D = np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz"); TR = D["tr"]; t_rec = float(D["t_rec"]); IDS = sorted(set(TR[:, 0].astype(int)))
LAM = 3e8 / 3.75e9; DF = 30e3; QP = np.array([1, 1j, -1, -1j])
HYP = [("truth", 0.0)] + [("dv%+.1f" % o, float(o)) for o in (-1.2, -0.6, 0.6, 1.2)] + [("rev", 99.0)]   # Doppler nulls: extra bistatic speed (m/s), or the phase history reversed
ACC = collections.defaultdict(float); NW = collections.Counter(); STAB = []
rng = np.random.default_rng(0)
def basis(bs): x = (np.asarray(bs) - 34) / 34.0; return np.stack([x ** d for d in range(DEG + 1)], 1)
def cofit(bs, G):
    A = basis(bs); best = None
    for st in range(5):
        phi = np.ones(len(bs), complex) if st == 0 else QP[rng.integers(0, 4, len(bs))]
        for _ in range(10):
            M = np.concatenate([A, phi[:, None] * A], 1); c, *_ = np.linalg.lstsq(M, G, rcond=None)
            B1 = A @ c[:DEG + 1]; B2 = A @ c[DEG + 1:]
            phi = QP[np.argmin(np.abs(G[:, None] - (B1[:, None] + QP[None, :] * B2[:, None])), 1)]
        M = np.concatenate([A, phi[:, None] * A], 1); c, *_ = np.linalg.lstsq(M, G, rcond=None); e = np.sum(np.abs(G - M @ c) ** 2)
        if best is None or e < best[0]: best = (e, phi)
    return best[1]
for k in KS:
    os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
    sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/cm.bin"]
    g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
        "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
    E, own, TH, S, Sc, Rbar = g["KEEP"]; osc = g["osc"]; TRS = g["TRS"]; KKD = g["KK_D"]; KKT = g["KK_T"]
    a0 = json.load(open("%s/b%d.json" % (GD, k)))["a0"]
    tiq = lambda a, l: k * 2.5 + ((a - a0) % 20480) * 5e-4 + l / 28e3
    # 1. phi per group: per transmission fit, aligned, voted
    PHS = []
    for i, (t, a, l, pm, Y, u) in enumerate(E[::3]):
        bs = np.flatnonzero(pm)
        if len(bs) < 24: continue
        G = np.array([np.vdot(Rbar[:, 12 * b:12 * b + 12], Y[:, 12 * b:12 * b + 12]) / np.vdot(Rbar[:, 12 * b:12 * b + 12], Rbar[:, 12 * b:12 * b + 12]).real for b in bs])
        PHS.append(dict(zip(bs, cofit(bs, G))))
    if not PHS: print("chunk %d: no transmission wide enough for the co-phasing fit, skipped" % k); continue
    ref = max(PHS, key=len); votes = collections.defaultdict(collections.Counter)
    for ph in PHS:
        com = [b for b in ph if b in ref]
        if not com: continue
        rot = QP[np.argmax([sum(np.isclose(ph[b] * q, ref[b]) for b in com) for q in QP])]
        for b, p in ph.items(): votes[b][int(np.round(np.angle(p * rot) / (np.pi / 2))) % 4] += 1
    PHI = {b: QP[c.most_common(1)[0][0]] for b, c in votes.items()}
    STAB += [c.most_common(1)[0][1] / sum(c.values()) for c in votes.values() if sum(c.values()) >= 10]
    # 2. streams
    trs = sorted([(tiq(a, l), R / osc(t) - Rbar) for t, a, l, R in TRS], key=lambda x: x[0])
    dm = collections.defaultdict(list)
    for (i, b), kk in own.items():
        if Sc[kk] < 10 or b not in PHI: continue
        t, a, l, pm, Y, u = E[i]; sl = slice(12 * b, 12 * b + 12); yv = Y[:, sl] * np.exp(-1j * TH[(i, b)])
        if np.sum(np.abs(yv - S[kk]) ** 2) > GATE * np.sum(np.abs(S[kk]) ** 2): continue
        dm[(kk, b)].append((tiq(a, l), yv / (np.linalg.norm(S[kk]) / np.linalg.norm(Rbar[:, sl])), b))
    def hp(series):
        T = np.array([x[0] for x in series]); V = np.array([x[1] for x in series]); C = np.concatenate([np.zeros((1,) + V.shape[1:], V.dtype), np.cumsum(V, 0)])
        lo = np.searchsorted(T, T - HP, "left"); hi = np.searchsorted(T, T + HP, "right"); M = (C[hi] - C[lo]) / (hi - lo)[(slice(None),) + (None,) * (V.ndim - 1)]
        return [(series[j][0], V[j] - M[j], *series[j][2:]) for j in range(len(series))]
    trs_h = hp(trs); dm_h = []
    for ser in dm.values():
        ser.sort(key=lambda x: x[0])
        if len(ser) >= 5: dm_h += hp(ser)
    def truth_at(te, off):
        out = []
        for tid in IDS:
            m = TR[:, 0] == tid; tk = TR[m, 1]
            if tk[0] <= te <= tk[-1]: out.append((tid, np.interp(te, tk, TR[m, 2]) + off))
        return out
    # 3. per window: normal equations for 1- and 2-unknown models, per antenna
    for st, rows in (("trs", trs_h), ("dmrs", dm_h)):
        W = collections.defaultdict(lambda: np.zeros((4, 2), complex))   # rhs per antenna: [sum conj(m) d, sum conj(phi m) d]
        WG = collections.defaultdict(lambda: np.zeros((2, 2), complex))  # gram: [[|m|^2, phi |m|^2], [conj(phi)|m|^2, |m|^2]]
        for r in rows:
            t = r[0]; d = r[1]
            if st == "trs": kk_ = KKT - 1638; ph = 1.0
            else: b = r[2]; kk_ = KKD[12 * b:12 * b + 12] - 1638; ph = PHI[b]
            te = t_rec + OFF + t
            for name, off in HYP:
                for tid, dr in truth_at(te, 0.0):
                    wc = (int(t // WIN) + 0.5) * WIN; tc = t_rec + OFF + wc
                    dc = np.interp(tc, TR[TR[:, 0] == tid, 1], TR[TR[:, 0] == tid, 2])
                    if name == "rev": ph_path = 2 * dc - dr                       # phase history mirrored about the window centre
                    else: ph_path = dr - 2 * off * (t - wc)                        # extra bistatic speed off (m/s): excess path changes by 2 off per s
                    m = np.exp(-2j * np.pi * (kk_ * DF * dr / 3e8 + ph_path / LAM))
                    key = (name, tid, int(t // WIN), (b if st == "dmrs" else -1)); n2 = np.sum(np.abs(m) ** 2)
                    W[key][:, 0] += d @ np.conj(m); W[key][:, 1] += d @ np.conj(ph * m)
                    WG[key] += np.array([[n2, ph * n2], [np.conj(ph) * n2, n2]])
        for key, rhs in W.items():
            Gm = WG[key]; name = key[0]
            e1 = np.sum(np.abs(rhs[:, 0]) ** 2) / Gm[0, 0].real                       # per group (DM-RS) or full band (TRS): a free factor
            ACC[(st, "1", name)] += e1
            if name == "truth": NW[st] += 1
    print("chunk %d done (%d groups with phi; %d DM-RS PRG-snapshots)" % (k, len(PHI), len(dm_h)), flush=True)
print("phi vote stability (share of the winning value per group): median %.2f" % np.median(STAB))
for st, mdl in (("trs", "1"), ("dmrs", "1")):
    tr_ = ACC[(st, mdl, "truth")]; nul = np.mean([ACC[(st, mdl, n)] for n, _ in HYP if n != "truth"])
    print("%-4s (%s): truth / Doppler nulls = %.2f   (truth %.3g, nulls %.3g; per null: %s)" % (st, "full band coherent" if st == "trs" else "free factor per 4-PRB group", tr_ / nul, tr_, nul, {n: round(ACC[(st, mdl, n)] / tr_, 2) for n, _ in HYP if n != "truth"}))
