"""joint_ports5 (eigen-clutter cancellation per group, r = 0,1,2,4,8): conditioning of the port fits + structured precoder models. Coherent TRS + DM-RS combination for a target through the CSI-RS port space (walk1 prototype; lidar = the track).
Port channels H (8-port CSI-RS, slot 0, receiver corrections as so_light), per PRB. Per 4-PRB group b:
  DM-RS beam (tablet, dominant per group):  S_b = H_b v_b          (v_b fitted)
  TRS:                                      R_b = H_b w_b          (w_b fitted per group; also one wideband w)
A target with port signature t (8) and receive steering a (4): its echo in a window w is
  DM-RS group b, antenna c:  kappa_w a_c (t . v_b) m      TRS group b, antenna c:  kappa_w a_c (t . w_b) m
(m: lidar range steering + Doppler phase). t, a learned on the OTHER windows of the chunk from both streams (weighted
by each stream's noise), then on the held-out window, per hypothesis (truth range, +12..+40 m):
  TRS alone / DM-RS alone / both combined: one complex amplitude fitted with the learned model, energy of the fit.
usage: joint_ports.py <chunks> <offset>"""
import sys, os, io, json, contextlib, collections, numpy as np
KS = [int(x) for x in sys.argv[1].split(",")]; OFF = float(sys.argv[2]); HP = 0.25; WIN = 0.2; GATE = 10 ** (-12 / 10)
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
GD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth/grid"
D = np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz"); TR = D["tr"]; t_rec = float(D["t_rec"]); IDS = sorted(set(TR[:, 0].astype(int)))
LAM = 3e8 / 3.75e9; DF = 30e3; HYP = [0.0, 12.0, 18.0, 25.0, 32.0, 40.0]; H_ = len(HYP)
TOT = collections.defaultdict(float)
def walker_at(te):
    for tid in IDS:
        m = TR[:, 0] == tid; tk = TR[m, 1]
        if tk[0] <= te <= tk[-1] and tid in (1, 4): return np.interp(te, tk, TR[m, 2])
    return None
def hp_series(T, V):
    C = np.concatenate([np.zeros((1,) + V.shape[1:], complex), np.cumsum(V, 0)]); lo = np.searchsorted(T, T - HP); hi = np.searchsorted(T, T + HP, "right")
    return V - (C[hi] - C[lo]) / (hi - lo)[(slice(None),) + (None,) * (V.ndim - 1)]
for k in KS:
    os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
    sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/jp.bin"]
    g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
        "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
    E, own, TH, S, Sc, Rbar = g["KEEP"]; KKD, KKT, TRS = g["KK_D"], g["KK_T"], g["TRS"]
    osc, files, shifts, drift, b0, PHs, Nn, NFq = g["osc"], g["files"], g["shifts"], g["drift"], g["b0"], g["PH"], g["N"], g["NFq"]
    a0 = json.load(open("%s/b%d.json" % (GD, k)))["a0"]; tiq = lambda a, l: k * 2.5 + ((a - a0) % 20480) * 5e-4 + l / 28e3
    # CSI-RS port channels
    import gs_trs_lib as L
    qc = 3 * np.arange(Nn) + 1; qd0 = 409.25
    Gc = [L.load_gsc(f) for f in files]; kc = sorted(set(Gc[0]).intersection(*[set(x) for x in Gc[1:]])); PP = []
    for aa, ll in kc:
        if ll != 13 or aa % 20 != 0: continue
        tt = aa * 0.5e-3 + 13 / 28e3; Pm = np.array([L.csi8_despread(*Gc[c][(aa, 13)], aa % 20)[:Nn] for c in range(4)])
        sh = np.array(shifts(aa), float); Pm = Pm * np.exp(-1j * PHs * sh[:, None, None] * (qc - qd0)[None, :, None])
        Pm = Pm * np.exp(-1j * np.array([np.polyval(drift[c], tt) for c in range(4)])[:, None, None] * (qc - qd0)[None, :, None])
        Pm = Pm * np.exp(2j * np.pi * np.array(b0)[:, None, None] * (qc - qd0)[None, :, None] / NFq); PP.append(Pm / osc(tt))
    ref = PP[len(PP) // 2]; Hc = np.mean([p * np.exp(-1j * np.angle(np.vdot(ref.ravel(), p.ravel()))) for p in PP], 0)   # [4, N, 8]
    # dominant tablet beam per group, its v_b; the TRS w_b
    cnt = collections.Counter(own.values()); best = {}
    for (i, b), kk in own.items():
        if Sc[kk] >= 50 and (b not in best or cnt[kk] > cnt[best[b]]): best[b] = kk
    B = sorted(best); bi = {b: j for j, b in enumerate(B)}; nB = len(B)
    def pfit(Xprb, b):
        A_ = Hc[:, 4 * b:4 * b + 4, :].reshape(16, 8); y_ = Xprb.reshape(16); v, *_ = np.linalg.lstsq(A_, y_, rcond=None)
        return v, np.sum(np.abs(y_ - A_ @ v) ** 2) / np.sum(np.abs(y_) ** 2)
    Vb = np.zeros((nB, 8), complex); Wb = np.zeros((nB, 8), complex); rv, rw = [], []
    for b in B:
        v, e = pfit(S[best[b]].reshape(4, 4, 3).mean(2), b); Vb[bi[b]] = v; rv.append(e)
        w, e = pfit(Rbar[:, 12 * b:12 * b + 12].reshape(4, 4, 3).mean(2), b); Wb[bi[b]] = w; rw.append(e)
    Aw = Hc.reshape(4 * Nn, 8)[np.repeat(np.arange(Nn), 1) * 0 + np.arange(4 * Nn) * 0 + np.arange(4 * Nn)] if False else Hc.transpose(1, 0, 2).reshape(-1, 8)
    wwide, *_ = np.linalg.lstsq(Aw, Rbar.reshape(4, Nn, 3).mean(2).T.reshape(-1), rcond=None)
    rwide = np.sum(np.abs(Rbar.reshape(4, Nn, 3).mean(2).T.reshape(-1) - Aw @ wwide) ** 2) / np.sum(np.abs(Rbar) ** 2 / 3)
    print("chunk %d: port fits -- DM-RS beams %.1f dB, TRS per group %.1f dB, TRS one wideband w %.1f dB" % (k, 10 * np.log10(np.median(rv)), 10 * np.log10(np.median(rw)), 10 * np.log10(rwide)), flush=True)
    nW = int(2.5 // WIN) + 1
    def pca_cut(V, r):
        if r == 0: return V
        X = V.reshape(len(V), -1); _, _, Vh = np.linalg.svd(X, full_matrices=False); P = Vh[:r]
        return (X - (X @ P.conj().T) @ P).reshape(V.shape)
    ser = collections.defaultdict(list)
    for (i, b), kk in own.items():
        if best.get(b) != kk: continue
        t, a, l, pm, Y, u = E[i]; sl = slice(12 * b, 12 * b + 12); yv = Y[:, sl] * np.exp(-1j * TH[(i, b)])
        if np.sum(np.abs(yv - S[kk]) ** 2) > GATE * np.sum(np.abs(S[kk]) ** 2): continue
        ser[b].append((tiq(a, l), yv))
    for b in ser: ser[b].sort(key=lambda x: x[0])
    TT = sorted([(tiq(a, l), R / osc(t)) for t, a, l, R in TRS], key=lambda x: x[0]); Tq = np.array([x[0] for x in TT]); Vt = np.array([x[1] for x in TT])
    WR = {}
    for j, t in enumerate(Tq):
        dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
        if dr is not None and 0 <= w < nW: WR[("t", j)] = (dr, w)
    def build(r):
        rd = np.zeros((H_, nW, nB, 4), complex); nd = np.zeros((H_, nW, nB)); rt = np.zeros((H_, nW, nB, 4), complex); nt = np.zeros((H_, nW, nB)); pd = [0.0, 0]; pt = [0.0, 0]
        for b, Ls in ser.items():
            T = np.array([x[0] for x in Ls]); Dd = hp_series(T, pca_cut(np.array([x[1] for x in Ls]), r)); kk_ = KKD[12 * b:12 * b + 12] - 1638
            pd[0] += float(np.sum(np.abs(Dd) ** 2)); pd[1] += Dd.size
            for j, t in enumerate(T):
                dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
                if dr is None or not 0 <= w < nW: continue
                for h, off in enumerate(HYP):
                    m = np.exp(-2j * np.pi * (kk_ * DF * (dr + off) / 3e8 + dr / LAM)); rd[h, w, bi[b]] += Dd[j] @ np.conj(m); nd[h, w, bi[b]] += np.sum(np.abs(m) ** 2)
        for b in B:
            sl = slice(12 * b, 12 * b + 12); Dt = hp_series(Tq, pca_cut(Vt[:, :, sl], r)); kt_ = KKT[sl] - 1638
            pt[0] += float(np.sum(np.abs(Dt) ** 2)); pt[1] += Dt.size
            for (_, j), (dr, w) in WR.items():
                for h, off in enumerate(HYP):
                    m = np.exp(-2j * np.pi * (kt_ * DF * (dr + off) / 3e8 + dr / LAM)); rt[h, w, bi[b]] += Dt[j] @ np.conj(m); nt[h, w, bi[b]] += np.sum(np.abs(m) ** 2)
        return rd, nd, rt, nt, pd[0] / max(pd[1], 1), pt[0] / max(pt[1], 1)
    Gst = np.array([np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), S[best[b]].ravel()) / np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), Rbar[:, 12 * b:12 * b + 12].ravel()).real for b in B])
    Z = lambda n: np.zeros((nB, n), complex); O = np.ones((nB, 1), complex)
    MODELS = {"free8": (Vb, Wb), "static": (np.hstack([Gst[:, None], Z(1)]), np.hstack([Z(1), O])), "flat": (np.hstack([O, Z(1)]), np.hstack([Z(1), O])),
              "pergroup": (np.hstack([np.eye(nB), Z(1)]), np.hstack([Z(nB), O]))}
    for r in (0, 1, 2, 4, 8):
        rd, nd, rt, nt, sd2, st2 = build(r)
        valid = [w for w in range(nW) if np.sum(nd[0, w] > 0) >= 20 and np.sum(nt[0, w]) > 0]
        def learn(train, Fd, Ft):
            P = Fd.shape[1]; kap = np.ones(len(train), complex); a = np.ones(4, complex)
            t, *_ = np.linalg.lstsq(Fd + Ft, np.ones(nB, complex), rcond=None)
            for _ in range(30):
                X = []; Y = []; Wt = []
                for wi, w in enumerate(train):
                    for (Pm_, rr, nn, s2) in ((Fd, rd, nd, sd2), (Ft, rt, nt, st2)):
                        gm = rr[0, w] / np.maximum(nn[0, w], 1e-9)[:, None]; wt = nn[0, w] / s2
                        X.append((kap[wi] * a[None, :, None] * Pm_[:, None, :]).reshape(-1, P)); Y.append(gm.reshape(-1)); Wt.append(np.repeat(wt, 4))
                X = np.concatenate(X); Y = np.concatenate(Y); Wt = np.sqrt(np.concatenate(Wt))
                t, *_ = np.linalg.lstsq(X * Wt[:, None], Y * Wt, rcond=None)
                cd = Fd @ t; ct = Ft @ t
                for wi, w in enumerate(train):
                    num = 0; den = 0
                    for (cc, rr, nn, s2) in ((cd, rd, nd, sd2), (ct, rt, nt, st2)):
                        M = cc[:, None] * a[None, :]; num += np.sum(np.conj(M) * rr[0, w]) / s2; den += np.sum(np.abs(M) ** 2 * nn[0, w][:, None]) / s2
                    kap[wi] = num / max(den, 1e-30)
                num = np.zeros(4, complex); den = np.zeros(4)
                for wi, w in enumerate(train):
                    for (cc, rr, nn, s2) in ((cd, rd, nd, sd2), (ct, rt, nt, st2)):
                        M = kap[wi] * cc[:, None]; num += np.sum(np.conj(M) * rr[0, w], 0) / s2; den += np.sum(np.abs(M) ** 2 * nn[0, w][:, None], 0) / s2
                a = num / np.maximum(den, 1e-30); a /= np.sqrt(np.mean(np.abs(a) ** 2)) + 1e-30
            return Fd @ t, Ft @ t, a
        for mn, (Fd, Ft) in MODELS.items():
            for orc in (0, 1):
                for wt_ in valid:
                    train = list(valid) if orc else [w for w in valid if w != wt_]
                    if len(train) < 4: continue
                    cd, ct, a = learn(train, Fd, Ft)
                    for h in range(H_):
                        for lab, pl in {"trs": [(ct, rt, nt, st2)], "dmrs": [(cd, rd, nd, sd2)], "both": [(ct, rt, nt, st2), (cd, rd, nd, sd2)]}.items():
                            num = 0; den = 0
                            for (cc, rr, nn, s2) in pl:
                                M = cc[:, None] * a[None, :]; num += np.sum(np.conj(M) * rr[h, wt_]) / s2; den += np.sum(np.abs(M) ** 2 * nn[h, wt_][:, None]) / s2
                            TOT[(r, mn, orc, lab, h)] += abs(num) ** 2 / max(den, 1e-30)
        print("chunk %d r=%d: noise per value DM-RS %.3g, TRS %.3g" % (k, r, sd2, st2), flush=True)
for r in (0, 1, 2, 4, 8):
    for orc in (0, 1):
        for mn in ("free8", "static", "flat", "pergroup"):
            out = []
            for lab in ("trs", "dmrs", "both"):
                tr_ = TOT[(r, mn, orc, lab, 0)]; nl = np.mean([TOT[(r, mn, orc, lab, h)] for h in range(1, H_)]); out.append("%s %.2f" % (lab, tr_ / max(nl, 1e-30)))
            print("r=%d %-7s %-9s truth/nulls: %s" % (r, "oracle" if orc else "heldout", mn, "   ".join(out)))
