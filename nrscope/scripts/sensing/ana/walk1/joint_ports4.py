"""joint_ports2: conditioning of the port fits + structured precoder models. Coherent TRS + DM-RS combination for a target through the CSI-RS port space (walk1 prototype; lidar = the track).
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
    # streams: DM-RS (dominant beam per group, beam frame) and TRS (per group)
    nW = int(2.5 // WIN) + 1
    rd = np.zeros((H_, nW, nB, 4), complex); nd = np.zeros((H_, nW, nB)); rt = np.zeros((H_, nW, nB, 4), complex); nt = np.zeros((H_, nW, nB))
    pwd = [0.0, 0]; pwt = [0.0, 0]
    ser = collections.defaultdict(list)
    for (i, b), kk in own.items():
        if best.get(b) != kk: continue
        t, a, l, pm, Y, u = E[i]; sl = slice(12 * b, 12 * b + 12); yv = Y[:, sl] * np.exp(-1j * TH[(i, b)])
        if np.sum(np.abs(yv - S[kk]) ** 2) > GATE * np.sum(np.abs(S[kk]) ** 2): continue
        ser[b].append((tiq(a, l), yv))
    for b, Ls in ser.items():
        Ls.sort(key=lambda x: x[0]); T = np.array([x[0] for x in Ls]); Dd = hp_series(T, np.array([x[1] for x in Ls])); kk_ = KKD[12 * b:12 * b + 12] - 1638
        pwd[0] += float(np.sum(np.abs(Dd) ** 2)); pwd[1] += Dd.size
        for j, t in enumerate(T):
            dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
            if dr is None or not 0 <= w < nW: continue
            for h, off in enumerate(HYP):
                m = np.exp(-2j * np.pi * (kk_ * DF * (dr + off) / 3e8 + dr / LAM)); rd[h, w, bi[b]] += Dd[j] @ np.conj(m); nd[h, w, bi[b]] += np.sum(np.abs(m) ** 2)
    TT = sorted([(tiq(a, l), R / osc(t)) for t, a, l, R in TRS], key=lambda x: x[0]); Tq = np.array([x[0] for x in TT]); Dt = hp_series(Tq, np.array([x[1] for x in TT]))
    pwt[0] += float(np.sum(np.abs(Dt) ** 2)); pwt[1] += Dt.size
    for j, t in enumerate(Tq):
        dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
        if dr is None or not 0 <= w < nW: continue
        for h, off in enumerate(HYP):
            m = np.exp(-2j * np.pi * ((KKT - 1638) * DF * (dr + off) / 3e8 + dr / LAM))
            for b in B:
                sl = slice(12 * b, 12 * b + 12); rt[h, w, bi[b]] += Dt[j][:, sl] @ np.conj(m[sl]); nt[h, w, bi[b]] += np.sum(np.abs(m[sl]) ** 2)
    sd2 = pwd[0] / max(pwd[1], 1); st2 = pwt[0] / max(pwt[1], 1)      # per-value residual power (noise proxy) of each stream
    valid = [w for w in range(nW) if np.sum(nd[0, w] > 0) >= 20 and np.sum(nt[0, w]) > 0]
    # ---- conditioning of H_b (16 x 8) ----
    sv = np.array([np.linalg.svd(Hc[:, 4 * b:4 * b + 4, :].reshape(16, 8), compute_uv=False) for b in B]); sv = sv / sv[:, :1]
    print("chunk %d: H_b singular values / largest, median over groups: %s" % (k, " ".join("%.3f" % x for x in np.median(sv, 0))), flush=True)
    # ---- Type I structured fit: v_b = [x_b u ; y_b u], u shared per beam set (kk) ----
    def cb_fit(perm):
        X = np.zeros((nB, 2), complex); res = []
        for kk in [0]:
            bs = list(B)
            Ha = {b: Hc[:, 4 * b:4 * b + 4, perm[:4]].reshape(16, 4) for b in bs}; Hb_ = {b: Hc[:, 4 * b:4 * b + 4, perm[4:]].reshape(16, 4) for b in bs}
            ys = {b: S[best[b]].reshape(4, 4, 3).mean(2).reshape(16) for b in bs}
            C = sum(np.outer(Vb[bi[b]][perm[:4]], Vb[bi[b]][perm[:4]].conj()) + np.outer(Vb[bi[b]][perm[4:]], Vb[bi[b]][perm[4:]].conj()) for b in bs)
            u = np.linalg.eigh(C)[1][:, -1]
            for _ in range(40):
                xy = {}
                for b in bs:
                    A2 = np.stack([Ha[b] @ u, Hb_[b] @ u], 1); xy[b], *_ = np.linalg.lstsq(A2, ys[b], rcond=None)
                Au = np.concatenate([xy[b][0] * Ha[b] + xy[b][1] * Hb_[b] for b in bs]); u, *_ = np.linalg.lstsq(Au, np.concatenate([ys[b] for b in bs]), rcond=None)
                u /= np.linalg.norm(u)
            for b in bs:
                A2 = np.stack([Ha[b] @ u, Hb_[b] @ u], 1); xy_, *_ = np.linalg.lstsq(A2, ys[b], rcond=None); X[bi[b]] = xy_
                res.append(np.sum(np.abs(ys[b] - A2 @ xy_) ** 2) / np.sum(np.abs(ys[b]) ** 2))
        return X, np.median(res)
    fits = {}
    for nm, perm in (("0-3|4-7", [0, 1, 2, 3, 4, 5, 6, 7]), ("even|odd", [0, 2, 4, 6, 1, 3, 5, 7])):
        fits[nm] = cb_fit(perm); print("chunk %d: Type I fit (pairing %s): residual %.1f dB" % (k, nm, 10 * np.log10(fits[nm][1])), flush=True)
    Xcb = min(fits.values(), key=lambda z: z[1])[0]
    QP = np.array([1, 1j, -1, -1j]); r_ = Xcb[:, 1] / (Xcb[:, 0] + 1e-30)
    qd = np.abs(np.angle(r_[:, None] / QP[None, :])).min(1)
    print("chunk %d: continuous co-phasing y/x: |ratio| median %.2f (p10 %.2f p90 %.2f); distance to QPSK median %.2f rad (uniform 0.39)"
          % (k, np.median(np.abs(r_)), np.percentile(np.abs(r_), 10), np.percentile(np.abs(r_), 90), np.median(qd)), flush=True)
    phq = QP[np.abs(np.angle(r_[:, None] / QP[None, :])).argmin(1)]
    # discrete codebook: u = DFT (N=4, O=4, 16 beams) or 2x2 (O1=O2=4, 64 beams); phi_b QPSK; one common scalar per beam set
    def dft_cands():
        out = []
        for m in range(16): out.append(("1D%d" % m, np.exp(2j * np.pi * m * np.arange(4) / 16) / 2))
        for m1 in range(8):
            for m2 in range(8):
                a1 = np.exp(2j * np.pi * m1 * np.arange(2) / 8); a2 = np.exp(2j * np.pi * m2 * np.arange(2) / 8)
                out.append(("2D%d,%d" % (m1, m2), np.kron(a1, a2) / 2))
        return out
    CANDS = dft_cands(); phd = np.ones(nB, complex); dres = []; Rnd = np.random.default_rng(1)
    for kk in [0]:
        bs = list(B); YS = {b: S[best[b]].reshape(4, 4, 3).mean(2).reshape(16) for b in B}
        bestc = None
        for perm in ([0, 1, 2, 3, 4, 5, 6, 7], [0, 2, 4, 6, 1, 3, 5, 7]):
            for nm, u in CANDS + [("rand", (Rnd.standard_normal(4) + 1j * Rnd.standard_normal(4)) / 2.8)]:
                A1 = {b: Hc[:, 4 * b:4 * b + 4, perm[:4]].reshape(16, 4) @ u for b in bs}; A2 = {b: Hc[:, 4 * b:4 * b + 4, perm[4:]].reshape(16, 4) @ u for b in bs}
                ph = {}; num = 0; den = 0
                for b in bs:
                    ys = YS[b]; e = [np.linalg.norm(ys - (A1[b] + q * A2[b]) * (np.vdot(A1[b] + q * A2[b], ys) / max(np.vdot(A1[b] + q * A2[b], A1[b] + q * A2[b]).real, 1e-30))) for q in QP]
                    ph[b] = QP[int(np.argmin(e))]; A = A1[b] + ph[b] * A2[b]; num += np.vdot(A, YS[b]); den += np.vdot(A, A).real
                sc = num / den; res = sum(np.sum(np.abs(YS[b] - sc * (A1[b] + ph[b] * A2[b])) ** 2) for b in bs) / sum(np.sum(np.abs(YS[b]) ** 2) for b in bs)
                if nm == "rand": print("   random-u control (pairing %s): %.1f dB" % (perm[:2], 10 * np.log10(res))); continue
                if bestc is None or res < bestc[0]: bestc = (res, nm, str(perm[:2]), ph)
        dres.append(bestc[0]); [phd.__setitem__(bi[b], bestc[3][b]) for b in bs]
        print("chunk %d beam set %s: best codebook beam %s pairing %s, residual %.1f dB" % (k, kk, bestc[1], bestc[2], 10 * np.log10(bestc[0])), flush=True)
    agree = np.mean(np.isclose(phd, phq)); print("chunk %d: QPSK co-phasing agreement (codebook search vs quantised continuous) %.2f" % (k, agree), flush=True)
    Gst = np.array([np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), S[best[b]].ravel()) / np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), Rbar[:, 12 * b:12 * b + 12].ravel()).real for b in B])
    np.savez("/mnt/mydrive/field/ana_walk1/jp_cache_b%d.npz" % k, rd=rd, nd=nd, rt=rt, nt=nt, Vb=Vb, Wb=Wb, Xcb=Xcb, phd=phd, phq=phq, r_=r_, Gst=Gst, sd2=sd2, st2=st2, valid=np.array(valid), B=np.array(B))
    Z = lambda n: np.zeros((nB, n), complex); O = np.ones((nB, 1), complex)
    MODELS = {"free8": (Vb, Wb), "typeI": (np.hstack([Xcb, Z(1)]), np.hstack([Z(2), O])),
              "static": (np.hstack([Gst[:, None], Z(1)]), np.hstack([Z(1), O])), "flat": (np.hstack([O, Z(1)]), np.hstack([Z(1), O])),
              "pergroup": (np.hstack([np.eye(nB), Z(1)]), np.hstack([Z(nB), O])),
              "cbq": (np.hstack([O, phq[:, None], Z(1)]), np.hstack([Z(2), O])), "cbdft": (np.hstack([O, phd[:, None], Z(1)]), np.hstack([Z(2), O])),
              "cbph": (np.hstack([O, np.exp(1j * np.angle(r_))[:, None], Z(1)]), np.hstack([Z(2), O]))}
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
        for wt_ in valid:
            train = list(valid) if os.environ.get("ORACLE") == "1" else [w for w in valid if w != wt_]
            if len(train) < 4: continue
            cd, ct, a = learn(train, Fd, Ft)
            for h in range(H_):
                for lab, pl in {"trs": [(ct, rt, nt, st2)], "dmrs": [(cd, rd, nd, sd2)], "both": [(ct, rt, nt, st2), (cd, rd, nd, sd2)]}.items():
                    num = 0; den = 0
                    for (cc, rr, nn, s2) in pl:
                        M = cc[:, None] * a[None, :]; num += np.sum(np.conj(M) * rr[h, wt_]) / s2; den += np.sum(np.abs(M) ** 2 * nn[h, wt_][:, None]) / s2
                    TOT[(mn, lab, h)] += abs(num) ** 2 / max(den, 1e-30)
    print("chunk %d: %d groups, %d test windows; noise per value DM-RS %.3g, TRS %.3g" % (k, nB, len(valid), sd2, st2), flush=True)
for mn in ("free8", "typeI", "static", "flat", "pergroup", "cbq", "cbph", "cbdft"):
    out = []
    for lab in ("trs", "dmrs", "both"):
        tr_ = TOT[(mn, lab, 0)]; nl = np.mean([TOT[(mn, lab, h)] for h in range(1, H_)]); out.append("%s %.2f" % (lab, tr_ / max(nl, 1e-30)))
    print("%-9s truth/nulls: %s" % (mn, "   ".join(out)))
