"""joint_ports6 (per-transmission precoder v_b(i) from its own static, truncated SVD of H_b; target factor t . v_b(i)): conditioning of the port fits + structured precoder models. Coherent TRS + DM-RS combination for a target through the CSI-RS port space (walk1 prototype; lidar = the track).
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
    nW = int(2.5 // WIN) + 1; RO = (4, "4s", "3s", "fixed")
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
    SVDb = {b: np.linalg.svd(Hc[:, 4 * b:4 * b + 4, :].reshape(16, 8), full_matrices=False) for b in B}
    def vraw(b, y, ro):
        U, s_, Vh = SVDb[b]; return Vh[:ro].conj().T @ ((U[:, :ro].conj().T @ y.reshape(4, 4, 3).mean(2).reshape(16)) / s_[:ro])
    SM = {}   # smoothed per-transmission precoder: mean of the raw estimates within +-SMW s on the same group
    SMW = float(os.environ.get("SMW", "0.05"))
    for b_, Ls_ in ser.items():
        if b_ not in SVDb: continue
        T_ = np.array([x[0] for x in Ls_])
        for ro_ in (3, 4):
            Vr = np.array([vraw(b_, x[1], ro_) for x in Ls_]); C_ = np.concatenate([np.zeros((1, 8), complex), np.cumsum(Vr, 0)])
            lo_ = np.searchsorted(T_, T_ - SMW); hi_ = np.searchsorted(T_, T_ + SMW, "right")
            for j_, x in enumerate(Ls_): SM[(b_, ro_, x[0])] = (C_[hi_[j_]] - C_[lo_[j_]]) / (hi_[j_] - lo_[j_])
    def vhat(b, y, ro, t=None):
        if ro == "fixed": return Vb[bi[b]]
        if isinstance(ro, str): return SM[(b, int(ro[0]), t)]
        return vraw(b, y, ro)
    TT = sorted([(tiq(a, l), R / osc(t)) for t, a, l, R in TRS], key=lambda x: x[0]); Tq = np.array([x[0] for x in TT]); Vt = np.array([x[1] for x in TT])
    def build(r):
        Rv = {ro: np.zeros((H_, nW, nB, 4, 8), complex) for ro in RO}; Gv = {ro: np.zeros((nW, nB, 8, 8), complex) for ro in RO}
        rt = np.zeros((H_, nW, nB, 4), complex); nt = np.zeros((H_, nW, nB)); cnt = np.zeros((nW, nB)); pd = [0.0, 0]; pt = [0.0, 0]
        for b, Ls in ser.items():
            T = np.array([x[0] for x in Ls]); Dd = hp_series(T, pca_cut(np.array([x[1] for x in Ls]), r)); kk_ = KKD[12 * b:12 * b + 12] - 1638
            pd[0] += float(np.sum(np.abs(Dd) ** 2)); pd[1] += Dd.size
            for j, t in enumerate(T):
                dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
                if dr is None or not 0 <= w < nW: continue
                vs = {ro: vhat(b, Ls[j][1], ro, t) for ro in RO}; cnt[w, bi[b]] += 1
                for ro in RO: Gv[ro][w, bi[b]] += 12 * np.outer(vs[ro].conj(), vs[ro])
                for h, off in enumerate(HYP):
                    m = np.exp(-2j * np.pi * (kk_ * DF * (dr + off) / 3e8 + dr / LAM)); pj = Dd[j] @ np.conj(m)
                    for ro in RO: Rv[ro][h, w, bi[b]] += pj[:, None] * vs[ro].conj()[None, :]
        for b in B:
            sl = slice(12 * b, 12 * b + 12); Dt = hp_series(Tq, pca_cut(Vt[:, :, sl], r)); kt_ = KKT[sl] - 1638
            pt[0] += float(np.sum(np.abs(Dt) ** 2)); pt[1] += Dt.size
            for j, t in enumerate(Tq):
                dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
                if dr is None or not 0 <= w < nW: continue
                for h, off in enumerate(HYP):
                    m = np.exp(-2j * np.pi * (kt_ * DF * (dr + off) / 3e8 + dr / LAM)); rt[h, w, bi[b]] += Dt[j] @ np.conj(m); nt[h, w, bi[b]] += 12
        return Rv, Gv, rt, nt, cnt, pd[0] / max(pd[1], 1), pt[0] / max(pt[1], 1)
    for r in (0,):
        Rv, Gv, rt, nt, cnt, sd2, st2 = build(r)
        valid = [w for w in range(nW) if np.sum(cnt[w] > 0) >= 20 and np.sum(nt[0, w]) > 0]
        for ro in RO:
            R_, G_ = Rv[ro], Gv[ro]
            def learn(train):
                kap = np.ones(len(train), complex); a = np.ones(4, complex); t = np.zeros(8, complex); gam = 1.0 + 0j
                t = np.linalg.lstsq(np.sum(G_[train], (0, 1)) + 1e-9 * np.eye(8), np.sum(R_[0, train], (0, 1, 2)), rcond=None)[0]
                for _ in range(30):
                    A = sum(abs(kap[wi]) ** 2 * np.sum(np.abs(a) ** 2) * np.sum(G_[w], 0) for wi, w in enumerate(train)) / sd2
                    bv = sum(np.conj(kap[wi]) * np.einsum("c,bcv->v", a.conj(), R_[0, w]) for wi, w in enumerate(train)) / sd2
                    t = np.linalg.solve(A + 1e-12 * np.trace(A).real * np.eye(8) + 1e-30 * np.eye(8), bv)
                    gam = sum(np.conj(kap[wi]) * np.sum(a.conj()[None, :] * rt[0, w]) for wi, w in enumerate(train)) / max(sum(abs(kap[wi]) ** 2 * np.sum(np.abs(a) ** 2) * np.sum(nt[0, w]) for wi, w in enumerate(train)), 1e-30)
                    for wi, w in enumerate(train):
                        num = np.sum(a.conj()[None, :] * (R_[0, w] @ t.conj()).conj()) / sd2 if False else np.sum(a.conj()[None, :] * np.einsum("bcv,v->bc", R_[0, w], t.conj())) / sd2
                        num += np.conj(gam) * np.sum(a.conj()[None, :] * rt[0, w]) / st2
                        den = np.sum(np.abs(a) ** 2) * (np.real(np.einsum("v,bvu,u->", t.conj(), G_[w], t)) / sd2 + abs(gam) ** 2 * np.sum(nt[0, w]) / st2)
                        kap[wi] = num / max(den, 1e-30)
                    num = np.zeros(4, complex); den = 0.0
                    for wi, w in enumerate(train):
                        num += np.conj(kap[wi]) * (np.einsum("bcv,v->c", R_[0, w], t.conj()) / sd2 + np.conj(gam) * np.sum(rt[0, w], 0) / st2)
                        den += abs(kap[wi]) ** 2 * (np.real(np.einsum("v,bvu,u->", t.conj(), G_[w], t)) / sd2 + abs(gam) ** 2 * np.sum(nt[0, w]) / st2)
                    a = num / max(den, 1e-30); a /= np.sqrt(np.mean(np.abs(a) ** 2)) + 1e-30
                return t, gam, a
            for orc in (0, 1):
                for wt_ in valid:
                    train = list(valid) if orc else [w for w in valid if w != wt_]
                    if len(train) < 4: continue
                    t, gam, a = learn(train)
                    for h in range(H_):
                        nd_ = np.einsum("c,bc->", a.conj(), np.einsum("bcv,v->bc", R_[h, wt_], t.conj())) / sd2; dd_ = np.sum(np.abs(a) ** 2) * np.real(np.einsum("v,bvu,u->", t.conj(), G_[wt_], t)) / sd2
                        nt_ = np.conj(gam) * np.sum(a.conj()[None, :] * rt[h, wt_]) / st2; dt_ = abs(gam) ** 2 * np.sum(np.abs(a) ** 2) * np.sum(nt[h, wt_]) / st2
                        TOT[(r, ro, orc, "trs", h)] += abs(nt_) ** 2 / max(dt_, 1e-30); TOT[(r, ro, orc, "dmrs", h)] += abs(nd_) ** 2 / max(dd_, 1e-30)
                        TOT[(r, ro, orc, "both", h)] += abs(nt_ + nd_) ** 2 / max(dt_ + dd_, 1e-30)
        print("chunk %d r=%d done (noise DM-RS %.3g TRS %.3g)" % (k, r, sd2, st2), flush=True)
for r in (0,):
    for orc in (0, 1):
        for ro in RO:
            out = []
            for lab in ("trs", "dmrs", "both"):
                tr_ = TOT[(r, ro, orc, lab, 0)]; nl = np.mean([TOT[(r, ro, orc, lab, h)] for h in range(1, H_)]); out.append("%s %.2f" % (lab, tr_ / max(nl, 1e-30)))
            print("clutter r=%d %-7s precoder %-5s truth/nulls: %s" % (r, "oracle" if orc else "heldout", ro, "   ".join(out)))
