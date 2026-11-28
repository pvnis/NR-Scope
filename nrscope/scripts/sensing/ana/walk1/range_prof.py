"""range_prof: walker range profiles (TRS, DM-RS anchored = divided by D_b, DM-RS not divided) and the delay profile of D_b in the physical frame. Base: joint_ports2: conditioning of the port fits + structured precoder models. Coherent TRS + DM-RS combination for a target through the CSI-RS port space (walk1 prototype; lidar = the track).
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
    DEL = np.arange(-60, 60.5, 1.0); Mdel = None
    # --- D_b per transmission in the physical frame (common phase per run removed; Rbar = TRS static) ---
    fbs = np.array(B, float) * 48 * 30e3; taus = np.arange(-1.5e-6, 1.5e-6, 5e-9); Ftau = np.exp(2j * np.pi * fbs[:, None] * taus[None, :])
    Dprof = np.zeros(len(taus)); nD = 0; byi = collections.defaultdict(dict)
    for (i, b), kk in own.items():
        if b in bi: byi[i][b] = kk
    for i, gb in byi.items():
        if len(gb) < 0.25 * nB: continue
        t, a, l, pm, Y, u = E[i]; D = np.zeros(nB, complex)
        for b, kk in gb.items():
            Rb = Rbar[:, 12 * b:12 * b + 12]; D[bi[b]] = np.vdot(Rb, Y[:, 12 * b:12 * b + 12] * np.exp(-1j * TH[(i, b)])) / np.vdot(Rb, Rb).real
        D /= np.sqrt(np.mean(np.abs(D) ** 2)) + 1e-30; Dprof += np.abs(D @ Ftau) ** 2; nD += 1
    Dprof /= max(nD, 1)
    # --- walker range profiles ---
    nW = int(2.5 // WIN) + 1; P = {lab: np.zeros((nW, len(DEL), 4), complex) for lab in ("trs", "dm_anch", "dm_raw", "trs_flip", "dm_anch_flip")}; N_ = {lab: np.zeros(nW) for lab in P}
    for b in B:
        Ls = [(tiq(E[i][1], E[i][2]), E[i][4][:, 12 * b:12 * b + 12] * np.exp(-1j * TH[(i, b)])) for (i, b_), kk in own.items() if b_ == b and kk == best[b]]
        Ls = [x for x in Ls if np.sum(np.abs(x[1] - S[best[b]]) ** 2) <= GATE * np.sum(np.abs(S[best[b]]) ** 2)]
        if len(Ls) < 20: continue
        Ls.sort(key=lambda x: x[0]); T = np.array([x[0] for x in Ls]); Dd = hp_series(T, np.array([x[1] for x in Ls])); kk_ = KKD[12 * b:12 * b + 12] - 1638
        G_ = np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), S[best[b]].ravel()) / np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), Rbar[:, 12 * b:12 * b + 12].ravel()).real
        M2 = np.exp(-2j * np.pi * kk_[None, :] * DF * DEL[:, None] / 3e8)
        for j, t in enumerate(T):
            dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
            if dr is None or not 0 <= w < nW: continue
            x = Dd[j] * np.exp(2j * np.pi * (kk_ * DF * dr / 3e8 + dr / LAM))[None, :]
            r = x @ M2.conj().T                                                  # [4, nDel]
            xf = Dd[j] * np.exp(2j * np.pi * (kk_ * DF * dr / 3e8 - dr / LAM))[None, :]; P["dm_anch_flip"][w] += (xf @ M2.conj().T).T / G_ * abs(G_)
            P["dm_raw"][w] += r.T; P["dm_anch"][w] += r.T / G_ * abs(G_); N_["dm_raw"][w] += 12; N_["dm_anch"][w] += 12
    TT = sorted([(tiq(a, l), R / osc(t)) for t, a, l, R in TRS], key=lambda x: x[0]); Tq = np.array([x[0] for x in TT]); Dt = hp_series(Tq, np.array([x[1] for x in TT]))
    M2t = np.exp(-2j * np.pi * (KKT - 1638)[None, :] * DF * DEL[:, None] / 3e8)
    for j, t in enumerate(Tq):
        dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
        if dr is None or not 0 <= w < nW: continue
        x = Dt[j] * np.exp(2j * np.pi * ((KKT - 1638) * DF * dr / 3e8 + dr / LAM))[None, :]; P["trs"][w] += (x @ M2t.conj().T).T; N_["trs"][w] += 1
        xf = Dt[j] * np.exp(2j * np.pi * ((KKT - 1638) * DF * dr / 3e8 - dr / LAM))[None, :]; P["trs_flip"][w] += (xf @ M2t.conj().T).T
    prof = {lab: np.sum(np.abs(P[lab]) ** 2, (0, 2)) for lab in P}
    for lab in prof: TOT[(lab, k)] = prof[lab]
    TOT[("D", k)] = Dprof
    print("chunk %d: %d transmissions for D_b; profiles done" % (k, nD), flush=True)
np.savez("/mnt/mydrive/field/ana_walk1/range_prof.npz", DEL=DEL, taus=taus, **{"%s_%d" % kk: v for kk, v in TOT.items()})
for lab in ("trs", "trs_flip", "dm_anch", "dm_anch_flip", "dm_raw"):
    p = sum(TOT[(lab, k)] for k in KS); p = p / sum(TOT[(lab.replace("_flip", ""), k)] for k in KS)[np.argmin(np.abs(DEL))]
    print("%-8s energy vs range offset (rel. to truth): %s" % (lab, " ".join("%+d:%.2f" % (d, p[np.argmin(np.abs(DEL - d))]) for d in (-40, -20, -10, -5, 0, 5, 10, 20, 40))))
Dp = sum(TOT[("D", k)] for k in KS); Dp /= Dp.max(); pk = [j for j in range(1, len(taus) - 1) if Dp[j] > Dp[j - 1] and Dp[j] >= Dp[j + 1] and Dp[j] > 0.1]
print("D_b delay profile peaks (> 0.1 of max): %s" % ", ".join("%.0f ns (%.1f m): %.2f" % (taus[j] * 1e9, taus[j] * 3e8, Dp[j]) for j in pk))
