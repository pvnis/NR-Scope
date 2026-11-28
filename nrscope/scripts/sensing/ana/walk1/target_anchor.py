"""Target anchoring, prototype on walk1 (lidar stands in for the TRS track). Tablet DM-RS, BA wide (one phase per
transmission), dominant beam per 4-PRB group, beam static removed (moving mean +-0.25 s), NO per-group division.
Per 0.2 s window w, group b, antenna c and hypothesis h (truth range, or +12..+40 m):
   rhs_h[w,b,c] = sum conj(m_h) d,   nrm[w,b] = sum |m|^2     (m: lidar range steering + Doppler phase)
Model of the walker echo: gamma[w,b,c] = kappa_w * c_b * a_c  (window factor absorbs the lidar phase drift; c_b the
target's per-group precoder factor; a_c the receive steering). c_b, a_c learned on the OTHER windows of the chunk
(alternating least squares, truth hypothesis); on the held-out window the energy of the projection onto c_b a_c m is
compared between the truth and the null ranges -- with the SAME c_b. Per-group corrections compared:
   target: learned c_b | static: c_b = the static's per-group factor (today's anchoring) | flat: c_b = 1.
TRS (no precoding, full band) alongside with its own per-window fit.  usage: target_anchor.py <chunks> <offset>"""
import sys, os, io, json, contextlib, collections, numpy as np
KS = [int(x) for x in sys.argv[1].split(",")]; OFF = float(sys.argv[2]); HP = 0.25; WIN = 0.2; GATE = 10 ** (-12 / 10)
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
GD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth/grid"
D = np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz"); TR = D["tr"]; t_rec = float(D["t_rec"]); IDS = sorted(set(TR[:, 0].astype(int)))
LAM = 3e8 / 3.75e9; DF = 30e3; HYP = [0.0, 12.0, 18.0, 25.0, 32.0, 40.0]
TOT = collections.defaultdict(float)
def walker_at(te):
    for tid in IDS:   # the labelled walker tracks; the first in view
        m = TR[:, 0] == tid; tk = TR[m, 1]
        if tk[0] <= te <= tk[-1] and tid in (1, 4): return np.interp(te, tk, TR[m, 2])
    return None
for k in KS:
    os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
    sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/ta.bin"]
    g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
        "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
    E, own, TH, S, Sc, Rbar = g["KEEP"]; KKD = g["KK_D"]; KKT = g["KK_T"]; TRS = g["TRS"]; osc = g["osc"]
    a0 = json.load(open("%s/b%d.json" % (GD, k)))["a0"]; tiq = lambda a, l: k * 2.5 + ((a - a0) % 20480) * 5e-4 + l / 28e3
    cnt = collections.Counter(own.values()); best = {}
    for (i, b), kk in own.items():
        if Sc[kk] >= 50 and (b not in best or cnt[kk] > cnt[best[b]]): best[b] = kk
    ser = collections.defaultdict(list)
    for (i, b), kk in own.items():
        if best.get(b) != kk: continue
        t, a, l, pm, Y, u = E[i]; sl = slice(12 * b, 12 * b + 12); yv = Y[:, sl] * np.exp(-1j * TH[(i, b)])
        if np.sum(np.abs(yv - S[kk]) ** 2) > GATE * np.sum(np.abs(S[kk]) ** 2): continue
        ser[b].append((tiq(a, l), yv / (np.linalg.norm(S[kk]) / np.linalg.norm(Rbar[:, sl]))))
    B = sorted(ser); bi = {b: j for j, b in enumerate(B)}
    Gst = np.array([np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), S[best[b]].ravel()) / np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), Rbar[:, 12 * b:12 * b + 12].ravel()).real for b in B])
    nW = int(2.5 // WIN) + 1; H = len(HYP)
    rhs = np.zeros((H, nW, len(B), 4), complex); nrm = np.zeros((H, nW, len(B)))
    for b, L in ser.items():
        L.sort(key=lambda x: x[0]); T = np.array([x[0] for x in L]); V = np.array([x[1] for x in L])
        C = np.concatenate([np.zeros((1, 4, 12), complex), np.cumsum(V, 0)]); lo = np.searchsorted(T, T - HP); hi = np.searchsorted(T, T + HP, "right")
        Dd = V - (C[hi] - C[lo]) / (hi - lo)[:, None, None]; kk_ = KKD[12 * b:12 * b + 12] - 1638
        for j, t in enumerate(T):
            dr = walker_at(t_rec + OFF + t)
            if dr is None: continue
            w = int((t - k * 2.5) // WIN)
            if not 0 <= w < nW: continue
            for h, off in enumerate(HYP):
                m = np.exp(-2j * np.pi * (kk_ * DF * (dr + off) / 3e8 + dr / LAM))
                rhs[h, w, bi[b]] += Dd[j] @ np.conj(m); nrm[h, w, bi[b]] += np.sum(np.abs(m) ** 2)
    # TRS: per window, full band, one factor per antenna
    trs_rhs = np.zeros((H, nW, 4), complex); trs_n = np.zeros((H, nW))
    TT = sorted([(tiq(a, l), R / osc(t)) for t, a, l, R in TRS], key=lambda x: x[0]); Tq = np.array([x[0] for x in TT]); V = np.array([x[1] for x in TT])
    C = np.concatenate([np.zeros((1,) + V.shape[1:], complex), np.cumsum(V, 0)]); lo = np.searchsorted(Tq, Tq - HP); hi = np.searchsorted(Tq, Tq + HP, "right")
    Dt = V - (C[hi] - C[lo]) / (hi - lo)[:, None, None]
    for j, t in enumerate(Tq):
        dr = walker_at(t_rec + OFF + t)
        if dr is None: continue
        w = int((t - k * 2.5) // WIN)
        if not 0 <= w < nW: continue
        for h, off in enumerate(HYP):
            m = np.exp(-2j * np.pi * ((KKT - 1638) * DF * (dr + off) / 3e8 + dr / LAM)); trs_rhs[h, w] += Dt[j] @ np.conj(m); trs_n[h, w] += np.sum(np.abs(m) ** 2)
    valid = [w for w in range(nW) if np.sum(nrm[0, w] > 0) >= 20]
    def learn(train):
        gam = rhs[0, train] / np.maximum(nrm[0, train], 1e-9)[:, :, None]; wt = nrm[0, train]      # [W, B, 4]
        c = np.ones(len(B), complex); a = np.ones(4, complex); kap = np.ones(len(train), complex)
        for _ in range(30):
            M = c[None, :, None] * a[None, None, :]
            kap = np.sum(wt[:, :, None] * np.conj(M) * gam, (1, 2)) / np.maximum(np.sum(wt[:, :, None] * np.abs(M) ** 2, (1, 2)), 1e-30)
            M = kap[:, None, None] * a[None, None, :]
            c = np.sum(wt[:, :, None] * np.conj(M) * gam, (0, 2)) / np.maximum(np.sum(wt[:, :, None] * np.abs(M) ** 2, (0, 2)), 1e-30)
            c /= np.sqrt(np.mean(np.abs(c) ** 2)) + 1e-30
            M = kap[:, None, None] * c[None, :, None]
            a = np.sum(wt[:, :, None] * np.conj(M) * gam, (0, 1)) / np.maximum(np.sum(wt[:, :, None] * np.abs(M) ** 2, (0, 1)), 1e-30)
        return c, a
    for wt_ in valid:
        train = [w for w in valid if w != wt_]
        if len(train) < 4: continue
        c, a = learn(train)
        for lab, cb in (("target", c), ("static", Gst / np.sqrt(np.mean(np.abs(Gst) ** 2))), ("flat", np.ones(len(B), complex))):
            Mbc = cb[:, None] * a[None, :]
            for h in range(H):
                e = abs(np.sum(np.conj(Mbc) * rhs[h, wt_])) ** 2 / max(np.sum(np.abs(Mbc) ** 2 * nrm[h, wt_][:, None]), 1e-30)
                TOT[(lab, h)] += e
        for h in range(H):
            TOT[("trs", h)] += np.sum(np.abs(trs_rhs[h, wt_]) ** 2) / max(trs_n[h, wt_], 1e-30)
    print("chunk %d: %d groups, %d test windows" % (k, len(B), len(valid)), flush=True)
for lab in ("trs", "target", "static", "flat"):
    tr_ = TOT[(lab, 0)]; nl = np.mean([TOT[(lab, h)] for h in range(1, len(HYP))])
    print("%-7s truth / range nulls = %.2f   (per null: %s)" % (lab, tr_ / nl, " ".join("+%dm:%.2f" % (HYP[h], TOT[(lab, h)] / tr_) for h in range(1, len(HYP)))))
