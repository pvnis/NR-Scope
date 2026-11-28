"""wtest: is the per-transmission precoder change real, SRS-periodic, or tied to MU-MIMO co-scheduling? (built on conditioning of the port fits + structured precoder models. Coherent TRS + DM-RS combination for a target through the CSI-RS port space (walk1 prototype; lidar = the track).
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
    import gs_all
    SVDb = {b: np.linalg.svd(Hc[:, 4 * b:4 * b + 4, :].reshape(16, 8), full_matrices=False) for b in B}
    def vh(b, y, ro):
        U, s_, Vh = SVDb[b]; return Vh[:ro].conj().T @ ((U[:, :ro].conj().T @ y.reshape(4, 4, 3).mean(2).reshape(16)) / s_[:ro])
    shp = lambda v, w: 1 - abs(np.vdot(v, w)) ** 2 / (np.vdot(v, v).real * np.vdot(w, w).real + 1e-30)   # shape change, scalar-free
    # raw records: port 1001 power on the tablet's PRBs, and other UEs' DM-RS in the slot
    tab = g["tab"]; P0 = {}; P1 = {}; OTH = collections.Counter()
    for c in range(1):
        os.environ.pop("P1001", None); gs_all.PRBM.clear(); d0 = gs_all.load(files[c])
        os.environ["P1001"] = "1"; gs_all.PRBM.clear(); d1 = gs_all.load(files[c]); os.environ.pop("P1001", None)
        for key, (m, v) in d0.items():
            a_, rn, l_ = key
            if rn == tab and key in d1:
                P0[(a_, l_)] = P0.get((a_, l_), 0) + float(np.sum(np.abs(v) ** 2)); P1[(a_, l_)] = P1.get((a_, l_), 0) + float(np.sum(np.abs(d1[key][1]) ** 2))
            if rn == 0 and c == 0: OTH[a_] += int(m.sum())
    rat = {key: P1[key] / max(P0[key], 1e-30) for key in P0}
    if False: print("chunk %d: port-1001 / port-1000 power on the tablet PRBs: median %.1f dB, p90 %.1f dB (%d symbols)" % (k, 10 * np.log10(np.median(list(rat.values()))), 10 * np.log10(np.percentile(list(rat.values()), 90)), len(rat)), flush=True)
    for ro in (2, 4):
        rec = collections.defaultdict(list)   # (b, kk) -> [(t, a, l, v)]
        for (i, b), kk in own.items():
            if b not in SVDb: continue
            t, a, l, pm, Y, u = E[i]; rec[(b, kk)].append((t, a, l, vh(b, Y[:, 12 * b:12 * b + 12] * np.exp(-1j * TH[(i, b)]), ro)))
        same, dev, lagd, jumps = [], [], collections.defaultdict(list), []; devkey = []
        for (b, kk), L_ in rec.items():
            if len(L_) < 30: continue
            L_.sort(key=lambda x: (x[1], x[2])); mu = np.mean([x[3] for x in L_], 0)
            for x in L_: dev.append(shp(x[3], mu)); devkey.append((x[1], x[2], kk))
            for x, y in zip(L_[:-1], L_[1:]):
                if x[1] == y[1]: same.append(shp(x[3], y[3]))
            for p in range(len(L_)):
                for q in range(p + 1, min(p + 60, len(L_))):
                    dt = (L_[q][1] - L_[p][1]) * 0.5
                    if dt > 0: lagd[int(np.clip(np.log2(dt), 0, 10))].append(shp(L_[p][3], L_[q][3]))
            # change times: consecutive transmissions (different slots) whose shape change exceeds 3x the same-slot noise
            sl_ = [(y[1], shp(x[3], y[3])) for x, y in zip(L_[:-1], L_[1:]) if y[1] != x[1]]
            jumps += [(a_, d_) for a_, d_ in sl_]
        ns = np.median(same) if same else np.nan
        print("  rank %d: shape change between the 2 DM-RS symbols of one slot (noise) median %.3f | vs the beam mean %.3f | %d groups" % (ro, ns, np.median(dev), len(rec)), flush=True)
        print("     vs time lag (ms): " + "  ".join("%d-%d: %.3f" % (2 ** e if e else 0, 2 ** (e + 1), np.median(v)) for e, v in sorted(lagd.items()) if len(v) > 50), flush=True)
        if ro == 4 and jumps:
            th = 3 * ns; ja = np.array([a_ for a_, d_ in jumps if d_ > th]); allA = np.array([a_ for a_, d_ in jumps])
            out = []
            for P in (10, 20, 40, 80, 160, 320):   # slots of 0.5 ms: P ms = 2P slots; Rayleigh concentration of change slots mod period vs all pairs
                ph = lambda A: abs(np.mean(np.exp(2j * np.pi * (A % (2 * P)) / (2 * P)))) if len(A) else np.nan
                out.append("%d ms: %.2f (all %.2f)" % (P, ph(ja), ph(allA)))
            print("     changes > 3x noise: %d of %d consecutive pairs; phase-lock of change slots to period: %s" % (len(ja), len(allA), "  ".join(out)), flush=True)
            dv = np.array(dev); r_ = np.array([rat.get((a_, l_), np.nan) for a_, l_, kk in devkey]); o_ = np.array([OTH.get(a_, 0) for a_, l_, kk in devkey])
            ok = np.isfinite(r_)
            if False:
                lo, hi = np.percentile(r_[ok], [33, 67])
                print("     MU-MIMO proxy: shape deviation when port 1001 is quiet (<= p33) %.3f vs busy (>= p67) %.3f; Spearman %.2f"
                      % (np.median(dv[ok & (r_ <= lo)]), np.median(dv[ok & (r_ >= hi)]), np.corrcoef(np.argsort(np.argsort(dv[ok])), np.argsort(np.argsort(r_[ok])))[0, 1]), flush=True)
            print("     other UEs in the slot: deviation with none %.3f vs with others %.3f (share of tablet symbols with others %.2f)"
                  % (np.median(dv[o_ == 0]) if np.any(o_ == 0) else np.nan, np.median(dv[o_ > 0]) if np.any(o_ > 0) else np.nan, np.mean(o_ > 0)), flush=True)
