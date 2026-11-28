"""How strongly does each stream illuminate the walker? Per chunk (tablet decoded), so_light with beam anchoring (BA wide,
tablet only) is run in-process; then:
  TRS:    d = R e^{-j osc} - Rbar                 (TRS static removed)
  DM-RS:  d = (Y e^{-j theta} - S_beam) / G_beam  (beam static removed, mapped onto the TRS static: TRS units)
both high-passed (moving mean over +-HP s within the stream: TRS; each beam and PRG), then per 0.2 s window and antenna
the walker echo amplitude is fitted: alpha = sum conj(m) d / sum |m|^2, m = exp(-j2pi (k df tau + d_ex / lambda)) with
the lidar excess path d_ex(t). Noise: the same at wrong ranges (+12..+40 m) and with the Doppler history reversed.
Echo power = mean |alpha|^2 (truth) - mean over the nulls. DM-RS / TRS ratio = the tablet beam's illumination of the
walker relative to the TRS's (both normalised to their own direct-path level).
usage: beamgain.py <chunk list> <offset s>"""
import sys, os, io, json, contextlib, collections, numpy as np
KS = [int(x) for x in sys.argv[1].split(",")]; OFF = float(sys.argv[2]); HP = 0.25; WIN = 0.2; GATE = 10 ** (float(os.environ.get("BG_GATE", "100")) / 10); ngated = 0
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
GD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth/grid"
D = np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz"); TR = D["tr"]; t_rec = float(D["t_rec"]); IDS = sorted(set(TR[:, 0].astype(int)))
LAM = 3e8 / 3.75e9; DF = 30e3
HYP = [("truth", 0.0, 1)] + [("r+%d" % o, float(o), 1) for o in (12, 18, 25, 32, 40)] + [("rev", 0.0, -1), ("rev+18", 18.0, -1), ("rev+32", 32.0, -1)]
acc = {s: collections.defaultdict(float) for s in ("trs", "dmrs")}; nwin = collections.Counter(); PREV = {"trs": {}, "dmrs": {}}
PRED = {"trs": [], "dmrs": []}; NV = {"trs": [], "dmrs": []}; PW = {"trs": [], "dmrs": []}
for k in KS:
    os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
    sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/bg.bin"]
    g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
        "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), dict(G), Rbar); E = []\n")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
    E, own, TH, S, Sc, Gb, Rbar = g["KEEP"]; osc = g["osc"]; TRS = g["TRS"]; KKD = g["KK_D"]; KKT = g["KK_T"]
    a0 = json.load(open("%s/b%d.json" % (GD, k)))["a0"]
    tiq = lambda a, l: k * 2.5 + ((a - a0) % 20480) * 5e-4 + l / 28e3
    # streams: list of (t_iq, value index array, d 4 x n)
    trs = [(tiq(a, l), np.arange(816), R / osc(t) - Rbar, KKT) for t, a, l, R in TRS]
    dm = collections.defaultdict(list)   # (beam, prg) -> [(t_iq, d 4x12)]
    for (i, b), kk in own.items():
        if kk not in Gb or Sc[kk] < 10: continue
        t, a, l, pm, Y, u = E[i]; sl = slice(12 * b, 12 * b + 12)
        yv = Y[:, sl] * np.exp(-1j * TH[(i, b)])
        if np.sum(np.abs(yv - S[kk]) ** 2) > GATE * np.sum(np.abs(S[kk]) ** 2): ngated += 1; continue   # heavy tail, as BA_GATE
        dm[(kk, b)].append((tiq(a, l), (yv - S[kk]) / Gb[kk], b))
    def hp(series):   # series [(t, d, ...)] in time -> d minus the moving mean over +-HP s (cumulative sums)
        T = np.array([x[0] for x in series]); V = np.array([x[1] for x in series]); C = np.concatenate([np.zeros((1,) + V.shape[1:], V.dtype), np.cumsum(V, 0)])
        lo = np.searchsorted(T, T - HP, "left"); hi = np.searchsorted(T, T + HP, "right")
        M = (C[hi] - C[lo]) / (hi - lo)[(slice(None),) + (None,) * (V.ndim - 1)]
        return [(series[j][0], V[j] - M[j], *series[j][2:]) for j in range(len(series))]
    trs.sort(key=lambda x: x[0]); trs_h = hp([(t, d, idx, kk_) for t, idx, d, kk_ in trs])
    dm_h = []
    for key, ser in dm.items():
        ser.sort(key=lambda x: x[0])
        if len(ser) >= 5: dm_h += hp(ser)
    def truth_at(te, off, sg):
        r = []
        for tid in IDS:
            m = TR[:, 0] == tid; tk = TR[m, 1]
            if tk[0] <= te <= tk[-1]:
                tt = te if sg > 0 else 2 * tk[0] - te
                if tk[0] <= tt <= tk[-1]: r.append((tid, np.interp(te, tk, TR[m, 2]) + off, np.interp(tt, tk, TR[m, 2]) + off))
        return r
    for st, rows in (("trs", trs_h), ("dmrs", dm_h)):
        W = collections.defaultdict(lambda: [np.zeros(4, complex), 0.0])
        WN = collections.defaultdict(lambda: [0.0, 0])   # window -> [sum |d|^2 over antennas and values, values]
        for r in rows:
            t, d = r[0], r[1]
            wn = WN[int(t // WIN)]; wn[0] += float(np.sum(np.abs(d) ** 2)) / 4; wn[1] += d.shape[1]
            kk_ = (r[3] if st == "trs" else KKD[3 * 4 * r[2]:3 * 4 * r[2] + 12]) - 1638
            te = t_rec + OFF + t
            for name, off, sg in HYP:
                for tid, dr, dp in truth_at(te, off, sg):
                    m = np.exp(-2j * np.pi * (kk_ * DF * dr / 3e8 + dp / LAM))
                    w = W[(name, tid, int(t // WIN))]; w[0] += d @ np.conj(m); w[1] += np.sum(np.abs(m) ** 2)
        for wi_, (sp_, n_) in WN.items():
            if n_ > 0: PRED[st].append((sp_ / n_) / n_); NV[st].append(n_); PW[st].append(sp_ / n_)
        for (name, tid, wi), (num, den) in W.items():
            acc[st][name] += float(np.mean(np.abs(num / den) ** 2));
            if name == "truth": nwin[st] += 1
    snap = {st: dict(acc[st]) for st in acc}
    if PREV is not None:
        row = []
        for st in ("trs", "dmrs"):
            dtr = snap[st].get("truth", 0) - PREV[st].get("truth", 0); dn = np.mean([snap[st].get(n, 0) - PREV[st].get(n, 0) for n, _, _ in HYP if n != "truth"])
            row.append("%s truth %.3g nulls %.3g (x%.2f)" % (st, dtr, dn, dtr / dn if dn else float("nan")))
        print("chunk %d: TRS %d snaps, DM-RS %d PRG-snaps | %s" % (k, len(trs_h), len(dm_h), " | ".join(row)), flush=True)
    PREV = snap
res = {}
for st in ("trs", "dmrs"):
    nul = np.array([acc[st][n] for n, _, _ in HYP if n != "truth"]); tr_ = acc[st]["truth"]
    res[st] = (tr_ - nul.mean()) / max(nwin[st], 1)
    print("%-4s windows %d: truth %.3g, nulls mean %.3g (spread %.3g) -> echo power per window %.3g  (truth/nulls %.2f)"
          % (st, nwin[st], tr_, nul.mean(), nul.std(), res[st], tr_ / nul.mean()))
if res["trs"] > 0 and res["dmrs"] > 0:
    print("tablet-beam illumination of the walker relative to the TRS: %.1f dB" % (10 * np.log10(res["dmrs"] / res["trs"])))
else:
    print("one stream has no echo above its nulls: TRS %.3g, DM-RS %.3g" % (res["trs"], res["dmrs"]))

for st in ("trs", "dmrs"):
    nul = {n: acc[st][n] / max(nwin[st], 1) for n, _, _ in HYP}
    print("%-4s per window: values median %d, per-value residual power median %.3g; white-noise prediction %.3g vs measured nulls %.3g (x%.1f)"
          % (st, np.median(NV[st]), np.median(PW[st]), np.median(PRED[st]), np.mean([v for k, v in nul.items() if k != "truth"]), np.mean([v for k, v in nul.items() if k != "truth"]) / np.median(PRED[st])))
    print("     nulls per hypothesis (per window):", {k: round(v, 4) for k, v in nul.items()})
