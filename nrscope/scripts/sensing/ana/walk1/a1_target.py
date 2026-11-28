"""Does assumption A1 (the precoder scales all significant paths of a PRG by the same factor) hold for a TARGET path?
Per chunk, tablet DM-RS (BA wide, tablet only), dynamic part high-passed per beam and PRG, NO per-PRG division:
  d_b(t) = Y_b e^{-j theta} - S_local          (beam frame: carries the precoder factor of every path)
For the walker (lidar range and Doppler), per 4-PRB group b and antenna c, the echo factor over the whole chunk:
  gamma_bc = sum_t conj(m) d / sum |m|^2,  m = exp(-j2pi (k df tau(t) + d_ex(t)/lambda))
and, for the static, the per-PRG factor G_bc = <R_bc, S_bc>/|R_bc|^2 (beam static over TRS static, per antenna).
A1 for the target means gamma_bc / G_bc is the same complex number for every group (per antenna). Reported:
  - coherence of gamma/G across groups |mean(u)| with u the unit phasors (1 = A1 holds), per antenna, vs
  - the same for gamma alone (no static factor applied) and for a null (wrong range +25 m) as the chance level,
  - and the phase change of gamma/G at group boundaries vs inside the chunk's adjacent groups.
usage: a1_target.py <chunk list> <offset s>"""
import sys, os, io, json, contextlib, collections, numpy as np
KS = [int(x) for x in sys.argv[1].split(",")]; OFF = float(sys.argv[2]); HP = 0.25; GATE = 10 ** (-12 / 10)
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
GD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth/grid"
D = np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz"); TR = D["tr"]; t_rec = float(D["t_rec"]); IDS = sorted(set(TR[:, 0].astype(int)))
LAM = 3e8 / 3.75e9; DF = 30e3
for k in KS:
    os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
    sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/a1.bin"]
    g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
        "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
    E, own, TH, S, Sc, Rbar = g["KEEP"]; KKD = g["KK_D"]
    a0 = json.load(open("%s/b%d.json" % (GD, k)))["a0"]
    tiq = lambda a, l: k * 2.5 + ((a - a0) % 20480) * 5e-4 + l / 28e3
    # dominant tablet beam per group (most members): A1 is tested per (beam, group) with that beam's static
    cnt = collections.Counter(own.values()); best = {}
    for (i, b), kk in own.items():
        if Sc[kk] >= 50 and (b not in best or cnt[kk] > cnt[best[b]]): best[b] = kk
    ser = collections.defaultdict(list)
    for (i, b), kk in own.items():
        if best.get(b) != kk: continue
        t, a, l, pm, Y, u = E[i]; sl = slice(12 * b, 12 * b + 12); yv = Y[:, sl] * np.exp(-1j * TH[(i, b)])
        if np.sum(np.abs(yv - S[kk]) ** 2) > GATE * np.sum(np.abs(S[kk]) ** 2): continue
        ser[b].append((tiq(a, l), yv))
    WIN = 0.2; res = {}
    for name, roff in (("truth", 0.0), ("null+25m", 25.0)):
        gam = collections.defaultdict(dict)   # window -> group -> 4 complex
        for b, L in ser.items():
            L.sort(key=lambda x: x[0]); T = np.array([x[0] for x in L]); V = np.array([x[1] for x in L])
            C = np.concatenate([np.zeros((1, 4, 12), complex), np.cumsum(V, 0)]); lo = np.searchsorted(T, T - HP); hi = np.searchsorted(T, T + HP, "right")
            Dd = V - (C[hi] - C[lo]) / (hi - lo)[:, None, None]
            kk_ = KKD[12 * b:12 * b + 12] - 1638; num = collections.defaultdict(lambda: np.zeros(4, complex)); den = collections.Counter()
            for j, t in enumerate(T):
                te = t_rec + OFF + t
                for tid in IDS:
                    m_ = TR[:, 0] == tid; tk = TR[m_, 1]
                    if not (tk[0] <= te <= tk[-1]): continue
                    dr = np.interp(te, tk, TR[m_, 2]); m = np.exp(-2j * np.pi * (kk_ * DF * (dr + roff) / 3e8 + dr / LAM)); w = int(t // WIN)
                    num[w] += Dd[j] @ np.conj(m); den[w] += np.sum(np.abs(m) ** 2)
            for w in num:
                if den[w] > 0: gam[w][b] = num[w] / den[w]
        Gs = {b: np.array([np.vdot(Rbar[c, 12 * b:12 * b + 12], S[best[b]][c]) / np.vdot(Rbar[c, 12 * b:12 * b + 12], Rbar[c, 12 * b:12 * b + 12]).real for c in range(4)]) for b in ser}
        cr, ca, pw = [], [], []
        for w, gg in gam.items():
            bs = sorted(gg)
            if len(bs) < 20: continue
            Gam = np.array([gg[b] for b in bs]); Gst = np.array([Gs[b] for b in bs])
            R = Gam / Gst; cr.append(np.mean([abs(np.sum(R[:, c] * np.abs(Gst[:, c]) ** 2)) / np.sum(np.abs(R[:, c]) * np.abs(Gst[:, c]) ** 2) for c in range(4)]))
            ca.append(np.mean([abs(np.sum(Gam[:, c])) / np.sum(np.abs(Gam[:, c])) for c in range(4)])); pw.append(np.mean(np.abs(Gam) ** 2))
        res[name] = (np.mean(cr), np.mean(ca), np.mean(pw), len(cr))
    print("chunk %d: windows %d | across-group coherence of echo/static factor: truth %.2f, null %.2f | of echo factor alone: truth %.2f, null %.2f | per-window echo power truth/null %.2f"
          % (k, res["truth"][3], res["truth"][0], res["null+25m"][0], res["truth"][1], res["null+25m"][1], res["truth"][2] / res["null+25m"][2]), flush=True)
