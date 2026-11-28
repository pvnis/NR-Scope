"""Injected-target sanity test on walk1 (so_light INJ: the whole raw signal copied at range R(t)=r0-2v(t-T0), level lv dB,
i.e. an echo that sees the same precoder as the static = a target inside the beam). Same matched filter as range_prof
(HP +-0.25 s, tablet dominant beams, gate, anchored by the beam's G phase): injected cell / floor (+30..+60 m beyond it) for TRS and DM-RS.
usage: inj_test.py <chunk> <r0> <v> <lv_db>"""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); r0, v, lv = (float(x) for x in sys.argv[2:5]); HP = 0.25; WIN = 0.2; GATE = 10 ** (-12 / 10); LAM = 3e8 / 3.75e9; DF = 30e3
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1", INJ="%g,%g,%g,0" % (r0, v, lv))
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/inj.bin"]; g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
    "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, TH, S, Sc, Rbar = g["KEEP"]; KKD, KKT, TRS, osc, T0 = g["KK_D"], g["KK_T"], g["TRS"], g["osc"], g["T0"][0]
Rt = lambda t: r0 - 2 * v * (t - T0)
def hp(T, V):
    C = np.concatenate([np.zeros((1,) + V.shape[1:], complex), np.cumsum(V, 0)]); lo = np.searchsorted(T, T - HP); hi = np.searchsorted(T, T + HP, "right")
    return V - (C[hi] - C[lo]) / (hi - lo)[:, None, None]
DL = np.r_[0.0, np.arange(30, 61, 3.0)]; t0 = min(x[0] for x in TRS); nW = int(2.6 // WIN) + 1
def prof(series, kk_, rel):
    acc = np.zeros((nW, len(DL), 4), complex); M = np.exp(-2j * np.pi * kk_[None, :] * DF * DL[:, None] / 3e8)
    for t, x in series:
        R = Rt(t); w = int((t - t0) // WIN)
        if not 0 <= w < nW: continue
        acc[w] += ((x * rel * np.exp(2j * np.pi * (kk_ * DF * R / 3e8 + R / LAM))[None, :]) @ M.conj().T).T
    return acc
cnt = collections.Counter(own.values()); best = {}
for (i, b), kk in own.items():
    if Sc[kk] >= 50 and (b not in best or cnt[kk] > cnt[best[b]]): best[b] = kk
TT = sorted([(x[0], x[3] / osc(x[0])) for x in TRS], key=lambda z: z[0]); Tq = np.array([x[0] for x in TT]); Dt = hp(Tq, np.array([x[1] for x in TT]))
Et = np.sum(np.abs(prof(list(zip(Tq, Dt)), KKT - 1638, 1.0)) ** 2, (0, 2)); print("TRS cell/floor %.1f" % (Et[0] / Et[1:].mean()))
def run(lab, use_th=True, gate=True, use_hp=True, anch=True, allbeams=False, mode="shared", synth=False):
    Ed = np.zeros((nW, len(DL), 4), complex); n = 0
    groups = collections.defaultdict(list)
    for (i, b), kk in own.items():
        if not allbeams and best.get(b) != kk: continue
        groups[(b, kk)].append(i)
    for (b, kk), ii in groups.items():
        Ls = [(E[i][0], E[i][4][:, 12 * b:12 * b + 12] * (np.exp(-1j * TH[(i, b)]) if use_th else 1)) for i in ii]
        if gate: Ls = [x for x in Ls if np.sum(np.abs(x[1] * (1 if use_th else 1) - S[kk]) ** 2) <= GATE * np.sum(np.abs(S[kk]) ** 2)]
        if len(Ls) < 20: continue
        Ls.sort(key=lambda x: x[0]); T = np.array([x[0] for x in Ls]); V = np.array([x[1] for x in Ls]); Dd = hp(T, V) if use_hp else V
        if synth:
            kq = KKD[12 * b:12 * b + 12] - 1638
            Dd = np.array([S[kk] * (10 ** (lv / 20) * np.exp(-2j * np.pi * Rt(t) / LAM) * np.exp(-2j * np.pi * kq * DF * Rt(t) / 3e8))[None, :] for t in T])
        G_ = np.vdot(Rbar[:, 12 * b:12 * b + 12].ravel(), S[kk].ravel()); n += len(T)
        Rb = Rbar[:, 12 * b:12 * b + 12]; Sb = S[kk]
        if mode == "shared": rel = abs(G_) / G_ if anch else 1.0
        elif mode == "perant":
            Gc = np.sum(Rb.conj() * Sb, 1); rel = (np.abs(Gc) / (Gc + 1e-30))[:, None]
        elif mode == "perant_ramp":   # per antenna: one number + a linear phase across the PRG
            q = np.arange(12) - 5.5; rel = np.ones((4, 12), complex)
            for c in range(4):
                r_ = Sb[c] / (Rb[c] + 1e-30); A_ = np.stack([np.ones(12), q], 1); co, *_ = np.linalg.lstsq(A_, np.unwrap(np.angle(r_)), rcond=None); rel[c] = np.exp(-1j * (A_ @ co))
        elif mode == "equalise":      # per antenna and value: TRS static / DM-RS static (phase only)
            r_ = Rb.conj() * Sb; rel = np.abs(r_) / (r_ + 1e-30)
        Ed += prof(list(zip(T, Dd)), KKD[12 * b:12 * b + 12] - 1638, rel)
    Ed = np.sum(np.abs(Ed) ** 2, (0, 2)); print("%-38s DM-RS cell/floor %.2f (%d snapshots-PRG)" % (lab, Ed[0] / Ed[1:].mean(), n), flush=True)
run("shared anchoring (as before)")
run("per-antenna anchoring", mode="perant")
run("per-antenna + in-PRG phase ramp", mode="perant_ramp")
run("per-antenna per-value equalisation", mode="equalise")
run("SYNTHETIC echo only, shared", synth=True)
run("SYNTHETIC echo only, equalise", synth=True, mode="equalise")
print("t range DM-RS %.3f..%.3f, TRS %.3f..%.3f, T0 %.3f, nW %d" % (min(e[0] for e in E), max(e[0] for e in E), Tq.min(), Tq.max(), T0, nW))
# raw check: one PRG, correlate data with the expected injected echo shape (static x g(t))
b = sorted(best)[len(best) // 2]; kk = best[b]; ii = sorted([i for (i, b_), k_ in own.items() if b_ == b and k_ == kk], key=lambda i: E[i][0])
kk_ = KKD[12 * b:12 * b + 12] - 1638; num = 0; den = 0; num0 = 0
for i in ii:
    t = E[i][0]; R = Rt(t); y = E[i][4][:, 12 * b:12 * b + 12] * np.exp(-1j * TH[(i, b)])
    gexp = 10 ** (lv / 20) * np.exp(-2j * np.pi * R / LAM) * np.exp(-2j * np.pi * kk_ * DF * R / 3e8)
    e = S[kk] * gexp[None, :]; d = y - S[kk]; num += np.vdot(e, d); den += np.vdot(e, e).real; num0 += abs(np.vdot(e, d)) ** 2 / max(np.vdot(e, e).real, 1e-30)
print("PRG %d, %d snapshots: LS amplitude of the expected echo in (y - static) = %.3f (1 = fully present), phase %.2f rad" % (b, len(ii), abs(num / den), np.angle(num / den)))
