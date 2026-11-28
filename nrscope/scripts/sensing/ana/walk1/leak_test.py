"""Is the tablet DM-RS residual a per-transmission component (layer leakage with its own random phase) rather than noise?
Wideband-anchored (BA_PH=wide) values y = Y e^{-j theta}: compare, per PRG, symbol 2 vs symbol 11 of the SAME slot (same
transmission: noise only) and symbol 2 of consecutive slots of the same beam (different transmissions), relative to
the beam's static power. If across-slot >> within-slot, a per-transmission component dominates the residual."""
import sys, os, io, contextlib, numpy as np
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/lk.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
                                         "        KEEP = (E, dict(own), dict(THETA), {k: Ssum[k] / Scnt[k] for k in Ssum}); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, TH, S = g["KEEP"]
idx = {(E[i][1], E[i][2]): i for i in range(len(E))}
def y(i, b): return E[i][4][:, 12 * b:12 * b + 12] * np.exp(-1j * TH[(i, b)])
w, x, r = [], [], []
for (a, l), i in idx.items():
    if l != 2: continue
    for (a2, l2), out in (((a, 11), w), ((a + 1, 2), x)):
        j = idx.get((a2, l2))
        if j is None: continue
        for b in np.flatnonzero(E[i][3] & E[j][3]):
            if (i, b) in own and (j, b) in own and own[(i, b)] == own[(j, b)]:
                Sb = S[own[(i, b)]]; P = np.sum(np.abs(Sb) ** 2)
                out.append(np.sum(np.abs(y(j, b) - y(i, b)) ** 2) / P / 2)     # /2: per snapshot
                if out is w: r.append(np.sum(np.abs(y(i, b) - Sb) ** 2) / P)
f = lambda v: 10 * np.log10(np.median(v))
print("tablet, per snapshot, relative to beam static power (median):")
print("  residual about the beam static          %.1f dB" % f(r))
print("  difference symbol 2 vs 11, SAME slot    %.1f dB  (noise only)" % f(w))
print("  difference symbol 2, CONSECUTIVE slots  %.1f dB  (noise + whatever changes per transmission)" % f(x))
