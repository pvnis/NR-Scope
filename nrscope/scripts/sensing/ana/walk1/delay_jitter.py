"""Per-snapshot delay jitter: a random timing error per snapshot puts a random linear phase slope across the band; the
static survives an average, but every echo is smeared in range while keeping its Doppler. Per tablet DM-RS
transmission: r_b = <S_b, Y_b e^{-j theta}> / |S_b|^2 per 4-PRB group (S_b its beam static), linear fit of angle(r_b)
over the group index -> delay. Same for every TRS symbol against the TRS static. Delay = slope / (2 pi 1.44 MHz).
usage: delay_jitter.py <chunk>"""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/dj.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
    "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, TH, S, Sc, Rbar = g["KEEP"]; TRS = g["TRS"]; osc = g["osc"]
BW = 4 * 12 * 30e3
def slope(bs, r):
    ph = np.unwrap(np.angle(r)); p = np.polyfit(bs, ph, 1); return p[0]
dm = []
for i, (t, a, l, pm, Y, u) in enumerate(E):
    bs = [b for b in np.flatnonzero(pm) if (i, b) in own and Sc[own[(i, b)]] >= 20]
    if len(bs) < 20: continue
    r = np.array([np.vdot(S[own[(i, b)]], Y[:, 12 * b:12 * b + 12] * np.exp(-1j * TH[(i, b)])) / np.vdot(S[own[(i, b)]], S[own[(i, b)]]).real for b in bs])
    dm.append(slope(np.array(bs, float), r))
tr = []
for (t, a, l, R) in TRS:
    X = R / osc(t); r = np.array([np.vdot(Rbar[:, 12 * b:12 * b + 12], X[:, 12 * b:12 * b + 12]) / np.vdot(Rbar[:, 12 * b:12 * b + 12], Rbar[:, 12 * b:12 * b + 12]).real for b in range(68)])
    tr.append(slope(np.arange(68.0), r))
for lab, v in (("tablet DM-RS vs beam static", np.array(dm)), ("TRS vs TRS static", np.array(tr))):
    ns = 1e9 * v / (2 * np.pi * BW)
    print("chunk %d, %-28s: %5d snapshots, delay per snapshot: std %.2f ns, MAD %.2f ns, p1..p99 %.1f..%.1f ns (2.44 m range bin = %.1f ns)"
          % (k, lab, len(ns), ns.std(), 1.4826 * np.median(np.abs(ns - np.median(ns))), *np.percentile(ns, [1, 99]), 2.44 / 3e8 * 1e9))
