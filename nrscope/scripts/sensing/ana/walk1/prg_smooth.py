"""Is the tablet's precoder wideband or per PRG? For each beam template set built in one run (consistent phase across
PRGs), G_b = <R_b, S_b> / |R_b|^2 per PRG (beam static over TRS static). Wideband precoding: G varies smoothly across
adjacent PRGs (both statics share the direct path); per-PRG precoders: G jumps. Reported: phase change of G between
adjacent PRGs of the same snapshot's beams, versus the same for the TRS static against itself shifted (baseline).
usage: prg_smooth.py <chunk>"""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/ps.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
    "        KEEP = (E, dict(own), dict(THETA), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, TH, S, Sc, Rbar = g["KEEP"]
def Gof(Sb, b):
    Rb = Rbar[:, 12 * b:12 * b + 12]; return np.vdot(Rb, Sb) / np.vdot(Rb, Rb).real
# per snapshot: the beams on adjacent PRGs (both with >= 20 members)
d_adj = []; d_far = []
for i, (t, a, l, pm, Y, u) in enumerate(E[::5]):
    i *= 5; bs = [b for b in np.flatnonzero(E[i][3]) if (i, b) in own and Sc[own[(i, b)]] >= 20]
    Gs = {b: Gof(S[own[(i, b)]], b) for b in bs}
    for b in bs:
        if b + 1 in Gs: d_adj.append(np.angle(Gs[b + 1] / Gs[b]))
        if b + 8 in Gs: d_far.append(np.angle(Gs[b + 8] / Gs[b]))
d_adj = np.array(d_adj); d_far = np.array(d_far)
mad = lambda v: 1.4826 * np.median(np.abs(np.angle(np.exp(1j * (v - np.angle(np.mean(np.exp(1j * v))))))))
# baseline: how much the TRS static's own phase pattern changes between adjacent PRGs (frequency selectivity of the channel)
r = [np.angle(np.vdot(Rbar[:, 12 * b:12 * b + 12], Rbar[:, 12 * (b + 1):12 * (b + 1) + 12])) for b in range(67)]
print("chunk %d: beam/TRS factor G, phase change between adjacent PRGs: MAD %.2f rad (|mean phasor| %.2f, n %d); 8 PRGs apart: MAD %.2f rad"
      % (k, mad(d_adj), abs(np.mean(np.exp(1j * d_adj))), len(d_adj), mad(d_far)))
print("          |G| spread across PRGs (dB, p10..p90): %s" % np.round(np.percentile([20 * np.log10(abs(Gof(S[own[(i, b)]], b))) for (i, b) in list(own)[::50] if Sc[own[(i, b)]] >= 20], [10, 50, 90]), 1).tolist())
# template-free: per snapshot, r_b = <R_b, Y_b>/|R_b|^2 on each PRG of the snapshot; adjacent-PRG phase change of r
dr = []; dr8 = []; per_snap_var = []
for i, (t, a, l, pm, Y, u) in enumerate(E[::5]):
    bs = np.flatnonzero(pm)
    if len(bs) < 10: continue
    r = {b: np.vdot(Rbar[:, 12 * b:12 * b + 12], Y[:, 12 * b:12 * b + 12]) / np.vdot(Rbar[:, 12 * b:12 * b + 12], Rbar[:, 12 * b:12 * b + 12]).real for b in bs}
    for b in bs:
        if b + 1 in r: dr.append(np.angle(r[b + 1] / r[b]))
        if b + 8 in r: dr8.append(np.angle(r[b + 8] / r[b]))
dr = np.array(dr); dr8 = np.array(dr8)
print("          template-free, one transmission: DM-RS/TRS ratio phase change between adjacent PRGs MAD %.2f rad (|mean phasor| %.2f); 8 apart MAD %.2f" % (mad(dr), abs(np.mean(np.exp(1j * dr))), mad(dr8)))
# and is that pattern the SAME from one transmission to the next (a fixed per-PRG precoder) or does it change?
pat = collections.defaultdict(list)
for i, (t, a, l, pm, Y, u) in enumerate(E[::3]):
    for b in np.flatnonzero(pm):
        if b + 1 < 68 and pm[b + 1]:
            r0 = np.vdot(Rbar[:, 12 * b:12 * b + 12], Y[:, 12 * b:12 * b + 12]); r1 = np.vdot(Rbar[:, 12 * (b + 1):12 * (b + 1) + 12], Y[:, 12 * (b + 1):12 * (b + 1) + 12])
            pat[b].append(np.angle(r1 / r0))
st = [mad(np.array(v)) for v in pat.values() if len(v) > 30]
print("          same adjacent pair over time: phase change spread (per PRG pair, median) %.2f rad over %d pairs" % (np.median(st), len(st)))
