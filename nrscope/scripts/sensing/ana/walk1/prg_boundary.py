"""Per-PRG precoding test: within ONE tablet DM-RS transmission, the ratio DM-RS / TRS static per PRB (3 values per PRB,
4 antennas) and its phase change between adjacent PRBs, grouped by the PRB index mod 4 (and mod 2). A wideband
precoder (or multipath) changes smoothly: the same at every position. A precoder chosen per PRG of 2 or 4 PRBs jumps
at the PRG boundaries (positions 1->2 / 3->0). Baseline: the TRS against its own mean, same metric."""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/pb.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n", "        KEEP = (E, Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, Rbar = g["KEEP"]; TRS = g["TRS"]; osc = g["osc"]
def prb_ratio(Y, mask_prb):
    r = {}
    for p in np.flatnonzero(mask_prb):
        y = Y[:, 3 * p:3 * p + 3].ravel(); R = Rbar[:, 3 * p:3 * p + 3].ravel(); r[p] = np.vdot(R, y) / np.vdot(R, R).real
    return r
def jumps(rows):
    J = collections.defaultdict(list)
    for r in rows:
        for p in r:
            if p + 1 in r: J[p % 4].append(np.angle(r[p + 1] / r[p]))
    return {m: 1.4826 * np.median(np.abs(np.angle(np.exp(1j * (np.array(v) - np.angle(np.mean(np.exp(1j * np.array(v))))))))) for m, v in sorted(J.items())}, \
           {m: abs(np.mean(np.exp(1j * np.array(v)))) for m, v in sorted(J.items())}
rows = []
for (t, a, l, pm, Y, u) in E[::4]:
    rows.append(prb_ratio(Y, np.repeat(pm, 4)))
mj, coh = jumps(rows)
print("chunk %d, tablet DM-RS / TRS static: phase change PRB p -> p+1, by p mod 4 (MAD rad): %s" % (k, {m: round(v, 2) for m, v in mj.items()}))
print("            consistency over transmissions (|mean phasor|): %s" % {m: round(v, 2) for m, v in coh.items()})
tr = [prb_ratio(R / osc(t), np.ones(272, bool)) for (t, a, l, R) in TRS[::4]]
mj2, coh2 = jumps(tr)
print("            baseline TRS / TRS static, same metric (MAD rad): %s" % {m: round(v, 2) for m, v in mj2.items()})
