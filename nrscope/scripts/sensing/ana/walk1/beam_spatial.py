"""Does the tablet's beam suppress the directions other than the cart's? Per PRG: 4x4 spatial covariance of the static
channel over its 12 values, share of power in the dominant eigenvector (the direct path's direction) -- for the TRS
static (wide beam) and for the tablet's beam statics (BA wide). Off-dominant power = what other directions (reflections,
a walker) get relative to the direct path: a beam that suppresses them shows it lower than the TRS."""
import sys, os, io, contextlib, numpy as np
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
CH = sys.argv[1] if len(sys.argv) > 1 else "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3"
sys.argv = ["so_light.py", CH, "/tmp/e2e/bs.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
                                         "        KEEP = (E, dict(own), {k: Ssum[k] / Scnt[k] for k in Ssum}, Scnt.copy(), Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, S, Scnt, Rbar = g["KEEP"]
def offdom(X):   # X: 4 x 12 -> power share outside the dominant spatial direction, dB
    C = X @ X.conj().T; w = np.linalg.eigvalsh(C); return 10 * np.log10(max(1 - w[-1] / w.sum(), 1e-6))
tr = [offdom(Rbar[:, 12 * b:12 * b + 12]) for b in range(68)]
bm = {}
for (i, b), k in own.items():
    if Scnt[k] >= 20 and k not in bm: bm[k] = offdom(S[k])
bv = np.array(list(bm.values()))
print("power outside the dominant spatial direction (per PRG, 4 antennas): TRS static median %.1f dB (p25 %.1f, p75 %.1f); tablet beams median %.1f dB (p25 %.1f, p75 %.1f), %d beams"
      % (np.median(tr), np.percentile(tr, 25), np.percentile(tr, 75), np.median(bv), np.percentile(bv, 25), np.percentile(bv, 75), len(bv)))
