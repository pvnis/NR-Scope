import sys, os, numpy as np, collections as CC
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="prg", BA_W="0", MINPRB="48")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/phd.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): exec(compile(open(__file__).read(), "so_light", "exec"))
res_lin = []; slopes = []; tot = []; ws = []
for i, (t, a, l, pm, Y, u) in enumerate(E):
    bs = [b for b in np.flatnonzero(pm) if Scnt[own[(i, b)]] >= nmin]
    if len(bs) < 6: continue
    ph = np.unwrap(np.array([PHI[(i, b)] for b in bs])); bb = np.array(bs, float)
    p = np.polyfit(bb, ph, 1); r = ph - np.polyval(p, bb)
    res_lin.append(np.std(r)); tot.append(np.std(ph)); slopes.append(p[0]); ws.append(shifts(a)[0] if hasattr(shifts(a), "__len__") else shifts(a))
slopes = np.array(slopes)
print("snapshots with >= 6 PRGs: %d" % len(slopes))
print("per-PRG phase spread inside a snapshot: %.2f rad raw, %.2f rad after a linear (delay) fit" % (np.median(tot), np.median(res_lin)))
print("slope per PRG: median %.3f rad, MAD %.3f rad -> delay %.1f ns (MAD %.1f ns)" % (np.median(slopes), 1.4826 * np.median(np.abs(slopes - np.median(slopes))), np.median(slopes) / (2 * np.pi * 48 * 30e3) * 1e9, 1.4826 * np.median(np.abs(slopes - np.median(slopes))) / (2 * np.pi * 48 * 30e3) * 1e9))
# phases of the TRS-model residual: does the slope cluster (a few window positions) ?
h = np.histogram(slopes, bins=np.linspace(-1, 1, 21))[0]; print("slope histogram (-1..1 rad/PRG, 0.1 steps):", h.tolist())
