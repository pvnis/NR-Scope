"""How much of a target echo does a per-PRG-symbol complex gain on the beam template absorb? For real templates S
(4 antennas x 12 pilots, carry3 chunk), a synthetic echo t (excess range r, angle phi, ULA half-wavelength, pilot
spacing 4 subcarriers x 30 kHz): absorbed share |<S,t>|^2 / (|S|^2 |t|^2), median over PRGs and beams."""
import sys, os, io, contextlib, numpy as np
CH = sys.argv[1]
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="prg", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", CH, "/tmp/abs_tmp.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(g["__file__"]).read(), "so_light", "exec"), g)
if os.path.exists("/tmp/abs_tmp.bin"): os.remove("/tmp/abs_tmp.bin")
Ssum, Scnt, own = g["Ssum"], g["Scnt"], g["own"]
kprg = {k: b for (i, b), k in own.items()}
T = [(Ssum[k] / Scnt[k], kprg[k]) for k in Ssum if Scnt[k] >= 10]
df = 4 * 30e3; c = 3e8; n = np.arange(12)
print("templates:", len(T))
for r in (0, 5, 10, 20, 30):
    row = []
    for phi in (0, 15, 30, 45):
        a = np.exp(1j * np.pi * np.arange(4) * np.sin(np.radians(phi)))[:, None]
        f = []
        for S, b in T:
            t = a * np.exp(-2j * np.pi * (n + 12 * b) * df * (r / c))[None, :]
            f.append(abs(np.vdot(S, t)) ** 2 / (np.vdot(S, S).real * np.vdot(t, t).real))
        row.append("%5.1f%%" % (100 * np.median(f)))
    print("excess %2d m | angle 0/15/30/45 deg: %s" % (r, " ".join(row)))
print("(range 0 and angle of the direct path = the static itself: absorbed fully by construction; 1/48 = 2.1 %)")
