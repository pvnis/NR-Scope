import sys, os, numpy as np
os.environ.update(MODE="tab", USE_DMRS="0", USE_CSI="0")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/x.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
exec(compile(open(__file__).read().split("# DM-RS\nif USE:")[0], "head", "exec"))
t = np.array([x[0] for x in TRS]); R0 = TRS[0][3]
ph = np.unwrap(np.array([np.angle(np.vdot(R0, x[3])) for x in TRS]))
d = np.diff(ph); dt = np.diff(t)
w = dt < 0.002; a = dt > 0.015
print("TRS symbols %d over %.2f s" % (len(t), t[-1] - t[0]))
print("within an occasion: phase rate median %.1f rad/s (%.1f Hz), spread (MAD) %.1f rad/s" % (np.median(d[w] / dt[w]), np.median(d[w] / dt[w]) / 2 / np.pi, 1.4826 * np.median(np.abs(d[w] / dt[w] - np.median(d[w] / dt[w])))))
print("across occasions (20 ms): phase step median %.2f rad, MAD %.2f rad" % (np.median(d[a]), 1.4826 * np.median(np.abs(d[a] - np.median(d[a])))))
p = np.polyfit(t, ph, 1); r = ph - np.polyval(p, t)
print("linear fit over the chunk: %.2f Hz; residual std %.3f rad, max %.2f rad" % (p[0] / 2 / np.pi, r.std(), np.abs(r).max()))
p2 = np.polyfit(t, ph, 3); r2 = ph - np.polyval(p2, t); print("cubic fit residual std %.3f rad" % r2.std())
