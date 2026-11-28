import sys, os, numpy as np, collections as CC
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", MINPRB="48")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/stab.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
exec(compile(open(__file__).read().split("# DM-RS\nif USE:")[0], "head", "exec"))
Dd = [load(f) for f in files]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
tab = CC.Counter(k[1] for k in keys if k[1] != 0).most_common(1)[0][0]
X = {}
for k in sorted(keys):
    if k[1] != tab or k[2] != 2: continue
    pm = np.logical_and.reduce([d[k][0] for d in Dd]); t = k[0] * 0.5e-3 + k[2] / 28e3
    Y = np.array([d[k][1] for d in Dd]).astype(complex); X[k[0]] = (pm, dedrift(t, undo_wsh(Y, shifts(k[0]), -1.0)) * ramp)
A = sorted(X); print(len(A), "tablet slots (symbol 2)")
def fr(X1, X0, b):
    sl = slice(12 * b, 12 * b + 12); x = X1[:, sl].ravel(); r = X0[:, sl].ravel(); g = np.vdot(r, x) / np.vdot(r, r)
    return np.sum(np.abs(x - g * r) ** 2) / np.sum(np.abs(g * r) ** 2)
for lag in (1, 2, 5, 10, 20, 40, 100, 200, 1000, 2000):
    r = []
    for a in A[::3]:
        if a + lag in X:
            m = X[a][0] & X[a + lag][0]
            r += [fr(X[a + lag][1], X[a][1], b) for b in np.flatnonzero(m)[::4]]
    if r: print("lag %5d slots (%6.1f ms): residual median %6.1f dB, p10 %6.1f dB (%d)" % (lag, lag * 0.5, 10 * np.log10(np.median(r)), 10 * np.log10(np.percentile(r, 10)), len(r)))
# TRS for comparison: same thing at 20 ms / 100 ms / 1 s
