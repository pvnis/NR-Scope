import sys, os, io, contextlib, numpy as np, collections as CC
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="cplx", BA_W="0", MINPRB="48", TAB_ONLY="1")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/jp.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(open(g["__file__"]).read(), "so_light", "exec"), g)
E, ALPHA, own = g["E"], g["ALPHA"], g["own"]
idx = {(E[i][1], E[i][2]): i for i in range(len(E))}
J = {}
for (a, l), i in idx.items():
    if l != 2 or (a + 1, 2) not in idx: continue
    j = idx[(a + 1, 2)]
    d = [np.angle(ALPHA[(j, b)] / ALPHA[(i, b)]) for b in np.flatnonzero(E[i][3] & E[j][3]) if (i, b) in ALPHA and (j, b) in ALPHA and own[(i, b)] == own[(j, b)]]
    if len(d) >= 8: J[a] = np.angle(np.mean(np.exp(1j * np.array(d))))
A = np.array(sorted(J)); v = np.array([J[a] for a in A])
big = A[np.abs(v - np.median(v)) > 0.6]
print("slot pairs %d, large jumps %d (%.0f %%)" % (len(A), len(big), 100 * len(big) / len(A)))
print("slot index mod 20 of large jumps:", dict(sorted(CC.Counter(big % 20).items())))
print("slot index mod 20 of all pairs   :", dict(sorted(CC.Counter(A % 20).items())))
print("gaps between large jumps (slots):", CC.Counter(np.diff(big)).most_common(8))
print("jump per slot-boundary position (mod 20): circular mean, spread (MAD), count")
for p in sorted(set(A % 20)):
    vv = v[A % 20 == p]; m = np.angle(np.mean(np.exp(1j * vv))); sp = 1.4826 * np.median(np.abs(np.angle(np.exp(1j * (vv - m)))))
    print("  %2d -> %2d: %+.2f rad, MAD %.2f rad, n %d" % (p, p + 1, m, sp, len(vv)))
