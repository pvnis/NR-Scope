import sys, os, io, contextlib, numpy as np, collections as CC
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="cplx", BA_W="0", MINPRB="48", TAB_ONLY="1")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/js.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(open(__file__).read(), "so_light", "exec"))
idx = {(E[i][1], E[i][2]): i for i in range(len(E))}
within = []; common = []; nb_ = []
for (a, l), i in idx.items():
    if l != 2 or (a + 1, 2) not in idx: continue
    j = idx[(a + 1, 2)]
    d = np.array([np.angle(ALPHA[(j, b)] / ALPHA[(i, b)]) for b in np.flatnonzero(E[i][3] & E[j][3]) if (i, b) in ALPHA and (j, b) in ALPHA and own[(i, b)] == own[(j, b)]])
    if len(d) < 8: continue
    c = np.angle(np.mean(np.exp(1j * d))); common.append(c); within += list(np.angle(np.exp(1j * (d - c)))); nb_.append(len(d))
common = np.array(common); within = np.array(within)
print("slot pairs %d (median %d PRGs each): common jump MAD %.2f rad; per-PRG deviation from the common jump MAD %.2f rad" % (len(common), np.median(nb_), 1.4826 * np.median(np.abs(common - np.median(common))), 1.4826 * np.median(np.abs(within))))
print("common jump histogram (-pi..pi, 12 bins):", np.histogram(common, bins=np.linspace(-np.pi, np.pi, 13))[0].tolist())
