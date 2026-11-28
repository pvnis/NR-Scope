import sys, os, io, contextlib, numpy as np
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/tl.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(__file__).read().replace("        E = []   # the per-PRG path below then writes nothing\n", "        E_keep, own_k, G_k, Ss_k, Sc_k, TH_k = E, dict(own), dict(G), dict(Ssum), Scnt.copy(), dict(THETA); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"))
E, own, G, Ssum, Scnt, THETA = E_keep, own_k, G_k, Ss_k, Sc_k, TH_k; r = []
for (i, b), kk_ in own.items():
    if kk_ not in G: continue
    Sb = Ssum[kk_] / Scnt[kk_]; y = E[i][4][:, 12 * b:12 * b + 12] * np.exp(-1j * THETA[(i, b)])
    r.append(np.sum(np.abs(y - Sb) ** 2) / np.sum(np.abs(Sb) ** 2))
r = np.array(r); q = np.percentile(r, [50, 90, 99, 99.9])
print("tablet PRG-snapshot residual vs beam static: median %.1f dB, p90 %.1f, p99 %.1f, p99.9 %.1f; MEAN %.1f dB" % tuple(list(10 * np.log10(q)) + [10 * np.log10(r.mean())]))
for th in (-15, -12, -10, -6):
    k = r < 10 ** (th / 10); print("  gate %d dB: keeps %.1f %% of PRG-snapshots, mean residual of kept %.1f dB" % (th, 100 * k.mean(), 10 * np.log10(r[k].mean())))
