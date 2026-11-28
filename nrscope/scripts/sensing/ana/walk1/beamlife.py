import sys, os, io, contextlib, numpy as np, collections as CC
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="cplx", BA_W="0", MINPRB="48", TAB_ONLY="1")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/bl.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(open(__file__).read(), "so_light", "exec"))
span = CC.defaultdict(lambda: [9e9, -9e9]); n = CC.Counter()
for (i, b), k in own.items():
    t = E[i][0]; s = span[k]; s[0] = min(s[0], t); s[1] = max(s[1], t); n[k] += 1
life = np.array([span[k][1] - span[k][0] for k in span]); cnt = np.array([n[k] for k in span])
w = cnt / cnt.sum()
print("tablet beams: %d; PRG-snapshot-weighted lifetime: median %.3f s, p25 %.3f, p75 %.3f; beams per PRG over the chunk: %.1f" % (len(life), np.median(np.repeat(life, cnt)), np.percentile(np.repeat(life, cnt), 25), np.percentile(np.repeat(life, cnt), 75), len(life) / 68))
