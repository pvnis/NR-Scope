"""Real data, tablet DM-RS, wideband beam anchoring: where does an injected echo get lost? so_light run twice (with and
without INJ), the BA_PH=wide internals kept: per PRG-snapshot, the raw echo (Y_inj - Y_0), after the per-run phase
(theta_inj vs theta_0), after the beam static (S_inj vs S_0); and whether beam assignments differ."""
import sys, os, io, contextlib, numpy as np
CHK = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3"
def run(inj):
    os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
    if inj: os.environ["INJ"] = inj
    else: os.environ.pop("INJ", None)
    sys.argv = ["so_light.py", CHK, "/tmp/e2e/sw.bin"]
    g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
                                             "        KEEP = (E, dict(own), dict(THETA), {k: Ssum[k] / Scnt[k] for k in Ssum}); E = []\n")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
    return g["KEEP"]
E1, own1, th1, S1 = run("30,1.0,-30,45"); E0, own0, th0, S0 = run(None)
same = sum(1 for k in own0 if own1.get(k) == own0[k]) / len(own0)
print("PRG-snapshots %d; same beam assignment with/without the echo: %.3f" % (len(own0), same))
raw = []; aft_th = []; aft_s = []
for (i, b), k0 in list(own0.items())[::7]:
    if (i, b) not in own1: continue
    sl = slice(12 * b, 12 * b + 12); Y1 = E1[i][4][:, sl]; Y0 = E0[i][4][:, sl]
    e_raw = (Y1 - Y0) * np.exp(-1j * th0[(i, b)])                         # the echo as injected, de-phased like the clean one
    e_th = Y1 * np.exp(-1j * th1[(i, b)]) - Y0 * np.exp(-1j * th0[(i, b)])  # after each its own phase
    e_s = (Y1 * np.exp(-1j * th1[(i, b)]) - S1[own1[(i, b)]]) - (Y0 * np.exp(-1j * th0[(i, b)]) - S0[k0])   # after its beam static
    n = np.linalg.norm(e_raw) ** 2 + 1e-30
    raw.append(1.0); aft_th.append(np.vdot(e_raw, e_th) / n); aft_s.append(np.vdot(e_raw, e_s) / n)
for lab, v in (("after the per-run phase", aft_th), ("after the beam static", aft_s)):
    v = np.array(v); print("%-26s echo kept: median |gamma| %.2f, coherent |mean| %.2f" % (lab, np.median(abs(v)), abs(np.mean(v))))
dth = np.array([np.angle(np.exp(1j * (th1[k] - th0[k]))) for k in own0 if k in th1])
print("theta change caused by the echo: MAD %.3f rad (an echo at -30 dB should move it by ~0.03 at most)" % (1.4826 * np.median(np.abs(dth))))
