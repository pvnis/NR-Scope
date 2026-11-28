"""A1 on real static paths: within a 4-PRB group, is the beam's static channel one complex number per antenna times the
TRS static (the paper's Observation 3 / App. B)? Per tablet beam (>= 50 members) and group: per antenna, fit
S_c = g_c R_c over the 12 values; residual |S - gR|^2 / |S|^2. Noise floor: the TRS static of the first half of the chunk
fitted to that of the second half the same way. Also: one complex number for ALL 4 antennas (stricter).
usage: a1_static.py <chunk>"""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/as.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
    "        KEEP = (E, dict(own), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, S, Sc, Rbar = g["KEEP"]; TRS = g["TRS"]; osc = g["osc"]
def res_per_ant(X, R):
    out = 0.0
    for c in range(4):
        gc = np.vdot(R[c], X[c]) / np.vdot(R[c], R[c]).real; out += np.sum(np.abs(X[c] - gc * R[c]) ** 2)
    return out / np.sum(np.abs(X) ** 2)
def res_all(X, R):
    gg = np.vdot(R.ravel(), X.ravel()) / np.vdot(R.ravel(), R.ravel()).real; return np.sum(np.abs(X - gg * R) ** 2) / np.sum(np.abs(X) ** 2)
bg = collections.defaultdict(set)
for (i, b), kk in own.items():
    if Sc[kk] >= 50: bg[kk].add(b)
rp, ra = [], []
for kk, bs in bg.items():
    for b in bs:
        sl = slice(12 * b, 12 * b + 12); rp.append(res_per_ant(S[kk], Rbar[:, sl])); ra.append(res_all(S[kk], Rbar[:, sl]))
X = [x[3] / osc(x[0]) for x in TRS]; R0 = X[0]; X = [x * np.exp(-1j * np.angle(np.vdot(R0, x))) for x in X]; h = len(X) // 2
T1 = np.mean(X[:h], 0); T2 = np.mean(X[h:], 0)
np_ = [res_per_ant(T2[:, 12 * b:12 * b + 12], T1[:, 12 * b:12 * b + 12]) for b in range(68)]
d = lambda v: "%.1f dB (p25 %.1f, p75 %.1f)" % tuple(10 * np.log10(np.percentile(v, [50, 25, 75])))
print("chunk %d: tablet beam static vs TRS static, per 4-PRB group: one factor per antenna -> residual %s; one factor for all antennas -> %s; %d beam-groups"
      % (k, d(rp), d(ra), len(rp)))
print("          noise floor (TRS static, half vs half, one factor per antenna): %s" % d(np_))
