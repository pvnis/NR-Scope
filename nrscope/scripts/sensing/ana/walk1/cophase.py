"""Type I codebook test: per transmission, the tablet DM-RS / TRS-static ratio per 4-PRB group G_b. Model
G_b = B1(b) + phi_b B2(b), B1, B2 smooth across the band (complex polynomials of degree DEG in the group index),
phi_b in {1, j, -1, -j} (subband co-phasing between polarisations). Fitted by alternation (assign phi_b, refit B1, B2).
Residual |G - model|^2 / |G|^2 compared with: no co-phasing (G = B(b) smooth), and phi_b drawn at random (same
alternation, phi_b fixed random): the co-phasing model must beat both clearly. Also reports how stable phi_b is
from one transmission to the next on the same PRG group.  usage: cophase.py <chunk> [deg]"""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); DEG = int(sys.argv[2]) if len(sys.argv) > 2 else 2
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/cp.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n", "        KEEP = (E, Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, Rbar = g["KEEP"]
QP = np.array([1, 1j, -1, -1j])
def Gvec(Y, pm):
    bs = np.flatnonzero(pm); G = np.array([np.vdot(Rbar[:, 12 * b:12 * b + 12], Y[:, 12 * b:12 * b + 12]) / np.vdot(Rbar[:, 12 * b:12 * b + 12], Rbar[:, 12 * b:12 * b + 12]).real for b in bs])
    return bs, G
def basis(bs): x = (bs - 34) / 34.0; return np.stack([x ** d for d in range(DEG + 1)], 1)
def fit_smooth(bs, G):
    A = basis(bs); c, *_ = np.linalg.lstsq(A, G, rcond=None); return np.sum(np.abs(G - A @ c) ** 2) / np.sum(np.abs(G) ** 2)
def fit_cophase(bs, G, phi0=None, fixed=False, it=8):
    A = basis(bs); phi = phi0 if phi0 is not None else np.ones(len(bs), complex)
    for _ in range(it):
        M = np.concatenate([A, phi[:, None] * A], 1); c, *_ = np.linalg.lstsq(M, G, rcond=None)
        B1 = A @ c[:DEG + 1]; B2 = A @ c[DEG + 1:]
        if fixed: break
        phi = QP[np.argmin(np.abs(G[:, None] - (B1[:, None] + QP[None, :] * B2[:, None])), 1)]
    M = np.concatenate([A, phi[:, None] * A], 1); c, *_ = np.linalg.lstsq(M, G, rcond=None)
    return np.sum(np.abs(G - M @ c) ** 2) / np.sum(np.abs(G) ** 2), phi
rng = np.random.default_rng(0); r0, r1, rr = [], [], []; PH = collections.defaultdict(list)
for i, (t, a, l, pm, Y, u) in enumerate(E[::6]):
    bs, G = Gvec(Y, pm)
    if len(bs) < 24: continue
    r0.append(fit_smooth(bs, G))
    best = None
    for start in range(4):   # a few starts of the alternation
        e, phi = fit_cophase(bs, G, QP[rng.integers(0, 4, len(bs))] if start else None)
        if best is None or e < best[0]: best = (e, phi)
    r1.append(best[0]); rr.append(fit_cophase(bs, G, QP[rng.integers(0, 4, len(bs))], fixed=True)[0])
    ph = best[1] / best[1][0]   # relative to the first group (a global phi is absorbed by B1/B2)
    for b, p in zip(bs, ph): PH[b].append(p)
d = lambda v: 10 * np.log10(np.median(v))
print("chunk %d (DEG %d): G per 4-PRB group fitted by: smooth only %.1f dB | co-phasing {1,j,-1,-j} %.1f dB | random co-phasing %.1f dB  (%d transmissions)"
      % (k, DEG, d(r0), d(r1), d(rr), len(r1)))
stab = [max(collections.Counter(np.round(np.angle(v) / (np.pi / 2)).astype(int) % 4).values()) / len(v) for v in PH.values() if len(v) > 20]
print("          phi_b stability over transmissions (share of the most frequent value per group): median %.2f" % np.median(stab))
