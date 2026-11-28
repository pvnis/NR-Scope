"""8-port CSI-RS with the receiver corrections (window shifts, delay drift, range ramp, TRS oscillator) as so_light
applies to the TRS/DM-RS; slots 0 and 10 (slot 2 has unstable port powers).
 1. per port: shape coherence over occasions (should approach 1 for fixed antenna ports);
 2. static port channels H[prb, ant, port] per slot type (occasions aligned by their common phase, averaged);
 3. tablet beam statics (BA wide) per 4-PRB group fitted as S = H v (8 weights; 4 PRBs x 4 antennas = 16 equations):
    residual vs a random 8-dim control;
 4. co-phasing read from v: ports 0-3 vs 4-7 as the two polarisations, phi = <v_lo, v_hi> / |v_lo|^2 (and its QPSK
    distance), and the same through the other pairing; stability of phi over beams on the same group.
usage: port_fit2.py <chunk>"""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/pf2.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
    "        KEEP = (E, dict(own), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, S, Sc, Rbar = g["KEEP"]
osc, files, shifts, drift, b0, PH, N, NFq = g["osc"], g["files"], g["shifts"], g["drift"], g["b0"], g["PH"], g["N"], g["NFq"]
import gs_trs_lib as L
qc = 3 * np.arange(N) + 1; qd0 = float(os.environ.get("WSH_Q0", "409.25"))
Gc = [L.load_gsc(f) for f in files]; keys = sorted(set(Gc[0]).intersection(*[set(x) for x in Gc[1:]]))
by = collections.defaultdict(list)
for a, l in keys:
    if l == 13 and a % 20 in (0, 10): by[a % 20].append(a)
def ports(a):
    t = a * 0.5e-3 + 13 / 28e3
    Pm = np.array([L.csi8_despread(*Gc[c][(a, 13)], a % 20)[:N] for c in range(4)])             # 4 x N x 8
    sh = np.array(shifts(a), float); Pm = Pm * np.exp(-1j * PH * sh[:, None, None] * (qc - qd0)[None, :, None])
    Pm = Pm * np.exp(-1j * np.array([np.polyval(drift[c], t) for c in range(4)])[:, None, None] * (qc - qd0)[None, :, None])
    Pm = Pm * np.exp(2j * np.pi * np.array(b0)[:, None, None] * (qc - qd0)[None, :, None] / NFq)
    return Pm / osc(t)
Hs = {}
for slot, A in sorted(by.items()):
    P = np.array([ports(a) for a in A]); ref = P[len(P) // 2]
    P = np.array([p * np.exp(-1j * np.angle(np.vdot(ref.ravel(), p.ravel()))) for p in P])
    coh = [[abs(np.vdot(ref[:, :, q].ravel(), P[i, :, :, q].ravel())) / (np.linalg.norm(ref[:, :, q]) * np.linalg.norm(P[i, :, :, q]) + 1e-30) for q in range(8)] for i in range(len(P))]
    H = P.mean(0); Hs[slot] = H
    flu = np.mean(np.abs(P - H) ** 2) / np.mean(np.abs(H) ** 2)
    print("slot %2d: %d occasions; per-port shape coherence median %s; fluctuation about the mean %.1f dB"
          % (slot, len(A), " ".join("%.2f" % v for v in np.median(coh, 0)), 10 * np.log10(flu)))
QP = np.array([1, 1j, -1, -1j])
for slot, H in Hs.items():
    res, rres, phis = [], [], collections.defaultdict(list); seen = set(); rng = np.random.default_rng(0)
    for (i, b), kk in own.items():
        if Sc[kk] < 20 or (kk, b) in seen: continue
        seen.add((kk, b))
        Sprb = S[kk].reshape(4, 4, 3).mean(2)                         # antenna x PRB
        A_ = H[:, 4 * b:4 * b + 4, :].reshape(16, 8); y = Sprb.reshape(16)
        v, *_ = np.linalg.lstsq(A_, y, rcond=None); res.append(np.sum(np.abs(y - A_ @ v) ** 2) / np.sum(np.abs(y) ** 2))
        Ar = (rng.standard_normal((16, 8)) + 1j * rng.standard_normal((16, 8))); vr, *_ = np.linalg.lstsq(Ar, y, rcond=None)
        rres.append(np.sum(np.abs(y - Ar @ vr) ** 2) / np.sum(np.abs(y) ** 2))
        lo, hi = v[:4], v[4:]; ph = np.vdot(lo, hi) / (np.vdot(lo, lo).real + 1e-30); phis[b].append(ph)
    allp = np.array([p for v_ in phis.values() for p in v_])
    qd = np.abs(np.angle(allp[:, None] / QP[None, :])).min(1)
    print("slot %2d ports: beam static = H v fit residual median %.1f dB (p25 %.1f, p75 %.1f), random 8-dim control %.1f dB; %d beam-groups"
          % (slot, 10 * np.log10(np.median(res)), 10 * np.log10(np.percentile(res, 25)), 10 * np.log10(np.percentile(res, 75)), 10 * np.log10(np.median(rres)), len(res)))
    print("          co-phasing from v (ports 4-7 vs 0-3): |phi| median %.2f, angle distance to the nearest of {1,j,-1,-j}: median %.2f rad (uniform would be 0.39)"
          % (np.median(np.abs(allp)), np.median(qd)))
