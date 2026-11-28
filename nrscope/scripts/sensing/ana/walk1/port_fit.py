"""Are the tablet's DM-RS precoders combinations of the 8 CSI-RS ports? Per CSI-RS slot type, the static 8-port channel
per PRB and antenna, H[prb, ant, port] (CSI-RS symbols de-rotated by the TRS oscillator model, averaged over the chunk).
For each tablet beam on each 4-PRB group: static S (4 antennas x 4 PRBs, the mean of each PRB's 3 values) fitted as
H v with 8 weights v (16 complex equations): residual |S - H v|^2 / |S|^2. Codebook precoding in the CSI-RS port space
gives a small residual; precoding outside it (e.g. from uplink sounding on all antennas) does not. Baseline: the TRS
static fitted the same way (the TRS is not precoded through the ports, but a wide beam may still be in their span).
usage: port_fit.py <chunk>"""
import sys, os, io, contextlib, collections, numpy as np
k = int(sys.argv[1]); CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="wide", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/pf.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open(g["__file__"]).read().replace("        E = []   # the per-PRG path below then writes nothing\n",
    "        KEEP = (E, dict(own), {kk: Ssum[kk] / Scnt[kk] for kk in Ssum}, Scnt.copy(), Rbar); E = []\n")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src, "so_light", "exec"), g)
E, own, S, Sc, Rbar = g["KEEP"]; osc = g["osc"]; files = g["files"]; dedrift = g["dedrift"]; undo_wsh = g["undo_wsh"]; shifts = g["shifts"]; ramp = g["ramp"]
import gs_trs_lib as L
Gc = [L.load_gsc(f) for f in files]
keys = sorted(set(Gc[0]).intersection(*[set(x) for x in Gc[1:]]))
by = collections.defaultdict(list)
for a, l in keys:
    if l == 13: by[a % 20].append(a)
N = 272
def csi_ports(a):   # -> H[prb, ant, port] for one CSI-RS occasion, oscillator out (as the DM-RS in E: dedrift, window, ramp)
    t = a * 0.5e-3 + 13 / 28e3; out = np.zeros((N, 4, 8), complex)
    for c in range(4):
        P = L.csi8_despread(*Gc[c][(a, 13)], a % 20)[:N]           # N x 8
        out[:, c, :] = P
    return t, out
res = {}
for slot, A in sorted(by.items()):
    if len(A) < 20: continue
    Hs = []
    for a in A:
        t, H = csi_ports(a); Hs.append(H / osc(t))
    Hs = np.array(Hs); ref = Hs[0]
    Hs = np.array([h * np.exp(-1j * np.angle(np.vdot(ref.ravel(), h.ravel()))) for h in Hs])   # common phase out
    res[slot] = Hs.mean(0)
    tv = np.mean(np.abs(Hs - Hs.mean(0)) ** 2) / np.mean(np.abs(Hs.mean(0)) ** 2)
    print("CSI-RS slot %d: %d occasions, port-channel fluctuation about its mean %.1f dB" % (slot, len(A), 10 * np.log10(tv)))
def fit(Sprb, H):   # Sprb: 4 PRBs x 4 antennas; H: 4 PRBs x 4 antennas x 8 ports -> residual dB
    A = H.reshape(-1, 8); y = Sprb.reshape(-1); v, *_ = np.linalg.lstsq(A, y, rcond=None)
    return 10 * np.log10(np.sum(np.abs(y - A @ v) ** 2) / np.sum(np.abs(y) ** 2)), v
for slot, H in res.items():
    rb, rt = [], []
    seen = set()
    for (i, b), kk in own.items():
        if Sc[kk] < 20 or kk in seen: continue
        seen.add(kk); Sb = S[kk]                                   # 4 x 12: antenna x (4 PRBs x 3)
        Sprb = Sb.reshape(4, 4, 3).mean(2).T                       # 4 PRBs x 4 antennas
        rb.append(fit(Sprb, H[4 * b:4 * b + 4])[0])
    for b in range(68):
        Rprb = Rbar[:, 12 * b:12 * b + 12].reshape(4, 4, 3).mean(2).T
        rt.append(fit(Rprb, H[4 * b:4 * b + 4])[0])
    # control: a random 8-dim subspace of the same shape (what a fit of 16 values with 8 free weights gives by chance)
    rng = np.random.default_rng(0); rr = []
    for b in range(68):
        Rprb = Rbar[:, 12 * b:12 * b + 12].reshape(4, 4, 3).mean(2).T; Hr = rng.standard_normal((4, 4, 8)) + 1j * rng.standard_normal((4, 4, 8)); rr.append(fit(Rprb, Hr)[0])
    print("ports of slot %d: tablet beam statics fit residual median %.1f dB (p25 %.1f, p75 %.1f, %d beams); TRS static %.1f dB; random 8-dim control %.1f dB"
          % (slot, np.median(rb), np.percentile(rb, 25), np.percentile(rb, 75), len(rb), np.median(rt), np.median(rr)))
