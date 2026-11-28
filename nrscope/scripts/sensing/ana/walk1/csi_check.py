"""Does the 8-port CSI-RS of a cell despread coherently? For symbol 13 of the CSI-RS slots: 8 consecutive subcarriers
starting at k0 (0..4) per PRB, FD-CDM2 pairs (ports 2i, 2i+1 = (z_2i +- z_2i+1)/2), sequence m = 2n + (k & 1) with
scrambling ID = PCI (and 30 random IDs as the noise reference). Coherence per port across neighbouring PRBs.
usage: csi_check.py <ue_dmrs.bin> <PCI>"""
import sys, collections, numpy as np
sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
import gs_trs_lib as L
f, PCI = sys.argv[1], int(sys.argv[2])
G = L.load_gsc(f)
keys = sorted(k for k in G if k[1] == 13)
by = collections.defaultdict(list)
for a, l in keys: by[a % 20].append(a)
def despread(fc, z, slot, nid, k0):
    L.PCI = nid; L._cache.clear()
    npb = len(z) // 12; r = L.seq(slot, 13, 2 * (fc + npb) + 4); n = fc + np.arange(npb)
    zz = np.empty((npb, 8), complex)
    for k in range(8): zz[:, k] = z[12 * np.arange(npb) + k0 + k] * np.conj(r[2 * n + (k & 1)]) * 2
    P = np.empty((npb, 8), complex)
    for i in range(4): P[:, 2 * i] = (zz[:, 2 * i] + zz[:, 2 * i + 1]) / 2; P[:, 2 * i + 1] = (zz[:, 2 * i] - zz[:, 2 * i + 1]) / 2
    return P
def coh(P):
    return float(np.median([abs(np.sum(P[1:, p] * np.conj(P[:-1, p]))) / (np.sum(abs(P[:, p]) ** 2) + 1e-30) for p in range(8)]))
rng = np.random.default_rng(1); RID = rng.choice(1008, 30, replace=False)
print("CSI-RS slots (mod 20):", {s: len(v) for s, v in sorted(by.items())})
for s, A in sorted(by.items()):
    if len(A) < 5: continue
    for k0 in range(5):
        c = np.median([coh(despread(*G[(a, 13)], a % 20, PCI, k0)) for a in A[:12]])
        cr = np.percentile([coh(despread(*G[(a, 13)], a % 20, int(nid), k0)) for a in A[:4] for nid in RID[:8]], 95)
        print("slot %2d, subcarriers %d-%d: port coherence median %.2f (random IDs p95 %.2f)%s" % (s, k0, k0 + 7, c, cr, "   <==" if c > 0.6 and c > 2 * cr else ""))
L.PCI = PCI
