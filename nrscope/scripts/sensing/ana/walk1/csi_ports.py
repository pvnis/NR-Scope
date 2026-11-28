"""Are the 8 CSI-RS ports antenna ports (fixed) or beams that change between occasions? Per CSI-RS slot type, per
occasion: port channels H[prb, ant, port] (despread on symbol 13, subcarriers 0-7, FD-CDM2). For each port p:
  shape coherence over occasions = |<H_p(occ), H_p(ref)>| / (|H_p(occ)| |H_p(ref)|) over all PRBs and antennas
  (1: the same spatial/frequency shape every time, up to a scalar); port power relative to its median;
and the relative phase between ports p and 0 (complex ratio of their scalar fits to the reference) over occasions.
usage: csi_ports.py <chunk>"""
import sys, collections, numpy as np
sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing"); import gs_trs_lib as L
k = int(sys.argv[1]); B = "/mnt/mydrive/field/20261124_140050_walk1_chunks_v2/b%d" % k
L.PCI = 630; L._cache.clear()
G = [L.load_gsc("%s/replay_blind_ch%d/ue_dmrs.bin" % (B, c)) for c in range(4)]
keys = sorted(set(G[0]).intersection(*[set(x) for x in G[1:]]))
by = collections.defaultdict(list)
for a, l in keys:
    if l == 13: by[a % 20].append(a)
for slot, A in sorted(by.items()):
    if len(A) < 20: continue
    H = np.array([[L.csi8_despread(*G[c][(a, 13)], a % 20)[:272] for c in range(4)] for a in A])   # [occ, ant, prb, port]
    ref = H[len(H) // 2]
    coh = np.zeros((len(A), 8)); pw = np.zeros((len(A), 8)); rel = np.zeros((len(A), 8), complex)
    for i in range(len(A)):
        for p in range(8):
            x = H[i, :, :, p].ravel(); r = ref[:, :, p].ravel()
            coh[i, p] = abs(np.vdot(r, x)) / (np.linalg.norm(r) * np.linalg.norm(x) + 1e-30)
            pw[i, p] = np.sum(np.abs(x) ** 2); rel[i, p] = np.vdot(r, x) / np.vdot(r, r).real
    rph = np.angle(rel / rel[:, :1])                                    # port p relative to port 0, per occasion
    mad = lambda v: 1.4826 * np.median(np.abs(np.angle(np.exp(1j * (v - np.angle(np.mean(np.exp(1j * v))))))))
    print("CSI-RS slot %2d (%d occasions): per-port shape coherence vs a reference occasion, median per port: %s"
          % (slot, len(A), " ".join("%.2f" % np.median(coh[:, p]) for p in range(8))))
    print("        port power spread (dB, p10..p90, per port): %s" % " ".join("%.1f..%.1f" % tuple(10 * np.log10(np.percentile(pw[:, p] / np.median(pw[:, p]), [10, 90]))) for p in range(8)))
    print("        relative phase to port 0 over occasions, MAD (rad): %s" % " ".join("%.2f" % mad(rph[:, p]) for p in range(1, 8)))
    # do ports come in polarisation pairs? correlation of port shapes 0-3 vs 4-7 (same column, other polarisation)
    pc = [abs(np.vdot(ref[:, :, p].ravel(), ref[:, :, p + 4].ravel())) / (np.linalg.norm(ref[:, :, p]) * np.linalg.norm(ref[:, :, p + 4])) for p in range(4)]
    print("        shape similarity port p vs p+4 (candidate polarisation pairs): %s" % " ".join("%.2f" % v for v in pc))
