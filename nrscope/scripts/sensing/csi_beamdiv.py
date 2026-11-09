"""How different are the real precoders across PRGs and over time? Per tablet snapshot, every PRG's weights w_b from
the nearest 8-port CSI-RS (least squares), then: within the snapshot, the share of each w_b along the snapshot's main
beam (|<w_b, u>|^2 / |w_b|^2, u = principal direction of its PRGs); across time, the same share along the burst's main
beam. 1 = the precoder only changes by a complex number (what stitching assumes), 1/8 = random direction.
usage: csi_beamdiv.py <burst dir>"""
import sys, os, numpy as np, collections as C
BD = sys.argv[1].rstrip("/") + "/"; saved = list(sys.argv); os.environ["CH"] = "0,1,2,3"
src = open("/home/sens/field/gs_all.py").read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
sys.path.insert(0, "/home/sens/field"); from gs_trs_lib import load_gsc, csi8_despread
Dd = [load(BD + "replay_ch%d/ue_dmrs.bin" % c) for c in range(4)]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
U = C.Counter(k[1] for k in keys if k[1] != 0).most_common(1)[0][0]
A = sorted((k[0] * 0.5e-3 + k[2] / 28e3, k[0], np.logical_and.reduce([d[k][0] for d in Dd]), np.array([d[k][1] for d in Dd]).astype(complex)) for k in keys if k[1] == U)
Gs = [load_gsc(BD + "replay_ch%d/ue_dmrs.bin" % c) for c in range(4)]; gk = sorted(set(Gs[0]).intersection(*[set(g) for g in Gs[1:]]))
CS = [(a * 0.5e-3 + 13 / 28e3, np.array([csi8_despread(*Gs[c][(a, l)], a % 20)[:N] for c in range(4)])) for a, l in gk if l == 13 and a % 20 in (2, 10)]
TC = np.array([c[0] for c in CS]); within = []; allw = []; nprg = []
for t, a, pm, Y in A[::6]:
    j = int(np.argmin(np.abs(TC - t)))
    if abs(TC[j] - t) > 0.003 or pm.sum() < 4: continue
    W = []
    for b in np.flatnonzero(pm):
        H = np.repeat(CS[j][1][:, 4 * b:4 * b + 4, :][:, :, None, :], 3, 2).reshape(-1, 8); y = Y[:, 12 * b:12 * b + 12].reshape(-1)
        w, *_ = np.linalg.lstsq(H, y, rcond=None); W.append(w / max(np.linalg.norm(w), 1e-30))
    W = np.array(W); ev, V = np.linalg.eigh(W.T @ W.conj()); u = V[:, -1]
    within += list(np.abs(W.conj() @ u) ** 2); allw += list(W); nprg.append(len(W))
W = np.array(allw); ev, V = np.linalg.eigh(W.T @ W.conj()); u = V[:, -1]; across = np.abs(W.conj() @ u) ** 2
print("%d snapshots, %d PRGs each (median)" % (len(nprg), np.median(nprg)))
print("share of a PRG's precoder along its snapshot's main beam: median %.2f, 10th pct %.2f" % (np.median(within), np.percentile(within, 10)))
print("share along the burst's main beam: median %.2f, 10th pct %.2f" % (np.median(across), np.percentile(across, 10)))
print("eigenvalues of the burst's precoder directions (share): %s" % np.round(ev[::-1] / ev.sum(), 2))
