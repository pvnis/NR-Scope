"""Data for the DM-RS sparsity figure (18 Nov recording, PCI 632), from the twin's grid masks
(/mnt/mydrive/field/precsim_masks18/b<k>.npz, written by masks_export.py) of all 48 chunks:
  per 146 ms window (292 slots, non-overlapping): our UE's (0x7688) DM-RS symbols = sensing snapshots;
  the slot x PRG grid of one window whose count is the median: 0 UL, 1 DL without DM-RS,
  2 DM-RS not decoded (blind), 3 DM-RS of another decoded UE, 4 our UE's DM-RS; and the TRS slots.
usage: extract.py OUT.npz"""
import sys, glob, numpy as np
D = "/mnt/mydrive/field/precsim_masks18"; TAB = 0x7688; W = 292; NG = 68
counts, wins = [], []
files = [f for f in glob.glob(D + "/b[0-9]*.npz") if "_tab" not in f]
for f in sorted(files, key=lambda s: int(s.split("/b")[-1].split(".")[0])):
    z = np.load(f); a, r = z["a"], z["rnti"]; k = int(f.split("/b")[-1].split(".")[0])
    for w0 in range(0, a.max() + 1 - W, W):
        m = (a >= w0) & (a < w0 + W)
        counts.append(int(((r == TAB) & m).sum())); wins.append((k, w0))
counts = np.array(counts); med = np.median(counts)
i = int(np.argmin(np.abs(counts - med))); k, w0 = wins[i]
z = np.load(D + "/b%d.npz" % k); a, r, pm, trs = z["a"], z["rnti"], z["pm"], z["trs"]
G = np.zeros((W, NG), np.uint8)
ul = np.array([(w0 + s) % 5 == 4 for s in range(W)])          # DDDSU: slot 4 of 5 is uplink
G[~ul] = 1
m = (a >= w0) & (a < w0 + W)
for pri, sel in ((2, r == 0), (3, (r != 0) & (r != TAB)), (4, r == TAB)):
    for ai, p in zip(a[m & sel], pm[m & sel]):
        g = G[ai - w0]; g[p & (g < pri)] = pri
tr = np.unique(trs[(trs[:, 0] >= w0) & (trs[:, 0] < w0 + W), 0]) - w0
dl = int((~ul).sum())
print("files %d, windows %d; our-UE snapshots per window: median %.0f, p25 %.0f, p75 %.0f, min %d, max %d; chosen b%d slot %d (%d)" % (
      len(files), len(counts), med, *np.percentile(counts, [25, 75]), counts.min(), counts.max(), k, w0, counts[i]))
print("chosen window: DL slots %d, slots with our DM-RS %d; share of DL slot x PRG cells: ours %.1f%%, any DM-RS %.1f%%" % (
      dl, int((G == 4).any(1).sum()), 100 * (G == 4).sum() / (dl * NG), 100 * (G >= 2).sum() / (dl * NG)))
np.savez_compressed(sys.argv[1], counts=counts, grid=G, trs=tr, chunk=k, w0=w0, dense=2 * dl, slot_ms=0.5)
