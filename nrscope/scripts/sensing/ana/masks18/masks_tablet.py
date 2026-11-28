"""Tablet-only scheduling from an exported 5 Nov burst (masks_export.py): only the tablet's records (0x7304), plus a
blind record per symbol covering the same PRGs (what dump_blind would detect of it). mode real: its real PRGs; mode
full: the same slots and symbols with the whole band (68 PRGs). TRS occasions kept.
usage: masks_tablet.py <in.npz> <out.npz> <real|full> [tablet RNTI, default 0x7304]"""
import sys, numpy as np
M = np.load(sys.argv[1]); mode = sys.argv[3]; TAB = int(sys.argv[4], 0) if len(sys.argv) > 4 else 0x7304   # the tablet RNTI (5 Nov 0x7304, 18 Nov 0x7688)
k = M["rnti"] == TAB
a, l, pm = M["a"][k], M["l"][k], M["pm"][k].copy()
if mode == "full": pm[:] = True
A = np.r_[a, a]; L = np.r_[l, l]; RN = np.r_[np.full(len(a), TAB), np.zeros(len(a), int)].astype(np.uint16); PM = np.r_[pm, pm]
np.savez_compressed(sys.argv[2], a=A, l=L, rnti=RN, pm=PM, trs=M["trs"], a0=M["a0"])
print("%s: %d tablet records (%d slots), PRGs per record median %d" % (sys.argv[2], len(a), len(set(a.tolist())), int(np.median(pm.sum(1)))))
