"""Export the real scheduling of a 5 Nov burst for the precoded simulator: per DM-RS symbol (slot a, symbol l) every
record's RNTI (0 = blind, the whole slot's PCI-scrambled DM-RS) and its PRG mask (68 groups of 4 PRBs), plus the TRS
occasions (slot, symbol) seen in the raw-symbol records. From antenna 0's dump (the masks are the same on all four).
usage: masks_export.py <burst dir> <out.npz>"""
import sys, struct, numpy as np
BD, OUT = sys.argv[1].rstrip("/"), sys.argv[2]
N = 272; NG = 68
raw = open(BD + "/replay_blind_ch0/ue_dmrs.bin", "rb").read()
i = 0; base = None; A, L, RN, PM, TRS = [], [], [], [], set()
def unwrap(a):
    global base
    if base is None: base = a
    return a + 20480 if a < base - 1000 else a
while i + 12 <= len(raw):
    tag = raw[i:i + 4]
    if tag in (b"GSD1", b"GSD2"):
        two = tag == b"GSD2"
        if i + 20 > len(raw): break
        sfn, sl, ns, rn, npb, occ, coh = struct.unpack_from("<HBBHHff", raw, i + 4); j = i + 20; per = 3 * npb * 8 * (2 if two else 1)
        if j + npb + ns * (1 + per) > len(raw): break
        m = np.frombuffer(raw, np.uint8, npb, j).astype(bool); j += npb; a = unwrap(sfn * 20 + sl)
        for _ in range(ns):
            l = raw[j]; j += 1
            if npb >= N and occ < 4 and coh >= 0.5:      # the records gs_all.load keeps
                A.append(a); L.append(l); RN.append(rn); PM.append(m[:N].reshape(NG, 4).all(1))
            j += per
        i = j
    elif tag == b"GSC1":
        sfn, sl, l, fc, npb = struct.unpack_from("<HBBHH", raw, i + 4); a = unwrap(sfn * 20 + sl)
        if l in (4, 8) and a % 20 in (5, 6): TRS.add((a, l))
        i += 12 + 12 * npb * 8
    else:
        i += 1
A = np.array(A); a0 = min(A.min(), min(x[0] for x in TRS)) if TRS else A.min()
T = np.array(sorted(TRS))
np.savez_compressed(OUT, a=A - a0, l=np.array(L, np.uint8), rnti=np.array(RN, np.uint16), pm=np.array(PM, bool),
                    trs=T - [a0, 0] if len(T) else np.zeros((0, 2), int), a0=a0)
u, c = np.unique(RN, return_counts=True)
print("%s: %d DM-RS records over %d slots, %d TRS symbols; RNTIs (records): %s" % (
    OUT, len(A), len(set(A.tolist())), len(T), ", ".join("%s:%d" % ("blind" if r == 0 else hex(r), n) for r, n in sorted(zip(u, c), key=lambda x: -x[1])[:8])))
