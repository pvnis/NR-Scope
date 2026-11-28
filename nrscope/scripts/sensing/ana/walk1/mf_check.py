"""Transparent matched filter on a so_light output (NRSEST1 .bin) with an injected target (<bin>.truth.csv): per snapshot
(slot, symbol, chain) the relative channel D(k) = X(k) / M(k) - 1 (M: the time-mean static channel per chain after
taking out each snapshot's common phase), then per 120 ms window the coherent sum at the target's (range, Doppler) and
at empty ranges (noise). Reports SNR per window, split into TRS-only snapshots and all snapshots, with the snapshots
weighted 1 or by inverse residual power.  usage: mf_check.py <bin> [<bin>...]"""
import sys, struct, numpy as np, csv, collections
NSC = 3276; LAM = 3e8 / 3.75e9; DF = 30e3
def read(path):
    f = open(path, "rb"); h = f.read(struct.calcsize("<8sIIQdIIII")); S = []
    hs = struct.calcsize("<IHHBBBBHHIIQ")
    while True:
        b = f.read(hs)
        if len(b) < hs: break
        mg, c, _, _, _, _, _, mask, _, fr, sl, a = struct.unpack("<IHHBBBBHHIIQ", b)
        for l in [l for l in range(14) if mask >> l & 1]:
            n = struct.unpack("<I", f.read(4))[0]; vb = np.unpackbits(np.frombuffer(f.read((NSC + 7) // 8), np.uint8), bitorder="little")[:NSC].astype(bool)
            v = np.frombuffer(f.read(8 * n), np.complex64); x = np.zeros(NSC, complex); x[vb] = v
            S.append((a, l, c, vb, x))
    return S
for path in sys.argv[1:]:
    S = read(path); tr = {}
    for r in csv.DictReader(open(path + ".truth.csv")): tr[int(r["slot_abs"])] = float(r["range_m"])
    A = sorted(tr);
    # per chain: common phase per snapshot vs the chain's mean over all full-band snapshots, then static mean
    by = collections.defaultdict(list)
    for a, l, c, vb, x in S: by[c].append((a, l, vb, x))
    rows = []
    nmax = max(vb.sum() for a, l, c, vb, x in S); full = lambda vb: vb.sum() >= 0.95 * nmax
    for c, L in by.items():
        ref = np.mean([x for a, l, vb, x in L if full(vb)], 0)
        Ls = []
        for a, l, vb, x in L:
            ph = np.vdot(ref[vb], x[vb]); x = x * np.exp(-1j * np.angle(ph)); Ls.append((a, l, vb, x))
        M = np.mean([x for a, l, vb, x in Ls if full(vb)], 0); Mabs = np.abs(M) + 1e-12
        for a, l, vb, x in Ls:
            d = np.zeros(NSC, complex); d[vb] = x[vb] / M[vb] - 1; rows.append((a, l, c, vb, d, l in (4, 8)))
    t0 = min(r[0] for r in rows); kk = np.arange(NSC) - NSC // 2
    out = collections.defaultdict(list)
    wins = sorted(set((r[0] - t0) // 240 for r in rows))
    for w in wins:
        R = [r for r in rows if (r[0] - t0) // 240 == w]
        for lab, sel in (("TRS only", [r for r in R if r[5]]), ("all", R)):
            if len(sel) < 8: continue
            # residual power per snapshot (for inverse-variance weights), measured away from the target: the d values themselves
            def mf(rng_off, vsgn=1.0):
                z = 0; zw = 0; wsum = 0
                for a, l, c, vb, d, ist in sel:
                    rt = np.interp(a, A, [tr[x] for x in A]) + rng_off
                    # injected: R(t) = r0 - 2 v (t - T0): echo phase exp(-2j pi R/lam) and delay R/c
                    st = np.exp(-2j * np.pi * rt / LAM) * np.exp(-2j * np.pi * kk[vb] * DF * rt / 3e8)
                    val = np.vdot(st, d[vb]) / max(vb.sum(), 1)
                    p = np.mean(np.abs(d[vb]) ** 2) + 1e-12
                    z += val; zw += val / p; wsum += 1 / p
                return z / len(sel), zw / wsum
            s0, s0w = mf(0.0)
            nz = [mf(o) for o in (37.0, 43.0, 51.0, 57.0, 63.0, 71.0)]
            n0 = np.sqrt(np.mean([abs(x[0]) ** 2 for x in nz])); n0w = np.sqrt(np.mean([abs(x[1]) ** 2 for x in nz]))
            out[lab].append((20 * np.log10(abs(s0) / n0), 20 * np.log10(abs(s0w) / n0w), len(sel)))
    for lab, v in out.items():
        v = np.array(v); print("%-30s %-9s %2d windows, %4.0f snapshot-chains/window: SNR equal weights %5.1f dB, inverse-var weights %5.1f dB"
                               % (path.split("/")[-1], lab, len(v), v[:, 2].mean(), np.median(v[:, 0]), np.median(v[:, 1])))
