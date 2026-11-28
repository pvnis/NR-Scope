"""Echo SNR of an injected target from TRS snapshots alone and from DM-RS snapshots alone (each on its own band), per
120 ms window, equal weights. Uses mf_check.read.  usage: mf_split.py <bin>..."""
import sys, collections, csv, numpy as np
src = open("/tmp/inj/mf_check.py").read().split("for path in sys.argv[1:]:")[0]
exec(src)
OFFS = (37.0, 43.0, 51.0, 57.0, 63.0, 71.0)
for path in sys.argv[1:]:
    S = read(path); tr = {int(r["slot_abs"]): float(r["range_m"]) for r in csv.DictReader(open(path + ".truth.csv"))}; A = sorted(tr)
    rA = np.array([tr[q] for q in A]); nmax = max(vb.sum() for a, l, c, vb, x in S)
    by = collections.defaultdict(list)
    for a, l, c, vb, x in S: by[c].append((a, l, vb, x))
    kk = np.arange(NSC) - NSC // 2; res = collections.defaultdict(list)
    for c, L in by.items():
        T_ = [x for a, l, vb, x in L if l in (4, 8)]; r0 = T_[0]
        ref = np.mean([x * np.exp(-1j * np.angle(np.vdot(r0, x))) for x in T_], 0)   # TRS static channel, de-rotated to the first TRS
        for a, l, vb, x in L:
            ph = np.vdot(ref[vb], x[vb]); x = x * np.exp(-1j * np.angle(ph)); r_ = ref[vb]; d = x[vb] - r_
            rt = np.interp(a, A, rA); kind = "trs" if l in (4, 8) else "dmrs"
            for off in (0.0,) + OFFS:
                st = np.exp(-2j * np.pi * (rt + off) / LAM) * np.exp(-2j * np.pi * kk[vb] * DF * (rt + off) / 3e8)
                e = r_ * st; res[(kind, off)].append((a, np.vdot(e, d) / np.vdot(e, e).real, np.sum(np.abs(d) ** 2) / np.sum(np.abs(r_) ** 2), vb.sum()))
    for kind in ("trs", "dmrs"):
        if not res[(kind, 0.0)]: continue
        t0 = min(r[0] for r in res[(kind, 0.0)]); W = collections.defaultdict(lambda: collections.defaultdict(complex))
        for off in (0.0,) + OFFS:
            for a, v, p, n in res[(kind, off)]: W[(a - t0) // 240][off] += v
        snr = [20 * np.log10(abs(w[0.0]) / np.sqrt(np.mean([abs(w[o]) ** 2 for o in OFFS]))) for w in W.values()]
        pw = np.median([p for a, v, p, n in res[(kind, 0.0)]]); bw = np.median([n for a, v, p, n in res[(kind, 0.0)]])
        print("%-13s %-4s: %5d snapshot-chains, median band %4d pilots, |x-static|^2/|static|^2 median %6.1f dB, window SNR median %5.1f dB"
              % (path, kind, len(res[(kind, 0.0)]), bw, 10 * np.log10(pw), np.median(snr)))
