"""The walker as probe: matched filter on so_light outputs (.bin per chunk) at the lidar walkers' excess path and Doppler,
TRS snapshots and DM-RS snapshots separately, per antenna, coherent over 120 ms windows, power summed over antennas and
windows. Null hypotheses: the same filter at wrong ranges (+12..+40 m) and with the Doppler reversed (time-reversed
phase). Output per stream: power at the truth over the nulls' mean, and z against the nulls' spread.
usage: walker_mf.py <sync npz (lidar truth, t_rec)> <offset s> <grid dir> <bin dir> <chunk list "0,1,2"> <tag>"""
import sys, json, os, collections, numpy as np
exec(open("/tmp/inj/mf_check.py").read().split("for path in sys.argv[1:]:")[0])
NPZ, OFF, GD, BD, KS, TAG = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], [int(x) for x in sys.argv[5].split(",")], sys.argv[6]
D = np.load(NPZ); TR = D["tr"]; t_rec = float(D["t_rec"]); ids = sorted(set(TR[:, 0].astype(int)))
kk = np.arange(NSC) - NSC // 2
HYP = [("truth", 0.0, 1)] + [("r+%d" % o, float(o), 1) for o in (12, 18, 25, 32, 40)] + [("rev", 0.0, -1)] + [("rev+%d" % o, float(o), -1) for o in (12, 25, 40)]
acc = {s: collections.defaultdict(float) for s in ("trs", "dmrs")}; nwin = collections.Counter()
for k in KS:
    p = os.path.join(BD, "b%d.bin" % k)
    if not os.path.exists(p): continue
    a0 = json.load(open(os.path.join(GD, "b%d.json" % k)))["a0"]; S = read(p)
    by = collections.defaultdict(list)
    for a, l, c, vb, x in S: by[c].append((a, l, vb, x))
    for c, L in by.items():
        T_ = [x for a, l, vb, x in L if l in (4, 8)]; r0 = T_[0]
        ref = np.mean([x * np.exp(-1j * np.angle(np.vdot(r0, x))) for x in T_], 0)
        W = collections.defaultdict(lambda: collections.defaultdict(complex))
        for a, l, vb, x in L:
            ph = np.vdot(ref[vb], x[vb]); x = x * np.exp(-1j * np.angle(ph)); d = x[vb] - ref[vb]
            ti = k * 2.5 + ((a - a0) % 20480) * 5e-4 + l / 28e3; te = t_rec + OFF + ti; st_ = "trs" if l in (4, 8) else "dmrs"
            w = int(ti // 0.12)
            for tid in ids:
                m = TR[:, 0] == tid; tk = TR[m, 1]
                if not (tk[0] <= te <= tk[-1]): continue
                dp = np.interp(te, tk, TR[m, 2])                     # excess path (m)
                for name, ro, sg in HYP:
                    rr = dp + ro; tt = te if sg > 0 else (2 * tk[0] - te)  # reversed Doppler: phase history run backwards
                    rph = np.interp(tt, tk, TR[m, 2]) + ro if sg < 0 else rr
                    e = ref[vb] * np.exp(-2j * np.pi * rph / LAM) * np.exp(-2j * np.pi * kk[vb] * DF * rr / 3e8)
                    W[(st_, w, tid)][name] += np.vdot(e, d) / (np.vdot(e, e).real + 1e-30)
        for (st_, w, tid), hz in W.items():
            for name, z in hz.items(): acc[st_][name] += abs(z) ** 2
            nwin[st_] += 1
for st_ in ("trs", "dmrs"):
    if not acc[st_]: continue
    nul = np.array([acc[st_][n] for n, _, _ in HYP if n != "truth"])
    print("%-12s %-4s %5d window-antenna-walker sums: truth / nulls = %.2f (%.1f dB), z = %.1f   [nulls %s]"
          % (TAG, st_, nwin[st_], acc[st_]["truth"] / nul.mean(), 10 * np.log10(acc[st_]["truth"] / nul.mean()), (acc[st_]["truth"] - nul.mean()) / (nul.std() + 1e-30),
             " ".join("%.2f" % (x / nul.mean()) for x in nul)))
