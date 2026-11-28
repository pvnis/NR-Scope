"""Walker visibility of one or more sync_walk.py npz files at a fixed radio<->lidar offset: mean cell-normalised map
power at the lidar walkers' (range, Doppler) cells, against the same at far offsets (chance), over all chunks and
split by whether the tablet's DCIs were decoded in the chunk.  usage: score_split.py <offset s> <npz>..."""
import sys, numpy as np
LOST = {4, 6, 8, 10, 11, 16, 19, 21, 22, 24, 25, 27, 29, 30, 31}
OFF = float(sys.argv[1])
for name in sys.argv[2:]:
    D = np.load(name); T, Z, mpb, t_rec, TR = D["T"], D["Z"].astype(np.float32), float(D["mpb"]), float(D["t_rec"]), D["tr"]
    nf = Z.shape[2]; v = np.linspace(-1, 1, nf) * 650.173071 * 0.08 / 2; dv = v[1] - v[0]; ch = (T // 2.5).astype(int)
    ids = sorted(set(TR[:, 0].astype(int))); row = []
    for lab, sel in (("all", np.ones(len(T), bool)), ("decoded", ~np.isin(ch, list(LOST))), ("lost", np.isin(ch, list(LOST)))):
        Zs = Z[sel]; med = np.median(Zs, 0); mad = 1.4826 * np.median(np.abs(Zs - med), 0) + 1e-6; Zn = (Zs - med) / mad; Ts = T[sel]
        def sc(off):
            s = []
            for i, t in enumerate(Ts):
                tt = t_rec + off + t
                for k in ids:
                    m = TR[:, 0] == k; tk = TR[m, 1]
                    if tk[0] <= tt <= tk[-1]:
                        d = np.interp(tt, tk, TR[m, 2]); vv = np.interp(tt, tk, TR[m, 3])
                        if abs(vv) < 0.2: continue
                        b = int(round(d / mpb)); j = int(round((vv - v[0]) / dv))
                        if 1 <= b < Z.shape[1] - 1 and 1 <= j < nf - 1: s.append(Zn[i, b - 1:b + 2, j - 1:j + 2].max())
            return np.mean(s) if s else np.nan
        bg = [sc(o) for o in np.r_[np.arange(-5, 12, 1.0), np.arange(30, 40, 1.0)]]; x = sc(OFF)
        row.append("%s %.1f" % (lab, (x - np.nanmean(bg)) / np.nanstd(bg)))
    print("%-28s sigma above chance: %s" % (name.split("/")[-1], " | ".join(row)))
