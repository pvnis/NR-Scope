"""Walker visibility split by the walker's angle off the gNB -> receiver line as seen from the gNB (a DM-RS beam aimed at
the cart illuminates the walker only near that line). Same score as score_split (cell-normalised map power at the
lidar walkers' cells vs far offsets), per angle class.  usage: score_angle.py <offset s> <lidar dir> <npz>..."""
import sys, os, json, numpy as np
OFF, LD = float(sys.argv[1]), sys.argv[2]
S = json.load(open(os.path.join(os.path.dirname(LD.rstrip("/")), "session.json"))); lat0, lon0 = S["lidar"]["lat"], S["lidar"]["lon"]; hd = np.radians(S["lidar"]["heading"])
enu = lambda la, lo: np.array([(lo - lon0) * 111320 * np.cos(np.radians(lat0)), (la - lat0) * 111320])
TX = enu(S["gnb"]["lat"], S["gnb"]["lon"]); u0 = -TX / np.linalg.norm(TX)          # gNB -> receiver direction
ANG = {}
for tr in json.load(open(os.path.join(LD, "movers.json")))["tracks"]:
    P = np.array(tr["path"])
    if len(P) < 50: continue
    e, n = P[:, 1] * np.cos(hd) + P[:, 2] * np.sin(hd), -P[:, 1] * np.sin(hd) + P[:, 2] * np.cos(hd)
    w = np.c_[e, n] - TX; ang = np.degrees(np.arccos(np.clip((w @ u0) / np.linalg.norm(w, axis=1), -1, 1)))
    ANG[tr["id"]] = (P[:, 0], ang)
CL = (("< 4 deg", 0, 4), ("4-8 deg", 4, 8), ("> 8 deg", 8, 99))
for name in sys.argv[3:]:
    D = np.load(name); T, Z, mpb, t_rec, TR = D["T"], D["Z"].astype(np.float32), float(D["mpb"]), float(D["t_rec"]), D["tr"]
    nf = Z.shape[2]; v = np.linspace(-1, 1, nf) * 650.173071 * 0.04; dv = v[1] - v[0]; ids = sorted(set(TR[:, 0].astype(int)))
    med = np.median(Z, 0); mad = 1.4826 * np.median(np.abs(Z - med), 0) + 1e-6; Zn = (Z - med) / mad
    def sc(off, lo, hi):
        s = []
        for i, t in enumerate(T):
            tt = t_rec + off + t
            for k in ids:
                m = TR[:, 0] == k; tk = TR[m, 1]
                if not (tk[0] <= tt <= tk[-1]) or k not in ANG: continue
                a_ = np.interp(t_rec + OFF + t, *ANG[k])            # the class is fixed by the TRUE time (same cells for every offset)
                if not (lo <= a_ < hi): continue
                d = np.interp(tt, tk, TR[m, 2]); vv = np.interp(tt, tk, TR[m, 3])
                if abs(vv) < 0.2: continue
                b = int(round(d / mpb)); j = int(round((vv - v[0]) / dv))
                if 1 <= b < Z.shape[1] - 1 and 1 <= j < nf - 1: s.append(Zn[i, b - 1:b + 2, j - 1:j + 2].max())
        return (np.mean(s) if s else np.nan), len(s)
    row = []
    for lab, lo, hi in CL:
        x, n = sc(OFF, lo, hi); bg = [sc(o, lo, hi)[0] for o in np.r_[np.arange(-5, 12, 1.0), np.arange(30, 40, 1.0)]]
        row.append("%s: %.1f (%d cells)" % (lab, (x - np.nanmean(bg)) / np.nanstd(bg), n))
    print("%-26s %s" % (name.split("/")[-1], " | ".join(row)))
