"""Radio <-> lidar time offset of an outdoor recording from the radar maps themselves: for each candidate IQ start (and
Doppler sign convention), the lidar movers predict a (range, Doppler) cell per map; score = mean normalised map power
there (best of +-1 range / +-1 Doppler bin). The right offset lines the walkers up with their echoes in range AND
Doppler sign at once, which a wrong offset does not. Also writes the per-map arrays for plotting (npz).
usage: sync_walk.py <map dir> <lidar rec dir> <recorder start, UTC 'YYYY-mm-dd HH:MM:SS.fff'> <out.npz>"""
import sys, os, json, glob, re, datetime, numpy as np
MD, LD, T0S, OUT = sys.argv[1:5]
t_rec = datetime.datetime.strptime(T0S, "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=datetime.timezone.utc).timestamp()
S = json.load(open(os.path.join(os.path.dirname(LD.rstrip("/")), "session.json"))); lat0, lon0 = S["lidar"]["lat"], S["lidar"]["lon"]
hd = np.radians(S["lidar"]["heading"]); enu = lambda la, lo: np.array([(lo - lon0) * 111320 * np.cos(np.radians(lat0)), (la - lat0) * 111320])
TX = np.r_[enu(S["gnb"]["lat"], S["gnb"]["lon"]), S["gnb"]["h"]]; RX = np.array([0, 0, 1.5]); LTR = np.linalg.norm(TX - RX); LAM = 3e8 / 3.75e9
TR = []
for tr in json.load(open(os.path.join(LD, "movers.json")))["tracks"]:
    P = np.array(tr["path"])
    if len(P) < 50: continue
    e, n = P[:, 1] * np.cos(hd) + P[:, 2] * np.sin(hd), -P[:, 1] * np.sin(hd) + P[:, 2] * np.cos(hd); Q = np.c_[e, n, np.ones(len(e))]
    d = np.linalg.norm(Q - TX, axis=1) + np.linalg.norm(Q - RX, axis=1) - LTR; ds = np.convolve(d, np.ones(5) / 5, "same")
    TR.append(dict(id=tr["id"], t=P[2:-2, 0], d=d[2:-2], v=(-np.gradient(ds, P[:, 0]) / 2)[2:-2]))
# maps
T, Z, mpb, fmx = [], [], None, None
for f in sorted(glob.glob(os.path.join(MD, "b*.map2d.csv")), key=lambda p: int(re.match(r"b(\d+)", os.path.basename(p)).group(1))):
    k = int(re.match(r"b(\d+)", os.path.basename(f)).group(1)); a0 = json.load(open(os.path.join(LD, "grid", "b%d.json" % k)))["a0"]
    for line in open(f):
        if line.startswith("#"): continue
        c = line.rstrip("\n").split(","); a = int(c[0]) * 20 + int(c[1]); nb, nf = int(c[12]), int(c[13]); mpb, fmx = float(c[9]), float(c[10])
        p = np.array(c[21:21 + nb * nf], float).reshape(nb, nf)
        lp = np.log10(np.maximum(p, 1e-12)); z = (lp - np.median(lp)) / (1.4826 * np.median(np.abs(lp - np.median(lp))) + 1e-9)
        T.append(k * 2.5 + ((a - a0) % 20480) * 5e-4 - float(c[8]) / 2); Z.append(z.astype(np.float32))
T = np.array(T); nb = Z[0].shape[0]; NF = np.array([z.shape[1] for z in Z]); VB = 2 * fmx / NF * LAM / 2   # m/s per Doppler bin, per map
vbin = float(np.median(VB))
print("%d maps over %.1f s; %d range bins x %.2f m, Doppler bins %d-%d (%.3f m/s)" % (len(T), np.ptp(T), nb, mpb, NF.min(), NF.max(), vbin))
def cells(off, sgn):
    rows = []
    for i, t in enumerate(T):
        tt = t_rec + off + t
        for R in TR:
            if R["t"][0] <= tt <= R["t"][-1]:
                d = np.interp(tt, R["t"], R["d"]); v = sgn * np.interp(tt, R["t"], R["v"])
                b = int(round(d / mpb)); j = int(round(v / VB[i])) + NF[i] // 2
                if 1 <= b < nb - 1 and 1 <= j < NF[i] - 1 and abs(v) > 0.2: rows.append((i, b, j))
    return rows
def score(off, sgn):
    r = cells(off, sgn)
    if len(r) < 50: return np.nan, len(r)
    return float(np.mean([Z[i][b - 1:b + 2, j - 1:j + 2].max() for i, b, j in r])), len(r)
offs = np.arange(-5, 40, 0.1); res = {}
for sgn in (1, -1):
    res[sgn] = np.array([score(o, sgn) for o in offs])
    s = res[sgn][:, 0]; i = int(np.nanargmax(s)); mu, sd = np.nanmedian(s), np.nanstd(s)
    print("sign %+d: best IQ start = recorder start %+.1f s, score %.2f (%d cells); median %.2f, std %.2f -> peak %.1f sigma" % (sgn, offs[i], s[i], res[sgn][i, 1], mu, sd, (s[i] - mu) / sd))
np.savez_compressed(OUT, T=T, Z=np.array([np.array([np.interp(np.linspace(-1, 1, 313), np.linspace(-1, 1, z.shape[1]), row) for row in z], np.float32) for z in Z]), offs=offs, sp=res[1], sm=res[-1], mpb=mpb, vbin=vbin, t_rec=t_rec,
                    tr=np.array([(R["id"], t, d, v) for R in TR for t, d, v in zip(R["t"], R["d"], R["v"])]))
print("wrote", OUT)
