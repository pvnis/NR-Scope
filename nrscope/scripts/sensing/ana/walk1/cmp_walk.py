"""Radar map vs lidar truth for one outdoor recording: per map, its detections (range = excess path over the direct path,
speed, SNR, verdict) against the bistatic excess path and its rate of every lidar mover labelled person/car.
Geometry: gNB and lidar (= receiver cart) from session.json, lidar frame turned by its heading (azimuth 0 = +Y).
Time: map slot -> IQ time through the chunk's first decoded slot (grid/b<k>.json a0, k * 2.5 s); the IQ start in
lidar time is fitted (offset scan): the one that puts the most target detections within 3 m of a truth.
usage: cmp_walk.py <map dir (b<k>.map2d.csv)> <lidar rec dir> <IQ start guess, epoch s> <out.png> [label of k list]"""
import sys, os, json, glob, re, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
MD, LD, T0G, OUT = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
S = json.load(open(os.path.join(os.path.dirname(LD.rstrip("/")), "session.json")))
lat0, lon0 = S["lidar"]["lat"], S["lidar"]["lon"]; hd = np.radians(S["lidar"]["heading"])
def enu(lat, lon): return np.array([(lon - lon0) * 111320 * np.cos(np.radians(lat0)), (lat - lat0) * 111320])
TX = np.r_[enu(S["gnb"]["lat"], S["gnb"]["lon"]), S["gnb"]["h"]]; RX = np.array([0.0, 0.0, 1.5]); LTR = np.linalg.norm(TX - RX)
M = json.load(open(os.path.join(LD, "movers.json"))); LB = json.load(open(os.path.join(LD, "labels.json"))).get("tracks", {})
truth = []
for tr in M["tracks"]:
    cls = LB.get(str(tr["id"]), {}).get("cls", tr.get("cls"))
    if cls in ("ignore", None) or len(tr["path"]) < 10: continue
    P = np.array(tr["path"]); x, y = P[:, 1], P[:, 2]
    e, n = x * np.cos(hd) + y * np.sin(hd), -x * np.sin(hd) + y * np.cos(hd)
    Q = np.c_[e, n, np.full(len(e), 1.0)]
    d = np.linalg.norm(Q - TX, axis=1) + np.linalg.norm(Q - RX, axis=1) - LTR
    truth.append(dict(id=tr["id"], cls=cls, name=LB.get(str(tr["id"]), {}).get("name", ""), t=P[:, 0], d=d, v=np.gradient(d, P[:, 0])))
# detections, in IQ time
det = []
for f in sorted(glob.glob(os.path.join(MD, "b*.map2d.csv")), key=lambda p: int(re.match(r"b(\d+)", os.path.basename(p)).group(1))):
    k = int(re.match(r"b(\d+)", os.path.basename(f)).group(1))
    g = os.path.join(LD, "grid", "b%d.json" % k)
    if not os.path.exists(g): continue
    a0 = json.load(open(g))["a0"]
    for line in open(f):
        if line.startswith("#"): continue
        c = line.rstrip("\n").split(","); a = int(c[0]) * 20 + int(c[1])
        ti = k * 2.5 + ((a - a0) % 20480) * 0.5e-3 - float(c[8]) / 2   # the map's centre
        nb, nf, nm = int(c[12]), int(c[13]), int(c[14]); p0 = 21 + nb * nf
        for i in range(nm):
            r, s, snr, vd = (float(v) for v in c[p0 + 4 * i: p0 + 4 * i + 4]); det.append((ti, r, s, snr, int(vd)))
det = np.array(det)
print("%d detections, %d accepted targets, over %.1f s of IQ" % (len(det), (det[:, 4] == 0).sum(), np.ptp(det[:, 0]) if len(det) else 0))
def score(off, sel):
    n = 0
    for t, r in det[sel][:, :2]:
        for T in truth:
            tt = T0G + off + t
            if T["t"][0] <= tt <= T["t"][-1] and abs(np.interp(tt, T["t"], T["d"]) - r) < 3: n += 1; break
    return n
sel = det[:, 4] == 0
offs = np.arange(-15, 15.01, 0.25); sc = [score(o, sel) for o in offs]; ob = offs[int(np.argmax(sc))]
print("IQ start offset vs guess: %+.2f s (%d of %d targets within 3 m of a truth; at 0: %d)" % (ob, max(sc), sel.sum(), sc[len(offs) // 2]))
fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
for T in truth:
    tt = T["t"] - T0G - ob; ax[0].plot(tt, T["d"], lw=2, alpha=.6, label="%s %d %s" % (T["cls"], T["id"], T["name"])); ax[1].plot(tt, -T["v"] / 2, lw=2, alpha=.6)
for vd, mk, cl, lb in ((3, "x", "0.6", "uncertain"), (1, ".", "0.75", "replica"), (0, "o", "C3", "target")):
    s = det[:, 4] == vd
    ax[0].scatter(det[s, 0], det[s, 1], s=14, marker=mk, c=cl, label=lb, zorder=3); ax[1].scatter(det[s, 0], det[s, 2], s=14, marker=mk, c=cl, zorder=3)
ax[0].set_ylabel("excess path (m)"); ax[0].set_ylim(0, 120); ax[0].legend(fontsize=7, ncol=3)
ax[1].set_ylabel("speed (m/s); truth: -d(path)/dt / 2"); ax[1].set_ylim(-30, 30); ax[1].set_xlabel("s into the IQ (lidar shifted by the fitted start)")
ax[0].set_title("%s: map detections vs lidar truth (offset %+.2f s)" % (os.path.basename(LD.rstrip("/")), ob))
fig.tight_layout(); fig.savefig(OUT, dpi=110); print("wrote", OUT)
