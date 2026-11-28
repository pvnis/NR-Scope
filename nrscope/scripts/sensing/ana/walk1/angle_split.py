"""DM-RS vs TRS walker score split by the walker's departure angle from the gNB, relative to the gNB->receiver direction.
Per window (jp_cache): energy at the lidar range and at +12..+40 m, antennas combined freely (sum over c of |.|^2), PRG factor =
the static's anchoring factor (Gst) for DM-RS, flat for TRS. Bins of |angle offset|; score = sum truth / sum mean-null.
usage: angle_split.py <chunks> <offset>"""
import sys, json, numpy as np, collections
KS = [int(x) for x in sys.argv[1].split(",")]; OFF = float(sys.argv[2])
LD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth"; S = json.load(open(LD + "/../session.json"))
lat0, lon0 = S["lidar"]["lat"], S["lidar"]["lon"]; hd = np.radians(S["lidar"]["heading"])
enu = lambda la, lo: np.array([(lo - lon0) * 111320 * np.cos(np.radians(lat0)), (la - lat0) * 111320])
TX = enu(S["gnb"]["lat"], S["gnb"]["lon"]); RX = np.zeros(2); az = lambda p: np.degrees(np.arctan2(*(p - TX)[::-1]))
TRK = []
for tr in json.load(open(LD + "/movers.json"))["tracks"]:
    if tr["id"] not in (1, 4): continue
    P = np.array(tr["path"]); e, n = P[:, 1] * np.cos(hd) + P[:, 2] * np.sin(hd), -P[:, 1] * np.sin(hd) + P[:, 2] * np.cos(hd); TRK.append((P[:, 0], e, n))
t_rec = float(np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz")["t_rec"])
def pos(te):
    for t, e, n in TRK:
        if t[0] <= te <= t[-1]: return np.array([np.interp(te, t, e), np.interp(te, t, n)])
a_rx = az(RX); print("gNB->receiver distance %.0f m, azimuth %.1f deg" % (np.linalg.norm(TX), a_rx))
acc = collections.defaultdict(lambda: np.zeros(4)); offs = []
for k in KS:
    Z = np.load("/mnt/mydrive/field/ana_walk1/jp_cache_b%d.npz" % k); rd, nd, rt, nt, G = Z["rd"], Z["nd"], Z["rt"], Z["nt"], Z["Gst"]
    c = G / np.sqrt(np.mean(np.abs(G) ** 2))
    for w in Z["valid"]:
        p = pos(t_rec + OFF + k * 2.5 + (w + 0.5) * 0.2)
        if p is None: continue
        da = (az(p) - a_rx + 180) % 360 - 180; offs.append(da)
        ed = [np.sum(np.abs(np.einsum("b,bc->c", c.conj(), rd[h, w])) ** 2) / max(np.sum(np.abs(c) ** 2 * nd[h, w]), 1e-30) for h in range(rd.shape[0])]
        et = [np.sum(np.abs(rt[h, w].sum(0)) ** 2) / max(np.sum(nt[h, w]), 1e-30) for h in range(rt.shape[0])]
        bn = "|off| < 2 deg" if abs(da) < 2 else ("2-4 deg" if abs(da) < 4 else ">= 4 deg")
        acc[bn] += [ed[0], np.mean(ed[1:]), et[0], np.mean(et[1:])]; acc[bn + " n"] += [1, 0, 0, 0]
print("walker azimuth offset from the gNB->receiver direction: %.1f .. %.1f deg (median |off| %.1f)" % (min(offs), max(offs), np.median(np.abs(offs))))
for bn in ("|off| < 2 deg", "2-4 deg", ">= 4 deg"):
    if bn in acc: v = acc[bn]; print("%-14s %3d windows: DM-RS truth/nulls %.2f | TRS %.2f" % (bn, acc[bn + " n"][0], v[0] / v[1], v[2] / v[3]))
