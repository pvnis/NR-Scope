"""Per-chunk injection tables following the lidar walker (ids 1, 4): for every slot of the chunk, his bistatic excess e and
the sine of his angle on the array axis minus the direct path's (so_light's copy of the received signal comes from the gNB
direction). Same geometry as track_walk.py. usage: mk_injtrk.py <out dir>"""
import sys, os, json, glob, re, numpy as np
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True); OFF = 20.2
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"; LD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth"
S = json.load(open(LD + "/../session.json")); Li = S["lidar"]; KY = 6378137 * np.pi / 180; KX = KY * np.cos(np.radians(Li["lat"])); hd = np.radians(Li["heading"])
L2E = lambda x, y: x * np.array([np.cos(hd), -np.sin(hd)]) + y * np.array([np.sin(hd), np.cos(hd)])
RX = np.r_[L2E(-0.445, -0.018), Li["h"] - 0.07]; TX = np.array([(S["gnb"]["lon"] - Li["lon"]) * KX, (S["gnb"]["lat"] - Li["lat"]) * KY, S["gnb"]["h"]])
b = np.radians(90 - Li["heading"]); AX = np.array([-np.sin(b), np.cos(b), 0.0]); LTR = np.linalg.norm(TX - RX); s_g = float((TX - RX) @ AX / LTR)
MV = [np.array(t["path"]) for t in json.load(open(LD + "/movers.json"))["tracks"] if t["id"] in (1, 4)]
t_rec = float(np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz")["t_rec"])
for d in glob.glob(CHD + "/b*"):
    k = int(re.match(r"b(\d+)$", os.path.basename(d)).group(1)); a0 = json.load(open(LD + "/grid/b%d.json" % k))["a0"]
    A, E, SS = [], [], []
    for j in range(5200):
        te = t_rec + OFF + k * 2.5 + j * 5e-4
        for P in MV:
            if P[0, 0] <= te <= P[-1, 0]:
                q = np.r_[L2E(np.interp(te, P[:, 0], P[:, 1]), np.interp(te, P[:, 0], P[:, 2])), 1.0]; u = q - RX
                A.append((a0 + j) % 20480); E.append(np.linalg.norm(q - TX) + np.linalg.norm(u) - LTR); SS.append(u @ AX / np.linalg.norm(u) - s_g); break
    np.savez(os.path.join(OUT, "b%d.npz" % k), amod=np.array(A), e=np.array(E), s=np.array(SS))
print("tables written for", len(glob.glob(OUT + "/b*.npz")), "chunks")
