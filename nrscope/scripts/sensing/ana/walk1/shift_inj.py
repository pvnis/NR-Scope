"""Time-shifted injection of the lidar walker (robustness): at time t the injected target is where the walker was at t - SHIFT
(where that exists). Writes the per-chunk injection tables (as mk_injtrk.py) and a truth CSV for import_capture: id 1 = the
injected (shifted) walker, ids 101.. = the real lidar movers. usage: shift_inj.py <shift s> <tables dir> <truth.csv>"""
import sys, os, json, glob, re, numpy as np
SH = float(sys.argv[1]); OUT, TRU = sys.argv[2], sys.argv[3]; os.makedirs(OUT, exist_ok=True); OFF = 20.2
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"; LD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth"
S = json.load(open(LD + "/../session.json")); Li = S["lidar"]; KY = 6378137 * np.pi / 180; KX = KY * np.cos(np.radians(Li["lat"])); hd = np.radians(Li["heading"])
L2E = lambda x, y: x * np.array([np.cos(hd), -np.sin(hd)]) + y * np.array([np.sin(hd), np.cos(hd)])
RX = np.r_[L2E(-0.445, -0.018), Li["h"] - 0.07]; TX = np.array([(S["gnb"]["lon"] - Li["lon"]) * KX, (S["gnb"]["lat"] - Li["lat"]) * KY, S["gnb"]["h"]])
b = np.radians(90 - Li["heading"]); AX = np.array([-np.sin(b), np.cos(b), 0.0]); LTR = np.linalg.norm(TX - RX); s_g = float((TX - RX) @ AX / LTR)
TRK = json.load(open(LD + "/movers.json"))["tracks"]; W = [np.array(t["path"]) for t in TRK if t["id"] in (1, 4)]
t_rec = float(np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz")["t_rec"])
def walker(te):
    for P in W:
        if P[0, 0] <= te <= P[-1, 0]: return np.interp(te, P[:, 0], P[:, 1]), np.interp(te, P[:, 0], P[:, 2])
    return None
for d in glob.glob(CHD + "/b*"):
    k = int(re.match(r"b(\d+)$", os.path.basename(d)).group(1)); a0 = json.load(open(LD + "/grid/b%d.json" % k))["a0"]; A, E, SS = [], [], []
    for j in range(5200):
        p = walker(t_rec + OFF + k * 2.5 + j * 5e-4 - SH)
        if p is None: continue
        q = np.r_[L2E(*p), 1.0]; u = q - RX
        A.append((a0 + j) % 20480); E.append(np.linalg.norm(q - TX) + np.linalg.norm(u) - LTR); SS.append(u @ AX / np.linalg.norm(u) - s_g)
    np.savez(os.path.join(OUT, "b%d.npz" % k), amod=np.array(A, int), e=np.array(E), s=np.array(SS))
with open(TRU, "w") as f:
    f.write("unix_s,id,lat,lon,cls\n")
    for P in W:
        for t, x, y in P:
            e, n = L2E(x, y); f.write("%.3f,1,%.8f,%.8f,person\n" % (t + SH, Li["lat"] + n / KY, Li["lon"] + e / KX))
    for tr in TRK:
        for t, x, y in tr["path"]:
            e, n = L2E(x, y); f.write("%.3f,%d,%.8f,%.8f,%s\n" % (t, 100 + tr["id"], Li["lat"] + n / KY, Li["lon"] + e / KX, tr.get("cls") or "person"))
print("shift %g s: tables -> %s, truth -> %s" % (SH, OUT, TRU))
