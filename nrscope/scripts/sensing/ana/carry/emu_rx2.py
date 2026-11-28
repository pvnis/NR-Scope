#!/usr/bin/env python3
"""EMULATED second receiver (not a measurement): a run whose receiver 1 is a real stream (the base run's DM-RS maps and
truth) and whose receiver 2 is virtual, at another place. Receiver 2's maps keep the real map times; each takes its
clutter from a random other real map (>= 20 s away: real structure, no walker there), plus the walker, detected with
probability PD (|bistatic speed| >= 0.6 m/s inside its +-60 deg view; PD_SLOW below), with Gaussian range / speed / angle
errors. Its truth is recomputed for its place (3D bistatic range, Doppler -1/lam d/dt from the lidar velocity, angle from
its boresight); its line-of-sight grid = its field of view only (buildings ignored: optimistic).
usage: emu_rx2.py BASE_RUN NEW_RUN x,y[,h] boresight_deg [pd pd_slow sig_ang_deg sig_r_m sig_v_ms seed]"""
import csv, json, os, shutil, sys, numpy as np
L = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/tracker/logs"
base, new = sys.argv[1], sys.argv[2]
p2 = [float(v) for v in sys.argv[3].split(",")]; R2 = np.array(p2 + [0.79] * (3 - len(p2))); bore2 = float(sys.argv[4])
PD, PD_SLOW, SA, SR, SV, SEED = [float(v) for v in (sys.argv[5:11] + ["0.57", "0.25", "14", "1.0", "0.2", "1"][len(sys.argv[5:11]):])]
rng = np.random.default_rng(int(SEED)); lam = 299792458 / 3.75e9; FOV = 60.0
B, N = os.path.join(L, base), os.path.join(L, new); shutil.rmtree(N, ignore_errors=True); shutil.copytree(B, N)
st = json.load(open(N + "/settings.json")); G = np.array(st["gnb"], float)
# receiver 1 := the real DM-RS stream
for a, b in (("map2d_rx2.csv", "map2d.csv"), ("map2d_rx2.csv.time.csv", "map2d.csv.time.csv"), ("ground_truth_rx2.csv", "ground_truth.csv"), ("ground_truth_rx2_cars.csv", "ground_truth_cars.csv")):
    if os.path.exists(os.path.join(B, a)): shutil.copy(os.path.join(B, a), os.path.join(N, b))
# receiver 2's truth
rows = list(csv.DictReader(open(os.path.join(B, "ground_truth_rx2.csv")))); cols = list(rows[0].keys())
def bist(P): return np.linalg.norm(P - G) + np.linalg.norm(P - R2)
out = []
for r in rows:
    r = dict(r)
    if r["path"] == "0":
        d = np.linalg.norm(G - R2); r["bistatic_range_m"] = "%.3f" % d; r["delay_s"] = "%.6e" % (d / 299792458); r["doppler_hz"] = "0.000"; r["aoa_deg"] = "%.2f" % ((np.degrees(np.arctan2(G[1] - R2[1], G[0] - R2[0])) - bore2 + 180) % 360 - 180)
    elif r["via_car"] == "1" and r["car_x"] != "nan":
        P = np.array([float(r["car_x"]), float(r["car_y"]), float(r["car_z"])]); V = np.array([float(r[k]) if r[k] != "nan" else 0.0 for k in ("car_vx", "car_vy", "car_vz")])
        b = bist(P); ug = (P - G) / np.linalg.norm(P - G); ur = (P - R2) / np.linalg.norm(P - R2)
        r["bistatic_range_m"] = "%.3f" % b; r["delay_s"] = "%.6e" % (b / 299792458); r["doppler_hz"] = "%.3f" % (-(V @ (ug + ur)) / lam)
        r["aoa_deg"] = "%.2f" % ((np.degrees(np.arctan2(P[1] - R2[1], P[0] - R2[0])) - bore2 + 180) % 360 - 180)
    out.append(r)
with open(N + "/ground_truth_rx2.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
los2 = float([r for r in out if r["path"] == "0"][0]["bistatic_range_m"])
tw = [r for r in out if r["via_car"] == "1" and r["car_x"] != "nan"]
tt = np.array([float(r["t_s"]) for r in tw]); ex = np.array([float(r["bistatic_range_m"]) - los2 for r in tw]); sp = np.array([float(r["doppler_hz"]) * lam / 2 for r in tw]); an = np.array([float(r["aoa_deg"]) for r in tw])
# receiver 2's maps: real times and headers, clutter from another real map, the walker added
hdr = [l for l in open(os.path.join(B, "map2d_rx2.csv")) if l.startswith("#")]
body = [l.strip().split(",") for l in open(os.path.join(B, "map2d_rx2.csv")) if not l.startswith("#")]
unix = [float(q["unix_s"]) for q in csv.DictReader(open(os.path.join(B, "map2d_rx2.csv.time.csv")))]
t0 = float(st["t0_unix"]); n_add = n_in = 0; lines = []
for i, v in enumerate(body):
    nb, nf = int(v[12]), int(v[13]); P = v[21:21 + nb * nf]; t = unix[i] - t0 - float(v[8]) / 2
    far = [j for j in range(len(body)) if abs(unix[j] - unix[i]) >= 20]; j = int(rng.choice(far)); u = body[j]
    unb, unf, unm, una = int(u[12]), int(u[13]), int(u[14]), int(u[15]); off = 21 + unb * unf
    mk = [u[off + 4 * q: off + 4 * q + 4] for q in range(unm)]; aoa = [u[off + 4 * unm + 5 * q: off + 4 * unm + 5 * q + 5] for q in range(una)]
    k = int(np.argmin(abs(tt - t)))
    if abs(tt[k] - t) < 0.5 and abs(an[k]) <= FOV:
        n_in += 1
        if rng.random() < (PD if abs(sp[k]) >= 0.6 else PD_SLOW):
            e = ex[k] + rng.normal(0, SR); s = sp[k] + rng.normal(0, SV); a = an[k] + rng.normal(0, SA)
            mk.append(["%.3f" % e, "%.3f" % s, "20.0", mk[0][3] if mk else "1"]); aoa.append(["%.3f" % e, "%.3f" % s, "%.2f" % a, "20.0", "10.0"]); n_add += 1
    head = v[:21]; head[14] = str(len(mk)); head[15] = str(len(aoa))
    lines.append(",".join(head + P + [x for m in mk for x in m] + [x for a_ in aoa for x in a_]))
open(N + "/map2d_rx2.csv", "w").write("".join(hdr) + "\n".join(lines) + "\n")
shutil.copy(os.path.join(B, "map2d_rx2.csv.time.csv"), N + "/map2d_rx2.csv.time.csv")
st["rx2"] = list(map(float, R2)); st["boresight2_deg"] = bore2; json.dump(st, open(N + "/settings.json", "w"), indent=1)
# receiver 2's line-of-sight grid: its field of view only (buildings ignored)
z = np.load(N + "/pd_grids.npz"); g0 = z["g0"]; x0, y0 = z["o0"]
Y, X = np.mgrid[0:g0.shape[0], 0:g0.shape[1]]; X = x0 + X; Y = y0 + Y
a = (np.degrees(np.arctan2(Y - R2[1], X - R2[0])) - bore2 + 180) % 360 - 180
np.savez(N + "/pd_grids_rx2.npz", key=np.array([*G, *R2, bore2], float), o0=z["o0"], g0=np.abs(a) <= FOV)
print("emulated receiver 2 at (%.1f, %.1f, %.2f), boresight %.0f deg: walker in its view in %d of %d maps, detection added in %d (pd %.2f / %.2f, sig angle %.0f deg)" % (*R2, bore2, n_in, len(body), n_add, PD, PD_SLOW, SA))
