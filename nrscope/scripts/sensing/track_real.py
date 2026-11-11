"""Track NR-Scope's real-data maps (map2d.csv + injected truth) with the tracker of ~/GIOVANNI/sim/track.py and score
the tracks: confirmed tracks never within 6 m of a target (false tracks per burst), share of target-maps with a confirmed
track within 6 m (coverage), GOSPA (c = 5 m).
Geometry (stand-in, exact positions not surveyed): receiver at the origin 15 m up (BC 4th floor), gNB 1 km behind the
array (Cubotron, 25 m up), array facing away from the gNB (boresight -y). Injected targets have a fixed angle and a
range that moves with their speed, i.e. they move along their ray, which this geometry localises consistently.
usage: track_real.py <glob of map dirs> [label]"""
import sys, glob, os, csv, numpy as np, collections as C
sys.path.insert(0, "/home/sens/GIOVANNI/sim"); import track as TK
T = np.array([0.0, 1000.0, 25.0]); R = np.array([0.0, 0.0, 15.0]); BORE_DEG = -90.0; LTR = float(np.linalg.norm(R - T))
geo = (T, R, LTR, np.deg2rad(BORE_DEG))
def read_maps(p):
    out = []
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split(",")
        fr, sl = int(f[0]), int(f[1]); tsp = float(f[8]); mpb = float(f[9]); fc = float(f[11])
        nb, nf, nm, na = int(f[12]), int(f[13]), int(f[14]), int(f[15]); p0 = 21 + nb * nf
        mk = [tuple(float(f[p0 + 4 * i + j]) for j in range(4)) for i in range(nm)]; q0 = p0 + 4 * nm
        ao = [tuple(float(f[q0 + 5 * i + j]) for j in range(5)) for i in range(na)]
        out.append(dict(fr=fr, sl=sl, mpb=mpb, res=max(3e8 / fc / 2 / max(tsp, 1e-3), TK.SPEED_FLOOR_MS), mk=mk, aoa=ao))
    return out
def xy_of(e, ang):
    p = TK.localize(e, ang, T, R, BORE_DEG); return None if p is None else np.array(p[:2])
res = C.defaultdict(list); lab = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1]
for d in sorted(glob.glob(sys.argv[1])):
    p = d + "/map2d.csv"
    if not os.path.exists(p) or not os.path.exists(p + ".truth.csv"): continue
    TR = C.defaultdict(list)
    for r in csv.DictReader(open(p + ".truth.csv")): TR[(int(r["frame"]), int(r["slot"]))].append((float(r["range_m"]), float(r["angle_deg"])))
    frames, truths = [], []; a_prev = None; off = 0
    for m in read_maps(p):
        a = m["fr"] * 20 + m["sl"]
        if a_prev is not None and a + off < a_prev - 10000: off += 20480
        a += off; a_prev = a
        meas = []
        for rr, vv, snr, ver in m["mk"]:
            if ver != 0: continue
            cand = [x for x in m["aoa"] if abs(x[0] - rr) <= m["mpb"] and abs(x[1] - vv) <= 1.5 * m["res"]]
            ang = max(cand, key=lambda x: x[3])[2] if cand else None
            xy = xy_of(rr, ang) if ang is not None else None
            meas.append(dict(e=rr, v=vv, theta=np.deg2rad(ang) if (ang is not None and xy is not None) else None, xy=xy, p_db=snr))
        frames.append((a * 0.5e-3, meas, m["mpb"], m["res"]))
        truths.append([q for q in (xy_of(e, ang) for e, ang in TR.get((m["fr"], m["sl"]), [])) if q is not None])
    if len(frames) < 4: continue
    out = TK.run_ekf_tracker(frames, geo)
    gos = []; near = C.defaultdict(bool); seen = set(); cn = cd = 0
    for trks, truth in zip(out, truths):
        est = [(tq[1], tq[2]) for tq in trks]; gos.append(TK.xy_gospa(est, [tuple(t) for t in truth], 5.0)[0])
        for tq in trks:
            seen.add(tq[0])
            if any(np.hypot(tq[1] - x, tq[2] - y) <= 6.0 for x, y in truth): near[tq[0]] = True
        for x, y in truth: cd += 1; cn += any(np.hypot(tq[1] - x, tq[2] - y) <= 6.0 for tq in trks)
    res["false"].append(sum(1 for i in seen if not near[i])); res["cov"].append(cn / max(cd, 1)); res["gospa"].append(np.mean(gos))
    res["dets"].append(np.mean([len(f[1]) for f in frames]))
print("%-34s bursts %2d | detections/map %.1f | false tracks/burst %.1f | coverage %3.0f%% | GOSPA %.1f m" % (
    lab, len(res["false"]), np.mean(res["dets"]), np.mean(res["false"]), 100 * np.mean(res["cov"]), np.mean(res["gospa"])))
