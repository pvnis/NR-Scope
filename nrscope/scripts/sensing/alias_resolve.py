"""Resolve the TRS/CSI-RS speed ambiguity (copies every ALIAS m/s) from the range drift across maps: detections are
chained across successive maps by range (within GATE m), the chain's range-vs-time slope gives the range rate
(= -2 v in the map convention), and each detection's speed is replaced by the copy v + k*ALIAS closest to the
range-rate speed. Scored against the injected truth (right speed and within 10 m), before and after.
usage: alias_resolve.py <map2d.csv> [...]"""
import sys, csv, os, numpy as np
ALIAS = 2.0; GATE = 4.0; RESMAX = float(os.environ.get("RESMAX", "1.5"))
sys_argv = list(sys.argv); sys.argv = ["x"]; exec(open("/home/sens/field/so_score.py").read().split("import os")[0]); sys.argv = sys_argv
def score(M, T, fix):
    f = nt = 0
    for i, m in enumerate(M):
        acc = fix[i]
        for rr, ss in T.get((m["frame"], m["slot"]), []):
            nt += 1; f += any(abs(v - ss) <= 2 * m["res"] and abs(r - rr) <= 10 for r, v in acc)
    return f, nt
for p in sys.argv[1:]:
    M = maps(p); T = {}; tm = {}
    for r in csv.DictReader(open(p + ".truth.csv")):
        T.setdefault((int(r["frame"]), int(r["slot"])), []).append((float(r["range_m"]), float(r["speed_ms"]))); tm[(int(r["frame"]), int(r["slot"]))] = int(r["slot_abs"])
    t = np.array([tm.get((m["frame"], m["slot"]), np.nan) * 0.5e-3 for m in M])
    raw = [[(x[0], x[1]) for x in m["markers"] if x[3] == 0] for m in M]
    fixed = []
    for i, dets in enumerate(raw):
        out = []
        for r, v in dets:
            pts = [(t[i], r)]
            for j in range(len(M)):                 # nearest detection in every other map within GATE*(maps apart)
                if j == i or not np.isfinite(t[j]): continue
                c = [rr for rr, vv in raw[j] if abs(rr - r) <= GATE * (1 + abs(j - i) / 3)]
                if c: pts.append((t[j], min(c, key=lambda x: abs(x - r))))
            if len(pts) >= 4 and np.isfinite(t[i]):
                tt, rr_ = np.array(pts).T; pf = np.polyfit(tt, rr_, 1); slope = pf[0]; vr = -slope / 2
                res = np.sqrt(np.mean((rr_ - np.polyval(pf, tt)) ** 2))
                # only when the chain is clean and clearly disagrees with the measured speed
                k = np.round((vr - v) / ALIAS) if (res <= RESMAX and abs(vr - v) > 1.0) else 0
                out.append((r, v + k * ALIAS))
            else: out.append((r, v))
        fixed.append(out)
    f0, n = score(M, T, raw); f1, _ = score(M, T, fixed)
    print("%-30s found before %3d/%-3d after range-drift alias resolution %3d/%-3d" % (os.path.dirname(p), f0, n, f1, n))
