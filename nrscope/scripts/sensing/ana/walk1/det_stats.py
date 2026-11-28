"""Per-map detections against the lidar: a map 'hits' when a detection (verdict target, or also uncertain) lies within
DR m of a walker's excess path and DV m/s of its bistatic speed (sign included); every other detection within the
plausible region (5-80 m, 0.15-3 m/s) is a false one. Per variant: hit rate over maps with a walker in view, false
detections per map and per minute.  usage: det_stats.py <offset s> <lidar dir> <map dir>..."""
import sys, os, json, glob, re, numpy as np
sys.argv[1:1] = []
OFF, LD = float(sys.argv[1]), sys.argv[2]
exec(open("/tmp/sync_walk.py").read().split("# maps")[0].replace("MD, LD, T0S, OUT = sys.argv[1:5]", "MD = None; T0S = '2026-11-24 14:00:56.515'"))
DR, DV = 3.0, 0.5
for MD in sys.argv[3:]:
    hits = {0: 0, 3: 0}; nmaps = 0; fa = {0: 0, 3: 0}; dur = 0
    for f in sorted(glob.glob(os.path.join(MD, "b*.map2d.csv"))):
        k = int(re.match(r"b(\d+)", os.path.basename(f)).group(1)); a0 = json.load(open(os.path.join(LD, "grid", "b%d.json" % k)))["a0"]
        for line in open(f):
            if line.startswith("#"): continue
            c = line.rstrip("\n").split(","); a = int(c[0]) * 20 + int(c[1]); ts = float(c[8])
            ti = k * 2.5 + ((a - a0) % 20480) * 5e-4 - ts / 2; te = t_rec + OFF + ti
            nb, nf, nm = int(c[12]), int(c[13]), int(c[14]); p0 = 21 + nb * nf
            mk = [tuple(float(x) for x in c[p0 + 4 * i:p0 + 4 * i + 4]) for i in range(nm)]
            tw = [(np.interp(te, R["t"], R["d"]), np.interp(te, R["t"], R["v"])) for R in TR if R["t"][0] <= te <= R["t"][-1]]
            tw = [(d, v) for d, v in tw if abs(v) >= 0.2]
            if tw: nmaps += 1
            dur += ts
            for want in (0, 3):
                sel = [m for m in mk if m[3] == 0 or (want == 3 and m[3] == 3)]
                hit = any(abs(m[0] - d) < DR and abs(m[1] - v) < DV for m in sel for d, v in tw)
                if tw and hit: hits[want] += 1
                fa[want] += sum(1 for m in sel if 5 < m[0] < 80 and 0.15 < abs(m[1]) < 3 and not any(abs(m[0] - d) < DR and abs(m[1] - v) < DV for d, v in tw))
    nm_ = max(len(glob.glob(os.path.join(MD, "b*.map2d.csv"))), 1)
    print("%-34s maps with a moving walker: %d | hit rate: targets %.0f %%, targets+uncertain %.0f %% | false detections per map: %.2f / %.2f"
          % (os.path.basename(MD.rstrip("/")), nmaps, 100 * hits[0] / max(nmaps, 1), 100 * hits[3] / max(nmaps, 1), fa[0] / max(1, sum(1 for _ in [0]) * 1) / max(nmaps, 1), fa[3] / max(nmaps, 1)))
