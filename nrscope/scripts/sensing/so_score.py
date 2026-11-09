"""Score NR-Scope sensing maps (map2d.csv from sensing_offline) against the injected targets (map2d.csv.truth.csv).
A target counts as found in a map when an accepted detection (marker verdict 0) lies within max(2 range bins, 3 m) and
2 speed cells of its truth at the map's slot. Other accepted detections are "other": false detections or real movers
of the recorded scene (compare with a run without injection).
usage: so_score.py <map2d.csv> [...]"""
import sys, csv, collections as C, numpy as np

def maps(path):
    out = []
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split(",")
        frame, slot = int(f[0]), int(f[1]); t_span = float(f[8]); mpb = float(f[9]); fc = float(f[11])
        nb, nf, nm, na = int(f[12]), int(f[13]), int(f[14]), int(f[15]); p0 = 21 + nb * nf
        mk = [(float(f[p0 + 4 * i]), float(f[p0 + 4 * i + 1]), float(f[p0 + 4 * i + 2]), int(f[p0 + 4 * i + 3])) for i in range(nm)]
        out.append(dict(frame=frame, slot=slot, t_span=t_span, mpb=mpb, res=3e8 / fc / 2 / max(t_span, 1e-3), markers=mk))
    return out

import os
for path in sys.argv[1:]:
    if not os.path.exists(path): print("%-34s maps  0 (no map produced)" % path); continue
    M = maps(path); T = C.defaultdict(list)
    try:
        for r in csv.DictReader(open(path + ".truth.csv")):
            if int(r["active"]): T[(int(r["frame"]), int(r["slot"]))].append((int(r["target"]), float(r["range_m"]), float(r["speed_ms"])))
    except FileNotFoundError:
        pass
    found = C.Counter(); n_t = C.Counter(); other = []; snr = []
    for m in M:
        acc = [x for x in m["markers"] if x[3] == 0]; used = set()
        for tid, r, v in T.get((m["frame"], m["slot"]), []):
            n_t[tid] += 1
            ok = [i for i, x in enumerate(acc) if abs(x[0] - r) <= max(2 * m["mpb"], 3.0) and abs(x[1] - v) <= 2 * m["res"]]
            if ok: found[tid] += 1; used |= set(ok); snr.append(max(acc[i][2] for i in ok))
        other.append(len(acc) - len(used))
    tg = " ".join("target %d %d/%d" % (k, found[k], n_t[k]) for k in sorted(n_t))
    print("%-34s maps %2d | %s | other detections per map %.1f | SNR of found %s" % (
        path, len(M), tg or "no truth", np.mean(other) if other else 0, ("%.1f dB" % np.median(snr)) if snr else "-"))
