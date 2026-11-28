"""P3: the association of CellSense (Kumar, Jain & Shah, arXiv 2606.07900, 2026): every detection
localised in closed form on its bistatic ellipse and AoA ray (their eq. 8; ours, track.localize),
then bipartite matching (Hungarian, Euclidean distance) to the tracks of a sliding repository of
the last `window` snapshots (theirs 10). No motion model or filter: a track is its latest matched
position. Unmatched detections start tracks; tracks not matched within the window are dropped.
Not given in the paper (assumed): the matching gate (gate_m) and which tracks are reported (those
matched in the last `report` scans). The velocity reported is the last displacement."""
import numpy as np
from scipy.optimize import linear_sum_assignment


def run(frames, window=10, gate_m=6.0, report=2):
    tracks, out, nid = [], [], 1
    for k, (t, meas, rb, sc, geo, rx) in enumerate(frames):
        pts = [m["xy"] for m in meas if m["xy"] is not None]
        tracks = [tr for tr in tracks if k - tr["k"] <= window]
        used = set()
        if tracks and pts:
            d = np.array([[np.hypot(*(tr["xy"] - p)) for p in pts] for tr in tracks])
            for i, j in zip(*linear_sum_assignment(d)):
                if d[i, j] <= gate_m:
                    tr = tracks[i]
                    tr["v"] = (pts[j] - tr["xy"]) / max(t - tr["t"], 1e-3)
                    tr["xy"], tr["t"], tr["k"] = pts[j], t, k
                    used.add(j)
        for j, p in enumerate(pts):
            if j not in used:
                tracks.append(dict(id=nid, xy=p, v=np.zeros(2), t=t, k=k))
                nid += 1
        out.append([(tr["id"], *map(float, tr["xy"]), *map(float, tr["v"])) for tr in tracks if k - tr["k"] < report])
    return out
