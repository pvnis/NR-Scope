"""Place the lidar on the site map by fitting its scan's walls to the OSM building outlines.

The pose is the lidar's origin (x, y, site frame: x east, y north) and its heading (compass
bearing of its +Y, azimuth 0, clockwise from north). A lidar point (X, Y) lands on the map at
origin + X * (cos h, -sin h) + Y * (sin h, cos h).

Steps: the scan's walls (cells of a top-view grid where points stack up over WALL_MIN_M
between the ground and 4 m: people, cars and the ground drop out), a coarse search over
shifts and turns around the starting guess (GNSS and compass) scored on a distance image of
the outlines, then a point-to-line ICP. The ICP's normal matrix says which directions the
walls fix: one straight facade fixes the distance to it and the heading, not the position
along it, and that direction is left at the guess and reported."""
import numpy as np
from scipy.ndimage import distance_transform_edt

CELL_M = 0.2          # top-view grid of the scan
WALL_MIN_M = 1.5      # a cell is a wall when its points span this much height
WALL_BAND_M = (0.3, 4.0)   # heights above the ground looked at
MAX_RANGE_M = 150.0   # scan points used (the OT128 reaches 200 m; the far ones are sparse)
N_PTS = 3000          # wall points kept for the fit
SEARCH_M, SEARCH_STEP_M = 10.0, 0.5
SEARCH_DEG, SEARCH_STEP_DEG = 15.0, 1.0
RASTER_M = 0.2        # distance image of the outlines
INLIER_M = 0.5        # a wall point "sits on" an outline within this


def wall_points(cloud, h, rng=None):
    """Wall cells of a scan, [n, 2] in the lidar frame. cloud: [n, >=3] x, y, z (lidar
    frame), h: the lidar's height above the ground (flat ground at z = -h)."""
    p = cloud[:, :3]
    z = p[:, 2] + h
    p = p[(z > WALL_BAND_M[0]) & (z < WALL_BAND_M[1]) & (np.hypot(p[:, 0], p[:, 1]) < MAX_RANGE_M)]
    if not len(p):
        return np.zeros((0, 2))
    ij = np.floor(p[:, :2] / CELL_M).astype(np.int64)
    _, inv, cnt = np.unique(ij, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    lo, hi = np.full(len(cnt), np.inf), np.full(len(cnt), -np.inf)
    np.minimum.at(lo, inv, p[:, 2])
    np.maximum.at(hi, inv, p[:, 2])
    keep = (hi - lo >= WALL_MIN_M) & (cnt >= 3)
    xy = np.zeros((len(cnt), 2))
    np.add.at(xy, inv, p[:, :2])
    xy = (xy / cnt[:, None])[keep]
    if len(xy) > N_PTS:
        xy = xy[(rng or np.random.default_rng(0)).choice(len(xy), N_PTS, replace=False)]
    return xy


def outline_segments(buildings, centre, reach):
    """The building outlines' edges within reach of centre: [n, 2, 2] (site frame)."""
    seg = []
    for b in buildings:
        f = np.asarray(b["footprint"], float)[:, :2]
        if len(f) < 2:
            continue
        e = np.stack([f, np.roll(f, -1, 0)], 1)
        if np.allclose(f[0], f[-1]):
            e = e[:-1]
        d = np.hypot(*(np.clip(centre, e.min(1), e.max(1)) - centre).T)   # bounding-box distance
        seg.append(e[d < reach])
    seg = np.concatenate(seg) if seg else np.zeros((0, 2, 2))
    return seg[np.hypot(*(seg[:, 1] - seg[:, 0]).T) > 1e-6]


def to_site(xy, pose):
    x, y, h = pose
    c, s = np.cos(np.radians(h)), np.sin(np.radians(h))
    return np.c_[x + xy[:, 0] * c + xy[:, 1] * s, y - xy[:, 0] * s + xy[:, 1] * c]


def nearest_on(seg, q):
    """For each point q, the nearest point on any segment and that segment's unit normal."""
    a, d = seg[:, 0], seg[:, 1] - seg[:, 0]
    t = np.clip(((q[:, None] - a) * d).sum(2) / (d * d).sum(1), 0, 1)
    c = a + t[..., None] * d
    dist = np.hypot(*(q[:, None] - c).transpose(2, 0, 1))
    k = dist.argmin(1)
    n = d[k][:, ::-1] * [1, -1]
    return c[np.arange(len(q)), k], n / np.hypot(*n.T)[:, None], dist[np.arange(len(q)), k]


def align(cloud, h, guess, buildings):
    """Fit the scan to the outlines. guess: (x, y, heading) from GNSS and compass. Returns a
    dict: pose (x, y, heading), shift (m) and turn (deg) from the guess, inliers (share of
    wall points within INLIER_M of an outline), rms_m, sigma (x, y m, heading deg),
    free_dir (a unit vector the walls do not fix, or None), ambiguous, n_walls, warnings."""
    W = wall_points(np.asarray(cloud, float), h)
    out = dict(n_walls=len(W), warnings=[])
    if len(W) < 50:
        return dict(out, ok=False, warnings=[f"only {len(W)} wall points in the scan: nothing to fit"])
    reach = MAX_RANGE_M + SEARCH_M
    seg = outline_segments(buildings, np.asarray(guess[:2], float), reach)
    if len(seg) < 2:
        return dict(out, ok=False, warnings=["no OSM building within reach of the guess"])

    # distance image of the outlines around the guess
    x0, y0 = guess[0] - reach, guess[1] - reach
    n = int(2 * reach / RASTER_M) + 1
    occ = np.ones((n, n), bool)
    for a, b in seg:
        k = int(np.hypot(*(b - a)) / (RASTER_M / 2)) + 2
        t = np.linspace(0, 1, k)[:, None]
        ij = np.round((a + t * (b - a) - [x0, y0]) / RASTER_M).astype(int)
        ij = ij[((ij >= 0) & (ij < n)).all(1)]
        occ[ij[:, 1], ij[:, 0]] = False
    D = distance_transform_edt(occ) * RASTER_M

    def score(P, tau):   # mean of 1 - d/tau over the points, 0 beyond tau or off the image
        ij = np.round((P - [x0, y0]) / RASTER_M).astype(int)
        ok = ((ij >= 0) & (ij < n)).all(-1)
        d = np.where(ok, D[np.clip(ij[..., 1], 0, n - 1), np.clip(ij[..., 0], 0, n - 1)], tau)
        return np.clip(1 - d / tau, 0, None).mean(-1)

    # coarse search: every turn, every shift
    sh = np.arange(-SEARCH_M, SEARCH_M + 1e-9, SEARCH_STEP_M)
    SX, SY = np.meshgrid(sh, sh)
    S = np.c_[SX.ravel(), SY.ravel()]
    best = []
    for dh in np.arange(-SEARCH_DEG, SEARCH_DEG + 1e-9, SEARCH_STEP_DEG):
        P = to_site(W, (guess[0], guess[1], guess[2] + dh))
        sc = np.concatenate([score(P[None] + S[i:i + 64, None], 1.0) for i in range(0, len(S), 64)])
        for i in np.argsort(-sc)[:5]:
            best.append((sc[i], guess[0] + S[i, 0], guess[1] + S[i, 1], guess[2] + dh))
    best.sort(key=lambda b: -b[0])
    pose = np.array(best[0][1:])
    rival = next((b for b in best[1:] if np.hypot(b[1] - pose[0], b[2] - pose[1]) > 2 or abs(b[3] - pose[2]) > 3), None)
    out["ambiguous"] = bool(rival and rival[0] > 0.9 * best[0][0])

    # point-to-line ICP (Gauss-Newton on x, y, heading), trimmed
    for it in range(30):
        P = to_site(W, pose)
        c, nrm, d = nearest_on(seg, P)
        keep = d < (1.0 if it < 10 else INLIER_M)
        if keep.sum() < 20:
            break
        r = ((P - c) * nrm).sum(1)[keep]
        hr = np.radians(pose[2])
        dP = np.c_[-W[:, 0] * np.sin(hr) + W[:, 1] * np.cos(hr), -W[:, 0] * np.cos(hr) - W[:, 1] * np.sin(hr)]   # dP/dheading (per rad)
        J = np.c_[nrm, (dP * nrm).sum(1)][keep]
        H = J.T @ J
        step = -np.linalg.lstsq(H, J.T @ r, rcond=1e-6)[0]
        pose = pose + [step[0], step[1], np.degrees(step[2])]
        if np.hypot(step[0], step[1]) < 1e-3 and abs(np.degrees(step[2])) < 1e-3:
            break

    P = to_site(W, pose)
    c, nrm, d = nearest_on(seg, P)
    keep = d < INLIER_M
    r = ((P - c) * nrm).sum(1)[keep]
    hr = np.radians(pose[2])
    dP = np.c_[-W[:, 0] * np.sin(hr) + W[:, 1] * np.cos(hr), -W[:, 0] * np.cos(hr) - W[:, 1] * np.sin(hr)]
    J = np.c_[nrm, (dP * nrm).sum(1)][keep]
    H = J.T @ J if keep.sum() >= 3 else np.eye(3) * 1e-9
    rms = float(np.sqrt(np.mean(r * r))) if len(r) else None

    # what the walls fix: the position block's weak direction (one straight facade: along it)
    ev, evec = np.linalg.eigh(H[:2, :2])
    free = None
    if ev[0] < 0.02 * max(ev[1], 1e-9):
        free = evec[:, 0]
        along = (np.asarray(guess[:2]) - pose[:2]) @ free   # keep the guess along that direction
        pose[:2] += along * free
        out["warnings"].append("the scan sees walls in one direction only: the position along them stays the GNSS guess")
    cov = np.linalg.pinv(H) * (max(rms or 0.1, 0.05) ** 2)
    sig = np.sqrt(np.clip(np.diag(cov), 0, None))
    if free is not None:
        sig[:2] = np.nan
    inl = float(keep.mean())
    if free is not None:
        out["ambiguous"] = False   # the free direction is the rival, already said
    if out["ambiguous"]:
        out["warnings"].append(f"another pose {rival[1] - pose[0]:+.1f} m, {rival[2] - pose[1]:+.1f} m, {rival[3] - pose[2]:+.1f}° fits almost as well: check the overlay")
    if inl < 0.3:
        out["warnings"].append(f"only {inl:.0%} of the wall points sit on an outline: wrong site, bad guess, or walls OSM lacks")
    return dict(out, ok=True, pose=[float(pose[0]), float(pose[1]), float(pose[2] % 360)],
                shift_m=float(np.hypot(*(pose[:2] - guess[:2]))), turn_deg=float((pose[2] - guess[2] + 180) % 360 - 180),
                inliers=inl, rms_m=rms, sigma=[None if np.isnan(v) else float(v) for v in (sig[0], sig[1], np.degrees(sig[2]))],
                free_dir=None if free is None else [float(free[0]), float(free[1])])
