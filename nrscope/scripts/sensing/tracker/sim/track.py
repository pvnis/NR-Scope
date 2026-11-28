"""Evaluation-side tracking: does NR-Scope's output follow the targets from map to map?

Each map's detections are grouped into one measurement per object (group_measurements),
localised with their angle of arrival, and tracked in x-y with an EKF on the bistatic
measurements (run_ekf_tracker), then scored against the ground truth (localize_and_score).

Extended objects. An object longer than a range bin gives several detections in one map:
a car's echo spreads over about 8 m of bistatic range (front and back, seen from the gNB
and from the receiver), so NR-Scope often reports it twice, about 2.7 bins apart at the same
speed. The detections of a map are grouped (group()): two detections belong to one object
when their speeds differ by at most GROUP_SPEED_CELLS speed cells and their ranges by at
most 2 * object_m (the bistatic extent of an object of length object_m), transitively. Each
group is one measurement, at its SNR-weighted centre.
"""
import numpy as np
from scipy.optimize import brentq, linear_sum_assignment
from scipy.stats import chi2
from pd_map import FOV_DEG, seen_at

CONFIRM = 3
MAX_MISS = 3
XY_DUP_M = 6.0          # a track within this distance of a target is near it (a duplicate, not a false track)
OBJECT_M = 4.5          # default object size: a car
SMALL_M = 1.0           # the size used instead for detections slower than SMALL_BELOW_MS:
SMALL_BELOW_MS = 2.5    # people and drones walk/hover slowly, vehicles don't (0 = off)
GROUP_SPEED_CELLS = 1.5
LIMB_SPREAD_MS = 1.5   # a walking person's legs and arms reach this far from its body speed (m/s); 0 = off (see limb_merge)
LIMB_ANGLE_DEG = 5.0
# the speed tolerances (grouping, gates, noise) are in speed cells, but never finer than this:
# an object's echoes spread over about +-0.5 m/s (a car's parts, a walker's limbs) whatever
# the map's resolution, so a finer map must not make one object look like several
SPEED_FLOOR_MS = 0.32


def group(dets, span, speed_tol, key=lambda d: (d[0], d[1]), weight=lambda d: 1.0, scale=None, scale_tol=False):
    """Connected groups of detections closer than span in range (first coordinate) and
    speed_tol in speed (second); returns each group's weighted centre. scale(di, dj), if
    given, multiplies span for that pair (and speed_tol if scale_tol)."""
    n = len(dets)
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            a, b = key(dets[i]), key(dets[j])
            f = scale(dets[i], dets[j]) if scale else 1.0
            if abs(a[0] - b[0]) <= span * f and abs(a[1] - b[1]) <= speed_tol * (f if scale_tol else 1.0):
                parent[root(j)] = root(i)
    out = {}
    for i in range(n):
        out.setdefault(root(i), []).append(dets[i])
    cent = []
    for g in out.values():
        w = np.array([weight(d) for d in g])
        cent.append(tuple(np.average(np.array([key(d) for d in g]), axis=0, weights=w)))
    return cent


def limb_merge(det, aoa, rb):
    """Micro-Doppler: a walking person is one body (torso) plus legs and arms swinging at other
    speeds, so it gives several detections at the same range and angle, spread in speed by up
    to about LIMB_SPREAD_MS (Chen's walking model). Merge them into one detection at the
    strongest one (the torso), before grouping and tracking: then the tracker sees one
    measurement per person at its body speed. det: [[e, v, power dB]]; aoa: [[e, v, angle deg]]
    (NR-Scope's AoA cells); rb: range bin (m). Only slow pairs (a person, not a car)."""
    if not LIMB_SPREAD_MS or len(det) < 2:
        return det
    ang = []
    for e, v, _ in det:   # each detection's angle: the AoA cell on it, if any
        c = [a for a in aoa if abs(a[0] - e) <= rb and abs(a[1] - v) <= 0.5]
        ang.append(min(c, key=lambda a: abs(a[0] - e) + abs(a[1] - v))[2] if c else None)
    order = sorted(range(len(det)), key=lambda i: -det[i][2])   # strongest first: the torsos
    taken, out = set(), []
    for i in order:
        if i in taken:
            continue
        taken.add(i)
        e, v, _ = det[i]
        if abs(v) < SMALL_BELOW_MS:
            for j in order:
                if j not in taken and abs(det[j][0] - e) <= rb and abs(det[j][1] - v) <= LIMB_SPREAD_MS \
                        and abs(det[j][1]) < SMALL_BELOW_MS + LIMB_SPREAD_MS \
                        and (ang[i] is None or ang[j] is None or abs(ang[i] - ang[j]) <= LIMB_ANGLE_DEG):
                    taken.add(j)
        out.append(det[i])
    return out


def size_scale(object_m, speed_of):
    """group()'s scale: SMALL_M / object_m for a pair of detections both slower than
    SMALL_BELOW_MS (a person, not a car), else 1."""
    if not SMALL_BELOW_MS or object_m <= SMALL_M:
        return None
    return lambda a, b: SMALL_M / object_m if max(abs(speed_of(a)), abs(speed_of(b))) < SMALL_BELOW_MS else 1.0


def group_measurements(rows, maps, object_m=OBJECT_M):
    """rows, maps: analyze.summary's rows and the maps they came from. Adds to every row
    r['groups'], the map's detections grouped into one measurement per object (the 2D
    tracker's input), and r['res'], its range bin and speed cell."""
    for r, m in zip(rows, maps):
        sc = max(abs(m["speed_step"]) * 2, SPEED_FLOOR_MS)
        z = group(limb_merge(r["det"], r.get("aoa", []), m["m_per_bin"]), 2 * object_m, GROUP_SPEED_CELLS * sc, weight=lambda d: 10 ** (d[2] / 10),
                  scale=size_scale(object_m, lambda d: d[1]))
        r["groups"] = [[float(e), float(v)] for e, v in z]
        r["res"] = (m["m_per_bin"], sc)


# ---- 2D positions (steps 4-5) ---------------------------------------------------------------

# The target's height is not known and not assumed: anywhere between the ground and 3 m.
# Heights are above the ground (the radios' z as set in the run, the ground at 0).
TARGET_Z_M = (0.0, 3.0)
TARGET_Z_MID = sum(TARGET_Z_M) / 2
TARGET_Z_VAR = (TARGET_Z_M[1] - TARGET_Z_M[0]) ** 2 / 12   # uniform over the band


def as3(p):
    p = [float(v) for v in p]
    return np.array(p + [0.0] * (3 - len(p)))


def localize_at(excess, angle_deg, gnb, rx, boresight_deg, z):
    """The target point at height z on both the bistatic ellipsoid (foci gnb, rx in 3D, total
    path |rx-gnb| + excess) and the cone NR-Scope's linear array measures (sin angle =
    u . axis, the axis horizontal), in front of the array. None when there is none."""
    T, R = as3(gnb), as3(rx)
    b = np.deg2rad(boresight_deg)
    n, ax = np.array([np.cos(b), np.sin(b)]), np.array([-np.sin(b), np.cos(b)])
    s = np.sin(np.deg2rad(angle_deg))
    D = np.linalg.norm(R - T) + excess
    dz = z - R[2]
    if excess <= 0 or abs(s) >= 1:
        return None

    def point(rho):   # range rho from rx: horizontal part along the axis s*rho, the rest along n
        h2 = rho * rho - dz * dz
        a = s * rho
        return np.r_[R[:2] + a * ax + np.sqrt(max(h2 - a * a, 0.0)) * n, z]

    f = lambda rho: np.linalg.norm(point(rho) - T) + rho - D
    lo = abs(dz) / np.sqrt(1 - s * s) + 1e-6
    if f(lo) > 0 or f(D) < 0:
        return None
    return point(brentq(f, lo, D))


def localize(excess, angle_deg, gnb, rx, boresight_deg):
    """NR-Scope's bistatic localisation in 3D without a target height: the curve of points
    (localize_at) for heights over TARGET_Z_M. Returns (x, y, spread): the point at mid
    height and the horizontal distance between the ends of the curve (the position
    uncertainty the unknown height leaves), or None."""
    pts = [localize_at(excess, angle_deg, gnb, rx, boresight_deg, z) for z in (TARGET_Z_M[0], TARGET_Z_MID, TARGET_Z_M[1])]
    ok = [p for p in pts if p is not None]
    if not ok:
        return None
    mid = pts[1] if pts[1] is not None else ok[0]
    spread = float(np.hypot(*(ok[0][:2] - ok[-1][:2]))) if len(ok) > 1 else 0.0
    return float(mid[0]), float(mid[1]), spread


XY_GATE_M = 6.0      # association gate around a predicted position
XY_HIT_M = 3.0       # a track or position estimate "is on" a car within this distance
# Targets NR-Scope cannot tell apart (in one peak of the map): within res_bins range bins (EKF_DEFAULTS)
# and RES_SPEED_CELLS speed cells of each other, chained. Scored as one group (localize_and_score).
RES_SPEED_CELLS = 2.0
XY_R_M = 1.5         # position measurement noise (m)
XY_Q = 4.0           # acceleration noise, (m/s^2)^2 per second
XY_GOSPA_C = 5.0     # GOSPA cutoff, metres


class XYTrack:
    """Constant-velocity Kalman filter on [x, y, vx, vy]."""

    def __init__(self, tid, z, t):
        self.id, self.t = tid, t
        self.x = np.array([z[0], z[1], 0.0, 0.0])
        self.P = np.diag([XY_R_M ** 2, XY_R_M ** 2, 100.0, 100.0])  # speed unknown at birth
        self.hits, self.miss, self.confirmed = 1, 0, False

    def predict(self, t):
        dt = t - self.t
        F = np.eye(4)
        F[0, 2] = F[1, 3] = dt
        q = XY_Q * np.array([[dt ** 3 / 3, dt ** 2 / 2], [dt ** 2 / 2, dt]])
        Q = np.zeros((4, 4))
        Q[np.ix_([0, 2], [0, 2])] = q
        Q[np.ix_([1, 3], [1, 3])] = q
        self.x, self.P, self.t = F @ self.x, F @ self.P @ F.T + Q, t

    def update(self, z):
        H = np.eye(2, 4)
        S = H @ self.P @ H.T + np.eye(2) * XY_R_M ** 2
        K = self.P @ H.T @ np.linalg.inv(S)
        self.x = self.x + K @ (np.asarray(z) - H @ self.x)
        self.P = (np.eye(4) - K @ H) @ self.P
        self.hits += 1
        self.miss = 0
        self.confirmed |= self.hits >= CONFIRM


def run_xy_tracker(frames):
    """frames: [(t, [(x, y), ...])]; returns per frame the confirmed tracks [(id, x, y, vx, vy)]."""
    tracks, out, next_id = [], [], 1
    for t, pts in frames:
        for tr in tracks:
            tr.predict(t)
        pairs = []
        if tracks and pts:
            d = np.array([[np.hypot(tr.x[0] - p[0], tr.x[1] - p[1]) for p in pts] for tr in tracks])
            rr, cc = linear_sum_assignment(np.where(d < XY_GATE_M, d, 1e6))
            pairs = [(i, j) for i, j in zip(rr, cc) if d[i, j] < XY_GATE_M]
        hit, used = {i for i, _ in pairs}, {j for _, j in pairs}
        for i, j in pairs:
            tracks[i].update(pts[j])
        for i, tr in enumerate(tracks):
            if i not in hit:
                tr.miss += 1
        tracks = [tr for tr in tracks if tr.miss < MAX_MISS and (tr.confirmed or tr.miss == 0)]
        for j, p in enumerate(pts):
            if j not in used:
                tracks.append(XYTrack(next_id, p, t))
                next_id += 1
        out.append([(tr.id, *map(float, tr.x)) for tr in tracks if tr.confirmed])
    return out


# ---- tracker 3: EKF in the plane on the raw bistatic measurements ---------------------------

EKF_SIGMA_THETA_DEG = 3.0   # angle noise
EKF_R_RANGE_BINS = 0.7      # excess-range noise, in range bins (as tracker 1)
EKF_R_SPEED_CELLS = 0.7     # speed noise, in speed cells
EKF_GATE2 = {2: 13.8, 3: 16.3}  # chi^2 gates, 99.9 % for 2 and 3 measured components

# Settings of the 2D tracker (GUI "Tracking settings"); any key left out takes this default.
EKF_DEFAULTS = dict(
    model="cv",            # motion model: "cv" constant velocity, "ct" coordinated turn
    q_acc=0.35,            # acceleration noise, (m/s^2)^2 per second: walkers and road traffic go straight (4 made tracks hop between neighbours)
    q_turn=0.05,           # CT: turn-rate noise, (rad/s)^2 per second
    sigma_theta_deg=EKF_SIGMA_THETA_DEG,
    r_range_bins=EKF_R_RANGE_BINS,
    r_speed_cells=EKF_R_SPEED_CELLS,
    gate_p=0.999,          # probability that a true measurement falls inside the gate
    # Track management by existence probability r (the Bernoulli of PMB/LMB filters, updated as
    # in IPDA): a detection raises r by its likelihood against clutter, a miss lowers it by
    # the receiver's P_D at the track (pd_map: low behind a building, so a miss there costs little).
    p_d=0.7,               # detection probability where the receiver sees the target (line of sight, in its field of view)
    p_d_hidden=0.1,        # elsewhere (a reflection can still reach it)
    view_gate="soft",      # several receivers, a track out of a receiver's view: off; "hard": none of its detections
                           # update it; "soft": only one with a position within view_gate_m of it
    view_gate_m=5.0,
    p_s=0.9,               # survival probability per second
    clutter=2.0,           # false detections per map per receiver, spread over CLUTTER_VOL
    r_birth=0.1,           # existence of a track started on an unexplained detection (times P_D / p_d there)
    r_confirm=0.95,        # shown from this existence on
    confirm_s=0.0,         # and only once the track is this old (s): a slower start hides short-lived false tracks (0 = off)
    dup_birth_m=4.0,       # no new track on a detection within this distance (m) of a confirmed track: a second track on the
                           # same object (another receiver's offset view, a limb, a split peak) is not started (0 = off;
                           # 4 m chosen on BC seed 1, Type I: false tracks -30 %, coverage -2..-4 points, PAPER_STATE 7.21)
    r_delete=0.02,         # deleted below this
    eot=1,                 # slow tracks (people) are extended: a GGIW extent, several detections per map (0 = point tracks)
    eot_tau=2.0,           # extent memory (s): how fast the extent forgets old detections
    merge_m=8.0,           # two vehicle tracks this close with the same velocity are one long vehicle (its front and back echoes): keep the older (0 = off)
    cross_m=2.5,           # crossing (m): two confirmed tracks this close, heading in opposite directions, coast straight
                           # (velocity kept, not deleted, no new track near them) until they separate (0 = off)
    res_bins=2.0,          # scoring only: targets this many range bins apart at the same speed are one peak (a group)
    ghost_map=1,           # multipath ghosts: no track starts on a wall-bounce echo predicted from walls (see ghost_echoes):
                           # 1 = the site map's walls, 2 = walls learned from this run's own tracks (wall_learn.py)
    ghost_weaker=1,        # 1: a predicted ghost blocks a detection only if its source is at least as strong (a bounce
                           # loses energy); 0: on position and speed alone (removes more ghosts, but can block people
                           # entering behind a strong false echo)
    replica_ms=0.0,        # Doppler replicas: the DM-RS sample slow time on the TDD pattern, so every peak repeats every
                           # 1/period in Doppler (5 ms: 200 Hz = 8 m/s at 3.75 GHz). A detection at the same range as a stronger
                           # one, k = 1..3 replica spacings away in speed, is dropped (0 = off)
    assoc="gnn",           # association: "gnn" one best pairing per map (Hungarian); "bp" soft, by belief propagation
                           # (Williams & Lau 2014, Meyer et al. 2018): each track updates with all its gated detections,
                           # weighted by their association probabilities, and a detection starts a track as far as it is nobody's
)
BP_ITERS = 30              # belief propagation iterations (converges in a few)
GROUP_M = 1.5              # bp: slow detections this close (same speed) are one extended target's (a person's limbs): one measurement
SECOND_ANGLE_DEG = 20.0   # angles of one detection this far apart are two targets (NR-Scope's second arrival)
GHOST_WALL_M = 40.0        # walls closer than this to a target can give it a ghost
GHOST_SEP_M = 5.0          # a predicted ghost less than this beyond its source lands on the source itself: ignored
GHOST_DEG = 8.0            # angle gate of a predicted ghost
CLUTTER_VOL = (100.0, 6.0, np.deg2rad(120))   # where the false detections fall: excess range (m), speed (m/s; near zero, like people), angle (rad)


EOT_X0_M = 0.3        # extent at birth: a person, a disc of about this standard deviation (m)
EOT_X_MAX_M = 1.5     # nor larger than a small group (semi-axis, 1 sigma): a wider one swallows its neighbours
EOT_ALPHA0 = 5.0      # its weight at birth, in detections (inverse-Wishart degrees of freedom - 6)
EOT_ALPHA_MIN = 2.0   # the extent never forgets below this weight
EOT_LIMB_MS = 1.0     # a person's legs and arms spread its speed: their echoes are its own (m/s, 1 sigma)
EOT_Z = 0.25          # detections spread over an ellipse: their covariance is EOT_Z x the extent's


def sqrtm2(A):
    """Symmetric square root of a 2x2 symmetric positive (semi)definite matrix."""
    w, V = np.linalg.eigh((A + A.T) / 2)
    return V @ np.diag(np.sqrt(np.maximum(w, 1e-12))) @ V.T


def ekf_cfg(cfg=None):
    c = dict(EKF_DEFAULTS)
    c.update({k: v for k, v in (cfg or {}).items() if k in c and v not in (None, "")})
    for k in c:
        if k not in ("model", "assoc", "view_gate"):
            c[k] = float(c[k])
    c["eot"] = int(c["eot"])
    c["ghost_map"] = int(c["ghost_map"])
    c["ghost_weaker"] = int(c["ghost_weaker"])
    c["model"] = "ct" if str(c["model"]).lower() == "ct" else "cv"
    c["assoc"] = "bp" if str(c["assoc"]).lower() == "bp" else "gnn"
    return c


def bp_assoc(beta, beta0, iters=BP_ITERS):
    """Association probabilities by belief propagation (Williams & Lau 2014). beta[i, j]: the
    likelihood ratio that track i made detection j (r pd g / kappa, 0 outside the gate);
    beta0[i]: that it made none (1 - r pd). Returns p[i, j], p0[i] (track i missed) and q0[j]
    (detection j is nobody's: clutter or a new target)."""
    nu = np.ones_like(beta)
    for _ in range(iters):
        bn = beta * nu
        mu = beta / (beta0[:, None] + bn.sum(1, keepdims=True) - bn)   # track -> detection
        nu = 1 / (1 + mu.sum(0, keepdims=True) - mu)                    # detection -> track
    bn = beta * nu
    tot = beta0 + bn.sum(1)
    return bn / tot[:, None], beta0 / tot, 1 / (1 + mu.sum(0))


def group_slow(meas):
    """bp: slow detections with positions closer than GROUP_M, at the same speed, merged into one
    measurement (mean range, speed, angle) that keeps their positions (pts) for the extent."""
    out, done = [], set()
    for j, m in enumerate(meas):
        if j in done:
            continue
        g = [j]
        if m["xy"] is not None and m["theta"] is not None and abs(m["v"]) < SMALL_BELOW_MS:
            for k in range(j + 1, len(meas)):
                n = meas[k]
                if k not in done and n["xy"] is not None and n["theta"] is not None and abs(n["v"] - m["v"]) <= EOT_LIMB_MS \
                        and min(np.hypot(*(np.asarray(meas[i]["xy"]) - n["xy"])) for i in g) <= GROUP_M:
                    g.append(k)
        done |= set(g)
        if len(g) == 1:
            out.append(dict(m, pts=[m["xy"]] if m["xy"] is not None else []))
            continue
        M = [meas[i] for i in g]
        p = [x["p_db"] for x in M if x["p_db"] is not None]
        out.append(dict(m, e=float(np.mean([x["e"] for x in M])), v=float(np.mean([x["v"] for x in M])),
                        theta=float(np.mean([x["theta"] for x in M])), xy=np.mean([x["xy"] for x in M], 0),
                        p_db=max(p) if p else None, pts=[x["xy"] for x in M]))
    return out


def ct_step(x, dt):
    """Coordinated turn: [px, py, vx, vy, w] moves dt seconds on a circle of turn rate w
    (the velocity rotates by w*dt); w -> 0 is constant velocity."""
    px, py, vx, vy, w = x
    wt = w * dt
    if abs(wt) < 1e-6:
        a, b = dt, 0.5 * w * dt * dt   # sin(wt)/w, (1-cos(wt))/w to first order
    else:
        a, b = np.sin(wt) / w, (1 - np.cos(wt)) / w
    c, s = np.cos(wt), np.sin(wt)
    return np.array([px + a * vx - b * vy, py + b * vx + a * vy, c * vx - s * vy, s * vx + c * vy, w])


def bistatic_meas(x, T, R, L, bore, z):
    """(e, v, theta) NR-Scope would measure from a target at (px, py, z) moving horizontally
    at (vx, vy): T, R the gNB and NR-Scope in 3D, L = |R - T|. e = |p-T| + |p-R| - L;
    v = -(1/2) (u_T + u_R) . vel (Doppler * lambda / 2); theta = the angle the linear array
    measures, asin(u_R . axis) (a cone: an echo from above or below reads nearer the normal)."""
    p = np.array([x[0], x[1], z])
    vel = np.array([x[2], x[3], 0.0])
    dT, dR = p - T, p - R
    nT, nR = max(np.linalg.norm(dT), 1e-6), max(np.linalg.norm(dR), 1e-6)
    ax = np.array([-np.sin(bore), np.cos(bore), 0.0])
    return np.array([nT + nR - L, -0.5 * (dT / nT + dR / nR) @ vel, np.arcsin(np.clip(dR @ ax / nR, -1, 1))])


def bistatic_h(x, T, R, L, bore):
    """Predicted (e, v, theta) of state x = [px, py, vx, vy] at mid height, its Jacobian on
    x, and the measurement covariance the unknown height adds (J_z J_z^T var_z): large near
    the radios, where the height matters, small far away. Numerical derivatives."""
    x = np.asarray(x[:4], float)
    h0 = bistatic_meas(x, T, R, L, bore, TARGET_Z_MID)
    H = np.empty((3, 4))
    for k in range(4):
        d = np.zeros(4)
        d[k] = 1e-4 * max(1.0, abs(x[k]))
        H[:, k] = (bistatic_meas(x + d, T, R, L, bore, TARGET_Z_MID) - bistatic_meas(x - d, T, R, L, bore, TARGET_Z_MID)) / (2 * d[k])
    jz = (bistatic_meas(x, T, R, L, bore, TARGET_Z_MID + 0.01) - bistatic_meas(x, T, R, L, bore, TARGET_Z_MID - 0.01)) / 0.02
    return h0, H, np.outer(jz, jz) * TARGET_Z_VAR


class EKFTrack(XYTrack):
    """Track in the plane, CV ([px, py, vx, vy]) or CT ([px, py, vx, vy, w]), updated with
    the bistatic measurements."""

    def __init__(self, tid, xy, vel, t, cfg):
        super().__init__(tid, xy, t)
        self.x[2:] = vel
        self.cfg = cfg
        self.p_db = None   # its latest echo's SNR (dB)
        self.r = cfg["r_birth"]   # existence probability
        self.t_birth = t
        self.X, self.alpha = np.eye(2) * EOT_X0_M ** 2, EOT_ALPHA0   # GGIW extent (see extended())
        if cfg["model"] == "ct":   # turn rate unknown at birth, a few degrees per second
            self.x = np.r_[self.x, 0.0]
            self.P = np.pad(self.P, ((0, 1), (0, 1)))
            self.P[4, 4] = 0.3 ** 2

    def extended(self):
        """A slow track (a person, or a group walking together) is an extended object: the
        random-matrix model of GGIW (Koch 2008, Feldmann et al. 2011). Its extent X, the 2x2
        covariance of where its detections come from, widens its gate and enters the update
        (a detection anywhere on the group is expected), and it can take several detections
        in one map. Vehicles stay point tracks (their front and back: merge_m)."""
        return self.cfg["eot"] and np.hypot(self.x[2], self.x[3]) < SMALL_BELOW_MS

    def predict(self, t):
        dt = t - self.t
        self.alpha = max(EOT_ALPHA_MIN, self.alpha * np.exp(-dt / self.cfg["eot_tau"]))
        q = self.cfg["q_acc"] * np.array([[dt ** 3 / 3, dt ** 2 / 2], [dt ** 2 / 2, dt]])
        n = len(self.x)
        Q = np.zeros((n, n))
        Q[np.ix_([0, 2], [0, 2])] = q
        Q[np.ix_([1, 3], [1, 3])] = q
        if n == 4:
            F = np.eye(4)
            F[0, 2] = F[1, 3] = dt
            x = F @ self.x
        else:   # CT: nonlinear, EKF prediction with the Jacobian of ct_step (numerical)
            x = ct_step(self.x, dt)
            F = np.empty((5, 5))
            for k in range(5):
                h = 1e-5 * max(1.0, abs(self.x[k]))
                d = np.zeros(5)
                d[k] = h
                F[:, k] = (ct_step(self.x + d, dt) - ct_step(self.x - d, dt)) / (2 * h)
            Q[4, 4] = self.cfg["q_turn"] * dt
        self.x, self.P, self.t = x, F @ self.P @ F.T + Q, t

    def innovation(self, z, Rz, geo):
        """z: [e, v] or [e, v, theta]; returns (nu, S, H) on the measured components."""
        zh, H4, Rh = bistatic_h(self.x[:4], *geo)
        k = len(z)
        nu = np.asarray(z) - zh[:k]
        if k == 3:
            nu[2] = (nu[2] + np.pi) % (2 * np.pi) - np.pi
        H = np.zeros((k, len(self.x)))
        H[:, :4] = H4[:k]   # the measurements do not depend on the turn rate
        S = H @ self.P @ H.T + Rz[:k, :k] + Rh[:k, :k]
        if self.extended():   # the detection comes from anywhere on the extent, and from any limb
            S = S + EOT_Z * H4[:k, :2] @ self.X @ H4[:k, :2].T
            S[1, 1] += EOT_LIMB_MS ** 2
        return nu, S, H

    @staticmethod
    def rxy(H4, Rz, Rh):
        """Position noise (x-y) of a localised detection: the (e, theta) noise through the
        inverse of their Jacobian on x-y."""
        J = H4[[0, 2], :2]
        Ji = np.linalg.inv(J) if abs(np.linalg.det(J)) > 1e-9 else np.linalg.pinv(J)
        Rs = Rz[np.ix_([0, 2], [0, 2])] + Rh[np.ix_([0, 2], [0, 2])]
        return Ji @ Rs @ Ji.T

    def extent_update(self, pts, Rxy):
        """Random-matrix update of the extent with this map's detections pts ([n, 2], x-y) of
        position noise Rxy (Feldmann et al. 2011), before the kinematic update."""
        pts = np.asarray(pts, float)
        n = len(pts)
        zbar = pts.mean(0)
        eps = zbar - self.x[:2]
        Y = EOT_Z * self.X + Rxy
        S = self.P[:2, :2] + Y / n
        Xs, Si, Yi = sqrtm2(self.X), np.linalg.inv(sqrtm2(S)), np.linalg.inv(sqrtm2(Y))
        N = Xs @ Si @ np.outer(eps, eps) @ Si.T @ Xs.T
        d = pts - zbar
        Z = Xs @ Yi @ (d.T @ d) @ Yi.T @ Xs.T
        X = (self.alpha * self.X + N + Z) / (self.alpha + n)
        w, V = np.linalg.eigh(X)   # never smaller than one person: with one detection a map it would shrink to nothing
        self.X = V @ np.diag(np.clip(w, EOT_X0_M ** 2, EOT_X_MAX_M ** 2)) @ V.T
        self.alpha += n

    def axes(self):
        """The extent's semi-axes (1 sigma, m), largest first, and the major axis' direction (rad)."""
        w, V = np.linalg.eigh(self.X)
        return float(np.sqrt(max(w[1], 0))), float(np.sqrt(max(w[0], 0))), float(np.arctan2(V[1, 1], V[0, 1]))

    def ekf_update(self, nu, S, H):
        K = self.P @ H.T @ np.linalg.inv(S)
        self.x = self.x + K @ nu
        self.P = (np.eye(len(self.x)) - K @ H) @ self.P

    def exist(self, pd, g=None, kappa=None):
        """Existence update: detected with likelihood g against clutter density kappa (IPDA,
        one measurement: it is the target's, or clutter while the target is missed or absent),
        or missed (g None)."""
        r = self.r
        if g is None:
            self.r = r * (1 - pd) / (1 - r * pd)
        else:
            self.r = (r * pd * g + kappa * r * (1 - pd)) / (r * pd * g + kappa * (1 - r * pd))
        self.confirmed |= self.r >= self.cfg["r_confirm"]


def ghost_echoes(x, T, R, L, bore, walls):
    """Map-based multipath ghost prediction (image method, as in automotive radar, e.g. Kamann
    et al. 2018): the (e, v, theta) of the echoes of a target at state x that also bounce once
    on a wall of the site map, gNB -> wall -> target -> NR-Scope and gNB -> target -> wall ->
    NR-Scope (blockage and wall heights ignored). Only ghosts at least GHOST_SEP_M beyond the
    target's own echo: closer ones can't be told from it."""
    from nonideal import bounce_points
    p = np.array([x[0], x[1], TARGET_Z_MID])
    vel = np.array([x[2], x[3], 0.0])
    unit = lambda d: d / max(np.linalg.norm(d), 1e-6)
    ax = np.array([-np.sin(bore), np.cos(bore), 0.0])
    uR, uT = unit(p - R), unit(p - T)
    e0 = np.linalg.norm(p - T) + np.linalg.norm(p - R) - L
    out = []
    for P, L1 in bounce_points(walls, T, p, GHOST_WALL_M):   # on the way in: arrives from the target
        out.append((L1 + np.linalg.norm(p - R) - L, -0.5 * (unit(p - P) + uR) @ vel, np.arcsin(np.clip(uR @ ax, -1, 1))))
    for P, L2 in bounce_points(walls, R, p, GHOST_WALL_M):   # on the way back: arrives from the wall
        out.append((np.linalg.norm(p - T) + L2 - L, -0.5 * (uT + unit(p - P)) @ vel, np.arcsin(np.clip(unit(P - R) @ ax, -1, 1))))
    return [g for g in out if g[0] - e0 > GHOST_SEP_M]


def gauss(nu, S):
    return float(np.exp(-0.5 * nu @ np.linalg.solve(S, nu)) / np.sqrt(np.linalg.det(2 * np.pi * S)))


def run_ekf_tracker(frames, geo, cfg=None, walls=None, pd=None, buildings=()):
    """frames: [(t, [dict(e, v, theta or None, xy or None)], range_bin, speed_cell[, geo, rx])].
    Every measurement can update a track (with or without its angle); only one with an
    angle, hence a position, and inside no track's gate can start a track. Association: Mahalanobis distance of the
    innovation, chi^2 gate, Hungarian. Returns per frame [(id, x, y, vx, vy)] of confirmed tracks.
    Several receivers (detection-level fusion): their frames in time order, each with its own
    geo (T, R, L, bore) and receiver number; every track takes each receiver's measurements
    through that receiver's bistatic model, so two ranges cross where the target is.
    pd: per receiver (1, 2, ...) its pd_map grid, or None: P_D from the field of view alone.
    buildings: footprints [[x, y], ...]: no track starts inside one (a wall ghost lands behind the wall)."""
    from scenario import inside
    cfg = ekf_cfg(cfg)
    gate2 = {k: float(chi2.ppf(cfg["gate_p"], k)) for k in (2, 3)}
    kappa = {k: cfg["clutter"] / np.prod(CLUTTER_VOL[:k]) for k in (2, 3)}
    tracks, out, next_id = [], [], 1
    multi = len({fr[5] for fr in frames if len(fr) > 5}) > 1   # detection-level fusion of several receivers
    for fr in frames:
        t, meas, rb, sc = fr[:4]
        geo_f = fr[4] if len(fr) > 4 else geo
        T, R, L, bore = geo_f
        grid = pd[fr[5] - 1] if pd and len(fr) > 5 else None

        def p_d(x):
            if grid is not None:
                seen = seen_at(grid, x[0], x[1])
            else:
                a = np.arctan2(x[1] - R[1], x[0] - R[0]) - bore
                seen = abs((a + np.pi) % (2 * np.pi) - np.pi) <= np.deg2rad(FOV_DEG)
            return cfg["p_d"] if seen else cfg["p_d_hidden"]
        Rz = np.diag([(cfg["r_range_bins"] * rb) ** 2, (cfg["r_speed_cells"] * sc) ** 2,
                      np.deg2rad(cfg["sigma_theta_deg"]) ** 2])
        for tr in tracks:
            tr.r *= cfg["p_s"] ** (t - tr.t)
            tr.predict(t)
        # crossing: confirmed tracks close together, heading in opposite directions
        crossing = set()
        if cfg["cross_m"] > 0:
            ct = [tr for tr in tracks if tr.confirmed]
            for a, ta in enumerate(ct):
                for tb in ct[a + 1:]:
                    if np.hypot(*(ta.x[:2] - tb.x[:2])) <= cfg["cross_m"] and ta.x[2:4] @ tb.x[2:4] < 0:
                        crossing |= {id(ta), id(tb)}
        bp = cfg["assoc"] == "bp"
        if bp:
            meas = group_slow(meas)
        zs = [[m["e"], m["v"]] + ([] if m["theta"] is None else [m["theta"]]) for m in meas]
        pairs, cache = [], {}
        if tracks and meas:
            cost = np.full((len(tracks), len(meas)), 1e6)
            for i, tr in enumerate(tracks):
                hidden = multi and cfg["view_gate"] in ("hard", "soft") and p_d(tr.x) < cfg["p_d"]
                if hidden and cfg["view_gate"] == "hard":
                    continue   # out of this receiver's view: its detections (other targets, clutter) would drag the track
                for j, z in enumerate(zs):
                    if hidden and (meas[j].get("xy") is None or np.hypot(*(np.asarray(meas[j]["xy"]) - tr.x[:2])) > cfg["view_gate_m"]):
                        continue   # soft: from a receiver that should not see it, only a detection placed on the track
                    nu, S, H = tr.innovation(z, Rz, geo_f)
                    d2 = float(nu @ np.linalg.solve(S, nu))
                    if d2 < gate2[len(z)]:
                        cost[i, j] = d2
                        cache[i, j] = (nu, S, H)
            rr, cc = linear_sum_assignment(cost)
            pairs = [(i, j) for i, j in zip(rr, cc) if (i, j) in cache]
        q0 = np.ones(len(meas))   # bp: how far each detection is nobody's
        if bp:
            pdv = np.array([p_d(tr.x) for tr in tracks])
            rv = np.array([tr.r for tr in tracks])
            beta = np.zeros((len(tracks), len(meas)))
            for (i, j), (nu, S, _) in cache.items():
                beta[i, j] = rv[i] * pdv[i] * gauss(nu, S) / kappa[len(zs[j])]
            pij, p0, q0 = bp_assoc(beta, 1 - rv * pdv) if len(meas) else (beta, np.ones(len(tracks)), q0)
            for i, tr in enumerate(tracks):
                js = [j for j in range(len(meas)) if (i, j) in cache]
                miss = tr.r * (1 - pdv[i]) / (1 - tr.r * pdv[i])   # it exists, given no detection of it
                if js:   # given it exists: the mixture of "missed" and each detection, moment-matched (PDA)
                    hyp = [(p0[i] * miss, tr.x, tr.P)]
                    for j in js:
                        nu, S, H = cache[i, j]
                        K = tr.P @ H.T @ np.linalg.inv(S)
                        hyp.append((pij[i, j], tr.x + K @ nu, (np.eye(len(tr.x)) - K @ H) @ tr.P))
                    w = np.array([h[0] for h in hyp]) / sum(h[0] for h in hyp)
                    x = sum(wk * h[1] for wk, h in zip(w, hyp))
                    v0 = tr.x[2:4].copy()
                    tr.P = sum(wk * (h[2] + np.outer(h[1] - x, h[1] - x)) for wk, h in zip(w, hyp))
                    tr.x = x
                    if id(tr) in crossing:   # coast straight
                        tr.x[2:4] = v0
                    jb = max(js, key=lambda j: pij[i, j])
                    if pij[i, jb] >= 0.5:
                        if meas[jb]["p_db"] is not None:
                            tr.p_db = meas[jb]["p_db"]
                        if tr.extended() and meas[jb]["pts"]:
                            _, H4, Rh = bistatic_h(tr.x[:4], T, R, L, bore)
                            tr.extent_update(meas[jb]["pts"], tr.rxy(H4, Rz, Rh))
                if js or id(tr) not in crossing:   # existence: detected with probability 1 - p0, else missed
                    tr.r = float(pij[i].sum() + p0[i] * miss)
                    tr.confirmed |= tr.r >= cfg["r_confirm"]
            used = {j for j in range(len(meas)) if q0[j] < 0.5}
        else:
            used = {j for _, j in pairs}
            # extended tracks also take the other detections that fall on them (a group: several)
            own = {i: [j] for i, j in pairs}
            lik = {i: (gauss(*cache[i, j][:2]), len(zs[j])) for i, j in pairs}   # each track's detection likelihood
            ext = [i for i, tr in enumerate(tracks) if tr.extended() and tr.confirmed]
            for j, m in enumerate(meas):
                if j in used or m["xy"] is None or not ext:
                    continue
                best = None
                for i in ext:   # the full (range, speed, angle) gate with the extent: a car passing a group has another speed
                    nu, S, _ = tracks[i].innovation(zs[j], Rz, geo_f)
                    d2 = float(nu @ np.linalg.solve(S, nu))
                    if d2 < gate2[len(zs[j])] and (best is None or d2 < best[0]):
                        best = (d2, i, gauss(nu, S))
                if best:
                    lik[best[1]] = max(lik.get(best[1], (0.0, 2)), (best[2], len(zs[j])))
                    own.setdefault(best[1], []).append(j)
                    used.add(j)
            for i, js in own.items():
                tr = tracks[i]
                pts = [meas[j]["xy"] for j in js if meas[j]["xy"] is not None]
                if tr.extended() and pts:
                    _, H4, Rh = bistatic_h(tr.x[:4], T, R, L, bore)
                    tr.extent_update(pts, tr.rxy(H4, Rz, Rh))
            for i, j in pairs:
                if meas[j]["p_db"] is not None:
                    tracks[i].p_db = meas[j]["p_db"]   # its echo's strength, for the ghost test
                v0 = tracks[i].x[2:4].copy()
                tracks[i].ekf_update(*cache[i, j])
                if id(tracks[i]) in crossing:   # coast straight: the detection corrects only the position
                    tracks[i].x[2:4] = v0
            for i, tr in enumerate(tracks):
                if i in lik:
                    tr.exist(p_d(tr.x), lik[i][0], kappa[lik[i][1]])
                elif id(tr) not in crossing:
                    tr.exist(p_d(tr.x))
        tracks = [tr for tr in tracks if tr.r >= cfg["r_delete"]]
        gated = set(used) if bp else {j for _, j in cache} | used  # explained by a track, so no new track
        if cfg["merge_m"] > 0:
            fast = sorted((tr for tr in tracks if np.hypot(tr.x[2], tr.x[3]) >= SMALL_BELOW_MS > 0), key=lambda tr: tr.id)
            drop = set()
            for a, ta in enumerate(fast):
                for tb in fast[a + 1:]:
                    if id(ta) not in drop and np.hypot(*(ta.x[:2] - tb.x[:2])) <= cfg["merge_m"] \
                            and np.hypot(*(ta.x[2:4] - tb.x[2:4])) <= 0.3 * np.hypot(ta.x[2], ta.x[3]):
                        wa, wb = 1 / np.trace(ta.P[:2, :2]), 1 / np.trace(tb.P[:2, :2])   # the older keeps its id, the better-known one's position counts more
                        ta.x[:4] = (wa * ta.x[:4] + wb * tb.x[:4]) / (wa + wb)
                        drop.add(id(tb))
            tracks = [tr for tr in tracks if id(tr) not in drop]
        if walls is not None and cfg["ghost_map"]:
            # multipath ghosts: no track starts on a predicted wall-bounce echo of any track (even a
            # new one) or of another detection of this map (a person and its ghost appear together)
            src = [(tr.x, None, tr.p_db) for tr in tracks]
            src += [(np.r_[m["xy"], 0.0, 0.0], m["v"], m["p_db"]) for m in meas if m["xy"] is not None and abs(m["v"]) > 0.3]
            G = [(g, g[1] if v0 is None else v0, p0) for x, v0, p0 in src for g in ghost_echoes(x, T, R, L, bore, walls)]
            weaker = lambda m, p0: not cfg["ghost_weaker"] or p0 is None or m["p_db"] is None or m["p_db"] <= p0
            for j, m in enumerate(meas):
                if m["theta"] is not None and any(abs(m["e"] - g[0]) <= 2 * rb and abs(m["v"] - gv) <= max(2 * sc, EOT_LIMB_MS)
                                                  and abs(m["theta"] - g[2]) <= np.deg2rad(GHOST_DEG) and weaker(m, p0)
                                                  for g, gv, p0 in G):
                    gated.add(j)
        for j, m in enumerate(meas):   # no new track inside a crossing: its detections belong to the two
            if m["xy"] is not None and any(id(tr) in crossing and np.hypot(*(tr.x[:2] - m["xy"])) <= cfg["cross_m"] for tr in tracks):
                gated.add(j)
        for j, m in enumerate(meas):
            if cfg.get("dup_birth_m", 0.0) > 0 and m["xy"] is not None and any(
                    tr.confirmed and np.hypot(*(tr.x[:2] - m["xy"])) <= cfg["dup_birth_m"] for tr in tracks):
                continue
            if j not in used and j not in gated and m["xy"] is not None and not any(inside(m["xy"], b) for b in buildings):
                # start at the localised position, the velocity along the bistatic bisector that
                # gives the measured speed (the component that is observed), the rest unknown
                x0 = np.r_[m["xy"], 0.0, 0.0]
                _, H, _ = bistatic_h(x0, T, R, L, bore)
                s = -2 * H[1, 2:]  # u_T + u_R
                vel = -2 * m["v"] * s / max(s @ s, 1e-6)
                tracks.append(EKFTrack(next_id, m["xy"], vel, t, cfg))
                tracks[-1].p_db = m["p_db"]
                tracks[-1].r *= p_d(m["xy"]) / cfg["p_d"] * q0[j]   # births where this receiver sees directly: a ghost lands behind walls
                next_id += 1
        out.append([(tr.id, *map(float, tr.x[:4]), *(tr.axes() if tr.extended() else (0.0, 0.0, 0.0)))
                    for tr in tracks if tr.confirmed and t - tr.t_birth >= cfg.get("confirm_s", 0.0)])
    return out


REPLICA_LAMBDA_M = 299792458.0 / 3.75e9


def drop_replicas(r, period_ms, tol_ms=1.0):
    """Remove Doppler replicas from one map row in place: a detection (and its AoA cells and groups) whose range is
    within one bin of a stronger detection and whose speed differs from it by k = 1..3 replica spacings
    (lambda / (2 period)), within tol_ms."""
    dets = r.get("det", [])
    if not dets or period_ms <= 0:
        return 0
    dv = REPLICA_LAMBDA_M / (2 * period_ms * 1e-3)
    rb = r.get("res", (2.44, 0.3))[0]
    def is_rep(e, v, p):
        for e2, v2, p2 in dets:
            if p2 > p and abs(e - e2) <= rb and any(abs(abs(v - v2) - k * dv) <= tol_ms for k in (1, 2, 3)):
                return True
        return False
    def power(e, v):
        pw = [p for e2, v2, p in dets if abs(e - e2) <= rb and abs(v - v2) <= tol_ms]
        return max(pw) if pw else -1e9
    keep = [d for d in dets if not is_rep(*d)]
    n = len(dets) - len(keep)
    r["det"] = keep
    if r.get("aoa"):
        r["aoa"] = [a for a in r["aoa"] if not is_rep(a[0], a[1], power(a[0], a[1]))]
    if r.get("groups"):
        r["groups"] = [g for g in r["groups"] if not is_rep(g[0], g[1], power(g[0], g[1]))]
    return n


def localize_and_score(rows, geo, object_m=OBJECT_M, tracker="ekf", cfg=None, half_len=None, walls=None, more=(), pd=None, buildings=()):
    """Steps 4-5. geo = (gnb, rx, boresight_deg), gnb and rx in 3D (heights above the ground). Localises every AoA cell of every row
    (r['pos'] = [[x, y, excess, speed, angle]]), scores the single-map positions against each
    car's true position (centre of the car), then tracks them in x-y and scores the tracks.
    Distances in metres. tracker: "ekf" (run_ekf_tracker, every detection group, with its
    angle when it has one) or "cv" (run_xy_tracker, positions only). cfg: the EKF's settings
    (EKF_DEFAULTS). half_len: {car: half its length}; a track is on a car within XY_HIT_M of
    its body, i.e. XY_HIT_M + half_len of its centre (the echo comes from its near side or
    its strongest part, not from its centre).
    more: other receivers' [(rows, geo)], for detection-level fusion: one tracker takes every
    receiver's maps in time order, each through its own geometry; the tracks are scored at
    every receiver's maps, a target counting as in view when either receiver has it in view.
    pd: per receiver (this one first, then more's) its pd_map grid; buildings: footprints (run_ekf_tracker)."""
    half_len, res_bins = half_len or {}, ekf_cfg(cfg)["res_bins"]
    sets = [(rows, geo)] + list(more)
    rep = ekf_cfg(cfg).get("replica_ms", 0.0)
    if rep > 0:
        for rs, _ in sets:
            for r in rs:
                drop_replicas(r, rep)
    order = sorted(((r["t"], k, r) for k, (rs, _) in enumerate(sets) for r in rs), key=lambda q: q[:2])
    frames, est_err, ekf_frames = [], {}, []
    for _, k, r in order:
        gnb, rx, bore = sets[k][1]
        T, R = as3(gnb), as3(rx)
        g3 = (T, R, float(np.linalg.norm(R - T)), np.deg2rad(bore))
        r["pos"] = []
        for e, v, ang in r.get("aoa", []):
            p = localize(e, ang, gnb, rx, bore)
            if p is not None:
                r["pos"].append([p[0], p[1], e, v, ang, p[2]])   # p[2]: spread from the unknown height
        # positions of one object: within its size of each other (no speed test: same object
        # ⇒ same place); usually one per object already, the AoA stage merges by angle
        frames.append((r["t"], group([(p[0], p[1], p[3]) for p in r["pos"]], object_m, object_m,
                                     key=lambda d: (d[0], d[1]), scale=size_scale(object_m, lambda d: d[2]), scale_tol=True)))
        # the EKF's measurements: tracker 1's detection groups, each with the angle of the AoA
        # cell nearest it (within 2 range bins and 2 speed cells), if any
        rb, sc = r.get("res", (2.44, 0.3))
        sc = max(sc, SPEED_FLOOR_MS)
        meas = []
        for e, v in r.get("groups", []):
            near = sorted((p for p in r["pos"] if abs(p[2] - e) <= 2 * rb and abs(p[3] - v) <= 2 * sc),
                          key=lambda p: abs(p[2] - e) / rb + abs(p[3] - v) / sc)
            # one measurement per direction: NR-Scope reports a second angle at the same range
            # and speed when a cell holds two arrivals (two people, or a person and a ghost)
            dirs = []
            for p in near:
                if all(abs(p[4] - q[4]) >= SECOND_ANGLE_DEG for q in dirs):
                    dirs.append(p)
            pw = [d[2] for d in r.get("det", []) if abs(d[0] - e) <= rb and abs(d[1] - v) <= 2 * sc]   # its detections' SNR (dB)
            for c in dirs or [None]:
                meas.append(dict(e=e, v=v, theta=None if c is None else np.deg2rad(c[4]),
                                 xy=None if c is None else np.array(c[:2]), p_db=max(pw) if pw else None))
        ekf_frames.append((r["t"], meas, rb, sc, g3, k + 1))
        # single-map localisation: each car in view takes its nearest estimate
        cars = r.get("cars") or [r["car"]]
        for g in r["targets"]:
            if not g["in_view"]:
                continue
            c = cars[g["car"] - 1]
            d = [np.hypot(p[0] - c[0], p[1] - c[1]) for p in r["pos"]]
            k = est_err.setdefault(g["car"], dict(maps=0, errs=[]))
            k["maps"] += 1
            if d and min(d) <= 2 * XY_GATE_M:
                k["errs"].append(min(d))
    confirmed = (run_ekf_tracker(ekf_frames, None, cfg, walls, pd, buildings) if tracker == "ekf" else run_xy_tracker(frames))
    if more:   # in view of any receiver: of this map's, or of the others' nearest maps
        for k, (rs, _) in enumerate(sets):
            t_k = np.array([r["t"] for r in rs])
            for k2, (rs2, _) in enumerate(sets):
                if k2 == k or not len(rs2):
                    continue
                t2 = np.array([r["t"] for r in rs2])
                for r in rs:
                    j = int(np.argmin(np.abs(t2 - r["t"])))
                    r["view_any"] = sorted(set(r.get("view_any", [])) | {g["car"] for g in rs2[j]["targets"] if g["in_view"]})
                    # the other receiver's own entry for a target this one does not list (its range
                    # and speed in that receiver's map, for the resolution groups below)
                    r["targets_any"] = {**{g["car"]: g for g in rs2[j]["targets"] if g["in_view"]}, **r.get("targets_any", {})}
    per_car, life, near, on_car, g_all, trk_err = {}, {}, {}, {}, [], {}
    for (_, _, r), tr in zip(order, confirmed):
        r["xy_tracks"] = [list(q) for q in tr]   # id, x, y, vx, vy[, extent semi-axes a, b (m), angle]
        r["xy_track_car"] = {}
        cars = r.get("cars") or [r["car"]]
        # in view of this receiver, or (fused) of any: a target only another receiver sees is not in
        # this receiver's list of targets, so the other receiver's list (view_any) adds it
        view = sorted(({g["car"] for g in r["targets"] if g["in_view"]} | set(r.get("view_any", ())))
                      & {k for k in range(1, len(cars) + 1) if np.isfinite(cars[k - 1][0])})
        for i, x, y, *_ in tr:
            life[i] = life.get(i, 0) + 1
            if any(np.hypot(x - c[0], y - c[1]) <= XY_DUP_M for c in cars):
                near[i] = near.get(i, 0) + 1
        if tr and view:
            d = np.array([[np.hypot(x - cars[k - 1][0], y - cars[k - 1][1]) for k in view] for _, x, y, *_ in tr])
            # gated assignment: pairs beyond the hit radius cannot be matched, so the assignment maximises the
            # matches within it (an ungated one can pair the track on a target with a far target to lower the total)
            hit = d <= XY_HIT_M + np.array([half_len.get(k, 0.0) for k in view])[None, :]
            rr, cc = linear_sum_assignment(np.where(hit, d, 1e6))
            for i, j in zip(rr, cc):
                if d[i, j] <= XY_HIT_M + half_len.get(view[j], 0.0):
                    r["xy_track_car"][view[j]] = tr[i][0]
                    on_car[tr[i][0]] = on_car.get(tr[i][0], 0) + 1
                    trk_err.setdefault(view[j], []).append(d[i, j])
        # resolution groups: the targets in view that fall in one peak of the map, chained
        rb, sc = r.get("res", (2.44, 0.3))
        tg = {**r.get("targets_any", {}), **{g["car"]: g for g in r["targets"]}}
        tg = {k: tg[k] for k in view if k in tg}
        comp = {k: {k} for k in view}
        for a in view:
            for b in view:
                if a < b and a in tg and b in tg and abs(tg[a]["excess"] - tg[b]["excess"]) <= res_bins * rb \
                        and abs(tg[a]["speed"] - tg[b]["speed"]) <= RES_SPEED_CELLS * max(sc, SPEED_FLOOR_MS):
                    u = comp[a] | comp[b]
                    for m in u:
                        comp[m] = u
        on = lambda q, k: np.hypot(q[1] - cars[k - 1][0], q[2] - cars[k - 1][1]) <= XY_HIT_M + half_len.get(k, 0.0) \
            or in_extent(q, np.asarray(cars[k - 1][:2], float))
        for k in view:
            pc = per_car.setdefault(k, dict(maps=0, covered=0, ids=[], grouped=0, alone=0, fair=0, fair_ids=[], merged=False))
            pc["maps"] += 1
            pc["grouped"] += any(on(q, k) for q in tr)
            tid = r["xy_track_car"].get(k)
            if tid is not None:
                pc["covered"] += 1
                if not pc["ids"] or pc["ids"][-1] != tid:
                    pc["ids"].append(tid)
            # fair score: alone, its own track; in a resolution group, any track on the group.
            # A change of track counts as a switch only between maps where it was alone, with
            # no merge in between (leaving a group under another track is unavoidable).
            if len(comp[k]) == 1:
                pc["alone"] += 1
                pc["fair"] += tid is not None
                fi = pc["fair_ids"]
                if tid is not None:
                    if pc["merged"] and fi:
                        fi.append(None)   # a merge in between: the next change is not a switch
                    if not fi or fi[-1] != tid:
                        fi.append(tid)
                pc["merged"] = False
            else:
                pc["fair"] += any(on(q, m) for q in tr for m in comp[k])
                pc["merged"] = True
        g_all.append(xy_gospa([(x, y) for _, x, y, *_ in tr], [cars[k - 1] for k in view]))
    rms = lambda e: float(np.sqrt(np.mean(np.square(e)))) if len(e) else None
    med = lambda e: float(np.median(e)) if len(e) else None
    cars_out = []
    for k in sorted(set(est_err) | set(per_car)):
        ee, pc = est_err.get(k, dict(maps=0, errs=[])), per_car.get(k, dict(maps=0, covered=0, ids=[], grouped=0, alone=0, fair=0, fair_ids=[]))
        fi = pc["fair_ids"]
        cars_out.append(dict(car=k, maps_in_view=ee["maps"], n_pos=len(ee["errs"]),
                             pos_rmse_m=rms(ee["errs"]), pos_median_m=med(ee["errs"]),
                             coverage=pc["covered"] / pc["maps"] if pc["maps"] else None,
                             group_coverage=pc["grouped"] / pc["maps"] if pc["maps"] else None,
                             switches=max(0, len(pc["ids"]) - 1),
                             alone_share=pc["alone"] / pc["maps"] if pc["maps"] else None,
                             fair_coverage=pc["fair"] / pc["maps"] if pc["maps"] else None,
                             fair_switches=sum(1 for a, b in zip(fi, fi[1:]) if a is not None and b is not None and a != b),
                             track_rmse_m=rms(trk_err.get(k, [])), track_median_m=med(trk_err.get(k, []))))
    g = np.array(g_all) if g_all else np.zeros((0, 4))
    off = [i for i, n in life.items() if on_car.get(i, 0) < n / 2]
    dup = [i for i in off if near.get(i, 0) >= life[i] / 2]
    return dict(cars=cars_out, n_tracks=len(life), false_tracks=len(off), duplicates=len(dup),   # false: off every target; of them, duplicates: near one
                gospa_m=float(g[:, 0].mean()) if len(g) else None,
                gospa_loc_m=float(g[:, 1].mean()) if len(g) else None,
                missed_per_map=float(g[:, 2].mean()) if len(g) else None,
                false_per_map=float(g[:, 3].mean()) if len(g) else None,
                params=dict(gate_m=XY_GATE_M, hit_m=XY_HIT_M, gospa_c_m=XY_GOSPA_C), ekf=ekf_cfg(cfg))


def in_extent(q, c, k2=9.0):
    """c inside track q's extent (id, x, y, vx, vy, a, b, angle) within 3 sigma."""
    if len(q) < 8 or q[5] <= 0:
        return False
    u = np.array([np.cos(q[7]), np.sin(q[7])])
    d = np.asarray(c, float) - [q[1], q[2]]
    return (d @ u / q[5]) ** 2 + ((d @ [-u[1], u[0]]) / max(q[6], 1e-3)) ** 2 <= k2


def xy_gospa(est, truth, c=XY_GOSPA_C):
    """GOSPA (p = 2, alpha = 2) in metres between estimated and true (x, y) points."""
    if not est and not truth:
        return 0.0, 0.0, 0, 0
    d = np.array([[min(np.hypot(e[0] - t[0], e[1] - t[1]), c) for t in truth] for e in est]).reshape(len(est), len(truth))
    loc, n = 0.0, 0
    if est and truth:
        for i, j in zip(*linear_sum_assignment(d)):
            if d[i, j] < c:
                loc += d[i, j] ** 2
                n += 1
    missed, false = len(truth) - n, len(est) - n
    return float(np.sqrt(loc + c * c / 2 * (missed + false))), float(np.sqrt(loc)), missed, false
