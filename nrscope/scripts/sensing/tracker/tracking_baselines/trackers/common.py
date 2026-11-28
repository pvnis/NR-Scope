"""What every tracker shares, so that only the tracking logic differs between them:
the bistatic measurement model and its noise (sim/track.py: bistatic_h, the run's range bin
and speed cell), the CV motion model, the birth state of a localised detection, the clutter
density, and the receiver's field of view (generic sensor knowledge, not site knowledge:
no walls, buildings or line-of-sight maps).

Input of a tracker (one call per run and setup): frames = [(t, meas, rb, sc, geo, rx)], in time
order, both receivers' maps interleaved when fused. meas: one dict per measurement group of the
map, e (excess range, m), v (bistatic speed, m/s), theta (AoA, rad, or None), xy (its localised
position, or None). geo = (T, R, L, bore). Output: per frame the reported tracks
[(id, x, y, vx, vy)]; localize_and_score scores them."""
import os
import sys

import numpy as np
from scipy.stats import chi2

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "sim"))
import track  # noqa: E402

SIGMA_E_BINS, SIGMA_V_CELLS, SIGMA_THETA_DEG = track.EKF_R_RANGE_BINS, track.EKF_R_SPEED_CELLS, track.EKF_SIGMA_THETA_DEG
P0_POS, P0_VEL = track.XY_R_M ** 2, 100.0          # birth covariance (as our tracker)
FOV = np.deg2rad(track.FOV_DEG)
CLUTTER_VOL = track.CLUTTER_VOL                   # (excess m, speed m/s, angle rad)


def Rz(rb, sc):
    return np.diag([(SIGMA_E_BINS * rb) ** 2, (SIGMA_V_CELLS * sc) ** 2, np.deg2rad(SIGMA_THETA_DEG) ** 2])


def zvec(m):
    """The measurement vector: (e, v, theta), or (e, v) without an angle."""
    return np.array([m["e"], m["v"]] + ([m["theta"]] if m["theta"] is not None else []))


def kappa(n_per_map, k):
    """Clutter density of k-component measurements: n_per_map false detections per map, uniform."""
    return n_per_map / np.prod(CLUTTER_VOL[:k])


def in_fov(x, geo):
    T, R, L, bore = geo
    a = np.arctan2(x[1] - R[1], x[0] - R[0]) - bore
    return abs((a + np.pi) % (2 * np.pi) - np.pi) <= FOV


def predict(x, P, dt, q):
    """CV, white acceleration of intensity q ((m/s^2)^2 per second)."""
    F = np.eye(4)
    F[0, 2] = F[1, 3] = dt
    qq = q * np.array([[dt ** 3 / 3, dt ** 2 / 2], [dt ** 2 / 2, dt]])
    Q = np.zeros((4, 4))
    Q[np.ix_([0, 2], [0, 2])] = qq
    Q[np.ix_([1, 3], [1, 3])] = qq
    return F @ x, F @ P @ F.T + Q


def innov(x, P, z, geo, R):
    """EKF innovation of measurement z (2 or 3 components): nu, S, H (on the 4-state)."""
    zh, H4, Rh = track.bistatic_h(x, *geo)
    k = len(z)
    nu = z - zh[:k]
    if k == 3:
        nu[2] = (nu[2] + np.pi) % (2 * np.pi) - np.pi
    H = H4[:k]
    return nu, H @ P @ H.T + R[:k, :k] + Rh[:k, :k], H


def update(x, P, nu, S, H):
    K = P @ H.T @ np.linalg.inv(S)
    return x + K @ nu, (np.eye(4) - K @ H) @ P


def gauss(nu, S):
    return float(np.exp(-0.5 * nu @ np.linalg.solve(S, nu)) / np.sqrt(np.linalg.det(2 * np.pi * S)))


def gate(k, p=0.999):
    return float(chi2.ppf(p, k))


def birth(m, geo):
    """Birth state of a localised measurement: its position, the velocity along the bistatic
    bisector that gives its measured speed (the observed component), the rest unknown."""
    x0 = np.r_[m["xy"], 0.0, 0.0]
    _, H, _ = track.bistatic_h(x0, *geo)
    s = -2 * H[1, 2:]
    x0[2:] = -2 * m["v"] * s / max(s @ s, 1e-6)
    return x0, np.diag([P0_POS, P0_POS, P0_VEL, P0_VEL])


def moment_match(ws, xs, Ps):
    ws = np.asarray(ws, float)
    ws = ws / ws.sum()
    x = sum(w * xi for w, xi in zip(ws, xs))
    P = sum(w * (Pi + np.outer(xi - x, xi - x)) for w, xi, Pi in zip(ws, xs, Ps))
    return x, P
