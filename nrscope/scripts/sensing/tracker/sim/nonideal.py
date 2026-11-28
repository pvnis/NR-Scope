"""Non-ideal receiver and targets for the Sionna channel (NR-Scope's copy only; the UE is untouched).

Off (enable false), the channel is the ideal one of before. On, it adds what a real capture has:

1. The array (NR-Scope's antennas)
   ant_gain_db     per-element gain error, standard deviation (dB)
   ant_phase_deg   per-element phase error, standard deviation (deg): the cables and front-ends
                   are not phase matched until calibrated
   coupling_db     mutual coupling to the adjacent elements (dB, random phase; -inf = none)
   An uncalibrated array: every path's steering vector is multiplied by M = C diag(g), C the
   coupling (identity plus the neighbours), g the element errors. The angle NR-Scope estimates
   is biased by a fixed, angle-dependent error.

2. The clocks (gNB and NR-Scope run on different oscillators)
   cfo_hz          carrier frequency offset of NR-Scope's oscillator against the gNB's (Hz; 0.1 ppm
                   at 3.4 GHz is 340 Hz). NR-Scope corrects it from the SSB, as on hardware: the run
                   turns its correction on (disable_cfo false) and the sensing sees what is left
   cfo_drift_hz_s  its drift (Hz/s), the oscillators warming up
   linewidth_hz    phase noise as a Wiener process of this 3-dB linewidth (Hz)
   timing_ppm      sample-clock offset (ppm): the delay of every path drifts by this, smoothly;
                   NR-Scope's synchronisation follows it from the SSB and re-aligns, as on hardware

4. The targets
   rcs_fluct_db    slow fluctuation of every target's echo (log-normal, standard deviation in dB)
   rcs_corr_s      its correlation time (s)
   micro_doppler   people swing their legs and arms: extra scatterers at the gait frequency
   vehicle_aspect  a vehicle reflects more broadside than end-on (bistatic aspect)

NR-Scope's side (set by the run, not here): carrier offset correction on, and with los_norm the
sensing divides each snapshot by its direct path (sensing.clutter_los_norm), as for hardware.

5. Target multipath
   target_multipath  echoes that also bounce once on a building wall, on the way to the target
                     or from it (image method on the footprints): ghosts at a longer range and
                     another angle, as on real captures
   wall_loss_db      the reflection loss of the wall (dB)
   wall_max_m        only walls closer than this to the target
"""
import numpy as np

NONIDEAL_DEFAULTS = dict(
    enable=False, seed=1,
    ant_gain_db=0.5, ant_phase_deg=10.0, coupling_db=-20.0,
    cfo_hz=300.0, cfo_drift_hz_s=1.0, linewidth_hz=50.0, timing_ppm=0.1,
    rcs_fluct_db=4.0, rcs_corr_s=0.5, micro_doppler=True, vehicle_aspect=True,
    target_multipath=True, wall_loss_db=6.0, wall_max_m=40.0,
    los_norm=True,
)


def settings(d=None):
    """The defaults, overridden by d (a JSON dict); None or enable false: off (returns None)."""
    c = dict(NONIDEAL_DEFAULTS)
    c.update({k: v for k, v in (d or {}).items() if k in c and v is not None and v != ""})
    for k, v in NONIDEAL_DEFAULTS.items():
        c[k] = type(v)(c[k]) if not isinstance(v, bool) else bool(c[k]) and str(c[k]).lower() not in ("0", "false")
    return c if c["enable"] else None


def array_matrix(c, n_ant):
    """M = C diag(g): what the uncalibrated array does to an ideal steering vector (s -> M s)."""
    rng = np.random.default_rng(c["seed"])
    g = 10 ** (rng.normal(0, c["ant_gain_db"], n_ant) / 20) * np.exp(1j * np.deg2rad(rng.normal(0, c["ant_phase_deg"], n_ant)))
    C = np.eye(n_ant, dtype=complex)
    if np.isfinite(c["coupling_db"]):
        k = 10 ** (c["coupling_db"] / 20)
        for i in range(n_ant - 1):
            C[i, i + 1] = C[i + 1, i] = k * np.exp(1j * rng.uniform(0, 2 * np.pi))
    return C @ np.diag(g)


class Clock:
    """Carrier phase (offset, drift, phase noise) and timing drift of NR-Scope against the gNB."""

    def __init__(self, c, fs):
        self.c, self.fs = c, fs
        self.rng = np.random.default_rng(c["seed"] + 1)
        self.phi = 0.0   # phase-noise phase carried from block to block

    def phase(self, n0, n):
        """Phase (rad) to apply to samples n0 .. n0+n-1."""
        t = (n0 + np.arange(n)) / self.fs
        p = 2 * np.pi * (self.c["cfo_hz"] * t + 0.5 * self.c["cfo_drift_hz_s"] * t ** 2)
        if self.c["linewidth_hz"] > 0:   # Wiener: increments of variance 2 pi linewidth / fs per sample
            w = np.cumsum(self.rng.normal(0, np.sqrt(2 * np.pi * self.c["linewidth_hz"] / self.fs), n)) + self.phi
            self.phi = float(w[-1])
            p = p + w
        return p

    def delay(self, t, max_s):
        """The common extra delay (s) at time t: the sample clock's drift, continuous. It is not
        wrapped here: NR-Scope's own synchronisation (SSB timing) re-aligns its window, as a
        real receiver does (a wrap here, one sample at once, put every echo at two ranges in a
        map). max_s: the most the channel's history can delay; the drift stops there."""
        return min(abs(self.c["timing_ppm"]) * 1e-6 * t, max_s) * np.sign(self.c["timing_ppm"] or 1)


class Fading:
    """Slow log-normal fluctuation of each target's echo (AR(1) in dB), one value per target."""

    def __init__(self, c):
        self.c = c
        self.rng = np.random.default_rng(c["seed"] + 2)
        self.x, self.t = {}, {}

    def amp(self, k, t):
        s, tau = self.c["rcs_fluct_db"], max(self.c["rcs_corr_s"], 1e-3)
        if s <= 0:
            return 1.0
        if k not in self.x:
            self.x[k] = self.rng.normal(0, s)
        else:
            r = np.exp(-max(t - self.t[k], 0.0) / tau)
            self.x[k] = r * self.x[k] + np.sqrt(1 - r * r) * self.rng.normal(0, s)
        self.t[k] = t
        return 10 ** (self.x[k] / 20)


def person_points(t, vel, phase0):
    """Body parts of a walking person, relative to the ground point under the torso: list of
    (offset (3,), velocity offset (3,), height for the visibility test, share of the RCS).
    Legs and arms swing along the walking direction at the gait frequency (about speed / 1.4 m),
    in opposition; the torso and head carry the rest (Chen's walking model, simplified)."""
    v = np.asarray(vel, float)[:2]
    sp = float(np.hypot(*v))
    h = v / sp if sp > 0.05 else np.array([1.0, 0.0])
    f = sp / 1.4   # strides per second
    w = 2 * np.pi * f
    out = [(np.array([0, 0, 1.0]), np.zeros(3), 1.0, 0.45), (np.array([0, 0, 1.6]), np.zeros(3), 1.6, 0.15)]
    for amp, z, share, z_vis, ph in ((0.35, 0.5, 0.15, 0.5, 0.0), (0.35, 0.5, 0.15, 0.5, np.pi),
                                     (0.25, 1.2, 0.05, 1.0, np.pi), (0.25, 1.2, 0.05, 1.0, 0.0)):
        s = np.sin(w * t + phase0 + ph) if sp > 0.05 else 0.0
        ds = w * np.cos(w * t + phase0 + ph) if sp > 0.05 else 0.0
        out.append((np.r_[amp * s * h, z], np.r_[amp * ds * h, 0.0], z_vis, share))
    return out


def aspect_gain(vel, dg, dr):
    """Power factor of a vehicle seen at a bistatic aspect: 1.75 broadside, 0.25 end-on (mean 1).
    dg = point - gNB, dr = point - receiver; the bisector of the two directions sees the side."""
    v = np.asarray(vel, float)[:2]
    if np.hypot(*v) < 0.1:
        return 1.0
    b = -(dg[:2] / np.linalg.norm(dg[:2]) + dr[:2] / np.linalg.norm(dr[:2]))
    c = abs(b @ v) / (np.linalg.norm(b) * np.linalg.norm(v) + 1e-9)   # 1 end-on, 0 broadside
    return 0.25 + 1.5 * (1 - c * c)


def walls(buildings):
    """Every footprint edge as a wall: arrays a (n,2), b (n,2), base z, top z."""
    A, B, z0, z1 = [], [], [], []
    for bd in buildings:
        q = np.asarray(bd["footprint"], float)
        if len(q) < 3:
            continue
        A.append(q); B.append(np.roll(q, -1, axis=0))
        base = float(bd.get("base") or 0.0)
        top = float(bd.get("top") or (base + float(bd.get("height") or 10.0)))
        z0.append(np.full(len(q), base)); z1.append(np.full(len(q), top))
    if not A:
        return None
    return np.concatenate(A), np.concatenate(B), np.concatenate(z0), np.concatenate(z1)


def bounce_points(W, src, tgt, wall_max):
    """One-bounce specular points on the walls between src and tgt (3D points): for each wall
    whose line mirrors src to an image from which the straight line to tgt crosses the wall
    segment on its face, the bounce point (3D) and the unfolded length. Returns [(P, length)]."""
    a, b, z0, z1 = W
    t2 = np.asarray(tgt, float)[:2]
    d_ab = b - a
    L2 = (d_ab ** 2).sum(1)
    # walls near the target (distance from the target to the segment)
    u = np.clip(((t2 - a) * d_ab).sum(1) / np.maximum(L2, 1e-9), 0, 1)
    near = np.hypot(*(t2 - (a + u[:, None] * d_ab)).T) < wall_max
    out = []
    s2 = np.asarray(src, float)[:2]
    for i in np.flatnonzero(near & (L2 > 1.0)):
        n = np.array([-d_ab[i, 1], d_ab[i, 0]]) / np.sqrt(L2[i])
        ds, dt = (s2 - a[i]) @ n, (t2 - a[i]) @ n
        if ds * dt <= 0:   # src and tgt on the same side of the wall, both off it
            continue
        img = s2 - 2 * ds * n
        # where the line image -> target crosses the wall's line
        f = (img - a[i]) @ n / ((img - t2) @ n)
        P2 = img + f * (t2 - img)
        v = ((P2 - a[i]) @ d_ab[i]) / L2[i]
        if not (0.02 < v < 0.98) or not (0 < f < 1):
            continue
        L3 = np.hypot(np.linalg.norm(t2 - img), np.asarray(tgt)[2] - np.asarray(src)[2])
        # the bounce height along the unfolded straight line
        zP = src[2] + (tgt[2] - src[2]) * f
        if not (z0[i] + 0.2 < zP < z1[i] - 0.2):
            continue
        out.append((np.r_[P2, zP], L3))
    return out
