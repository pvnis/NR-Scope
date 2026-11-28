"""Sionna RT channel for NR-Scope's copy of the downlink.

The scene is Sionna's built-in street canyon with one or more cars, each driving along its own
route (CarRoute):

  gNB (transmitter) ---- direct path, wall reflections -----> NR-Scope (receiver)
           \\                                                  /
            `---- reflections off the moving car (the target) '

Every `update_s` seconds of signal the car is moved to where it is at that time and the
paths are ray traced again (Sionna's PathSolver, on the GPU). Each path p has a complex gain
a_p, a delay tau_p and a Doppler shift f_p, and the received signal is

    y(t) = sum_p a_p * exp(j 2 pi f_p (t - t_update)) * x(t - tau_p) + noise

The delays are applied exactly (as a phase ramp in frequency, so fractions of a sample are
fine), and the Doppler term carries each path's phase between two ray traces. This runs on
the GPU with PyTorch.

Receive array (n_ant > 1): NR-Scope's antennas form a uniform linear array at `spacing`
wavelengths, centred on rx_pos is NOT assumed: element m sits at rx_pos + m * spacing * lambda
* axis, axis = the array normal (boresight_deg, counterclockwise from +x) turned +90 deg.
This is NR-Scope's convention (nr_ue_aoa.c, nr_ue_localize.c): element m receives a wave
arriving at angle theta off the normal (counterclockwise positive) with phase
exp(j 2 pi spacing m sin theta). Each path gets that phase per element: from its arrival
direction (plane wave) for Sionna's paths, from the exact distance to each element for the
point-scatterer cars. The delay across the array (under 0.5 ns) is neglected. Every element
gets its own noise. The output is then [n_ant, n] instead of [n].

Gains are scaled so the strongest path (normally the direct one) has gain 1, so NR-Scope
sees the same level as without the channel, and white noise is added `snr_db` below it.

Ground truth goes to `gt_path` (CSV, one row per path per update): what the receiver gets,
and which car each path touched (via_car: 0 = none, k = car k, numbered from 1; car_x.. are
car 1's state), and aoa_deg: its arrival angle off the array normal in NR-Scope's convention. Every car's position and velocity at each update go to <gt_path>_cars.csv
(t_s, car, x, y, z, vx, vy, vz). NR-Scope never reads them; they are only for evaluation.
"""
import csv
import os
import time

import mitsuba as mi
import numpy as np
import sionna.rt as rt
import torch

import nonideal

C = 299792458.0
NO_OBJECT = 0xFFFFFFFF  # Sionna's object id for "no interaction"


SCENES = ("simple_street_canyon", "etoile", "munich", "florence", "san_francisco")


def load_base(scene="simple_street_canyon"):
    """A Sionna built-in scene by name (SCENES), a real site "site:NAME" (osm_site.py), or a
    Mitsuba XML file (e.g. exported from Blender)."""
    if str(scene).endswith(".xml"):
        return rt.load_scene(scene)
    if str(scene).startswith("site:"):   # a real site built by osm_site.py
        return rt.load_scene(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sites", scene[5:], "scene.xml"))
    if scene not in SCENES:
        raise ValueError(f"unknown scene {scene}")
    return rt.load_scene(getattr(rt.scene, scene))


# How realistic the propagation is (the Scenarios page's "Propagation" settings):
#   pattern      the antennas of the gNB and of NR-Scope's elements, vertical: "iso" (the same in
#                every direction), "dipole" (short dipole, 1.5 sin^2) or "hw_dipole" (half-wave
#                dipole, the donut of a real omni, 2.15 dBi on the horizon, nothing straight up/down)
#   max_depth    reflections (and diffractions) per path
#   roughness    per material, Sionna's scattering coefficient S: the share of the field a surface
#                scatters in every direction instead of reflecting like a mirror (0 = mirror).
#                A material is matched by name ("concrete" also matches "itu_concrete").
#   diffraction  paths bending around building edges and corners (slower)
#   back_db      a shield behind NR-Scope's array (a metal plate, ~30 cm, that blocks, with a thin
#                absorber on its front so it does not reflect into the elements): waves arriving from
#                behind the array lose back_db (0 = no shield; such a shield 15-25, limited by what
#                bends around its edges). It fades in between 80 and 100 deg off the boresight.
PROPAGATION = dict(pattern="hw_dipole", max_depth=3, diffraction=False, back_db=0.0,
                   roughness=dict(concrete=0.3, brick=0.3, marble=0.3, plasterboard=0.3, wood=0.3, glass=0.1,
                                  metal=0.1, foliage=0.5, ground=0.2))
# the "mirror" setup of the first runs: isotropic antennas, 2 reflections, smooth surfaces
PROPAGATION_MIRROR = dict(pattern="iso", max_depth=2, diffraction=False, back_db=0.0, roughness={})


def propagation(p=None):
    """PROPAGATION updated by p (a dict, or JSON text); roughness merged material by material."""
    import json as _json
    p = _json.loads(p) if isinstance(p, str) else dict(p or {})
    out = dict(PROPAGATION, **{k: v for k, v in p.items() if k != "roughness"})
    out["roughness"] = dict(PROPAGATION["roughness"], **(p.get("roughness") or {}))
    return out


def pattern_gain(pattern, d):
    """Power gain of a vertical antenna towards directions d (n, 3): iso 1, dipole 1.5 sin^2(theta),
    hw_dipole 1.64 (cos(pi/2 cos theta) / sin theta)^2 (theta from the vertical)."""
    d = np.atleast_2d(np.asarray(d, float))
    ct = np.clip(d[:, 2] / np.maximum(np.linalg.norm(d, axis=1), 1e-12), -1.0, 1.0)
    st2 = np.maximum(1.0 - ct ** 2, 1e-12)
    if pattern == "dipole":
        return 1.5 * st2
    if pattern == "hw_dipole":
        return 1.64 * np.cos(np.pi / 2 * ct) ** 2 / st2
    return np.ones(len(d))


def roughness_of(material_name, rough):
    """The roughness (scattering coefficient) of a material by name: the longest key it contains."""
    n = material_name.lower()
    keys = [k for k in rough if k in n]
    return float(rough[max(keys, key=len)]) if keys else 0.0


def build_scene(gnb_pos, rx_pos, fc=3.41499e9, n_cars=1, scene="simple_street_canyon", prop=None):
    """The scene (load_base) with the gNB, NR-Scope and n_cars rough metal cars (returned
    separately, a list). prop: propagation() settings (antenna pattern, surface roughness); None:
    the mirror setup (PROPAGATION_MIRROR)."""
    prop = PROPAGATION_MIRROR if prop is None else prop
    scene = load_base(scene)
    scene.frequency = fc
    scene.tx_array = rt.PlanarArray(num_rows=1, num_cols=1, pattern=prop["pattern"], polarization="V")
    scene.rx_array = scene.tx_array
    scene.add(rt.Transmitter("gnb", position=[float(v) for v in gnb_pos]))
    scene.add(rt.Receiver("nrscope", position=[float(v) for v in rx_pos]))
    for obj in scene.objects.values():  # mirrors, or rough by material
        obj.radio_material.scattering_coefficient = roughness_of(obj.radio_material.name, prop["roughness"])
    # A small car only reflects towards the receiver by chance if it is a mirror;
    # a rough surface scatters in every direction, as real vehicles do.
    metal = rt.ITURadioMaterial("car-metal", "metal", thickness=0.01, scattering_coefficient=0.7)
    cars = [rt.SceneObject(fname=rt.scene.low_poly_car, name=f"car{i + 1}", radio_material=metal)
            for i in range(n_cars)]
    scene.edit(add=cars)
    return scene, cars


class CarRoute:
    """A target moving along a path at `speed` m/s. Path shapes (`path`):
      line    straight from `start` to `end`
      curve   from `start` to `end` along a quadratic Bezier bulging `bend` metres sideways at
              its middle (positive: to the left of the direction of travel)
      circle  around `end` (the centre) through `start`, counterclockwise; its length is the
              circumference and the target goes round and round whatever `mode` says
    At the end of a line or curve the target stops ("stop") or comes back the same way
    ("loop"). The path is sampled finely and walked by arc length; the velocity is along the
    path's tangent. On a line or curve the target moves like a vehicle: it accelerates from
    rest at `accel` m/s^2 up to `speed`, cruises, and brakes at `accel` to stop at the end, then
    ("loop") does the same on the way back. No instant reversal: speed and velocity are
    continuous. On a circle the speed is constant."""

    def __init__(self, start, end, speed, mode="stop", path="line", bend=0.0, accel=2.0):
        self.start = np.array(start, dtype=float)
        self.end = np.array(end, dtype=float)
        self.speed = float(speed)
        self.mode = mode
        self.path = path
        self.pts = route_points(self.start, self.end, path, bend)
        seg = np.linalg.norm(np.diff(self.pts, axis=0), axis=1)
        self.cum = np.r_[0.0, np.cumsum(seg)]
        self.length = float(self.cum[-1])
        self.accel = float(accel)
        # one leg (rest to rest over the path): trapezoidal speed profile, or triangular when
        # the path is too short to reach `speed`
        a, v, L = self.accel, self.speed, self.length
        self.v_top = min(v, np.sqrt(a * L)) if a > 0 and v > 0 else v
        self.t_ramp = self.v_top / a if a > 0 else 0.0
        self.t_leg = (L - self.v_top * self.t_ramp) / self.v_top + 2 * self.t_ramp if self.v_top > 0 else np.inf

    def _leg(self, tau):
        """Distance and speed tau seconds into a leg."""
        a, vt, tr, tl = self.accel, self.v_top, self.t_ramp, self.t_leg
        if tau <= 0:
            return 0.0, 0.0
        if tau >= tl:
            return self.length, 0.0
        if tau < tr:
            return 0.5 * a * tau * tau, a * tau
        if tau > tl - tr:
            r = tl - tau
            return self.length - 0.5 * a * r * r, a * r
        return 0.5 * a * tr * tr + vt * (tau - tr), vt

    def _along(self, d):
        """Position and unit tangent at arc length d (0 <= d <= length)."""
        i = int(np.clip(np.searchsorted(self.cum, d, side="right") - 1, 0, len(self.pts) - 2))
        seg = self.pts[i + 1] - self.pts[i]
        n = np.linalg.norm(seg)
        f = (d - self.cum[i]) / n if n > 0 else 0.0
        return self.pts[i] + f * seg, (seg / n if n > 0 else np.zeros(3))

    def at(self, t):
        """Position and velocity at time t (s)."""
        if self.length == 0 or self.speed == 0:
            return self.start.copy(), np.zeros(3)
        d = self.speed * t
        if self.path == "circle":
            p, u = self._along(d % self.length)
            return p, u * self.speed
        if self.accel <= 0:  # constant speed, instant reversal (the old behaviour)
            if self.mode == "loop":
                d = d % (2 * self.length)
                if d > self.length:
                    p, u = self._along(2 * self.length - d)
                    return p, -u * self.speed
            elif d >= self.length:
                return self.pts[-1].copy(), np.zeros(3)
            p, u = self._along(min(d, self.length))
            return p, u * self.speed
        k, tau = (int(t // self.t_leg), t % self.t_leg) if self.mode == "loop" else (0, min(t, self.t_leg))
        d, v = self._leg(tau)
        if k % 2:  # on the way back
            p, u = self._along(self.length - d)
            return p, -u * v
        p, u = self._along(d)
        return p, u * v


def route_points(start, end, path="line", bend=0.0, n=241):
    """The sampled path of CarRoute (the GUI draws the same shapes, routePoints in index.html)."""
    start, end = np.asarray(start, float), np.asarray(end, float)
    if path == "circle":
        r = np.hypot(*(start - end)[:2])
        a0 = np.arctan2(start[1] - end[1], start[0] - end[0])
        a = a0 + np.linspace(0, 2 * np.pi, n)
        return np.stack([end[0] + r * np.cos(a), end[1] + r * np.sin(a), np.full(n, start[2])], 1)
    if path == "curve" and bend:
        d = end - start
        left = np.array([-d[1], d[0], 0.0]) / max(np.hypot(d[0], d[1]), 1e-9)
        ctrl = (start + end) / 2 + 2 * bend * left  # a quadratic Bezier's middle sits halfway to its control point
        s = np.linspace(0, 1, n)[:, None]
        return (1 - s) ** 2 * start + 2 * s * (1 - s) * ctrl + s ** 2 * end
    return np.stack([start, end])


class CarBody:
    """A car as rigid point scatterers: bumper corners and centres, roof corners.

    Body frame: x forward (along the direction of travel), y to the left, z up, origin on the
    ground under the car's centre, as Sionna's car mesh. Each scatterer has a fixed complex
    weight (equal power, random phase drawn once); the weights' powers sum to one."""

    L, W, H_BUMPER, H_ROOF = 4.4, 1.8, 0.5, 1.4

    def __init__(self, seed=0, size=None):
        """size: (length, width, height) of another vehicle (scenario.CLASSES); default a car."""
        L, W, H = size if size is not None else (self.L, self.W, self.H_ROOF)
        hb, hr = (self.H_BUMPER, self.H_ROOF) if size is None else (0.36 * H, H)
        self.height = hr
        l, w = L / 2, W / 2
        self.body = np.array([[l, w, hb], [l, 0, hb], [l, -w, hb],
                              [-l, w, hb], [-l, 0, hb], [-l, -w, hb],
                              [0.6 * l, 0.8 * w, hr], [0.6 * l, -0.8 * w, hr],
                              [-0.6 * l, 0.8 * w, hr], [-0.6 * l, -0.8 * w, hr]])
        ph = np.random.default_rng(seed + 1).uniform(0, 2 * np.pi, len(self.body))
        self.weights = np.exp(1j * ph) / np.sqrt(len(self.body))
        self.heading = np.array([1.0, 0.0])

    def points(self, pos, vel):
        """World positions of the scatterers for a car at pos moving at vel."""
        v = np.asarray(vel, float)[:2]
        if np.linalg.norm(v) > 1e-6:  # keep the last heading while the car stands still
            self.heading = v / np.linalg.norm(v)
        fx, fy = self.heading
        rot = np.array([[fx, -fy, 0], [fy, fx, 0], [0, 0, 1]])
        return np.asarray(pos, float) + self.body @ rot.T


PERSON_RCS_M2 = 0.5   # a walking adult at 3.5 GHz: about 0.1-1 m^2 in the literature, fluctuating
# A person is three point scatterers, legs, torso and head-shoulders: heights above the ground (m)
# and the share of the RCS (power) each carries. Each is hidden or seen on its own, so a person
# behind a low wall still returns the part above it, weaker, as a real one does.
PERSON_POINTS = ((0.5, 0.3), (1.0, 0.5), (1.6, 0.2))


class SionnaChannel:
    def __init__(self, fs, fc=3.41499e9, update_s=0.005,
                 gnb_pos=(-40.0, 0.0, 8.0), rx_pos=(10.0, -6.0, 2.0),
                 car_start=(-20.0, 3.0, 0.0), car_end=(40.0, 3.0, 0.0), car_speed=5.0, car_mode="stop",
                 snr_db=30.0, max_depth=2, history=8192, gt_path=None, seed=0, stitch=True,
                 car_model="scatterers", cars=None, n_ant=1, boresight_deg=90.0, spacing=0.5,
                 scenario=None, scene="simple_street_canyon", above_ground=False, prop=None, nonideal=None):
        self.fs = fs
        # propagation settings (propagation()); None: the mirror setup of the first runs
        self.prop = dict(PROPAGATION_MIRROR, max_depth=max_depth) if prop is None else propagation(prop)
        max_depth = int(self.prop["max_depth"])
        self.lam = C / fc
        self.n_ant = int(n_ant)
        self.boresight = np.deg2rad(boresight_deg)
        b = self.boresight
        self.normal = np.array([np.cos(b), np.sin(b), 0.0])
        self.axis = np.array([-np.sin(b), np.cos(b), 0.0])  # element m at rx + m * spacing * lambda * axis
        self.elems = np.asarray(rx_pos, float) + np.outer(np.arange(self.n_ant) * spacing * self.lam, self.axis)
        self.car_model = car_model  # "scatterers" (rigid point scatterers) or "sionna" (diffuse paths)
        # cars: the targets, [dict(start=, end=, speed=, mode=, type="car"|"person", rcs=)]; without
        # it, the one car of car_start..car_mode. A car is a Sionna mesh (it also blocks and
        # reflects) whose echo is made of CarBody's scatterers; a person is one point scatterer
        # of radar cross-section rcs (m^2, default PERSON_RCS_M2), no mesh (see _person). A car
        # with an rcs gets its echo power from the radar equation too (see _scatterer_car);
        # without, from Sionna's paths off its mesh, which come out near 0.1 m^2.
        #
        # scenario: a file of scenario.py (many people and vehicles that come and go) instead of
        # cars. Its vehicles have no mesh (kind "vehicle"): CarBody scatterers of their class's
        # size, power from the radar equation; its people are _person's. Every agent is placed
        # on the ground under it and is seen only if the buildings block neither the gNB's ray
        # to it nor its ray to the receiver (_visible).
        self.agents, self.n_hidden = None, 0
        if scenario:
            from scenario import Agents
            self.agents = Agents(scenario)
            self.kinds = ["person" if c == "person" else "vehicle" for c in self.agents.cls]
            self.rcs = [float(r) for r in self.agents.rcs]
            self.routes = []
            self.car_bodies = [None if c == "person" else CarBody(seed + 7 * i, tuple(self.agents.size[i]))
                               for i, c in enumerate(self.agents.cls)]
        else:
            cars = cars or [dict(start=car_start, end=car_end, speed=car_speed, mode=car_mode)]
            self.kinds = [c.get("type", "car") for c in cars]
            self.rcs = [float(c["rcs"]) if c.get("rcs") else (PERSON_RCS_M2 if c.get("type") == "person" else None)
                        for c in cars]
            self.routes = [CarRoute(c["start"], c["end"], c["speed"], c.get("mode", "stop"), c.get("path", "line"),
                                    float(c.get("bend") or 0.0),
                                    accel=float(c["accel"]) if c.get("accel") is not None
                                    else (1.0 if c.get("type") == "person" else 2.0)) for c in cars]
            self.car_bodies = [CarBody(seed + 7 * i) for i in range(len(cars))]
        self.car_power = [None] * len(self.kinds)  # smoothed total power of Sionna's paths per car (scatterers model)
        self.update_s = update_s
        self.snr_db = snr_db
        self.stitch = stitch  # phase-continuous car echo across ray traces (see _stitch)
        self.history = history
        self.max_depth = max_depth
        self.dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.rng = torch.Generator(device=self.dev).manual_seed(seed)
        self.scene, meshes = build_scene(gnb_pos, rx_pos, fc, self.kinds.count("car"), scene, self.prop)
        self.top = float(np.array(self.scene.mi_scene.bbox().max)[2]) + 10.0  # above every roof
        self.bottom = float(np.array(self.scene.mi_scene.bbox().min)[2]) - 10.0  # below the ground
        if above_ground:  # gnb_pos and rx_pos z are heights above the ground under them
            pts = np.array([gnb_pos, rx_pos], float)
            pts[:, 2] += self._ground(pts[:, :2])
            self.scene.transmitters["gnb"].position = mi.Point3f(*map(float, pts[0]))
            self.scene.receivers["nrscope"].position = mi.Point3f(*map(float, pts[1]))
            self.elems = pts[1] + np.outer(np.arange(self.n_ant) * spacing * self.lam, self.axis)
            self.gnb_pos, self.rx_pos = pts[0], pts[1]
        it = iter(meshes)
        self.mesh = {k: next(it) for k, kind in enumerate(self.kinds, 1) if kind == "car"}  # target k -> its mesh
        self.cars = list(self.mesh.values())
        self.solver = rt.PathSolver()

        self.buf = torch.zeros(history, dtype=torch.complex64, device=self.dev)  # last input samples
        self.n = 0                    # samples processed so far (signal time = n / fs)
        self.next_update = 0          # sample index of the next ray trace
        self.paths = None
        self.sig_power = None         # running mean power of the gNB signal (non-silent blocks)
        self.trace_s = 0.0
        self.n_traces = 0

        self.gt_file = None
        if gt_path:
            self.gt_file = open(gt_path, "w", newline="")
            self.gt = csv.writer(self.gt_file)
            self.gt.writerow(["t_s", "car_x", "car_y", "car_z", "car_vx", "car_vy", "car_vz",
                              "path", "via_car", "delay_s", "bistatic_range_m", "gain_db", "doppler_hz", "aoa_deg"])
            self.cars_file = open(gt_path.rsplit(".", 1)[0] + "_cars.csv", "w", newline="")
            self.cars_gt = csv.writer(self.cars_file)
            self.cars_gt.writerow(["t_s", "car", "x", "y", "z", "vx", "vy", "vz"])

        # a real receiver and real targets (nonideal.py); None: the ideal channel
        self.ni = globals()["nonideal"].settings(nonideal)
        self.walls = None
        if self.ni:
            ni = globals()["nonideal"]
            self.ni_M = ni.array_matrix(self.ni, self.n_ant)
            self.clock = ni.Clock(self.ni, fs)
            self.fading = ni.Fading(self.ni)
            self.gait = np.random.default_rng(self.ni["seed"] + 3).uniform(0, 2 * np.pi, len(self.kinds))
            geo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sites", str(scene)[5:], "geometry.json")
            if self.ni["target_multipath"] and str(scene).startswith("site:") and os.path.exists(geo):
                import json
                self.walls = ni.walls(json.load(open(geo))["buildings"])

    # -- ray tracing ---------------------------------------------------------------------------
    def _rays(self, o, d, maxt, test=True):
        """Mitsuba rays from o along d (unit, (n, 3)): blocked before maxt (test) or the hit distance."""
        f = lambda a: mi.Float(np.ascontiguousarray(a, dtype=np.float32))
        ray = mi.Ray3f(mi.Point3f(f(o[:, 0]), f(o[:, 1]), f(o[:, 2])), mi.Vector3f(f(d[:, 0]), f(d[:, 1]), f(d[:, 2])))
        ray.maxt = f(maxt)
        if test:
            return np.array(self.scene.mi_scene.ray_test(ray)).astype(bool).ravel()
        si = self.scene.mi_scene.ray_intersect(ray)
        return np.where(np.array(si.is_valid()).astype(bool).ravel(), np.array(si.t).ravel(), np.inf)

    def _ground(self, xy):
        """The ground under each point: the lowest surface, hit going up from below the scene
        (so a roof overhang or a tree above is not taken for the ground)."""
        xy = np.atleast_2d(np.asarray(xy, float))
        n = len(xy)
        hit = self._rays(np.c_[xy, np.full(n, self.bottom)], np.tile([0.0, 0.0, 1.0], (n, 1)),
                         np.full(n, self.top - self.bottom), test=False)
        return np.where(np.isfinite(hit), self.bottom + hit, 0.0)

    def _agent_states(self, t):
        """Scenario agents at t: (pos on the ground, vel) or None (not there, or hidden from the gNB
        or the receiver by a building)."""
        st = self.agents.at(t)
        on = [i for i, x in enumerate(st) if x is not None]
        self.last_agents = (st, None, on)  # every agent's place, seen or not (ground truth)
        if not on:
            return st
        xy = np.array([st[i][0] for i in on])
        z = self._ground(xy[:, :2])  # the ground (or terrain) under each
        g = np.array(self.scene.transmitters["gnb"].position.numpy()).ravel()
        r = np.array(self.scene.receivers["nrscope"].position.numpy()).ravel()
        # the points tested: a vehicle's centre; a person's three points (PERSON_POINTS), seen if any is
        hs = [[q[0] for q in PERSON_POINTS] if self.kinds[i] == "person" else [0.5 * self.car_bodies[i].height] for i in on]
        own = np.repeat(np.arange(len(on)), [len(x) for x in hs])
        c = np.c_[xy[own, :2], z[own] + np.concatenate(hs)]
        vis = np.ones(len(c), bool)
        for end in (g, r):
            d = end - c
            n = np.linalg.norm(d, axis=1)
            vis &= ~self._rays(c, d / n[:, None], n - 0.05)
        seen = np.array([vis[own == j].any() for j in range(len(on))])
        self.point_vis = {i: vis[own == j] for j, i in enumerate(on)}   # per point, for _person
        self.n_hidden = int((~seen).sum())
        out = [None] * len(st)
        for j, i in enumerate(on):
            out[i] = (np.array([xy[j, 0], xy[j, 1], z[j]]), st[i][1]) if seen[j] else None
        return out

    def _trace(self, t):
        states = self._agent_states(t) if self.agents is not None else [r.at(t) for r in self.routes]
        for k, car in self.mesh.items():
            pos, vel = states[k - 1]
            car.position = mi.Point3f(*map(float, pos))
            car.velocity = mi.Vector3f(*map(float, vel))

        # A scenario's agents are not in the Sionna scene (they are added below as scatterers), so
        # the scene never changes: trace it once. Re-tracing would also re-draw the random rays of
        # the rough surfaces every update_s, a flicker of the whole static scene (false Doppler).
        if self.agents is not None and not self.mesh and getattr(self, "_static", None) is not None:
            a, tau, dop, u, steer, via_car = (v.copy() for v in self._static)
        else:
            t0 = time.monotonic()
            dif = bool(self.prop["diffraction"])
            p = self.solver(self.scene, max_depth=self.max_depth, diffuse_reflection=True,
                            samples_per_src=200000, diffraction=dif, edge_diffraction=dif,
                            seed=0 if getattr(self, 'fixed_seed', False) else self.n_traces)
            a, tau = p.cir(normalize_delays=False, out_type="numpy")
            self.trace_s += time.monotonic() - t0
            self.n_traces += 1

            valid = np.array(p.valid).ravel().astype(bool)
            a = a.ravel()[valid]
            tau = tau.ravel()[valid]
            dop = np.array(p.doppler).ravel()[valid]
            th, ph = np.array(p.theta_r).ravel()[valid], np.array(p.phi_r).ravel()[valid]
            u = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)  # rx -> where it came from
            tt, pt = np.array(p.theta_t).ravel()[valid], np.array(p.phi_t).ravel()[valid]
            self.ut_sionna = np.stack([np.sin(tt) * np.cos(pt), np.sin(tt) * np.sin(pt), np.cos(tt)], 1)  # gNB -> departure
            steer = np.exp(2j * np.pi / self.lam * (u @ (self.elems - self.elems[0]).T))  # [paths, n_ant]
            objs = np.array(p.objects).reshape(self.max_depth, -1)[:, valid]
            via_car = np.zeros(len(a), int)  # 0: no car, k: car k (the first car the path touched)
            for k, car in self.mesh.items():
                via_car[(via_car == 0) & (objs == car.object_id).any(axis=0)] = k
            if self.agents is not None and not self.mesh:
                self._static = tuple(v.copy() for v in (a, tau, dop, u, steer, via_car))

        if getattr(self, "no_car", False):  # debug: drop the paths that touch the car
            keep = via_car == 0
            a, tau, dop, via_car, u, steer = a[keep], tau[keep], dop[keep], via_car[keep], u[keep], steer[keep]
        if getattr(self, "only_los", False):  # debug: keep the direct (shortest) path only
            k = int(np.argmin(tau))
            a, tau, dop, via_car, u, steer = a[k:k + 1], tau[k:k + 1], dop[k:k + 1], via_car[k:k + 1], u[k:k + 1], steer[k:k + 1]
        if getattr(self, "zero_delay", False):  # debug: delays relative to the direct path, direct at 0
            tau = tau - tau.min()
        if self.paths is None and not len(a):
            raise RuntimeError("no propagation path from the gNB to NR-Scope: is one of them inside a building?")
        if self.paths is None:  # fix the scale once, on the strongest path of the first trace
            self.scale = 1.0 / np.abs(a).max()
        a = a * self.scale
        if self.agents is not None:
            if not getattr(self, "only_los", False) and not getattr(self, "no_car", False):
                for k, st in enumerate(states, 1):
                    if st is not None:
                        f = self._person if self.kinds[k - 1] == "person" else self._scatterer_car
                        a, tau, dop, via_car, u, steer = f(k, *st, a, tau, dop, via_car, u, steer)
        elif self.car_model == "scatterers" and not getattr(self, "no_car", False) \
                and not getattr(self, "only_los", False):
            for k, (pos, vel) in enumerate(states, 1):
                if self.kinds[k - 1] == "car":
                    a, tau, dop, via_car, u, steer = self._scatterer_car(k, pos, vel, a, tau, dop, via_car, u, steer)
        elif self.stitch:
            a = self._stitch(t, a, tau, dop, via_car)
        if self.agents is None and not getattr(self, "only_los", False):
            for k, (pos, vel) in enumerate(states, 1):
                if self.kinds[k - 1] == "person":
                    a, tau, dop, via_car, u, steer = self._person(k, pos, vel, a, tau, dop, via_car, u, steer)
        if self.agents is not None and self.ni and self.walls is not None and not getattr(self, "only_los", False):
            for k, st in enumerate(states, 1):
                if st is not None:
                    a, tau, dop, via_car, u, steer = self._wall_echoes(k, *st, t, a, tau, dop, via_car, u, steer)
        back_db = float(self.prop.get("back_db") or 0.0)
        if back_db > 0:   # the shield behind the array (back_db above)
            c = (u @ self.normal) / np.maximum(np.hypot(u @ self.normal, u @ self.axis), 1e-12)   # cos of the azimuth off boresight
            w = np.clip((np.cos(np.deg2rad(80)) - c) / (2 * np.cos(np.deg2rad(80))), 0.0, 1.0)
            a = a * 10 ** (-back_db * w / 20)
        if self.ni:   # the uncalibrated array: every steering vector s becomes M s
            steer = steer @ self.ni_M.T
        cone = np.rad2deg(np.arcsin(np.clip(u @ self.axis, -1, 1)))   # what the linear array measures: asin(cos(el) sin(az))
        aoa = np.where(u @ self.normal >= 0, cone, np.copysign(180.0, cone) - cone)  # behind the array: beyond +-90
        # Paths that do not move (no Doppler, no target) are summed once into one filter per
        # antenna (_static_filter); the others are applied one by one with their Doppler.
        # every path of this trace, for offline use (tools/two_ue/precode_offline.py): the first
        # len(self.ut_sionna) are Sionna's, the rest the targets' added above
        self.last = dict(t=t, a=a, tau=tau, dop=dop, u=u, via_car=via_car)
        still = (via_car == 0) & (np.abs(dop) < 1e-9)
        mv = ~still if still.sum() > getattr(self, "fir_min", 8) else np.ones(len(a), bool)
        h = self._static_filter(a[~mv], tau[~mv], steer[~mv]) if (~mv).any() else None
        self.paths = dict(t=t, a=torch.tensor(a[mv], dtype=torch.complex64, device=self.dev),
                          steer=torch.tensor(steer[mv].T, dtype=torch.complex64, device=self.dev),  # [n_ant, paths]
                          tau=torch.tensor(tau[mv], dtype=torch.float64, device=self.dev),
                          dop=torch.tensor(dop[mv], dtype=torch.float64, device=self.dev), h=h)
        on = via_car > 0
        self.prev_car = (t, a[on], tau[on], dop[on])

        if self.gt_file:
            if self.agents is not None:  # every agent there, seen or hidden; nan when not there
                every = self.last_agents[0]
                truth = [None if x is None else (np.r_[x[0][:2], 0.0], x[1]) for x in every]
            else:
                truth = states
            nan = (np.full(3, np.nan), np.full(3, np.nan))
            pos, vel = truth[0] if truth and truth[0] is not None else nan
            for k, pv in enumerate(truth, 1):
                p, v = pv if pv is not None else nan
                self.cars_gt.writerow([f"{t:.6f}", k, *np.round(p, 4), *np.round(v, 4)])
            self.cars_file.flush()
            # every target's path, and the 50 strongest others (a rough scene has thousands)
            rows = np.flatnonzero(via_car > 0).tolist()
            rest = np.flatnonzero(via_car == 0)
            rows += rest[np.argsort(-np.abs(a[rest]))][:50].tolist()
            for i in sorted(rows):
                self.gt.writerow([f"{t:.6f}", *np.round(pos, 4), *np.round(vel, 4), i, int(via_car[i]),
                                  f"{tau[i]:.4e}", f"{tau[i] * C:.3f}", f"{20 * np.log10(abs(a[i]) + 1e-30):.2f}",
                                  f"{dop[i]:.3f}", f"{aoa[i]:.2f}"])
            self.gt_file.flush()

    def _scatterer_car(self, k, pos, vel, a, tau, dop, via_car, u, steer):
        """Replace Sionna's paths off car k by the rigid point-scatterer car.

        Sionna's diffuse scattering is a new random speckle pattern at every trace: the car's
        echo level jumped by about 7 dB rms every 5 ms, which spread the car over the whole
        Doppler axis. Here the car is CarBody's fixed scatterers. Each has the exact path
        length gNB -> scatterer -> receiver (so its phase is continuous) and the Doppler of its
        own velocity, and the echo's total power follows Sionna's car paths smoothed over about
        0.25 s, so the level the ray tracer predicts is kept without its speckle."""
        mine = via_car == k
        p_sionna = float(np.sum(np.abs(a[mine]) ** 2))
        if p_sionna > 0:
            w = 1.0 - np.exp(-self.update_s / 0.25)
            old = self.car_power[k - 1]
            self.car_power[k - 1] = p_sionna if old is None else (1 - w) * old + w * p_sionna
        a, tau, dop, via_car, u, steer = a[~mine], tau[~mine], dop[~mine], via_car[~mine], u[~mine], steer[~mine]
        power = self.car_power[k - 1]
        if self.rcs[k - 1]:  # radar equation at the car's centre, as for a person (_person)
            g = np.array(self.scene.transmitters["gnb"].position.numpy()).ravel()
            r = np.array(self.scene.receivers["nrscope"].position.numpy()).ravel()
            c = np.asarray(pos, float) + [0.0, 0.0, 0.5 * self.car_bodies[k - 1].height]
            power = (self.lam * np.sqrt(self.rcs[k - 1]) / ((4 * np.pi) ** 1.5 * np.linalg.norm(c - g)
                                                            * np.linalg.norm(c - r)) * self.scale) ** 2
        if not power:
            return a, tau, dop, via_car, u, steer
        g = np.array(self.scene.transmitters["gnb"].position.numpy()).ravel()
        r = np.array(self.scene.receivers["nrscope"].position.numpy()).ravel()
        body = self.car_bodies[k - 1]
        pts = body.points(pos, vel)
        dg, dr = pts - g, pts - r
        lg, lr = np.linalg.norm(dg, axis=1), np.linalg.norm(dr, axis=1)
        length = lg + lr
        rate = (dg / lg[:, None] + dr / lr[:, None]) @ np.asarray(vel, float)   # d(length)/dt
        amp = body.weights * np.sqrt(power) * self._antennas(dg, dr)
        if self.ni:   # a real vehicle: its echo fluctuates, and is stronger seen from the side
            amp = amp * self.fading.amp(k, self.n / self.fs)
            if self.ni["vehicle_aspect"]:
                amp = amp * np.sqrt(nonideal.aspect_gain(vel, dg.mean(0), dr.mean(0)))
        a_car = amp * np.exp(-2j * np.pi * length / self.lam)
        # exact distance to every element, relative to element 0's (= r)
        de = np.linalg.norm(pts[:, None, :] - self.elems[None, :, :], axis=2)
        s_car = np.exp(-2j * np.pi * (de - de[:, :1]) / self.lam)
        return (np.concatenate([a, a_car]), np.concatenate([tau, length / C]),
                np.concatenate([dop, -rate / self.lam]), np.r_[via_car, np.full(len(pts), k)],
                np.concatenate([u, dr / lr[:, None]]), np.concatenate([steer, s_car]))

    def _antennas(self, dg, dr):
        """Amplitude factor of the two antennas for echoes off points: sqrt(G_gnb G_rx), the gNB's
        gain towards the point (dg = point - gNB) and NR-Scope's from it (dr = point - NR-Scope)."""
        pat = self.prop["pattern"]
        return np.sqrt(pattern_gain(pat, dg) * pattern_gain(pat, dr))

    def _person(self, k, pos, vel, a, tau, dop, via_car, u, steer):
        """Add target k, a person, as the point scatterers of PERSON_POINTS that both radios see
        (all of them outside a scenario), each with its share of the RCS.

        Its amplitude follows the bistatic radar equation, on the same absolute scale as
        Sionna's paths (a free-space path of length d has |a| = lambda / (4 pi d) with these
        isotropic antennas): received over transmitted power is
            lambda^2 sigma / ((4 pi)^3 R1^2 R2^2),
        R1, R2 the distances gNB - person and person - receiver, so
            |a| = lambda sqrt(sigma) / ((4 pi)^(3/2) R1 R2),
        times self.scale like every other path, sigma being the point's share. The phase is
        the exact path length; the Doppler that of the body (no limb micro-Doppler, constant RCS)."""
        g = np.array(self.scene.transmitters["gnb"].position.numpy()).ravel()
        r = np.array(self.scene.receivers["nrscope"].position.numpy()).ravel()
        vis = getattr(self, "point_vis", {}).get(k - 1, np.ones(len(PERSON_POINTS), bool))
        seen = {h: v for (h, _), v in zip(PERSON_POINTS, vis)}   # visibility by the tested heights
        if self.ni and self.ni["micro_doppler"]:   # legs and arms swing: each part its own position and velocity
            parts = [(off, dv, sh) for off, dv, zv, sh in nonideal.person_points(self.n / self.fs, vel, self.gait[k - 1]) if seen.get(zv, True)]
        else:
            parts = [(np.array([0.0, 0.0, h]), np.zeros(3), w) for h, w in PERSON_POINTS if seen[h]]
        if not parts:
            return a, tau, dop, via_car, u, steer
        pts = np.array([np.asarray(pos, float) + o for o, _, _ in parts])
        vels = np.array([np.asarray(vel, float) + d for _, d, _ in parts])
        share = np.array([w for _, _, w in parts])
        dg, dr = pts - g, pts - r
        lg, lr = np.linalg.norm(dg, axis=1), np.linalg.norm(dr, axis=1)
        length = lg + lr
        rate = ((dg / lg[:, None] + dr / lr[:, None]) * vels).sum(1)
        amp = self.lam * np.sqrt(self.rcs[k - 1] * share) / ((4 * np.pi) ** 1.5 * lg * lr) * self.scale
        amp *= self._antennas(dg, dr)
        if self.ni:
            amp = amp * self.fading.amp(k, self.n / self.fs)
        a_p = amp * np.exp(-2j * np.pi * length / self.lam)
        de = np.linalg.norm(pts[:, None, :] - self.elems[None, :, :], axis=2)
        s_p = np.exp(-2j * np.pi * (de - de[:, :1]) / self.lam)
        return (np.r_[a, a_p], np.r_[tau, length / C], np.r_[dop, -rate / self.lam], np.r_[via_car, np.full(len(pts), k)],
                np.vstack([u, dr / lr[:, None]]), np.vstack([steer, s_p]))

    def _wall_echoes(self, k, pos, vel, t, a, tau, dop, via_car, u, steer):
        """Target k's echoes that also bounce once on a building wall (nonideal target_multipath):
        gNB -> wall -> target -> NR-Scope and gNB -> target -> wall -> NR-Scope, image method on
        the footprints, each leg checked for blockage. Radar equation on the unfolded lengths, the
        wall's loss on top. They are not the target in the ground truth (via_car 0): NR-Scope
        reporting one is a ghost, as on a real capture."""
        g = np.array(self.scene.transmitters["gnb"].position.numpy()).ravel()
        r = np.array(self.scene.receivers["nrscope"].position.numpy()).ravel()
        person = self.kinds[k - 1] == "person"
        pt = np.asarray(pos, float) + [0.0, 0.0, 1.0 if person else 0.5 * self.car_bodies[k - 1].height]
        v = np.asarray(vel, float)
        rcs = self.rcs[k - 1] or 10.0
        gam = 10 ** (-self.ni["wall_loss_db"] / 20) * self.fading.amp(k, t)
        cand = []   # (bounce point, unfolded first leg, second leg, rate, arrival point at the receiver)
        for P, L1 in nonideal.bounce_points(self.walls, g, pt, self.ni["wall_max_m"]):   # on the way to the target
            img_dir = (pt - P) / max(np.linalg.norm(pt - P), 1e-6)
            L2 = np.linalg.norm(pt - r)
            cand.append(([(g, P), (P, pt), (pt, r)], L1, L2, img_dir @ v + (pt - r) / L2 @ v, pt))
        for P, L2 in nonideal.bounce_points(self.walls, r, pt, self.ni["wall_max_m"]):   # on the way back
            L1 = np.linalg.norm(pt - g)
            d_out = (pt - P) / max(np.linalg.norm(pt - P), 1e-6)
            cand.append(([(g, pt), (pt, P), (P, r)], L1, L2, (pt - g) / L1 @ v + d_out @ v, P))
        if not cand:
            return a, tau, dop, via_car, u, steer
        segs = np.array([s for c in cand for s in c[0]])   # [n * 3, 2, 3]: every leg of every candidate
        d = segs[:, 1] - segs[:, 0]
        n = np.linalg.norm(d, axis=1)
        blocked = self._rays(segs[:, 0] + 0.05 * d / n[:, None], d / n[:, None], n - 0.1).reshape(len(cand), 3).any(1)
        A, T, D, U, S = [], [], [], [], []
        for c, b in zip(cand, blocked):
            if b:
                continue
            _, L1, L2, rate, arr = c
            amp = self.lam * np.sqrt(rcs) / ((4 * np.pi) ** 1.5 * L1 * L2) * self.scale * gam
            du = arr - r
            A.append(amp * np.exp(-2j * np.pi * (L1 + L2) / self.lam)); T.append((L1 + L2) / C); D.append(-rate / self.lam)
            U.append(du / np.linalg.norm(du))
            de = np.linalg.norm(arr[None, :] - self.elems, axis=1)
            S.append(np.exp(-2j * np.pi * (de - de[0]) / self.lam))
        if not A:
            return a, tau, dop, via_car, u, steer
        return (np.r_[a, A], np.r_[tau, T], np.r_[dop, D], np.r_[via_car, np.zeros(len(A), int)],
                np.vstack([u, U]), np.vstack([steer, S]))

    def _stitch(self, t, a, tau, dop, via_car):
        """Keep the car's echo phase-continuous across ray traces.

        The car scatters diffusely and each trace samples it at slightly different points, so
        the phases Sionna gives the new car paths do not continue the rotation the old ones
        had: every trace boundary added a phase jump (about -0.66 rad every 5 ms in the street
        canyon), which pulled the echo's Doppler from -110 Hz to -89 Hz. Rotate all the new car
        paths by one common phase so that their frequency response matches the old paths'
        response carried forward to this instant; the Doppler then comes from Sionna alone."""
        prev = getattr(self, "prev_car", None)
        on = via_car > 0
        if prev is None or not on.any() or len(prev[1]) == 0:
            return a
        t0, a0, tau0, dop0 = prev
        f = np.linspace(-self.fs / 2.5, self.fs / 2.5, 257)  # inside the occupied band
        h_old = (a0 * np.exp(2j * np.pi * dop0 * (t - t0)))[None, :] * np.exp(-2j * np.pi * np.outer(f, tau0))
        h_new = a[on][None, :] * np.exp(-2j * np.pi * np.outer(f, tau[on]))
        phi = np.angle(np.vdot(h_new.sum(1), h_old.sum(1)))
        a = a.copy()
        a[on] = a[on] * np.exp(1j * phi)
        return a

    # -- applying the channel ------------------------------------------------------------------
    def _static_filter(self, a, tau, steer, half=32):
        """The summed impulse response of fixed paths, sampled at fs: h[ant, k] = sum_p
        steer[p, ant] a_p sinc(k - tau_p fs), the sinc Kaiser-windowed over +-half samples
        (the band-limited fractional delay the FFT ramp in _apply applies path by path). The
        same scene paths at every trace (a scenario): computed once."""
        key = (len(a), float(np.sum(tau)), complex(np.sum(a)))
        if getattr(self, "_h_key", None) == key:
            return self._h
        d = tau * self.fs
        L = int(np.ceil(d.max())) + half + 1
        if L > self.history // 2:
            raise RuntimeError(f"a path {d.max():.0f} samples long: longer than the channel's history")
        k = np.arange(L)
        H = np.zeros((self.n_ant, L), complex)
        for i in range(0, len(a), 2048):
            x = k[None, :] - d[i:i + 2048, None]
            w = np.where(np.abs(x) <= half, np.kaiser(2 * half + 1, 8.0)[np.clip(np.round(x).astype(int) + half, 0, 2 * half)], 0.0)
            H += (steer[i:i + 2048].T * a[i:i + 2048]) @ (np.sinc(x) * w)
        self._h_key, self._h = key, torch.tensor(H, dtype=torch.complex64, device=self.dev)
        return self._h

    def _apply(self, x, n0):
        """Channel over one stretch of samples that shares one ray trace."""
        n = len(x)
        full = torch.cat([self.buf, x])              # history + new samples
        m = len(full)
        f = torch.fft.fftfreq(m, d=1.0 / self.fs, device=self.dev, dtype=torch.float64)
        p = self.paths
        X = torch.fft.fft(full)
        t = (n0 + torch.arange(n, device=self.dev, dtype=torch.float64)) / self.fs - p["t"]
        y = torch.zeros((self.n_ant, n), dtype=torch.complex64, device=self.dev)
        if self.ni:   # the sample clock's drift: one more common delay for every path
            X = X * torch.exp(-2j * np.pi * f * self.clock.delay(n0 / self.fs, 0.4 * self.history / self.fs)).to(torch.complex64)
        if p.get("h") is not None:   # the fixed paths, summed: one filter per antenna
            y += torch.fft.ifft(X[None, :] * torch.fft.fft(p["h"], n=m))[:, -n:]
        B = 64  # paths at a time: many scenario agents would not fit in GPU memory at once
        for i in range(0, len(p["tau"]), B):
            j = slice(i, i + B)
            # Delay: X(f) * exp(-j 2 pi f tau_p), one row per path.
            ramp = torch.exp(-2j * np.pi * torch.outer(p["tau"][j], f)).to(torch.complex64)
            delayed = torch.fft.ifft(X[None, :] * ramp)[:, -n:]
            # Doppler of each path, from the time of the last ray trace.
            rot = torch.exp(2j * np.pi * torch.outer(p["dop"][j], t)).to(torch.complex64)
            y += p["steer"][:, j] @ (p["a"][j, None] * rot * delayed)  # [n_ant, n]
        self.buf = full[-self.history:]
        return y

    def __call__(self, block):
        """Apply the channel to one block from the gNB (numpy complex64), in order."""
        x = torch.from_numpy(np.array(block, dtype=np.complex64)).to(self.dev)  # copy: ZMQ buffers are read-only
        out = []
        start = 0
        while start < len(x):  # split the block where a new ray trace is due
            if self.n + start >= self.next_update:
                self._trace((self.n + start) / self.fs)
                self.next_update += int(round(self.update_s * self.fs))
            stop = min(len(x), self.next_update - self.n)
            out.append(self._apply(x[start:stop], self.n + start))
            start = stop
        y = torch.cat(out, dim=1)
        if self.ni:   # carrier offset, its drift and phase noise of NR-Scope's oscillator against the gNB's
            ph = torch.from_numpy(self.clock.phase(self.n, len(x))).to(self.dev)
            y = y * torch.exp(1j * ph).to(torch.complex64)[None, :]

        # Noise: snr_db below the signal power of the strongest (unit-gain) path.
        pw = float((x.abs() ** 2).mean())
        if pw > 0:
            self.sig_power = pw if self.sig_power is None else 0.99 * self.sig_power + 0.01 * pw
        if self.sig_power and not (getattr(self, "noise_only_active", False) and pw == 0):
            std = np.sqrt(self.sig_power / 10 ** (self.snr_db / 10) / 2)
            y = y + std * torch.complex(torch.randn(y.shape, generator=self.rng, device=self.dev),
                                        torch.randn(y.shape, generator=self.rng, device=self.dev))
        self.n += len(x)
        y = y.cpu().numpy()
        return y if self.n_ant > 1 else y[0]

    def close(self):
        if self.gt_file:
            self.gt_file.close()
            self.cars_file.close()
