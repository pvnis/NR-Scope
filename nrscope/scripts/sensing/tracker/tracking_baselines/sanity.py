"""Sanity test of the tracker implementations: an easy synthetic scene every tracker is designed
for, so every correct implementation must do well. Point targets (one detection each), constant
velocity, P_D 0.9 inside the field of view, Poisson clutter (1 per map, uniform), measurement noise
as the trackers assume; the run's two receivers and gNB, maps every 0.153 s alternating receivers.
Scored with the same GOSPA (c = 5 m) per map, after a 2 s warm-up.

python tracking_baselines/sanity.py [tracker ...] [--run RUN] [--setup 1|2|fused]"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "sim"))
sys.path.insert(0, HERE)
import track  # noqa: E402
from trackers import REGISTRY  # noqa: E402

RB, SC = 2.44, 0.32


def scene(run, setup, seed=0, dur=30.0, n=6):
    s = json.load(open(os.path.join(os.path.dirname(HERE), "logs", run, "settings.json")))
    geos = []
    for rx, bore in ((s["rx"], s["boresight_deg"]), (s["rx2"], s["boresight2_deg"])):
        T, R = track.as3(s["gnb"]), track.as3(rx)
        geos.append((T, R, float(np.linalg.norm(R - T)), np.deg2rad(bore)))
    use = {"1": [0], "2": [1], "fused": [0, 1]}[setup]
    rng = np.random.default_rng(seed)
    # targets: start in the area both receivers look at, walk or drive straight
    c = (np.array(s["rx"][:2]) + np.array(s["rx2"][:2])) / 2
    tg = []
    for k in range(n):
        p0 = c + rng.uniform(-12, 12, 2)
        sp = rng.choice([1.2, 1.2, 1.2, 6.0])
        a = rng.uniform(0, 2 * np.pi)
        tg.append((p0 - np.array([np.cos(a), np.sin(a)]) * sp * dur / 2, sp * np.array([np.cos(a), np.sin(a)])))
    R3 = np.diag([(track.EKF_R_RANGE_BINS * RB) ** 2, (track.EKF_R_SPEED_CELLS * SC) ** 2, np.deg2rad(3.0) ** 2])
    frames, truth = [], []
    t, k = 0.0, 0
    while t < dur:
        g = geos[use[k % len(use)]]
        T, R, L, bore = g
        meas, tru = [], []
        for p0, v in tg:
            x = np.r_[p0 + v * t, v]
            a = np.arctan2(x[1] - R[1], x[0] - R[0]) - bore
            inside = abs((a + np.pi) % (2 * np.pi) - np.pi) <= np.deg2rad(track.FOV_DEG)
            z = track.bistatic_meas(x, T, R, L, bore, track.TARGET_Z_MID)
            if not inside or z[0] > 100:
                continue
            tru.append(x[:2])
            if rng.random() < 0.9:
                meas.append(z + rng.multivariate_normal(np.zeros(3), R3))
        for _ in range(rng.poisson(1.0)):
            meas.append(np.array([rng.uniform(0, 100), rng.uniform(-3, 3), rng.uniform(-1, 1)]))
        md = []
        for z in meas:
            p = track.localize(z[0], np.rad2deg(z[2]), T, R, np.rad2deg(bore))
            md.append(dict(e=float(z[0]), v=float(z[1]), theta=float(z[2]), xy=None if p is None else np.array(p[:2]), p_db=None))
        frames.append((t, md, RB, SC, g, use[k % len(use)] + 1))
        truth.append(tru)
        t += 0.153 / len(use)
        k += 1
    return frames, truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trackers", nargs="*")
    ap.add_argument("--run", default="20261027_131056_bc_env_new_pos_2_receivers_dens6_seed1")
    ap.add_argument("--setup", default="fused")
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args()
    print(f"| tracker | GOSPA (m) | missed/map | false/map |  ({a.setup}, {a.seeds} seeds)\n|---|---|---|---|")
    for name in a.trackers or list(REGISTRY):
        res = []
        for seed in range(a.seeds):
            frames, truth = scene(a.run, a.setup, seed)
            fn = REGISTRY[name]
            if fn is None:
                cfg = track.ekf_cfg({})
                out = track.run_ekf_tracker(frames, None, cfg)
            else:
                if name.startswith("P1") and a.setup == "fused":
                    continue
                out = fn(frames)
            for (t, *_), est, tru in zip(frames, out, truth):
                if t >= 2.0:
                    res.append(track.xy_gospa([(q[1], q[2]) for q in est], tru))
        if res:
            r = np.array(res)
            print(f"| {name} | {r[:, 0].mean():.2f} | {r[:, 2].mean():.2f} | {r[:, 3].mean():.2f} |", flush=True)


if __name__ == "__main__":
    main()
