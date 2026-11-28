"""Same tuning budget for every tracker: a 3 x 3 grid over its two most influential parameters,
per setup; the best by mean GOSPA is kept (summary.py). Tuned on the scored run itself, so every
tracker's result is optimistic in the same way.

python tracking_baselines/tune.py <run> [tracker ...] [--jobs 6]
Writes results/<run>/<tracker>__grid<k>__<setup>.json (k: the grid point, params inside)."""
import argparse
import itertools
import json
import os
from concurrent.futures import ProcessPoolExecutor

import harness
from trackers import REGISTRY, SINGLE_SENSOR

GRIDS = {
    "OURS_full": dict(p_d=[0.5, 0.7, 0.9], r_confirm=[0.8, 0.9, 0.95]),
    "C1_GNN-EKF": dict(confirm=[2, 3, 4], max_miss=[2, 3, 6]),
    "C2_GM-PHD-VoMa2006": dict(p_d=[0.5, 0.7, 0.9], w_birth=[0.01, 0.05, 0.2]),
    "C3_JIPDA-Musicki2004": dict(p_d=[0.5, 0.7, 0.9], r_confirm=[0.5, 0.8, 0.95]),
    "C4_MHT-CoxHingorani1996": dict(p_d=[0.5, 0.7, 0.9], s_conf=[2.3, 4.6, 6.9]),
    "C5_BP-Meyer2017": dict(p_d=[0.5, 0.7, 0.9], mu_n=[0.01, 0.05, 0.2]),
    "C6_GGIW-PHD-Granstrom2012": dict(p_d=[0.5, 0.7, 0.9], w_birth=[0.01, 0.05, 0.2]),
    "P1_Bauhofer2026_GM-PHD-range-doppler": dict(p_d=[0.5, 0.8, 0.9], w_birth=[1e-6, 1e-3, 1e-1]),
    "P2_BSense2026_core": dict(g_s=[3.0, 6.0, 15.0], L=[3, 4, 6]),
    "P3_CellSense2026_window-Hungarian": dict(gate_m=[3.0, 6.0, 10.0], report=[1, 2, 3]),
    "P4_Gunes2025_plotfusion-GNN": dict(max_dt=[0.05, 0.08, 0.15], gate_p=[0.99, 0.999, 0.9999]),
}


def job(a):
    run, name, setup, k, params = a
    r = harness.score(run, name, setup, params)
    out = os.path.join(harness.HERE, "results", run)
    json.dump(r, open(os.path.join(out, f"{name}__grid{k}__{setup}.json"), "w"), indent=1, default=float)
    p = r["position"]
    return f"{name} {setup} {params} GOSPA {p.get('gospa_m')}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("trackers", nargs="*")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--setups", default="1,2,fused")   # a single-receiver run: --setups 1   # a long run's scoring takes ~4.6 GB: 4 fit in 30 GB
    a = ap.parse_args()
    os.makedirs(os.path.join(harness.HERE, "results", a.run), exist_ok=True)
    jobs = []
    for name in a.trackers or list(REGISTRY):
        g = GRIDS[name]
        for setup in [x for x in ([1, 2] if name in SINGLE_SENSOR else [1, 2, "fused"]) if str(x) in a.setups.split(",")]:
            for k, vals in enumerate(itertools.product(*g.values())):
                done = os.path.join(harness.HERE, "results", a.run, f"{name}__grid{k}__{setup}.json")
                if not os.path.exists(done):   # resume: skip finished grid points
                    jobs.append((a.run, name, setup, k, dict(zip(g, vals))))
    with ProcessPoolExecutor(a.jobs) as ex:
        for line in ex.map(job, jobs):
            print(line, flush=True)


if __name__ == "__main__":
    main()
