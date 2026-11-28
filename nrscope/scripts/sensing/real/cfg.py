"""Configuration of the real-data chain: defaults (the setup of record for real recordings, 2026-11-26), a JSON config
file merged on top, dotted --set overrides, and the receiver geometry from the recording's session.json + cart profile.

  python3 cfg.py [--config FILE] [--set key=value ...]      prints the merged configuration and the geometry"""
import argparse, copy, json, os
import numpy as np

H = os.path.expanduser("~")
DEFAULTS = {
    "recording": {   # set per recording: configs/<recording>.json (real_new.py writes it)
        "chunks": "",          # decoded chunks (the Real tab's Decode)
        "lidar": "",           # the lidar recording folder (movers.json, grid/b<k>.json)
        "session": "",         # the test's session.json (gNB, lidar pose)
        "cart": H + "/GIOVANNI/real/cart.json",   # NR-Scope relative to the lidar
        "t_rec": "",           # UTC time the IQ start is counted from ("YYYY-mm-dd HH:MM:SS.fff")
        "iq_offset_s": 0.0,    # IQ start after t_rec (0 with the wall-clock recorder)
        "truth": "",           # ground-truth CSV (unix_s,id,lat,lon,cls); empty: built from the lidar movers
        "site": "cubotron_epfl",   # sim/sites/<site>: buildings for the line-of-sight map and wall ghosts
    },
    "streams": {
        "trs": {"est": "", "so_light_env": {"USE_DMRS": "0"}, "omp_threshold_db": 14.0, "target_dv": 0.3},
        "dmrs": {"est": "", "so_light_env": {"NO_TRS_OUT": "1", "MODE": "tab", "BA": "1", "BA_PH": "wide"},
                 "omp_threshold_db": 18.0, "target_dv": 0.1},
    },
    "sensing": {
        "base_config": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so/c_omp_stk0_paper.yaml",
        "aoa_cal_deg": "auto",       # or [0, a1, a2, a3]
        "aoa_cal_chunks": "all",     # chunks the automatic calibration uses: "all" or "3,10,13"
        "aoa_cal_max_power_dev_db": 3.0,   # a chunk whose per-chain direct-path power pattern departs more is dropped
        "aoa_cal_max_spread_deg": 15.0,    # p90 spread over the kept chunks above which the calibration is refused
        "so_light_common_env": {"MINPRB": "48", "USE_CSI": "0"},
        "workers": 10,
    },
    "los": {"fov_deg": 60.0, "step_m": 1.0, "span_m": 80.0, "target_z_m": 1.0,
            "grid_file": ""},   # a precomputed line-of-sight map (pd_grids.npz with this geometry key) instead of computing one
    "tracker": {
        "params": {"pd_learn": 1, "pd_learn_tau": 30.0, "pd_learn_prior": 5.0, "pd_learn_deg": 10.0, "pd_learn_m": 10.0},
        "setups": ["1", "2", "fused"],
        "selfcal": False,
        "micro_doppler_ms": 0.0,
        "t_window_unix": None,   # [t0, t1] unix s: track and score only inside (e.g. while the walker is in the receiver FOV); None: everything
    },
    "tbd": {   # track-before-detect on range x speed x angle cubes (cube.py, tbd.py)
        "r_min_m": 4.0, "r_max_m": 80.0, "r_step_m": 1.5,           # excess path; below r_min the direct path leaks
        "v_min_ms": 0.1, "v_max_ms": 4.0, "v_step_ms": 0.17,        # speed, NR-Scope's convention (f = 2 v / lam)
        "angle_max_deg": 60.0, "angle_step_deg": 4.0,
        "whiten": "none",                                            # "map" / "chunk": angle stage whitened against the noise's spatial covariance
        "window_s": 0.0,                                             # coherent window; 0: the map's own (NR-Scope's) span
        "streams": ["trs", "dmrs"],                                  # evidence summed in one recursion
        "per_stream": {"trs": {"snr": 3.0, "weight": 1.0}, "dmrs": {"snr": 3.0, "weight": 1.0}},
        "bg_tau_s": 10.0, "survival": 0.95, "p_r": 0.1, "p_v": 0.15, "p_a": 0.1, "gap_s": 1.0,
        "output": "paths",          # "paths": tracks from an offline Viterbi; "detections": peaks for the paper tracker
        "threshold": 8.0, "map_threshold": 40.0, "max_dets": 6, "nms_r_cells": 5, "nms_a_cells": 5,
        "max_paths": 8, "path_v_max_ms": 2.5, "move_cost": [0.5, 0.5, 0.5], "on_cost": 0.0, "birth_cost": 10.0, "death_cost": 10.0,
        "supp": 2, "min_maps": 10, "min_score": 20.0,
    },
    "replay": {"setup": "fused"},   # which setup's tracks the Real tab replays
    "paths": {
        "paper": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/tracker", "python": H + "/GIOVANNI/.venv/bin/python",
        "sensing_offline": H + "/GIOVANNI/NR-Scope/build/nrscope/src/tests/sensing_offline",
        "so_light": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "aoa_cal": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/aoa_cal.py",
        "data": "/mnt/mydrive/field/real",   # per recording: estimates, captures, results
    },
}


def merge(a, b):
    out = copy.deepcopy(a)
    for k, v in b.items():
        out[k] = merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def set_key(cfg, dotted, value):
    keys = dotted.split("."); d = cfg
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    try:
        d[keys[-1]] = json.loads(value)
    except ValueError:
        d[keys[-1]] = value


def load(config=None, sets=()):
    cfg = merge(DEFAULTS, json.load(open(config))) if config else copy.deepcopy(DEFAULTS)
    for kv in sets:
        k, v = kv.split("=", 1); set_key(cfg, k, v)
    return cfg


def recording_dir(cfg):
    """per-recording data folder: <data>/<lidar recording name>"""
    return os.path.join(cfg["paths"]["data"], os.path.basename(cfg["recording"]["lidar"].rstrip("/")))


def est_dir(cfg, stream):
    return cfg["streams"][stream]["est"] or os.path.join(recording_dir(cfg), "est_" + stream)


def geometry(cfg):
    """origin = lidar; NR-Scope = lidar + cart offset in the lidar's axes (heights above the ground);
    boresight from east, counter-clockwise (the tracker's convention)."""
    rec = cfg["recording"]; S = json.load(open(rec["session"])); L = S["lidar"]; cart = json.load(open(rec["cart"]))
    KY = 6378137 * np.pi / 180; KX = KY * np.cos(np.radians(L["lat"])); h = np.radians(L["heading"])
    o = cart["off"]; e = o[0] * np.cos(h) + o[1] * np.sin(h); n = -o[0] * np.sin(h) + o[1] * np.cos(h)
    return dict(origin=(L["lat"], L["lon"]), gnb=(S["gnb"]["lat"], S["gnb"]["lon"]), gnb_h=float(S["gnb"]["h"]),
                rx=(L["lat"] + n / KY, L["lon"] + e / KX), rx_h=float(L["h"] + o[2]),
                boresight=float((90 - (L["heading"] + cart.get("yaw", 0))) % 360))


def local_xy(geo, lat, lon):
    """metres east / north of the origin (the import's frame)"""
    KY = 6378137 * np.pi / 180; KX = KY * np.cos(np.radians(geo["origin"][0]))
    return np.array([(lon - geo["origin"][1]) * KX, (lat - geo["origin"][0]) * KY])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[])
    a = ap.parse_args(); c = load(a.config, a.set)
    print(json.dumps(c, indent=1)); print(json.dumps(geometry(c), indent=1))
