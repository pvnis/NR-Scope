"""Angle calibration of the receive chains (aoa_cal_deg) from the TRS direct path of the recording (/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/aoa_cal.py:
chunks with an abnormal per-chain direct-path power pattern dropped, circular median of the rest, p90 spread check),
cached per recording and chunk selection. A configured list is used as is.

  python3 calib.py [--config FILE] [--set ...]"""
import argparse, glob, json, os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C


def get(cfg):
    v = cfg["sensing"]["aoa_cal_deg"]
    if v != "auto":
        return [float(x) for x in v], "configured"
    ch = cfg["recording"]["chunks"]; sel = cfg["sensing"]["aoa_cal_chunks"]
    if sel == "all":
        sel = ",".join(sorted((re.search(r"b(\d+)$", d).group(1) for d in glob.glob(ch + "/b*") if re.search(r"b\d+$", d)), key=int))
    d = C.recording_dir(cfg); os.makedirs(d, exist_ok=True)
    tag = sel.replace(",", "-") if len(sel) < 40 else "n%d" % len(sel.split(","))
    cache = os.path.join(d, "aoa_cal_%s_p%g_s%g.json" % (tag, cfg["sensing"]["aoa_cal_max_power_dev_db"], cfg["sensing"]["aoa_cal_max_spread_deg"]))
    if os.path.exists(cache):
        return json.load(open(cache))["aoa_cal_deg"], cache
    env = dict(os.environ, AOA_MAX_POWER_DEV_DB=str(cfg["sensing"]["aoa_cal_max_power_dev_db"]), AOA_MAX_SPREAD_DEG=str(cfg["sensing"]["aoa_cal_max_spread_deg"]))
    cart = ",".join(str(x) for x in json.load(open(cfg["recording"]["cart"]))["off"])
    r = subprocess.run(["python3", cfg["paths"]["aoa_cal"], ch, cfg["recording"]["session"], sel, cart], env=env, capture_output=True, text=True,
                       cwd=os.path.dirname(cfg["paths"]["so_light"]))
    line = [l for l in r.stdout.splitlines() if l.startswith("aoa_cal_deg:")]
    if r.returncode != 0 or not line:
        raise RuntimeError("angle calibration failed or unstable:\n" + (r.stdout + r.stderr)[-1500:])
    cal = json.loads(line[0].split(":", 1)[1])
    json.dump(dict(aoa_cal_deg=cal, chunks=sel, log=r.stdout[-4000:]), open(cache, "w"), indent=1)
    return cal, cache


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[])
    a = ap.parse_args(); cal, src = get(C.load(a.config, a.set)); print("aoa_cal_deg %s (%s)" % (cal, src))
