#!/usr/bin/env python3
"""Prepare a new recording for real_fuse.py and run it: one command after the Real tab's Store on drive, Decode and
Occupancy grid, the lidar placed (session.json with the gNB) and Find movers.

  python3 real_new.py LIDAR_RECORDING_DIR [--iq-offset S --t-rec "YYYY-mm-dd HH:MM:SS.fff"] [--set ...] [--name NAME] [--check]

Checks the inputs (decoded chunks, occupancy grid per chunk, lidar movers, session with gNB and lidar pose, cart).
IQ start time: t0_wall_utc_s from the recording's meta.txt (rec_direct since it logs the wall clock); without it
--t-rec and --iq-offset are required (walk1: the offset had to be fitted). Writes configs/<recording>.json and runs
real_fuse.py with it (unless --check)."""
import argparse, datetime, glob, json, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cfg as C


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lidar"); ap.add_argument("--iq-offset", type=float); ap.add_argument("--t-rec"); ap.add_argument("--site", default="cubotron_epfl")
    ap.add_argument("--set", nargs="+", action="extend", default=[]); ap.add_argument("--name", default="fused"); ap.add_argument("--check", action="store_true")
    a = ap.parse_args(); L = os.path.abspath(a.lidar.rstrip("/")); problems = []
    session = os.path.join(os.path.dirname(L), "session.json"); cart = os.path.join(C.H, "GIOVANNI", "real", "cart.json")
    cd = open(os.path.join(L, "chunks_dir")).read().strip() if os.path.exists(os.path.join(L, "chunks_dir")) else ""
    chunks = sorted(glob.glob(cd + "/b*")) if cd else []
    if not chunks:
        problems.append("no decoded chunks (Real tab: Store on drive, then Decode)")
    grids = [k for k in (re.search(r"b(\d+)$", c) for c in chunks) if k and not os.path.exists(os.path.join(L, "grid", "b%s.json" % k.group(1)))]
    if chunks and grids:
        problems.append("%d chunks without an occupancy grid (Real tab: Occupancy grid)" % len(grids))
    if not os.path.exists(os.path.join(L, "movers.json")):
        problems.append("no lidar movers (Real tab: Find movers)")
    S = json.load(open(session)) if os.path.exists(session) else {}
    if not (S.get("gnb", {}).get("lat") and S.get("lidar", {}).get("heading") is not None):
        problems.append("session.json lacks the gNB or the lidar pose (Real tab: place the lidar, set the gNB)")
    meta = os.path.join(re.sub(r"_chunks$", "", cd), "meta.txt") if cd else ""
    wall = None
    if meta and os.path.exists(meta):
        m = re.search(r"^t0_wall_utc_s\s+([\d.]+)", open(meta).read(), re.M)
        wall = float(m.group(1)) if m else None
    if wall is not None:
        t_rec = datetime.datetime.fromtimestamp(wall, datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]; off = 0.0; src = "meta.txt t0_wall_utc_s"
    elif a.t_rec and a.iq_offset is not None:
        t_rec, off, src = a.t_rec, a.iq_offset, "given"
    else:
        t_rec, off, src = None, None, None
        problems.append("no IQ start time: meta.txt has no t0_wall_utc_s (record with the wall-clock rec_direct) -> give --t-rec and --iq-offset")
    print("recording %s\n  chunks %s (%d)\n  IQ start %s + %s s (%s)" % (L, cd, len(chunks), t_rec, off, src))
    for p in problems:
        print("  MISSING:", p)
    if problems:
        sys.exit(1)
    name = os.path.basename(L); D = os.path.join(C.DEFAULTS["paths"]["data"], name)
    conf = dict(recording=dict(chunks=cd, lidar=L, session=session, cart=cart, t_rec=t_rec, iq_offset_s=off, truth="", site=a.site))
    os.makedirs(os.path.join(HERE, "configs"), exist_ok=True); cf = os.path.join(HERE, "configs", name + ".json")
    json.dump(conf, open(cf, "w"), indent=1); print("  config ->", cf, "(data in %s)" % D)
    if not a.check:
        subprocess.run([sys.executable, os.path.join(HERE, "real_fuse.py"), "--config", cf, "--name", a.name, "--set", *a.set] if a.set else
                       [sys.executable, os.path.join(HERE, "real_fuse.py"), "--config", cf, "--name", a.name])


if __name__ == "__main__":
    main()
