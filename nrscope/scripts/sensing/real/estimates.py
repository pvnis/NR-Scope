"""Per-stream channel estimates of a recording (so_light, one .bin per chunk), built once and kept:
the TRS stream and the DM-RS-only stream, each with its so_light options from the configuration.

  python3 estimates.py [--config FILE] [--set ...] trs|dmrs"""
import argparse, glob, os, re, subprocess, sys, concurrent.futures as CF
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C


def build(cfg, stream):
    out = C.est_dir(cfg, stream); os.makedirs(out, exist_ok=True)
    env = dict(os.environ, **cfg["sensing"]["so_light_common_env"], **cfg["streams"][stream]["so_light_env"])
    so = cfg["paths"]["so_light"]; ch = cfg["recording"]["chunks"]
    dirs = sorted((d for d in glob.glob(ch + "/b*") if re.search(r"b\d+$", d)), key=lambda p: int(re.search(r"b(\d+)$", p).group(1)))
    def one(d):
        b = os.path.join(out, os.path.basename(d) + ".bin")
        if os.path.exists(b) and os.path.getsize(b) > 0:
            return True
        r = subprocess.run(["nice", "python3", so, d, b], env=env, cwd=os.path.dirname(so), capture_output=True, text=True)
        if r.returncode != 0:
            open(b + ".err", "w").write(r.stdout[-2000:] + r.stderr[-2000:])
        return r.returncode == 0
    with CF.ThreadPoolExecutor(cfg["sensing"]["workers"]) as ex:
        ok = list(ex.map(one, dirs))
    return out, sum(ok), len(ok)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("stream", choices=["trs", "dmrs"]); ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[])
    a = ap.parse_args(); out, ok, n = build(C.load(a.config, a.set), a.stream); print("%s: %d/%d chunks" % (out, ok, n))
