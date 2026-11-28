"""Sensing of one stream: sensing_offline on every chunk's estimates with the base config + the stream's
omp_threshold_db / target_dv + the angle calibration. Maps per chunk (b<k>.map2d.csv) and the config used (sensing.yaml).

  python3 sense.py [--config FILE] [--set ...] trs|dmrs OUT_DIR"""
import argparse, glob, os, sys, concurrent.futures as CF
import subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C, calib


def config_yaml(cfg, stream, cal, out):
    st = cfg["streams"][stream]
    L = [l for l in open(cfg["sensing"]["base_config"]) if not l.strip().startswith(("target_dv:", "omp_threshold_db:", "aoa_cal_deg:"))]
    i = next(j for j, l in enumerate(L) if l.startswith("sensing:"))
    L[i + 1:i + 1] = ["  omp_threshold_db: %g\n" % st["omp_threshold_db"], "  target_dv: %g\n" % st["target_dv"],
                      "  aoa_cal_deg: [%s]\n" % ", ".join("%.2f" % x for x in cal)]
    open(out, "w").write("".join(L))
    return out


def run(cfg, stream, cal, out_dir):
    os.makedirs(out_dir, exist_ok=True); y = config_yaml(cfg, stream, cal, os.path.join(out_dir, "sensing.yaml"))
    so = cfg["paths"]["sensing_offline"]
    def one(b):
        d = os.path.join(out_dir, os.path.basename(b)[:-4] + ".map2d.csv")
        for _ in range(2):
            if subprocess.run([so, y, b, "--dump", d], capture_output=True).returncode == 0 and os.path.exists(d):
                return True
        return False
    with CF.ThreadPoolExecutor(cfg["sensing"]["workers"]) as ex:
        ok = list(ex.map(one, glob.glob(os.path.join(C.est_dir(cfg, stream), "b*.bin"))))
    for f in glob.glob(os.path.join(out_dir, "*.snap.csv")):
        os.remove(f)
    return y, sum(ok), len(ok)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("stream", choices=["trs", "dmrs"]); ap.add_argument("out"); ap.add_argument("--config"); ap.add_argument("--set", nargs="+", action="extend", default=[])
    a = ap.parse_args(); c = C.load(a.config, a.set); y, ok, n = run(c, a.stream, calib.get(c)[0], a.out); print("%s: %d/%d chunks sensed (%s)" % (a.out, ok, n, y))
