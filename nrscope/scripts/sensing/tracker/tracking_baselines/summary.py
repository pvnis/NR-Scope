"""Tables and a figure of the tracker comparison of one run.

python tracking_baselines/summary.py <run>
Per tracker and setup: its default settings, and the best grid point by GOSPA (tune.py).
Writes results/<run>/summary.md, summary.csv and summary.png."""
import csv
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from trackers import REGISTRY, SINGLE_SENSOR  # noqa: E402

SETUPS = ["1", "2", "fused"]
SETUP_NAME = {"1": "rx1 alone", "2": "rx2 alone", "fused": "rx1 + rx2 fused"}
COLORS = {"fused": "#2a78d6", "1": "#eb6834", "2": "#1baf7a"}   # as plots/common.py


def cov(r, cls):
    for c in r.get("per_class") or []:
        if c["cls"] == cls:
            return c["coverage"]
    return None


def row(r):
    p = r["position"]
    return dict(gospa=p.get("gospa_m"), loc=p.get("gospa_loc_m"), missed=p.get("missed_per_map"),
                false=p.get("false_per_map"), n_tracks=p.get("n_tracks"), false_tracks=p.get("false_tracks"),
                people=cov(r, "person"), cars=cov(r, "car"),
                switches=sum(c.get("switches") or 0 for c in r.get("per_class") or []), params=r["params"])


def main():
    run = sys.argv[1]
    d = os.path.join(HERE, "results", run)
    res = {}
    for f in glob.glob(os.path.join(d, "*.json")):
        r = json.load(open(f))
        rest = os.path.basename(f)[len(r["tracker"]) + 2:-5]   # "<setup>", "grid<k>__<setup>" or "<tag>__<setup>"
        kind = "default" if rest in SETUPS else "grid" if rest.startswith("grid") else None
        if kind is None:
            continue
        res.setdefault((r["tracker"], r["setup"]), {}).setdefault(kind, []).append(r)
    rows = []
    for name in REGISTRY:
        for s in SETUPS:
            e = res.get((name, s))
            if not e or (s == "fused" and name in SINGLE_SENSOR):
                continue
            dflt = row(e["default"][0]) if e.get("default") else None
            grid = [row(r) for r in e.get("grid", []) if r["position"].get("gospa_m") is not None]
            best = min(grid, key=lambda x: x["gospa"]) if grid else None
            rows.append(dict(tracker=name, setup=s, default=dflt, best=best))
    f = lambda v, n=2: "-" if v is None else (f"{100 * v:.0f}%" if n == "%" else f"{v:.{n}f}")
    lines = [f"# Tracker comparison: {run}", "",
             "Same detections, same scoring (sim/analyze.py + track.localize_and_score). Best: best GOSPA over",
             "the tracker's 3x3 grid (tune.py), tuned on this run. GOSPA: c = 5 m, p = 2, per map, mean.",
             "Coverage: share of maps where the target has a track on it.", ""]
    for s in SETUPS:
        lines += [f"## {SETUP_NAME[s]}", "",
                  "| tracker | GOSPA default (m) | GOSPA best (m) | missed/map | false/map | false tracks / tracks | people covered | cars covered | ID switches | best params |",
                  "|---|---|---|---|---|---|---|---|---|---|"]
        for r in sorted([r for r in rows if r["setup"] == s], key=lambda r: (r["best"] or r["default"])["gospa"]):
            b = r["best"] or r["default"]
            lines.append(f"| {r['tracker']} | {f(r['default'] and r['default']['gospa'])} | {f(b['gospa'])} | {f(b['missed'])} | "
                         f"{f(b['false'])} | {b['false_tracks']} / {b['n_tracks']} | {f(b['people'], '%')} | {f(b['cars'], '%')} | "
                         f"{b['switches']} | {json.dumps(b['params'])} |")
        lines.append("")
    open(os.path.join(d, "summary.md"), "w").write("\n".join(lines))
    with open(os.path.join(d, "summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tracker", "setup", "gospa_default", "gospa_best", "missed_per_map", "false_per_map",
                    "false_tracks", "n_tracks", "people_coverage", "car_coverage", "switches", "best_params"])
        for r in rows:
            b = r["best"] or r["default"]
            w.writerow([r["tracker"], r["setup"], r["default"] and r["default"]["gospa"], b["gospa"], b["missed"], b["false"],
                        b["false_tracks"], b["n_tracks"], b["people"], b["cars"], b["switches"], json.dumps(b["params"])])
    plot(rows, os.path.join(d, "summary.png"), run)
    print("\n".join(lines))


def plot(rows, out, run):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = [n for n in REGISTRY if any(r["tracker"] == n for r in rows)]
    names.sort(key=lambda n: min((r["best"] or r["default"])["gospa"] for r in rows if r["tracker"] == n))
    metrics = [("gospa", "GOSPA (m), lower is better"), ("missed", "missed targets per map"), ("false", "false tracks per map")]
    fig, axs = plt.subplots(1, 3, figsize=(11, 0.42 * len(names) + 1.4), sharey=True)
    y = np.arange(len(names))
    h = 0.26
    for ax, (k, title) in zip(axs, metrics):
        for i, s in enumerate(SETUPS):
            v = [next(((r["best"] or r["default"])[k] for r in rows if r["tracker"] == n and r["setup"] == s), np.nan) for n in names]
            ax.barh(y + (i - 1) * h, v, height=h - 0.03, color=COLORS[s], label=SETUP_NAME[s])
        ax.set_title(title, fontsize=9, loc="left")
        ax.grid(axis="x", color="#e4e3df", lw=0.6)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axs[0].set_yticks(y, [n.replace("_", " ") for n in names], fontsize=8)
    axs[0].invert_yaxis()
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, fontsize=8, frameon=False, loc="lower center", ncol=3)
    fig.suptitle(f"Trackers on the same detections (best of each tracker's 3x3 grid) - {run}", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out, dpi=150)


if __name__ == "__main__":
    main()
