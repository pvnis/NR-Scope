"""Tracker comparison averaged over all long runs: mean and std over runs of each tracker's best
grid point (summary.csv of each run, see summary.py).

python tracking_baselines/summary_all.py
Writes results/summary_all.png (GOSPA, missed/map, false/map), results/coverage_all.png
(people, cars covered) and results/summary_all.csv."""
import csv
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from summary import COLORS, SETUP_NAME, SETUPS  # noqa: E402
from trackers import REGISTRY  # noqa: E402

KEYS = dict(gospa="gospa_best", missed="missed_per_map", false="false_per_map",
            people="people_coverage", cars="car_coverage")


# runs whose rx1 angles are broken (known rx1 AoA bug: dens3_seed2 84 deg off, dens6_seed3 56% > 20 deg off)
EXCLUDE = ("_long_dens3_seed2", "_long_dens6_seed3")


def load():
    files = sorted(f for f in glob.glob(os.path.join(HERE, "results", "*_long_*", "summary.csv"))
                   if not any(x in f for x in EXCLUDE))
    vals = {}   # (tracker, setup, metric) -> values over runs
    for f in files:
        for r in csv.DictReader(open(f)):
            for k, col in KEYS.items():
                if r[col] not in ("", "None"):
                    vals.setdefault((r["tracker"], r["setup"], k), []).append(float(r[col]))
    return vals, len(files)


def plot(vals, n, metrics, out, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = [t for t in REGISTRY if any((t, s, "gospa") in vals for s in SETUPS)]
    names.sort(key=lambda t: np.mean(vals.get((t, "fused", "gospa"), vals[(t, "1", "gospa")])))
    fig, axs = plt.subplots(1, len(metrics), figsize=(3.7 * len(metrics), 0.42 * len(names) + 1.6), sharey=True)
    y = np.arange(len(names))
    h = 0.26
    for ax, (k, label, pct) in zip(axs, metrics):
        sc = 100 if pct else 1
        for i, s in enumerate(SETUPS):
            m = [sc * np.mean(vals[(t, s, k)]) if (t, s, k) in vals else np.nan for t in names]
            e = [sc * np.std(vals[(t, s, k)]) if (t, s, k) in vals else np.nan for t in names]
            ax.barh(y + (i - 1) * h, m, height=h - 0.03, xerr=e, color=COLORS[s], label=SETUP_NAME[s],
                    error_kw=dict(ecolor="#555", lw=0.8, capsize=2))
            for yy, mm, ee in zip(y + (i - 1) * h, m, e):
                if not np.isnan(mm):
                    ax.annotate(f"{mm:.0f}%" if pct else f"{mm:.2f}", (mm + ee, yy), xytext=(4, 0),
                                textcoords="offset points", va="center", fontsize=5.5, color="#333")
        ax.set_xlim(left=0)
        ax.set_title(label, fontsize=9, loc="left")
        ax.grid(axis="x", color="#e4e3df", lw=0.6)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axs[0].set_yticks(y, [t.replace("_", " ") for t in names], fontsize=8)
    axs[0].invert_yaxis()
    hh, ll = axs[0].get_legend_handles_labels()
    fig.legend(hh, ll, fontsize=8, frameon=False, loc="lower center", ncol=3)
    fig.suptitle(f"{title} - mean ± std over {n} runs (best of each tracker's 3x3 grid per run)",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out, dpi=150)


def main():
    vals, n = load()
    d = os.path.join(HERE, "results")
    plot(vals, n, [("gospa", "GOSPA (m), lower is better", False), ("missed", "missed targets per map", False),
                   ("false", "false tracks per map", False)], os.path.join(d, "summary_all.png"), "Trackers")
    plot(vals, n, [("people", "people covered (%), higher is better", True),
                   ("cars", "cars covered (%), higher is better", True)], os.path.join(d, "coverage_all.png"), "Coverage")
    with open(os.path.join(d, "summary_all.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tracker", "setup", "metric", "mean", "std", "n_runs"])
        for (t, s, k), v in sorted(vals.items()):
            w.writerow([t, s, k, np.mean(v), np.std(v), len(v)])
    print(f"{n} runs")


if __name__ == "__main__":
    main()
