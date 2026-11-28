# Tracking baselines

Is our tracker better than the established ones, on the same detections? Every tracker here
replaces only the tracker inside the normal scoring (`sim/analyze.py` -> `track.localize_and_score`):
same NR-Scope detections, same measurement groups, same bistatic measurement model and noise,
same birth state, same scoring (GOSPA, false tracks, coverage, switches). Nothing in `sim/` is changed.

## Trackers (`trackers/`)

| name | what | source |
|---|---|---|
| `OURS_full` | EKF on bistatic measurements, GNN, existence probability, GGIW extent for people, site rules (P_D map, wall ghosts, buildings) | `sim/track.py` |
| `C1_GNN-EKF` | GNN + EKF, M-of-N confirmation | Blackman & Popoli 1999 |
| `C2_GM-PHD-VoMa2006` | Gaussian-mixture PHD (EKF), measurement-driven birth, label tagging | Vo & Ma 2006; Ristic et al. 2012; Panta et al. 2009 |
| `C3_JIPDA-Musicki2004` | joint integrated PDA, exact joint events | Musicki & Evans 2004 |
| `C4_MHT-CoxHingorani1996` | hypothesis-oriented MHT, Murty k-best, LLR track scores | Reid 1979; Cox & Hingorani 1996; Blackman 2004 |
| `C5_BP-Meyer2017` | belief-propagation multi-sensor multitarget tracking | Meyer et al. 2017, 2018 |
| `C6_GGIW-PHD-Granstrom2012` | extended-target GGIW-PHD, distance partitioning | Granstrom & Orguner 2012 |
| `P1_Bauhofer2026_GM-PHD-range-doppler` | GM-PHD in range/range-rate, their parameters (single sensor) | Bauhofer et al., EuCNC 2026 |
| `P2_BSense2026_core` | confidence filtering + IMM, without their learned parts | Meng et al., MobiSys 2026 |
| `P3_CellSense2026_window-Hungarian` | Hungarian over a 10-snapshot window, no filter | Kumar et al., arXiv 2606.07900 |
| `P4_Gunes2025_plotfusion-GNN` | two-receiver plot fusion + deghosting, then GNN | Gunes & Caliskan, RadarConf 2025 |

Shared parts: `trackers/common.py`. Each file's docstring lists what follows the paper and what is
adapted or assumed. Only `OURS_full` uses site knowledge; the others know the receivers' field of view.

## Running

    python tracking_baselines/harness.py <run> [tracker ...]     # default settings -> results/<run>/<tracker>__<setup>.json
    python tracking_baselines/tune.py <run> [tracker ...]        # 3x3 grid each   -> results/<run>/<tracker>__grid<k>__<setup>.json
    python tracking_baselines/summary.py <run>                    # results/<run>/summary.md, .csv, .png

Setups: `1` (rx1 alone), `2` (rx2 alone), `fused` (both receivers' maps in one tracker).
Tuning is done on the scored run itself, for every tracker alike (optimistic for all).
