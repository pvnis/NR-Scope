#!/usr/bin/env python3
# SPDX-License-Identifier: LicenseRef-CSSL-1.0
"""Offline reprocessor for the sensing test_record_samples mode.

The C pipeline, with sensing.test_record_samples: true, writes one file per map,
<dump>.rec.NNNNN.csv, holding every (chain, layer) of the slow-time samples that map
was built from (nr_ue_sensing_dump_snapshots_tagged). This script rebuilds the
range-Doppler map from those files in pure Python, so any clutter parameter can be
swept offline without a rebuild or a new capture.

It is a faithful port of nr_ue_sensing_range_doppler() for the delay-profile-domain
stages: amplitude (gain) normalisation, clutter removal (none / mean / slow-trend /
kernel CLEAN), rank-1 comb removal, the non-uniform Doppler DFT, and the noise-floor
normalisation. Per-chain maps are averaged exactly as nr_ue_sensing_task_map() does.

A bare run (no --sweep) reproduces exactly the map normal mode would produce with the
same parameters, except that it does NOT apply the spatial null or the AoA: it is the
antenna- and layer-averaged clutter-removed range-Doppler map. To compare it against a
live map2d.csv, run the normal pipeline with spatial_null: false (the AoA only adds the
marker/angle payload, it does not change the power map).

The TDD detector and the +-v mirror rejection run as with tdd_detect and mirror_reject
on (tdd_detect.py, a port of nr_ue_tdd_detect.c): on the window's rx0 conditioned samples,
markers written into each map line and every marker to <label>.targets.csv. Disable with
--no-tdd-detect / --no-mirror-reject.

NOT ported (use the live pipeline, or extend here): the spatial null across chains,
AoA/MUSIC, and alignment (alignment is applied in the capture path
before these samples are stored, so it cannot be swept from a recording).

Each produced map is written one line per input file to <out>/<label>.csv (default
<recdir>/replay/) in the exact format nr_ue_sensing_dump_map() uses, so
scripts/sensing/plot_range_doppler.py plots it. The clutter parameters default to the
compile-time defaults (kernel mode, 4 paths, trend degree 2, comb on); pass --<field>
to match a different config, e.g. --clutter-mode 1 for mean mode.

Examples:
  # sweep the slow-trend degree 0..3 over a run's recordings
  replay_samples.py --rec '/tmp/map2d.csv.rec.*.csv' --sweep trend_degree --lo 0 --hi 3 --step 1

  # single reprocess at the defaults, plot the first map
  replay_samples.py --rec '/tmp/map2d.csv.rec.00001.csv'

  # sliding average of 4 consecutive maps (non-coherent), optionally velocity-compensated
  replay_samples.py --rec '~/runs/run_1/map2d.csv.rec.*.csv' --avg-maps 4 [--vcomp]

Residue vs slow movers (--residue = the three stages below; each also on its own):
  A --hygiene      before the maps: drop outlier snapshots, cut windows at scene steps
  C --clutter-map  after the maps, in place of the TDD detector: each range-Doppler cell
                   against its own history (stationary residue -> ~0 dB, a mover is new)
  E --track        after C: Kalman tracks in (range, range rate); a mover's range must
                   drift as its Doppler says, residue has Doppler but stays put
  Writes <label>.tracks.csv and <label>.time.png (range-time and Doppler-time of the
  run with the tracks: lime mover, cyan residue, grey unverified).
  --inject DB,R0,AMP,PERIOD adds a synthetic walker to the recording, to measure what
  survives on the real residue (indoor_1_record: found down to about -25 dB re the
  direct path, nothing labelled mover without it).

  replay_samples.py --rec '~/runs/indoor_1_record/map2d.csv.rec.*.csv' --residue
  replay_samples.py --rec '...' --residue --inject=-25,20,6,10 --sweep cm_pfa --lo 1e-5 --hi 1e-3 --step 1e-4
"""
import argparse
import glob
import math
import os
import sys

import numpy as np

import tdd_detect as tdd

C = 299792458.0
MAX_BINS_FREQ = 2048
MIN_SNAPSHOTS = 32
SLOW_TREND_MAX = 3
COMB_MAX_HARMONIC = 6
NR_MAP_NORM_SAMPLES = 4096

# Compile-time defaults from nr_ue_map.h; the sweep varies one of these.
DEFAULTS = dict(
    clutter_mode=2,       # 0 none, 1 mean, 2 kernel
    max_paths=4,
    trend_degree=2,
    kernel_half_span=48,
    static_min=0.8,
    snr_min=10.0,
    min_sep_bins=3,
    los_first_db=6.0,
    comb_remove=1,
    comb_harmonic=6,
    comb_tdd_multiple=2,
    comb_power_iters=24,
    sliding_la_ms=0.0,    # sliding ECA (ACRF-CS) averaging window, ms; 0 = off (slow trend used)
    sliding_ls_ms=0.5,    # its update block, ms; one slot pushes the block side peaks to 2 kHz
    # --kfactor analysis (fluctuation statistics per range bin over seconds)
    k_window_s=2.0,       # analysis window, s; a person standing still needs seconds to show
    k_hop_s=0.5,          # step between windows, s
    k_thresh_db=3.0,      # flag a bin when its (floor-corrected) Rician K is below this
    k_margin_db=6.0,      # ...and its slow fluctuation is this far above the window's floor (25th pct of bins)
    k_common_norm=1,      # remove a common complex gain per snapshot first (the slow common drift)
    k_floor_corr=1,       # use the slow (lagged-correlation) fluctuation power, which the floor barely reaches
    k_lag_ms=100.0,       # lag of that correlation; the floor decorrelates in ~100 ms, a swaying body
                          # does not (0: next snapshot)
    k_detrend=1,          # remove a linear drift per bin within the window
    # --hygiene (A): snapshots and scene steps the clutter removal cannot model
    a_ref_ms=100.0,       # span of the running local reference a snapshot is compared with
    a_outlier_x=5.0,      # drop a snapshot whose residue against it exceeds this x the run median
    a_block_ms=20.0,      # block for the scene-step test
    a_step_ms=100.0,      # the step test compares the mean of this span before and after
    a_step_x=4.0,         # a step is a change this x the run median of that comparison
    # --clutter-map (C): each cell against its own history
    cm_tau_s=5.0,         # memory of the per-cell background, s
    cm_learn_maps=10,     # maps that only build the background before anything is detected
    cm_pfa=1e-5,          # false-alarm rate per cell: power over its background mean > -ln(pfa)
    cm_min_db=6.0,        # ...and at least this many dB above it
    cm_vmin_ms=0.1,       # ignore |speed| below this (the static scene's own leakage)
    # --track (E): kinematic consistency of the detections
    tr_accel=1.0,         # random acceleration of a walker, m/s^2 (Kalman process noise)
    tr_gate=13.8,         # association gate, chi-square with 2 dof (13.8: 99.9%)
    tr_m=4,               # confirm a track after M hits ...
    tr_n=6,               # ... within its last N maps
    tr_max_miss=4,        # delete a track after this many maps in a row without a hit
    tr_check_s=1.5,       # history needed before range drift and Doppler are compared
    tr_tol_ms=0.5,        # allowed |range slope - Doppler range rate|, m/s (bistatic) ...
    tr_tol_frac=0.6,      # ... or this fraction of |Doppler range rate|, whichever is larger; residue
                          # has no range drift (slope ~0), so it misses by 1.0 x its rate
)

# field -> (lo, hi) hard range, for clamping a swept value (mirrors nr_sweep_apply)
RANGE = dict(
    clutter_mode=(0, 2), max_paths=(1, 4), trend_degree=(0, SLOW_TREND_MAX),
    kernel_half_span=(1, None), static_min=(0.0, None), snr_min=(0.0, None),
    min_sep_bins=(1, None), los_first_db=(0.0, None), comb_remove=(0, 1),
    comb_harmonic=(1, COMB_MAX_HARMONIC), comb_tdd_multiple=(1, None),
    sliding_la_ms=(0.0, None), sliding_ls_ms=(0.1, None),
    k_window_s=(0.05, None), k_hop_s=(0.01, None), k_thresh_db=(None, None), k_margin_db=(None, None),
    k_common_norm=(0, 1), k_floor_corr=(0, 1), k_detrend=(0, 1), k_lag_ms=(0.0, None),
    a_ref_ms=(1.0, None), a_outlier_x=(1.0, None), a_block_ms=(1.0, None), a_step_ms=(1.0, None),
    a_step_x=(1.0, None), cm_tau_s=(0.1, None), cm_learn_maps=(1, None), cm_pfa=(1e-12, 0.5),
    cm_min_db=(0.0, None), cm_vmin_ms=(0.0, None), tr_accel=(0.01, None), tr_gate=(0.1, None),
    tr_m=(1, None), tr_n=(1, None), tr_max_miss=(1, None), tr_check_s=(0.1, None), tr_tol_ms=(0.0, None), tr_tol_frac=(0.0, 0.95),
)
INT_FIELDS = {"k_common_norm", "k_floor_corr", "k_detrend", "clutter_mode", "max_paths", "trend_degree", "kernel_half_span",
              "min_sep_bins", "comb_remove", "comb_harmonic", "comb_tdd_multiple",
              "comb_power_iters", "cm_learn_maps", "tr_m", "tr_n", "tr_max_miss"}


# ---------------------------------------------------------------- primitives

def hann(L):
    if L < 8:
        return np.ones(L)
    i = np.arange(L)
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * i / (L - 1)))


def clutter_kernel(n_pilots, a_m, u0, idft_size, n_bins):
    """K[b] = sum_l w[l] exp(j 2pi (a_m + l)(b - u0) / idft_size), the delay response of
    one grant's Hann window. Port of nr_ue_sensing_clutter_kernel()."""
    w = hann(n_pilots)
    u = np.arange(n_bins) - u0
    l = np.arange(n_pilots)
    phase = 2.0 * np.pi * np.outer(u, a_m + l) / idft_size  # (n_bins, n_pilots)
    return (w[None, :] * np.exp(1j * phase)).sum(axis=1)


def slow_basis(t, degree):
    """Orthonormal real polynomial basis 1, tau, ... tau^degree on tau = (t-mid)/half.
    Port of nr_ue_sensing_slow_basis()."""
    n = len(t)
    if degree < 0:
        return np.zeros((0, n))
    degree = min(degree, SLOW_TREND_MAX)
    t_mid = t.mean()
    t_half = 0.5 * (t[-1] - t[0]) if t[-1] > t[0] else 1.0
    tau = (t - t_mid) / t_half
    Q = []
    for k in range(degree + 1):
        v = tau ** k
        for q in Q:
            v = v - (q @ v) * q
        nrm = np.linalg.norm(v)
        if nrm < 1e-9:
            continue
        Q.append(v / nrm)
    return np.array(Q) if Q else np.zeros((0, n))


def comb_basis(t, f0, n_harm, Q):
    """Orthonormal comb tones exp(j 2pi k f0 t), k=+-1..+-n_harm, made orthogonal to the
    real polynomial basis Q and to each other. Port of nr_ue_sensing_comb_basis()."""
    n = len(t)
    E = []
    for k in range(1, n_harm + 1):
        for sgn in (1, -1):
            f = sgn * k * f0
            turns = f * t
            ph = 2.0 * np.pi * (turns - np.trunc(turns))
            v = np.cos(ph) + 1j * np.sin(ph)
            for q in Q:                       # real polynomials
                v = v - (q @ v) * q
            for u in E:                       # tones already kept
                v = v - np.vdot(u, v) * u
            nrm2 = np.vdot(v, v).real
            if nrm2 <= 1e-3 * n:
                continue
            E.append(v / math.sqrt(nrm2))
    return np.array(E) if E else np.zeros((0, n), complex)


def peak_frac(p, b):
    """Sub-bin vertex of a parabola through p[b-1..b+1] in dB. Port of peak_frac()."""
    if b <= 0 or b >= len(p) - 1 or p[b - 1] <= 0 or p[b] <= 0 or p[b + 1] <= 0:
        return 0.0
    ym, y0, yp = (10.0 * math.log10(x) for x in (p[b - 1], p[b], p[b + 1]))
    den = ym - 2.0 * y0 + yp
    if den >= 0.0:
        return 0.0
    d = 0.5 * (ym - yp) / den
    return d if -0.5 < d < 0.5 else 0.0


def los_peak(e_prof, los_first_db):
    u0 = int(np.argmax(e_prof))
    e_min = e_prof[u0] * 10.0 ** (-los_first_db / 10.0)
    for b in range(u0):
        is_peak = (b == 0 or e_prof[b] >= e_prof[b - 1]) and e_prof[b] >= e_prof[b + 1]
        if is_peak and e_prof[b] >= e_min:
            return b
    return u0


def remove_path(res, snaps, u, half_span):
    """Fit one static path at bin u (single alpha over the window, per-snapshot kernel) and
    subtract it. Port of nr_ue_sensing_remove_path(); res is modified in place."""
    n_bins = res.shape[1]
    b_lo = max(0, int(math.floor(u)) - half_span)
    b_hi = min(n_bins - 1, int(math.ceil(u)) + half_span)
    if b_hi < b_lo:
        return
    span = b_hi - b_lo + 1
    cache, Ks = {}, []
    num, den = 0j, 0.0
    for i, s in enumerate(snaps):
        key = (s["n_pilots"], s["a_m"])
        if key not in cache:
            cache[key] = clutter_kernel(s["n_pilots"], s["a_m"], u - b_lo, s["idft"], span)
        K = cache[key]
        Ks.append(K)
        r = res[i, b_lo:b_hi + 1]
        num += np.vdot(K, r)          # sum r * conj(K)
        den += np.vdot(K, K).real
    if den > 0.0:
        alpha = num / den
        for i in range(len(snaps)):
            res[i, b_lo:b_hi + 1] -= alpha * Ks[i]


def sliding_filter(t, la_s, ls_s):
    """Sliding ECA in the channel domain (ACRF-CS, Liu et al. 2020; the ECA-S equivalent).

    The slow time is cut into blocks of ls_s seconds. Every sample of a block loses the mean
    of the samples within +-la_s/2 of the block centre, i.e. the static channel estimated
    around it, per delay bin. In this domain ECA's "project off the delayed copies of the
    reference" is exactly that subtraction, since the reference is already divided out.

    Unlike the slow trend, which fits one polynomial over the whole window, this follows a
    static path that varies during the window. The price is a Doppler notch about 1/la_s
    wide (la_s 40 ms: ~25 Hz, ~1.1 m/s at 3.45 GHz). The block structure puts side peaks at
    multiples of 1/ls_s; one slot (0.5 ms) moves them to 2 kHz, off the map. Blocks and
    windows are in time, not in sample counts, because TDD sampling is not uniform.

    Returns the filter as a function of a slow-time array [n] or [n, k] (axis 0 = time),
    so the TDD detector can model a target's tone through the same filter."""
    t = np.asarray(t, dtype=np.float64)
    ls_s = max(ls_s, 1e-6)
    la_s = max(la_s, ls_s)                          # the window must cover its own block
    blk = np.floor((t - t[0]) / ls_s).astype(np.int64)
    ub = np.unique(blk)
    centres = t[0] + (ub + 0.5) * ls_s
    lo = np.searchsorted(t, centres - 0.5 * la_s, side="left")
    hi = np.searchsorted(t, centres + 0.5 * la_s, side="right")
    pos = np.searchsorted(ub, blk)
    s_lo, s_hi = lo[pos], hi[pos]

    def f(v):
        v = np.asarray(v)
        cs = np.concatenate([np.zeros((1,) + v.shape[1:], dtype=v.dtype), np.cumsum(v, axis=0)])
        cnt = np.maximum(s_hi - s_lo, 1).reshape((-1,) + (1,) * (v.ndim - 1))
        return v - (cs[s_hi] - cs[s_lo]) / cnt
    return f


# ---------------------------------------------------------------- the map

def build_map(snaps, prm, n_freq_fixed=0):
    """Range-Doppler power map for one (chain, layer). Returns dict or None."""
    n = len(snaps)
    if n < MIN_SNAPSHOTS:
        return None
    h = np.stack([s["h"] for s in snaps])          # (n, n_bins) complex
    t = np.array([s["t"] for s in snaps])
    n_pil = np.array([s["n_pilots"] for s in snaps])
    n_bins = h.shape[1]
    s0 = snaps[0]
    scs, carrier, k_step, idft0 = s0["scs"], s0["carrier"], s0["k_step"], s0["idft"]

    tdd_period = prm["tdd_slots"] * (1e-3 / (scs / 15000.0))
    v = prm["max_speed"]
    f_max = (2.0 * v * carrier / C if carrier > 0 and v > 0 else 0.0) + 1.0 / tdd_period
    t_span = t[-1] - t[0]
    if t_span <= 0:
        return None
    m_per_bin = C / (idft0 * k_step * scs)

    n_freq = int(4.0 * f_max * t_span) | 1
    n_freq = max(n_freq, 3)
    if 0 < n_freq_fixed <= MAX_BINS_FREQ:
        n_freq = n_freq_fixed

    n_max = int(n_pil.max())
    if n_pil.min() <= 0:
        return None
    gain = n_max / n_pil                            # (n,)

    e_prof = (np.abs(h) ** 2).sum(axis=0)           # raw energy profile
    u0 = los_peak(e_prof, prm["los_first_db"])
    bin_los = u0 + peak_frac(e_prof, u0)

    mode = prm["clutter_mode"]
    res = h.astype(np.complex128).copy()            # LOS norm off by default

    # kernel CLEAN: direct path, then greedy static peaks (on the gain-weighted mean)
    if mode == 2:
        max_paths = int(np.clip(prm["max_paths"], 1, 4))
        remove_path(res, snaps, bin_los, prm["kernel_half_span"])
        u_fit = [bin_los]
        while len(u_fit) < max_paths:
            g = gain[:, None] * res
            coh = (np.abs(g.sum(axis=0)) ** 2) / (n * n)
            inc = (np.abs(g) ** 2).sum(axis=0) / n
            far = np.array([min(abs(b - uf) for uf in u_fit) >= prm["min_sep_bins"]
                            for b in range(n_bins)])
            if not far.any():
                break
            b_best = int(np.argmax(np.where(far, coh, -np.inf)))
            noise_mean = np.median(inc) / n
            s_static = coh[b_best] / inc[b_best] if inc[b_best] > 0 else 0.0
            if s_static < prm["static_min"] or coh[b_best] < prm["snr_min"] * noise_mean:
                break
            u = b_best + peak_frac(coh, b_best)
            remove_path(res, snaps, u, prm["kernel_half_span"])
            u_fit.append(u)

    # gain, then per-bin static removal: the slow trend (C default) or, when
    # sliding_la_ms > 0, the sliding ECA in its place (not in the C pipeline yet)
    slow_f = None
    z = gain[:, None] * res                         # (n, n_bins)
    if mode != 0 and prm.get("sliding_la_ms", 0) > 0:
        slow_f = sliding_filter(t, prm["sliding_la_ms"] * 1e-3, prm["sliding_ls_ms"] * 1e-3)
        z = slow_f(z)
        trend_degree = -1                           # nothing projected; the detector models slow_f
        Q = slow_basis(t, 0)                        # the comb is kept orthogonal to the mean only
    else:
        trend_degree = -1 if mode == 0 else prm["trend_degree"]
        Q = slow_basis(t, trend_degree)             # (n_q, n)
        for q in Q:
            z = z - np.outer(q, q @ z)              # remove projection on q, per bin
    res = z

    # rank-1 comb removal
    n_q = Q.shape[0]
    if prm["comb_remove"] and mode != 0 and n_q > 0:
        f0 = 1.0 / (prm["comb_tdd_multiple"] * tdd_period)
        E = comb_basis(t, f0, prm["comb_harmonic"], Q)   # (n_e, n)
        if E.shape[0] > 0:
            # P[b, j] = <e_j, z_b> = sum_i conj(e_j[i]) res[i, b]
            P = (E.conj() @ res).T                        # (n_bins, n_e)
            # dominant singular pair of P by power iteration
            n_e = E.shape[0]
            vvec = np.ones(n_e, complex) / math.sqrt(n_e)
            sigma = 0.0
            for _ in range(prm["comb_power_iters"]):
                u = P @ vvec
                nu = np.linalg.norm(u)
                if nu <= 0:
                    break
                u /= nu
                vvec = P.conj().T @ u
                sigma = np.linalg.norm(vvec)
                if sigma <= 0:
                    break
                vvec /= sigma
            if sigma > 0:
                # subtract sigma u[b] conj(v[j]) along e_j: res -= (P_rank1 @ E)
                rank1 = sigma * np.outer(u, vvec.conj())  # (n_bins, n_e)
                res -= (rank1 @ E).T                       # (n, n_bins)

    # non-uniform Doppler DFT, power, 1/n^2 normalisation
    df = (2.0 * f_max) / (n_freq - 1) if n_freq > 1 else 0.0
    freqs = -f_max + df * np.arange(n_freq)
    W = np.exp(-2j * np.pi * np.outer(t, freqs))    # (n, n_freq)
    A = res.T @ W                                    # (n_bins, n_freq)
    power = (np.abs(A) ** 2) / (n * n)

    # normalise to the noise floor, exactly as range_doppler does: median of a strided
    # (odd stride) sample of the row-major cells, not a full median.
    flat = power.reshape(-1)                          # row-major: b outer, f inner
    n_cells = flat.size
    stride = max(1, -(-n_cells // NR_MAP_NORM_SAMPLES)) | 1
    samp = np.sort(flat[0:n_cells:stride][:NR_MAP_NORM_SAMPLES])
    ref = samp[len(samp) // 2] if samp.size else 0.0  # sorted[n/2], as in C
    if ref > 0:
        power = power / ref

    return dict(power=power.astype(np.float32), n_bins=n_bins, n_freq=n_freq,
                f_max=f_max, m_per_bin=m_per_bin, t_span=t_span, carrier=carrier,
                bin_los=bin_los, n_snap=n, k_step=k_step, layer=s0["layer"],
                t_center=0.5 * (snaps[0]["t_abs"] + snaps[-1]["t_abs"]),
                # the conditioned slow-time samples and what the TDD detector needs with
                # them, as range_doppler hands them back in slow_out
                res=res, t=t, win_len=n_max, win_start=int(s0["a_m"]), idft=idft0,
                trend_degree=trend_degree, tdd_period=tdd_period, slow_filter=slow_f)


def build_averaged(groups, prm, antenna_avg, layer_avg, n_freq_fixed=0):
    """Average per-(chain,layer) maps on the first one's grid, as task_map does.
    n_freq_fixed > 0 puts the whole window on that Doppler grid (used to keep every
    window of a run on one grid, so maps can be averaged across windows)."""
    keys = sorted(groups)
    if not antenna_avg:
        keys = [k for k in keys if k[0] == min(a for a, _ in groups)]
    if not layer_avg:
        keys = [k for k in keys if k[1] == min(l for _, l in keys)]
    first = build_map(groups[keys[0]], prm, n_freq_fixed)
    if first is None:
        return None
    acc = first["power"].astype(np.float64).copy()
    n_comb = 1
    for k in keys[1:]:
        m = build_map(groups[k], prm, first["n_freq"])
        if m is None or m["n_bins"] != first["n_bins"] or m["n_freq"] != first["n_freq"]:
            continue
        acc += m["power"]
        n_comb += 1
    first["power"] = (acc / n_comb).astype(np.float32)
    first["n_combined"] = n_comb
    return first


def average_maps(ms, K, vcomp):
    """Non-coherent average of each map with the most recent older maps, K in all at most,
    exactly as nr_ue_sensing_avg_apply() in the C pipeline (sensing.avg_maps).

    Output i is window i: the first maps of a run average over 1, 2, ... maps until K are
    available. Older maps are picked by time, within (K-1)*1.5 window spans, so a dropped
    window shortens the average instead of reaching further back. Averaging power keeps the
    floor's mean and shrinks its fluctuation, so the noise speckle drops and a target that
    holds its cell gains margin.

    vcomp: velocity-compensated averaging. Every Doppler column f has a known bistatic
    path-length rate, dR/dt = -lambda*f (approaching = positive Doppler). An older map j,
    dt = t_i - t_j earlier, saw that content at path length R + lambda*f*dt, so column f of
    map j is read shifted by lambda*f*dt/m_per_bin bins before averaging. A mover then stays
    aligned across the K maps instead of smearing over range. Cells read from outside the
    map are left out of that cell's mean."""
    if K <= 1:
        return ms
    out = []
    for i, ref in enumerate(ms):
        nb, nf = ref["n_bins"], ref["n_freq"]
        horizon = (K - 1) * 1.5 * ref["t_span"]
        older = [j for j in range(i) if 0.0 < ref["t_center"] - ms[j]["t_center"] <= horizon]
        older = sorted(older, key=lambda j: -ms[j]["t_center"])[:K - 1]
        lam = C / ref["carrier"]
        freqs = -ref["f_max"] + np.arange(nf) * (2.0 * ref["f_max"] / (nf - 1))
        acc = ref["power"].astype(float).copy()
        cnt = np.ones((nb, nf))
        b = np.arange(nb, dtype=float)
        for j in older:
            P = ms[j]["power"].astype(float)
            if not vcomp:
                acc += P
                cnt += 1.0
                continue
            dt = ref["t_center"] - ms[j]["t_center"]
            shift = lam * freqs * dt / ref["m_per_bin"]          # bins, per column
            for f in range(nf):
                col = np.interp(b + shift[f], b, P[:, f], left=np.nan, right=np.nan)
                ok = ~np.isnan(col)
                acc[ok, f] += col[ok]
                cnt[ok, f] += 1.0
        avg = dict(ref)
        avg["power"] = (acc / cnt).astype(np.float32)
        avg["n_avg"] = 1 + len(older)
        out.append(avg)
    return out


# ---------------------------------------------------------------- io

def parse_rec(path):
    """Return {(aarx, layer): [snapshot dicts]} and tdd_slots/max_speed from one file."""
    groups, tdd_slots, max_speed = {}, 10, 0.0
    with open(path) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split(",")
            if len(parts) < 15:
                continue
            aarx = int(parts[0]); layer = int(parts[1])
            vals = parts
            n_pilots = int(vals[7]); idft = int(vals[8]); n_bins = int(vals[9])
            scs = int(vals[10]); carrier = float(vals[12])
            max_speed = float(vals[13]); tdd_slots = int(vals[14])
            k_step = int(vals[5]); k_first = int(vals[6])
            a_m = (k_first - k_first % k_step) // k_step
            iq = np.array(vals[15:15 + 2 * n_bins], dtype=np.float64)
            h = iq[0::2] + 1j * iq[1::2]
            groups.setdefault((aarx, layer), []).append(dict(
                t=float(vals[4]), t_abs=int(vals[3]) / float(vals[11]),
                n_pilots=n_pilots, idft=idft, scs=scs,
                carrier=carrier, k_step=k_step, a_m=a_m, layer=layer, h=h))
    return groups, tdd_slots, max_speed


def write_map_line(fh, idx, m):
    """One line in nr_ue_sensing_dump_map() format, for plot_range_doppler.py."""
    markers = m.get("markers", [])
    hdr = [idx, 0, 0, 0, m["layer"], m["k_step"], 0, m["n_snap"],
           m["t_span"], m["m_per_bin"], m["f_max"], m["carrier"],
           m["n_bins"], m["n_freq"], len(markers), 0, m["bin_los"]]
    power = m["power"].reshape(-1)  # row-major b outer, f inner
    # markers, 4 numbers each, as the C writes them: range_m, speed_ms, snr_dB, verdict
    mk = [x for t in markers for x in (t["range_m"], t["speed_ms"], t["snr_dB"], t["verdict"])]
    fh.write(",".join(f"{x:g}" for x in hdr) + "," +
             ",".join(f"{x:g}" for x in list(power) + mk) + "\n")


# ---------------------------------------------------------------- K-factor analysis

def load_run(files):
    """Every recorded snapshot of a run per (aarx, layer), time ordered, de-duplicated and
    normalised to the widest grant, so analysis windows can span consecutive recordings."""
    acc, m_per_bin = {}, None
    for path in files:
        g, _, _ = parse_rec(path)
        for k, snaps in g.items():
            acc.setdefault(k, []).extend(snaps)
            if m_per_bin is None and snaps:
                s0 = snaps[0]
                m_per_bin = C / (s0["idft"] * s0["k_step"] * s0["scs"])
    out = {}
    for k, snaps in acc.items():
        t = np.array([x["t_abs"] for x in snaps])
        order = np.argsort(t, kind="stable")
        keep = np.concatenate([[True], np.diff(t[order]) > 1e-9])
        idx = order[keep]
        npil = np.array([snaps[i]["n_pilots"] for i in idx], dtype=float)
        H = np.stack([snaps[i]["h"] for i in idx]) * (npil.max() / npil)[:, None]
        out[k] = dict(t=t[idx], h=H)
    return out, m_per_bin


def kfactor_window(H, t, prm):
    """Per-bin powers of one chain over one analysis window: static |mean|^2, total
    fluctuation, and its slow (lag-1 correlated) part.

    Common gain: over seconds every bin wanders together (gain and phase of the whole
    snapshot, several rad in 2 s on the lab runs), which makes walls look as alive as a
    person. A complex gain per snapshot, fitted against the static profile and refined a
    few times, takes it out; a single mover barely moves that fit since the static scene
    dominates it. Slow part: |mean(r(t+lag) conj(r(t)))|. The floor decorrelates over ~100 ms
    and mostly averages out of it, while a body's fluctuation, correlated over longer, stays."""
    if prm["k_common_norm"]:
        ref = H.mean(0)
        for _ in range(3):
            g = (H @ ref.conj()) / np.vdot(ref, ref).real
            g[np.abs(g) < 1e-12] = 1.0
            Hn = H / g[:, None]
            ref = Hn.mean(0)
        H = Hn
    mu = H.mean(0)
    r = H - mu
    if prm["k_detrend"]:
        A = np.vstack([np.ones_like(t), t - t.mean()]).T
        r = r - A @ np.linalg.lstsq(A, r, rcond=None)[0]
    p_tot = (np.abs(r) ** 2).mean(0)
    lag = prm.get("k_lag_ms", 0.0) * 1e-3
    if lag > 0:
        # pairs about `lag` apart; measured on the lab runs the floor keeps ~20% correlation
        # for tens of ms and ~5% at 100-200 ms, so a lag of that order strips far more of it
        # than the next snapshot does
        j = np.searchsorted(t, t + lag)
        i = np.nonzero(j < len(t))[0]
        j = j[i]
        ok = np.abs(t[j] - t[i] - lag) < max(0.25e-3, 0.1 * lag)
        i, j = i[ok], j[ok]
    else:
        i, j = np.arange(len(t) - 1), np.arange(1, len(t))
    p_slow = np.abs((r[j] * r[i].conj()).mean(0)) if len(i) else np.zeros_like(p_tot)
    return np.abs(mu) ** 2, p_tot, p_slow


def kfactor_run(data, m_per_bin, prm, label, out_dir):
    """Sliding K-factor analysis over a whole run: per window and range bin, chains summed in
    power. Writes <label>.kfactor.csv and a time x range picture, prints a summary."""
    W, hop = prm["k_window_s"], prm["k_hop_s"]
    t0 = min(d["t"][0] for d in data.values())
    t1 = max(d["t"][-1] for d in data.values())
    rows, times, kmap, fmap = [], [], [], []
    for start in np.arange(t0, t1 - W + 1e-9, hop):
        S = Pt = Ps = 0.0
        n_used = 0
        for d in data.values():
            sel = (d["t"] >= start) & (d["t"] < start + W)
            if sel.sum() < 200:
                continue
            a, b, c = kfactor_window(d["h"][sel], d["t"][sel], prm)
            S, Pt, Ps, n_used = S + a, Pt + b, Ps + c, n_used + 1
        if n_used == 0:
            continue
        los = S.max()
        fl = Ps if prm["k_floor_corr"] else Pt
        db = lambda x: 10.0 * np.log10(np.maximum(x, 1e-30))
        static_db, fl_db = db(S / los), db(fl / los)
        k_raw, k_db = db(S / Pt), db(S / fl)
        # what "significant" fluctuation means is relative to this window's floor: the white
        # floor keeps a small lag-1 correlation (~0.05 of its power) that an absolute level
        # would mistake for slow motion, and the floor itself drifts over a run
        floor_db = np.percentile(fl_db, 25)
        flag = (fl_db >= floor_db + prm["k_margin_db"]) & (k_db < prm["k_thresh_db"])
        tc = start + 0.5 * W - t0
        times.append(tc); kmap.append(k_db); fmap.append(flag)
        for b in range(len(S)):
            rows.append((tc, b, b * m_per_bin, static_db[b], k_raw[b], k_db[b], Ps[b] / Pt[b], fl_db[b],
                         fl_db[b] - floor_db, int(flag[b])))
    if not rows:
        print(f"{label}: no window long enough")
        return
    path = os.path.join(out_dir, f"{label}.kfactor.csv")
    with open(path, "w") as fh:
        fh.write("t_s,bin,range_m,static_db,k_raw_db,k_db,lag1_coh,fluct_db,fluct_over_floor_db,flag\n")
        for r in rows:
            fh.write(f"{r[0]:.3f},{r[1]},{r[2]:.2f},{r[3]:.2f},{r[4]:.2f},{r[5]:.2f},{r[6]:.3f},{r[7]:.2f},"
                     f"{r[8]:.2f},{r[9]}\n")
    K, F = np.array(kmap), np.array(fmap)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 5))
        nb = K.shape[1]
        im = ax.imshow(np.clip(K.T, -10, 20), aspect="auto", origin="lower", cmap="viridis",
                       extent=[times[0] - hop / 2, times[-1] + hop / 2, -0.5 * m_per_bin, (nb - 0.5) * m_per_bin])
        ti, bi = np.nonzero(F)
        ax.scatter(np.array(times)[ti], bi * m_per_bin, s=8, c="red", label=f"K < {prm['k_thresh_db']:g} dB")
        ax.set_xlabel("time in run (s)"); ax.set_ylabel("bistatic path length (m)")
        ax.set_title(f"{label}: Rician K per range bin, {W:g} s windows"); ax.legend(loc="upper right")
        fig.colorbar(im, ax=ax, label="K (dB)")
        fig.savefig(os.path.join(out_dir, f"{label}.kfactor.png"), dpi=110, bbox_inches="tight")
        plt.close(fig)
    except ImportError:
        pass
    print(f"{label:28s}  {len(times)} windows of {W:g} s | flagged cells {F.mean()*100:.1f}% | "
          f"windows with a flag {np.mean(F.any(1))*100:.0f}%")
    for tc, f in zip(times, F):
        if f.any():
            rng = ", ".join(f"{b*m_per_bin:.0f}" for b in np.nonzero(f)[0][:8])
            print(f"    t {tc:6.2f} s: {f.sum():2d} bin(s) at {rng}{' ...' if f.sum() > 8 else ''} m")


# ---------------------------------------------------------------- detection

def run_detector(m, mirror, cfg=None):
    """The TDD detector on the window's rx0 conditioned samples, on the map's own Doppler
    grid, then the mirror rejection on the (possibly averaged) map, as map_task does with
    tdd_detect and mirror_reject on. Attaches m["markers"]: targets first, then rejections."""
    obs = tdd.Obs(h=m["res"], t=m["t"], idft_size=m["idft"], win_len=m["win_len"],
                  win_start=m["win_start"], win=hann(m["win_len"]), m_per_bin=m["m_per_bin"],
                  n_freq=m["n_freq"], f_max_hz=m["f_max"], lambda_m=C / m["carrier"],
                  t_tdd_s=m["tdd_period"], trend_q=slow_basis(m["t"], m["trend_degree"]),
                  slow_filter=m.get("slow_filter"))
    targets, rejected = tdd.detect(obs, cfg)
    if mirror:
        targets, rejected = tdd.mirror_reject(targets, rejected, m["power"], m["f_max"])
    m["markers"] = targets + rejected


def write_targets(path, ms):
    """Every marker of every map, one per line, for persistence analysis."""
    names = {tdd.VERDICT_TARGET: "target", tdd.VERDICT_REPLICA: "replica",
             tdd.VERDICT_DUPLICATE: "duplicate", tdd.VERDICT_UNCERTAIN: "uncertain"}
    with open(path, "w") as fh:
        fh.write("map,t_s,verdict,bin,range_m,f_hz,speed_ms,snr_dB\n")
        t0 = ms[0]["t_center"]
        for i, m in enumerate(ms):
            for t in m["markers"]:
                fh.write(f"{i},{m['t_center'] - t0:.3f},{names[t['verdict']]},{t['bin']:.2f},"
                         f"{t['range_m']:.2f},{t['f_hz']:.1f},{t['speed_ms']:.2f},{t['snr_dB']:.1f}\n")


# ---------------------------------------------------------------- synthetic walker

def inject_walker(parsed, spec):
    """Add a synthetic point target to every recorded snapshot, to measure what A/C/E can
    find on the real residue. spec = "level_db,r0_m,amp_m,period_s": a walker whose
    bistatic excess range swings R(t) = r0 + amp sin(2 pi t / period), at level_db relative
    to the run's median direct-path peak. Its delay response is the grant's own Hann kernel
    at R/m_per_bin, its slow-time phase -2 pi R / lambda (so its Doppler is -Rdot/lambda,
    consistent with its range drift); chain a gets a fixed extra phase 0.7 a rad."""
    level_db, r0, amp, period = (float(x) for x in spec.split(","))
    peaks = [np.abs(sn["h"]).max() for g in parsed for ss in g.values() for sn in ss[::50]]
    a0 = 10.0 ** (level_db / 20.0) * float(np.median(peaks))
    t0 = min(sn["t_abs"] for g in parsed for ss in g.values() for sn in ss)
    cache = {}   # kernel on a 1/100-bin delay grid per grant shape: one per snapshot is too slow
    for g in parsed:
        for (aarx, _), ss in g.items():
            for sn in ss:
                tt = sn["t_abs"] - t0
                lam = C / sn["carrier"]
                m_per_bin = C / (sn["idft"] * sn["k_step"] * sn["scs"])
                R = r0 + amp * math.sin(2.0 * math.pi * tt / period)
                key = (sn["n_pilots"], sn["a_m"], sn["idft"], len(sn["h"]), int(round(R / m_per_bin * 100)))
                K = cache.get(key)
                if K is None:
                    K = clutter_kernel(key[0], key[1], key[4] / 100.0, key[2], key[3])
                    K = cache[key] = K / np.abs(K).max()
                sn["h"] = sn["h"] + a0 * np.exp(-2j * np.pi * R / lam + 0.7j * aarx) * K
    return a0


# ---------------------------------------------------------------- A: snapshot hygiene

def _running_mean(t, H, half_s, keep):
    """Mean of the kept rows of H within +-half_s of each row's time, per column."""
    w = keep.astype(float)
    cs = np.vstack([np.zeros((1, H.shape[1]), complex), np.cumsum(H * w[:, None], axis=0)])
    cw = np.concatenate([[0.0], np.cumsum(w)])
    lo = np.searchsorted(t, t - half_s, "left")
    hi = np.searchsorted(t, t + half_s, "right")
    n = np.maximum(cw[hi] - cw[lo], 1.0)
    return (cs[hi] - cs[lo]) / n[:, None]


def hygiene(parsed, prm):
    """A. Clean the run before any map is built, on the reference stream (lowest chain and
    layer), and apply the verdict to every (chain, layer) of the same symbol.

    1. Outlier snapshots. Each snapshot is compared with the mean of the run's snapshots
       within +-a_ref_ms/2: residue = |h - mean|^2 / |mean|^2. Indoors 2.5% of snapshots
       sit 10-500x above the median, in bursts, with 2-6x the power (a grant on another
       precoder or power level). The clutter removal cannot model them and each one lands
       on every range and Doppler cell of its map. Two passes, the second with the
       reference rebuilt without the first pass's outliers.
    2. Scene steps. The mean of a_step_ms after each a_block_ms boundary is compared with
       the mean before it; a local maximum above a_step_x x the run median is a step (a new
       effective channel: the static scene changes shape and power, x0.4-4 indoors). A
       window that straddles one smears the step over every Doppler cell, so it keeps only
       its larger side.

    parsed: one {(aarx, layer): [snaps]} per recorded window. Returns the cleaned list and
    a stats dict."""
    keys = sorted({k for g in parsed for k in g})
    ref = min(keys)
    snaps = [s for g in parsed for s in g.get(ref, [])]
    if len(snaps) < MIN_SNAPSHOTS:
        return parsed, dict(outliers=0, steps=[], steps_abs=[], dropped_side=0, total=len(snaps))
    t = np.array([s["t_abs"] for s in snaps])
    order = np.argsort(t, kind="stable")
    t = t[order]
    npil = np.array([snaps[i]["n_pilots"] for i in order], float)
    H = np.stack([snaps[i]["h"] for i in order]) * (npil.max() / npil)[:, None]

    keep = np.ones(len(t), bool)
    half = 0.5 * prm["a_ref_ms"] * 1e-3
    for _ in range(2):
        loc = _running_mean(t, H, half, keep)
        e = np.sum(np.abs(H - loc) ** 2, 1) / np.maximum(np.sum(np.abs(loc) ** 2, 1), 1e-30)
        keep = e < prm["a_outlier_x"] * np.median(e)
    bad_t = set(np.round(t[~keep] * 1e7).astype(np.int64).tolist())

    # scene steps on block means of the kept snapshots
    blk = prm["a_block_ms"] * 1e-3
    span = max(1, int(round(prm["a_step_ms"] / prm["a_block_ms"])))
    edges = np.arange(t[0], t[-1] + blk, blk)
    idx = np.clip(np.searchsorted(edges, t, "right") - 1, 0, len(edges) - 1)
    M = np.full((len(edges), H.shape[1]), np.nan + 0j)
    for i in np.unique(idx[keep]):
        M[i] = H[keep & (idx == i)].mean(0)
    d = np.zeros(len(edges))
    for i in range(span, len(edges) - span):
        a = M[i - span:i]; b = M[i:i + span]
        a = a[~np.isnan(a[:, 0].real)]; b = b[~np.isnan(b[:, 0].real)]
        if len(a) * 2 < span or len(b) * 2 < span:
            continue
        A = a.mean(0); B = b.mean(0)
        d[i] = np.linalg.norm(B - A) / max(np.linalg.norm(A), 1e-30)
    med = np.median(d[d > 0]) if np.any(d > 0) else 0.0
    steps = [edges[i] for i in range(1, len(d) - 1)
             if med > 0 and d[i] > prm["a_step_x"] * med and d[i] >= d[i - 1] and d[i] >= d[i + 1]]

    out, dropped_side = [], 0
    for g in parsed:
        ng = {}
        for k, ss in g.items():
            ss = [s for s in ss if int(round(s["t_abs"] * 1e7)) not in bad_t]
            if ss:
                ts = np.array([s["t_abs"] for s in ss])
                for st in steps:
                    if ts[0] < st < ts[-1]:
                        before = ts < st
                        n0 = int(before.sum())
                        side = before if n0 >= len(ss) - n0 else ~before
                        dropped_side += len(ss) - int(side.sum())
                        ss = [s for s, kp in zip(ss, side) if kp]
                        ts = ts[side]
            ng[k] = ss
        out.append(ng)
    return out, dict(outliers=int((~keep).sum()), total=len(t), steps=[s - t[0] for s in steps], steps_abs=steps,
                     dropped_side=dropped_side)


# ---------------------------------------------------------------- C: clutter-map CFAR

def clutter_map(ms, prm, steps=()):
    """C. Detect against each cell's own history instead of against its neighbours.

    The residue the clutter removal leaves is stationary: the same symmetric +-0.3-1.3 m/s
    skirt on the same range bins all run long. A CFAR across neighbouring cells reads it as
    targets; a cell's own past does not, while a person moving through the cell is new to
    it. Per (range, Doppler) cell the mean power is kept, exponentially over cm_tau_s (a
    classic clutter map); the first cm_learn_maps maps only build it (plain average). A cell
    is then a detection when its power exceeds its background mean by -ln(cm_pfa) (the
    exceedance of exponential power, which both noise and residue speckle follow) and by
    cm_min_db. Detected cells do not update the background, so a target standing in a cell
    is not learnt away while it is detected. Only |speed| in [cm_vmin_ms, max_speed]: past
    max_speed the map holds the TDD replicas.

    steps: absolute times of scene steps (hygiene()). A new effective channel leaves a new
    residue, so the background is learnt again from the first map after a step.

    Then local maxima only (+-1 bin, +-3 Doppler cells), and the TDD replica test: a
    detection with a stronger one at f +- k/T_TDD on the same range (+-1 bin) is a replica.
    Attaches m["markers"] (verdict 0 kept, 1 replica) and m["excess_db"]."""
    if not ms:
        return
    nb, nf = ms[0]["power"].shape
    mu = np.zeros((nb, nf)); n_seen = 0
    t_prev = None
    thr = -math.log(prm["cm_pfa"])
    thr = max(thr, 10.0 ** (prm["cm_min_db"] / 10.0))
    steps = sorted(steps)
    for m in ms:
        P = m["power"].astype(float)
        if P.shape != (nb, nf):
            m["markers"] = []
            continue
        if t_prev is not None and any(t_prev < st <= m["t_center"] for st in steps):
            n_seen = 0
        dt = m["t_span"] if t_prev is None else max(m["t_center"] - t_prev, 1e-3)
        t_prev = m["t_center"]
        alpha = max(1.0 - math.exp(-dt / prm["cm_tau_s"]), 1.0 / (n_seen + 1))
        freqs = -m["f_max"] + np.arange(nf) * (2.0 * m["f_max"] / (nf - 1))
        lam = C / m["carrier"]
        speed = freqs * lam / 2.0
        if n_seen == 0:
            mu[:] = P
        ratio = P / np.maximum(mu, 1e-30)
        exc = 10.0 * np.log10(np.maximum(ratio, 1e-12))
        det = np.zeros((nb, nf), bool)
        if n_seen >= prm["cm_learn_maps"]:
            vok = (np.abs(speed) >= prm["cm_vmin_ms"]) & (np.abs(speed) <= prm["max_speed"])
            det = (ratio > thr) & vok[None, :]
        m["excess_db"] = exc.astype(np.float32)
        upd = ~det
        mu[upd] += alpha * (P[upd] - mu[upd])
        n_seen += 1

        cand = []
        for b, f in zip(*np.nonzero(det)):
            b0, b1 = max(0, b - 1), min(nb, b + 2)
            f0, f1 = max(0, f - 3), min(nf, f + 4)
            if exc[b, f] >= exc[b0:b1, f0:f1].max():
                cand.append((exc[b, f], b, f))
        cand.sort(reverse=True)
        f_tdd = 1.0 / m["tdd_period"]
        df = freqs[1] - freqs[0]
        markers = []
        for z, b, f in cand:
            rep = False
            for (z2, b2, f2) in cand:
                if z2 <= z or abs(b2 - b) > 1:
                    continue
                k = (freqs[f] - freqs[f2]) / f_tdd
                if abs(k - round(k)) * f_tdd <= 2 * df and round(k) != 0:
                    rep = True
                    break
            frac = peak_frac(P[:, f], b) if 0 < b < nb - 1 else 0.0
            markers.append(dict(verdict=tdd.VERDICT_REPLICA if rep else tdd.VERDICT_TARGET,
                                bin=b + frac, range_m=(b + frac) * m["m_per_bin"], f_hz=freqs[f],
                                speed_ms=speed[f], snr_dB=float(z)))
        markers.sort(key=lambda x: x["verdict"])
        m["markers"] = markers


# ---------------------------------------------------------------- E: kinematic tracking

def track(ms, prm):
    """E. Track the kept detections over maps and test kinematic consistency.

    State [R, Rdot]: bistatic excess range and its rate, constant-velocity Kalman filter
    with random acceleration tr_accel. A detection measures both: R from its range bin and
    Rdot = -lambda f from its Doppler (approaching, positive Doppler, shortens the path).
    Association: greedy nearest, chi-square gate tr_gate (2 dof). A track is confirmed after
    tr_m hits in its last tr_n maps and deleted after tr_max_miss misses in a row.

    The discriminator: a real mover's range drifts as its Doppler says, a residue cell has
    a Doppler and never moves. Once a track spans tr_check_s, the slope of a straight line
    through its measured ranges is compared with its mean measured Rdot:
      mover      |slope - Rdot| <= tr_tol_ms, or tr_tol_frac |Rdot| if larger
      residue    it has Doppler (|Rdot| * span > one range bin) but no matching drift
      unverified too short, or too slow for its drift to show within a bin
    Returns the tracks; each marker gets "track" and its verdict is set: 0 for a mover's
    detection, 3 (uncertain) for everything else kept by C."""
    tracks, nid = [], 0
    for i, m in enumerate(ms):
        lam = C / m["carrier"]
        sR = m["m_per_bin"] / 2.5
        sV = lam * (1.0 / max(m["t_span"], 1e-3))
        Rm = np.diag([sR * sR, sV * sV])
        dets = [k for k in m.get("markers", []) if k["verdict"] == tdd.VERDICT_TARGET]
        z = [np.array([k["range_m"], -lam * k["f_hz"]]) for k in dets]
        t = m["t_center"]
        live = [tr for tr in tracks if not tr["dead"]]
        pred = []
        for tr in live:
            dt = t - tr["t"]
            F = np.array([[1.0, dt], [0.0, 1.0]])
            q = prm["tr_accel"] ** 2
            Q = q * np.array([[dt ** 3 / 3, dt ** 2 / 2], [dt ** 2 / 2, dt]])
            x = F @ tr["x"]; Pp = F @ tr["P"] @ F.T + Q
            pred.append((x, Pp))
        pairs = []
        for a, (x, Pp) in enumerate(pred):
            S = Pp + Rm; Si = np.linalg.inv(S)
            for j, zz in enumerate(z):
                r = zz - x
                d2 = float(r @ Si @ r)
                if d2 <= prm["tr_gate"]:
                    pairs.append((d2, a, j))
        pairs.sort()
        used_t, used_z = set(), set()
        for d2, a, j in pairs:
            if a in used_t or j in used_z:
                continue
            used_t.add(a); used_z.add(j)
            tr = live[a]; x, Pp = pred[a]
            S = Pp + Rm; K = Pp @ np.linalg.inv(S)
            tr["x"] = x + K @ (z[j] - x); tr["P"] = (np.eye(2) - K @ np.eye(2)) @ Pp
            tr["t"] = t; tr["miss"] = 0
            tr["hist"].append((i, t, z[j][0], z[j][1]))
            dets[j]["track"] = tr["id"]
        for a, tr in enumerate(live):
            if a not in used_t:
                tr["x"], tr["P"] = pred[a]; tr["t"] = t
                tr["miss"] += 1
                if tr["miss"] >= prm["tr_max_miss"]:
                    tr["dead"] = True
        for j, zz in enumerate(z):
            if j in used_z:
                continue
            tracks.append(dict(id=nid, x=zz.copy(), P=Rm * 4.0, t=t, miss=0, dead=False,
                               hist=[(i, t, zz[0], zz[1])]))
            dets[j]["track"] = nid
            nid += 1
        for tr in tracks:
            hits = [h[0] for h in tr["hist"] if h[0] > i - prm["tr_n"]]
            if len(hits) >= prm["tr_m"]:
                tr["confirmed"] = True

    for tr in tracks:
        h = np.array([(tt, R, V) for _, tt, R, V in tr["hist"]])
        tr["span"] = float(h[-1, 0] - h[0, 0]) if len(h) > 1 else 0.0
        tr["rdot"] = float(h[:, 2].mean())
        tr["range"] = float(h[:, 1].mean())
        tr["slope"] = float(np.polyfit(h[:, 0] - h[0, 0], h[:, 1], 1)[0]) if len(h) >= 3 and tr["span"] > 0 else 0.0
        bin_m = ms[0]["m_per_bin"]
        if not tr.get("confirmed") or tr["span"] < prm["tr_check_s"]:
            tr["class"] = "unverified"
        elif abs(tr["slope"] - tr["rdot"]) <= max(prm["tr_tol_ms"], prm["tr_tol_frac"] * abs(tr["rdot"])):
            tr["class"] = "mover" if abs(tr["rdot"]) * tr["span"] > 0.5 * bin_m else "unverified"
        elif abs(tr["rdot"]) * tr["span"] > bin_m:
            tr["class"] = "residue"
        else:
            tr["class"] = "unverified"
    cls = {tr["id"]: tr["class"] for tr in tracks}
    for m in ms:
        for k in m.get("markers", []):
            if k["verdict"] == tdd.VERDICT_TARGET:
                k["verdict"] = tdd.VERDICT_TARGET if cls.get(k.get("track")) == "mover" else tdd.VERDICT_UNCERTAIN
    return tracks


def write_tracks(path, tracks, t0):
    with open(path, "w") as fh:
        fh.write("track,class,t_start_s,t_end_s,span_s,hits,range_m,rdot_doppler_ms,range_slope_ms\n")
        for tr in tracks:
            if len(tr["hist"]) < 2:
                continue
            fh.write(f"{tr['id']},{tr['class']},{tr['hist'][0][1] - t0:.2f},{tr['hist'][-1][1] - t0:.2f},"
                     f"{tr['span']:.2f},{len(tr['hist'])},{tr['range']:.1f},{tr['rdot']:.2f},{tr['slope']:.2f}\n")


def plot_run(path, ms, tracks):
    """Range-time and Doppler-time of a whole run: the max excess over the background per
    range bin and per Doppler column, with C's detections and E's tracks on top."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    ok = [m for m in ms if "excess_db" in m]
    if not ok:
        return
    t0 = ms[0]["t_center"]
    tc = np.array([m["t_center"] - t0 for m in ok])
    nb, nf = ok[0]["excess_db"].shape
    lam = C / ok[0]["carrier"]
    sp = (-ok[0]["f_max"] + np.arange(nf) * (2.0 * ok[0]["f_max"] / (nf - 1))) * lam / 2.0
    vis = np.abs(sp) <= 3.0
    E = np.array([m["excess_db"] for m in ok])
    fig, ax = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    ax[0].pcolormesh(tc, np.arange(nb) * ok[0]["m_per_bin"], E[:, :, np.abs(sp) >= 0.1].max(2).T,
                     vmin=0, vmax=15, cmap="magma", shading="auto")
    ax[1].pcolormesh(tc, sp[vis], E[:, :, vis].max(1).T, vmin=0, vmax=15, cmap="magma", shading="auto")
    col = dict(mover="lime", residue="cyan", unverified="0.6")
    for tr in tracks:
        if len(tr["hist"]) < 2:
            continue
        h = np.array([(tt - t0, R) for _, tt, R, _ in tr["hist"]])
        spd = np.array([-V / 2.0 for _, _, _, V in tr["hist"]])   # Rdot = -lambda f, speed = f lambda / 2
        c = col[tr["class"]]
        lw = 2.0 if tr["class"] == "mover" else 1.0
        ax[0].plot(h[:, 0], h[:, 1], "-o", color=c, ms=2, lw=lw)
        ax[1].plot(h[:, 0], spd, "-o", color=c, ms=2, lw=lw)
    for m in ok:
        for k in m.get("markers", []):
            if k["verdict"] == tdd.VERDICT_REPLICA:
                continue
            ax[0].plot(m["t_center"] - t0, k["range_m"], ".", color="w", ms=2)
    ax[0].set_ylabel("bistatic excess range, m")
    ax[1].set_ylabel("speed, m/s (f lambda/2)")
    ax[1].set_xlabel("time, s")
    ax[0].set_title("excess over each cell's background (dB, 0-15); tracks: lime mover, cyan residue, grey unverified")
    plt.tight_layout()
    plt.savefig(path, dpi=80)
    plt.close(fig)


# ---------------------------------------------------------------- driver

def clamp(field, value):
    lo, hi = RANGE.get(field, (None, None))
    if lo is not None:
        value = max(lo, value)
    if hi is not None:
        value = min(hi, value)
    return int(round(value)) if field in INT_FIELDS else value


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rec", required=True, help="glob of <dump>.rec.*.csv files")
    ap.add_argument("--sweep", default=None, choices=sorted(RANGE),
                    help="field to vary; omit for a single run at the defaults/overrides")
    ap.add_argument("--lo", type=float, default=0.0)
    ap.add_argument("--hi", type=float, default=0.0)
    ap.add_argument("--step", type=float, default=1.0)
    ap.add_argument("--out-dir", default=None, help="write per-set map CSVs here")
    ap.add_argument("--no-antenna-avg", action="store_true")
    ap.add_argument("--no-layer-avg", action="store_true")
    ap.add_argument("--avg-maps", type=int, default=1,
                    help="sliding non-coherent average over K consecutive maps (default 1: off)")
    ap.add_argument("--no-tdd-detect", action="store_true",
                    help="skip the TDD detector (on by default, as tdd_detect: true)")
    ap.add_argument("--tdd-pfa", type=float, default=None,
                    help="detector CFAR false-alarm rate per cell (C default 1e-6)")
    ap.add_argument("--tdd-gamma", type=float, default=None,
                    help="detector replica power-drop requirement, Eq. (23) (C default 0.2)")
    ap.add_argument("--no-mirror-reject", action="store_true",
                    help="keep targets with an equally strong mirror at -v (mirror_reject: false)")
    ap.add_argument("--kfactor", action="store_true",
                    help="instead of maps: per-range-bin Rician K / slow fluctuation over sliding "
                         "windows of seconds (k_* parameters, sweepable), to <label>.kfactor.csv/.png")
    ap.add_argument("--vcomp", action="store_true",
                    help="with --avg-maps: shift each Doppler column by its path-length "
                         "rate before averaging, so movers stay aligned")
    ap.add_argument("--hygiene", action="store_true",
                    help="A: drop outlier snapshots and cut windows at scene steps before the maps (a_*)")
    ap.add_argument("--clutter-map", action="store_true",
                    help="C: detect each cell against its own history instead of the TDD detector (cm_*)")
    ap.add_argument("--track", action="store_true",
                    help="E: track C's detections and keep the kinematically consistent ones (tr_*); "
                         "implies --clutter-map")
    ap.add_argument("--residue", action="store_true", help="A + C + E together")
    ap.add_argument("--inject", default=None, metavar="DB,R0,AMP,PERIOD",
                    help="add a synthetic walker before everything: level re the direct path (dB), "
                         "range centre and swing (m, bistatic excess), swing period (s)")
    for fld, val in DEFAULTS.items():     # allow overriding any base parameter
        ap.add_argument("--" + fld.replace("_", "-"), type=float, default=None)
    args = ap.parse_args()
    if args.residue:
        args.hygiene = args.clutter_map = args.track = True
    if args.track:
        args.clutter_map = True

    # expand ~ and env vars ourselves: a quoted glob is not expanded by the shell, and
    # glob does not expand ~ either, so '~/runs/...' would otherwise match nothing.
    rec_pat = os.path.expanduser(os.path.expandvars(args.rec))
    files = sorted(glob.glob(rec_pat))
    if not files:
        print(f"no files match {rec_pat}", file=sys.stderr)
        return 1
    print(f"{len(files)} recorded maps")

    base = dict(DEFAULTS)
    for fld in DEFAULTS:
        ov = getattr(args, fld)
        if ov is not None:
            base[fld] = clamp(fld, ov)

    # parameter sets: a sweep of one field, or just the base
    if args.sweep:
        vals, v = [], args.lo
        while v <= args.hi + 1e-9:
            vals.append(v)
            v += args.step
        sets = []
        for v in vals:
            p = dict(base)
            p[args.sweep] = clamp(args.sweep, v)
            sets.append((f"{args.sweep}_{v:g}", p))
    else:
        sets = [("base", base)]

    # tdd_slots / max_speed come from the recordings (first file)
    _, tdd_slots, max_speed = parse_rec(files[0])
    for _, p in sets:
        p["tdd_slots"], p["max_speed"] = tdd_slots, max_speed

    # Default: write the maps next to the recordings, so a bare run reproduces the
    # normal-mode averaged map (one line per window, plottable with plot_range_doppler.py).
    if args.out_dir is None:
        args.out_dir = os.path.join(os.path.dirname(files[0]) or ".", "replay")
    args.out_dir = os.path.expanduser(os.path.expandvars(args.out_dir))
    os.makedirs(args.out_dir, exist_ok=True)
    if args.kfactor:
        data, m_per_bin = load_run(files)
        print(f"K-factor analysis of {len(files)} recordings, chains {sorted(data)}, to {args.out_dir}/")
        for label, prm in sets:
            kfactor_run(data, m_per_bin, prm, label, args.out_dir)
        return 0
    print(f"writing maps to {args.out_dir}/ (same pipeline as normal mode, without the "
          f"spatial null and AoA)")

    suffix = ""
    if args.avg_maps > 1:
        suffix = f"_avg{args.avg_maps}" + ("_vcomp" if args.vcomp else "")

    parsed_all = [parse_rec(path)[0] for path in files]
    if args.inject:
        inject_walker(parsed_all, args.inject)
        suffix += "_inj" + args.inject.replace(",", "_")
    if args.hygiene:
        suffix += "_A"
    if args.clutter_map:
        suffix += "_C"
    if args.track:
        suffix += "_E"
    for label, prm in sets:
        label += suffix
        parsed = parsed_all
        steps_abs = []
        if args.hygiene:
            parsed, st = hygiene(parsed_all, prm)
            steps_abs = st["steps_abs"]
            print(f"{label}: hygiene dropped {st['outliers']}/{st['total']} outlier snapshots "
                  f"({100.0 * st['outliers'] / max(st['total'], 1):.1f}%), scene steps at "
                  f"{', '.join(f'{x:.2f}' for x in st['steps']) or 'none'} s "
                  f"({st['dropped_side']} snapshots cut off on the short side)")
        # every window on the first window's Doppler grid, so maps can be averaged
        ms, nf_run = [], 0
        for groups in parsed:
            if not groups:
                continue
            m = build_averaged(groups, prm, not args.no_antenna_avg, not args.no_layer_avg, nf_run)
            if m is None:
                continue
            nf_run = nf_run or m["n_freq"]
            ms.append(m)
        ms = average_maps(ms, args.avg_maps, args.vcomp)
        if not ms:
            print(f"{label}: no maps built")
            continue
        tracks = None
        if args.clutter_map:
            clutter_map(ms, prm, steps_abs)
            if args.track:
                tracks = track(ms, prm)
        elif not args.no_tdd_detect:
            cfg = tdd.cfg_default()
            if args.tdd_pfa is not None:
                cfg["pfa"] = args.tdd_pfa
            if args.tdd_gamma is not None:
                cfg["gamma"] = args.tdd_gamma
            for m in ms:
                run_detector(m, not args.no_mirror_reject, cfg)
        maps = [m["power"] for m in ms]
        if args.out_dir:
            with open(os.path.join(args.out_dir, f"{label}.csv"), "w") as fh:
                for i, m in enumerate(ms):
                    write_map_line(fh, i, m)
            detected = args.clutter_map or not args.no_tdd_detect
            if detected:
                write_targets(os.path.join(args.out_dir, f"{label}.targets.csv"), ms)
            if tracks is not None:
                write_tracks(os.path.join(args.out_dir, f"{label}.tracks.csv"), tracks, ms[0]["t_center"])
            if args.clutter_map:
                plot_run(os.path.join(args.out_dir, f"{label}.time.png"), ms, tracks or [])
        peaks = np.array([10.0 * math.log10(max(p.max(), 1e-12)) for p in maps])
        # sliding outputs share K-1 of their K maps, so compare only outputs K apart
        rep = repeatability(maps[::max(1, args.avg_maps)])
        line = (f"{label:28s}  maps {len(maps):3d}  peak SNR {peaks.mean():6.1f} +- "
                f"{peaks.std():4.1f} dB   repeatability {rep:.3f}")
        if args.clutter_map or not args.no_tdd_detect:
            n_t = np.array([sum(t["verdict"] == tdd.VERDICT_TARGET for t in m["markers"]) for m in ms])
            n_u = sum(sum(t["verdict"] == tdd.VERDICT_UNCERTAIN for t in m["markers"]) for m in ms)
            n_r = sum(sum(t["verdict"] == tdd.VERDICT_REPLICA for t in m["markers"]) for m in ms)
            line += (f"   targets/map {n_t.mean():.2f} (maps with one: {np.mean(n_t > 0)*100:.0f}%)"
                     f"  replicas {n_r}  uncertain {n_u}")
        if tracks is not None:
            by = {}
            for tr in tracks:
                if len(tr["hist"]) >= 2:
                    by[tr["class"]] = by.get(tr["class"], 0) + 1
            line += "   tracks " + " ".join(f"{k} {v}" for k, v in sorted(by.items()))
        print(line)
    return 0


def repeatability(maps):
    """Mean correlation of consecutive maps in dB: higher is more repeatable. Only maps
    of the same shape are compared."""
    cors = []
    for a, b in zip(maps, maps[1:]):
        if a.shape != b.shape:
            continue
        x = 10.0 * np.log10(np.maximum(a, 1e-12)).ravel()
        y = 10.0 * np.log10(np.maximum(b, 1e-12)).ravel()
        x -= x.mean(); y -= y.mean()
        d = math.sqrt((x @ x) * (y @ y))
        if d > 0:
            cors.append((x @ y) / d)
    return float(np.mean(cors)) if cors else float("nan")


if __name__ == "__main__":
    sys.exit(main())
