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
)

# field -> (lo, hi) hard range, for clamping a swept value (mirrors nr_sweep_apply)
RANGE = dict(
    clutter_mode=(0, 2), max_paths=(1, 4), trend_degree=(0, SLOW_TREND_MAX),
    kernel_half_span=(1, None), static_min=(0.0, None), snr_min=(0.0, None),
    min_sep_bins=(1, None), los_first_db=(0.0, None), comb_remove=(0, 1),
    comb_harmonic=(1, COMB_MAX_HARMONIC), comb_tdd_multiple=(1, None),
)
INT_FIELDS = {"clutter_mode", "max_paths", "trend_degree", "kernel_half_span",
              "min_sep_bins", "comb_remove", "comb_harmonic", "comb_tdd_multiple",
              "comb_power_iters"}


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

    # gain + slow-trend removal, per bin (always runs; also applies the gain for mode 0/1)
    trend_degree = -1 if mode == 0 else prm["trend_degree"]
    Q = slow_basis(t, trend_degree)                 # (n_q, n)
    z = gain[:, None] * res                         # (n, n_bins)
    for q in Q:
        z = z - np.outer(q, q @ z)                  # remove projection on q, per bin
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
                trend_degree=trend_degree, tdd_period=tdd_period)


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


# ---------------------------------------------------------------- detection

def run_detector(m, mirror, cfg=None):
    """The TDD detector on the window's rx0 conditioned samples, on the map's own Doppler
    grid, then the mirror rejection on the (possibly averaged) map, as map_task does with
    tdd_detect and mirror_reject on. Attaches m["markers"]: targets first, then rejections."""
    obs = tdd.Obs(h=m["res"], t=m["t"], idft_size=m["idft"], win_len=m["win_len"],
                  win_start=m["win_start"], win=hann(m["win_len"]), m_per_bin=m["m_per_bin"],
                  n_freq=m["n_freq"], f_max_hz=m["f_max"], lambda_m=C / m["carrier"],
                  t_tdd_s=m["tdd_period"], trend_q=slow_basis(m["t"], m["trend_degree"]))
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
    ap.add_argument("--vcomp", action="store_true",
                    help="with --avg-maps: shift each Doppler column by its path-length "
                         "rate before averaging, so movers stay aligned")
    for fld, val in DEFAULTS.items():     # allow overriding any base parameter
        ap.add_argument("--" + fld.replace("_", "-"), type=float, default=None)
    args = ap.parse_args()

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
    print(f"writing maps to {args.out_dir}/ (same pipeline as normal mode, without the "
          f"spatial null and AoA)")

    suffix = ""
    if args.avg_maps > 1:
        suffix = f"_avg{args.avg_maps}" + ("_vcomp" if args.vcomp else "")

    for label, prm in sets:
        label += suffix
        # every window on the first window's Doppler grid, so maps can be averaged
        ms, nf_run = [], 0
        for path in files:
            groups, _, _ = parse_rec(path)
            m = build_averaged(groups, prm, not args.no_antenna_avg, not args.no_layer_avg, nf_run)
            if m is None:
                continue
            nf_run = nf_run or m["n_freq"]
            ms.append(m)
        ms = average_maps(ms, args.avg_maps, args.vcomp)
        if not ms:
            print(f"{label}: no maps built")
            continue
        if not args.no_tdd_detect:
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
            if not args.no_tdd_detect:
                write_targets(os.path.join(args.out_dir, f"{label}.targets.csv"), ms)
        peaks = np.array([10.0 * math.log10(max(p.max(), 1e-12)) for p in maps])
        # sliding outputs share K-1 of their K maps, so compare only outputs K apart
        rep = repeatability(maps[::max(1, args.avg_maps)])
        line = (f"{label:28s}  maps {len(maps):3d}  peak SNR {peaks.mean():6.1f} +- "
                f"{peaks.std():4.1f} dB   repeatability {rep:.3f}")
        if not args.no_tdd_detect:
            n_t = np.array([sum(t["verdict"] == tdd.VERDICT_TARGET for t in m["markers"]) for m in ms])
            n_u = sum(sum(t["verdict"] == tdd.VERDICT_UNCERTAIN for t in m["markers"]) for m in ms)
            n_r = sum(sum(t["verdict"] == tdd.VERDICT_REPLICA for t in m["markers"]) for m in ms)
            line += (f"   targets/map {n_t.mean():.2f} (maps with one: {np.mean(n_t > 0)*100:.0f}%)"
                     f"  replicas {n_r}  uncertain {n_u}")
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
