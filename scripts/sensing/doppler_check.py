#!/usr/bin/env python3
"""
Doppler profiles from the raw snapshot dump, computed here rather than in the UE.

nr_ue_sensing_range_doppler() turns the slow-time snapshots into a range-Doppler map
and keeps only |C|^2, so the map dump says what the pipeline decided but not why. The
snapshot dump written by nr_ue_sensing_dump_snapshots() is the input of that function,
before any of it ran, which is enough to redo the whole chain offline:

    gain        every snapshot scaled to the widest grant of the window, n_max / n_m
    clutter     none, the slow trend alone (mean), or the kernel fit then the trend
    Doppler     direct DFT at the recorded instants, the sampling being non-uniform

Each stage is a switch here, so the effect of one of them is a difference between two
curves on the same axes instead of two runs of the modem. --compare draws the
conditioning variants over each other for exactly that.

WHAT THE PROFILES SHOW

A profile is one row of the map: the Doppler spectrum at a fixed delay bin. Reading it
rather than the 2-D map is what makes the floor next to a peak visible, and that floor
is what decides whether a target a few dB up is detectable at all.

Three things put energy off zero Doppler that no target accounts for, and all three
are meant to be read here:

  TDD replicas    the reference symbols only arrive in the DL part of each period, so
                  the spectrum repeats every 1/T_TDD (200 Hz, 8.8 m/s at 3.41 GHz).
                  --replicas marks where the copies of each peak have to be.
  clutter residue whatever the clutter stage left at zero Doppler is copied to every
                  replica position as well, which is why the trend degree matters far
                  from 0 Hz and not only next to it.
  grant motion    a changing n_m or a_m modulates the amplitude and the phase of a
                  static path in slow time, and the transform reads that modulation as
                  motion. The bottom panel draws n_m and a_m against time next to the
                  profiles, so a feature can be checked against the grant that produced
                  it before it is taken for a target.

THE DOPPLER GRID

Built exactly as the C builds it, so a profile lands on the same points as the dumped
map: f_max is --max-speed converted to Hz plus one TDD period of margin, and the grid
holds 2 points per resolution cell of 1/t_span. The C then clamps the point count to
NR_SENSING_MAP_MAX_BINS_FREQ; nothing here does, so --oversample buys a smooth curve where
the map has a handful of samples across a mainlobe. Oversampling adds no resolution,
only points: 1/t_span is still the width of the narrowest feature.

Usage:
    ./doppler_check.py /tmp/map.csv.snap.csv
    ./doppler_check.py /tmp/map.csv.snap.csv --bin los --bin peak --oversample 4
    ./doppler_check.py /tmp/map.csv.snap.csv --compare none,mean,kernel
    ./doppler_check.py /tmp/map.csv.snap.csv --trend 0,1,2 --bin 30m --replicas
    ./doppler_check.py /tmp/map.csv.snap.csv --map --max-range 150

Bins are given as 'los' (the direct path), 'peak' (the strongest cell away from it),
a bin index, or a path length with an 'm' suffix.
"""

import argparse
import sys

import numpy as np

from model_check import A_L, C0, Snapshots

# nr_ue_sensing.h
TDD_PERIOD_SLOTS = 10

# nr_ue_map.h
CLUTTER_MAX_PATHS = 4
CLUTTER_STATIC_MIN = 0.8
CLUTTER_SNR_MIN = 10.0
CLUTTER_MIN_SEP_BINS = 3
CLUTTER_KERNEL_HALF_SPAN = 48
SLOW_TREND_DEGREE = 2


# ---------------------------------------------------------------------------
# the conditioning, as nr_ue_sensing_range_doppler() does it
# ---------------------------------------------------------------------------


def peak_frac(p, b):
    """Sub-bin offset of the peak at bin b, from a parabola through it in dB.

    nr_ue_sensing_peak_frac(). Fitted on the log of the profile because a Hann
    mainlobe is very nearly a parabola there and only roughly one on a linear scale.
    Zero at the edges, on a plateau, and when the vertex lands further than half a bin
    away, which a true local maximum cannot produce.
    """
    if b <= 0 or b >= len(p) - 1:
        return 0.0
    if p[b - 1] <= 0.0 or p[b] <= 0.0 or p[b + 1] <= 0.0:
        return 0.0
    ym, y0, yp = 10.0 * np.log10(p[b - 1:b + 2])
    den = ym - 2.0 * y0 + yp
    if den >= 0.0:
        return 0.0
    d = 0.5 * (ym - yp) / den
    return float(d) if -0.5 < d < 0.5 else 0.0


def los_bin(snap):
    """Fractional bin of the direct path: the bin holding the most energy.

    The direct path is static and tens of dB above every reflection, so the strongest
    bin of the incoherent profile is it. The sub-bin part is worth keeping: at 2.44
    m/bin the rounding alone is over a metre.
    """
    e = np.sum(np.abs(snap.h) ** 2, axis=0)
    b = int(np.argmax(e))
    return b + peak_frac(e, b)


def kernel_K(snap, u, b_lo, b_hi):
    """K_m(b - u) over bins b_lo..b_hi, one row per snapshot.

    K_m(x) = A_{n_m}(x) * exp(j2pi c_m x / N) with c_m = a_m + (n_m - 1)/2, which is
    nr_ue_sensing_clutter_kernel() written through the zero-phase transform instead of
    the pilot sum. Built once per distinct (n_m, a_m) because that pair is all it
    depends on, exactly as the C pools its kernels.
    """
    x = np.arange(b_lo, b_hi + 1, dtype=float) - u
    K = np.empty((snap.M, x.size), dtype=complex)
    grants = np.stack([snap.n_m, snap.a_m.astype(int)], axis=1)
    for n_val, a_val in np.unique(grants, axis=0):
        rows = (snap.n_m == n_val) & (snap.a_m.astype(int) == a_val)
        c_m = a_val + (n_val - 1) / 2.0
        K[rows] = A_L(int(n_val), x, snap.N) * np.exp(1j * 2 * np.pi * c_m * x / snap.N)
    return K


def remove_path(snap, res, u):
    """Fit one static path at the fractional bin u and subtract it, in place.

    alpha = sum_m sum_b r_m[b] conj(K_m(b-u)) / sum |K_m|^2: one amplitude for the
    whole window, not one per snapshot. That is what keeps this a clutter filter. K_m
    carries the grant variation while alpha stays constant, so a mover at the same
    delay survives, where a per-snapshot amplitude would delete everything there
    whatever its Doppler.

    Only bins within CLUTTER_KERNEL_HALF_SPAN of u are fitted and touched: beyond that
    the Hann skirt is some 100 dB down, below anything the c16_t response holds.
    """
    b_lo = max(0, int(np.floor(u)) - CLUTTER_KERNEL_HALF_SPAN)
    b_hi = min(snap.n_bins - 1, int(np.ceil(u)) + CLUTTER_KERNEL_HALF_SPAN)
    if b_hi < b_lo:
        return
    K = kernel_K(snap, u, b_lo, b_hi)
    r = res[:, b_lo:b_hi + 1]
    den = np.sum(np.abs(K) ** 2)
    if den <= 0.0:
        return
    alpha = np.sum(r * np.conj(K)) / den
    res[:, b_lo:b_hi + 1] = r - alpha * K


def slow_basis(t, degree):
    """Orthonormal basis of the slow trend on the sample times, as q_0..q_K.

    nr_ue_sensing_slow_basis(): the polynomials 1, tau, ..., tau^degree with tau the
    times centred and scaled to about [-1, 1] for conditioning, made orthonormal by
    Gram-Schmidt. Degree 0 alone gives q_0 = 1/sqrt(n), so the projection below is
    exactly the slow-time mean.
    """
    if degree < 0 or len(t) == 0:
        return np.zeros((0, len(t)))
    t_mid = t.mean()
    t_half = 0.5 * (t[-1] - t[0]) if t[-1] > t[0] else 1.0
    tau = (t - t_mid) / t_half

    q = []
    for k in range(degree + 1):
        v = tau ** k
        for qj in q:
            v = v - np.dot(qj, v) * qj
        norm = np.linalg.norm(v)
        if norm < 1e-9:
            continue  # not independent of the lower degrees on these sample times
        q.append(v / norm)
    return np.array(q) if q else np.zeros((0, len(t)))


def condition(snap, clutter="kernel", trend=SLOW_TREND_DEGREE, max_paths=CLUTTER_MAX_PATHS,
              gain=True, verbose=False):
    """The residual the Doppler transform works on, and what was removed to get it.

    Mirrors nr_ue_sensing_range_doppler() stage for stage:

      1. the pilot-count gain n_max / n_m, so a narrow grant is not read as a dip in
         slow time. Corrects the height of a peak, not its width.
      2. clutter 'kernel': the direct path is fitted with each snapshot's own point
         spread function and removed, then further static paths one at a time, while
         the strongest peak left in the slow-time mean is both static and above the
         noise.
      3. the slow trend of degree `trend` removed from every bin, on the normalised
         snapshots. This is the whole of clutter 'mean'.

    Returns (res, info). res is [M, n_bins] with the gain applied, as the C hands it
    to the transform.
    """
    res = snap.h.copy()
    g = (float(snap.n_m.max()) / snap.n_m) if gain else np.ones(snap.M)
    info = {"bin_los": los_bin(snap), "paths": [], "gain": g}

    if clutter == "kernel":
        # Path 0: the direct path, always, at the refined fractional bin.
        remove_path(snap, res, info["bin_los"])
        info["paths"].append((info["bin_los"], np.nan, np.nan))

        while len(info["paths"]) < max_paths:
            z = g[:, None] * res
            coh = np.abs(z.mean(axis=0)) ** 2                # |slow-time mean|^2
            inc = np.mean(np.abs(z) ** 2, axis=0)            # slow-time mean of |r|^2

            # strongest peak of the mean, away from the paths already fitted
            free = np.ones(snap.n_bins, dtype=bool)
            for u_l, _, _ in info["paths"]:
                lo = max(0, int(np.ceil(u_l - CLUTTER_MIN_SEP_BINS + 1)))
                hi = int(np.floor(u_l + CLUTTER_MIN_SEP_BINS - 1))
                free[lo:hi + 1] = False
            if not free.any():
                break
            b = int(np.argmax(np.where(free, coh, -np.inf)))

            # noise per snapshot from the median bin, which the few paths do not move;
            # averaging M snapshots divides it by M in the mean
            noise_mean = np.median(inc) / snap.M
            s_static = coh[b] / inc[b] if inc[b] > 0.0 else 0.0
            if s_static < CLUTTER_STATIC_MIN or coh[b] < CLUTTER_SNR_MIN * noise_mean:
                break

            u = b + peak_frac(coh, b)
            remove_path(snap, res, u)
            snr = 10.0 * np.log10(coh[b] / noise_mean)
            info["paths"].append((u, s_static, snr))
            if verbose:
                print(f"  static path {len(info['paths']) - 1} removed at bin {u:6.2f} "
                      f"({u * snap.m_per_bin:6.1f} m), S {s_static:.2f}, {snr:5.1f} dB "
                      f"above noise")
    elif clutter not in ("none", "mean"):
        raise ValueError(f"unknown clutter mode '{clutter}'")

    # Stage 2, and the whole of the mean mode. Normalised first: a trend fitted to
    # unnormalised snapshots carries the very amplitude modulation it is there to
    # cancel. 'none' keeps the gain but removes nothing, as trend_degree -1 does in C.
    res = g[:, None] * res
    q = slow_basis(snap.t, -1 if clutter == "none" else trend)
    info["n_q"] = len(q)
    for qj in q:
        res = res - np.outer(qj, qj @ res)
    return res, info


# ---------------------------------------------------------------------------
# the Doppler axis
# ---------------------------------------------------------------------------


def doppler_grid(snap, max_speed, oversample=1, f_max=None, n_freq=None):
    """The grid the C would build, optionally with more points on it.

    f_max reaches the fastest target asked for plus one TDD period of margin, which
    the TDD detector needs to see the replicas of a target at the edge of the axis.
    The point count is 2 per resolution cell of 1/t_span, times --oversample, and odd
    so that 0 Hz falls exactly on a point.
    """
    t_span = float(snap.t[-1] - snap.t[0])
    if f_max is None:
        f_max = snap.doppler_hz(max_speed) + 1.0 / snap.t_tdd
    if n_freq is None:
        n_freq = max(3, int(4.0 * f_max * t_span * oversample) | 1)
    return np.linspace(-f_max, f_max, n_freq), f_max, t_span


def doppler_power(res, t, freqs):
    """|sum_m res[m,b] exp(-j2pi f t_m)|^2 -> [n_bins, n_freq].

    A direct DFT at the recorded instants, not an FFT over a uniform grid: the
    reference symbols are not uniformly spaced, and pretending they are is what puts
    a target at the wrong speed.
    """
    e = np.exp(-2j * np.pi * np.outer(freqs, t))     # [F, M]
    return np.abs(e @ res).T ** 2                    # [B, F]


# ---------------------------------------------------------------------------
# picking the bins to profile
# ---------------------------------------------------------------------------


def parse_bin(spec):
    """'los', 'peak', a bin index, or a path length with an 'm' suffix."""
    s = str(spec).strip().lower()
    if s in ("los", "peak"):
        return s
    if s.endswith("m"):
        return ("m", float(s[:-1]))
    return ("b", int(float(s)))


def resolve_bins(specs, snap, power, bin_los):
    """Bin specs -> [(bin index, label)], in the order given."""
    out = []
    for spec in specs:
        if spec == "los":
            b = int(round(bin_los))
            out.append((b, f"LoS, bin {b}"))
        elif spec == "peak":
            # strongest cell away from the direct path: what the clutter stage left
            free = np.abs(np.arange(snap.n_bins) - bin_los) >= CLUTTER_MIN_SEP_BINS
            b = int(np.argmax(np.where(free, power.max(axis=1), -np.inf)))
            out.append((b, f"strongest off-LoS, bin {b}"))
        elif spec[0] == "m":
            b = int(np.clip(round(spec[1] / snap.m_per_bin), 0, snap.n_bins - 1))
            out.append((b, f"bin {b}"))
        else:
            b = int(np.clip(spec[1], 0, snap.n_bins - 1))
            out.append((b, f"bin {b}"))
    # keep the first mention of a bin, so 'los peak' does not draw one bin twice
    seen, uniq = set(), []
    for b, label in out:
        if b not in seen:
            seen.add(b)
            uniq.append((b, label))
    return uniq


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------


def plot_grants(ax, snap):
    """n_m and a_m against time, on the same axis as the profiles above.

    The panel the profiles have to be read against: the transform cannot tell a
    genuine Doppler from a grant that moved, so a feature off zero Doppler means one
    thing when these two traces are flat and another when they are not.
    """
    ax.step(snap.t * 1e3, snap.n_m, where="post", lw=1.2, color="tab:blue",
            label="n_m, pilots")
    ax.set_xlabel("time [ms]")
    ax.set_ylabel("n_m [pilots]", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")

    lo, hi = int(snap.n_m.min()), int(snap.n_m.max())
    ax.set_ylim(0, hi * 1.15)
    ax.grid(alpha=0.3)

    ax2 = ax.twinx()
    ax2.step(snap.t * 1e3, snap.a_m, where="post", lw=1.2, color="tab:red",
             label="a_m, lattice start")
    ax2.set_ylabel("a_m [lattice index]", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    a_lo, a_hi = snap.a_m.min(), snap.a_m.max()
    ax2.set_ylim(a_lo - 1 if a_lo == a_hi else a_lo, a_hi + 1 if a_lo == a_hi else a_hi)

    n_grants = len(np.unique(np.stack([snap.n_m, snap.a_m.astype(int)], axis=1), axis=0))
    fixed = "fixed" if n_grants == 1 else f"{n_grants} distinct (n_m, a_m)"
    ax.set_title(f"grant over the window: {fixed},  n_m {lo}..{hi},  "
                 f"a_m {a_lo:.0f}..{a_hi:.0f},  {snap.M} snapshots", fontsize=9)

    # The sampling pattern itself, at the bottom: DMRS symbols land in bursts inside a
    # slot and only in the DL part of the TDD period, which is what puts the replicas
    # on the profiles above.
    ax.plot(snap.t * 1e3, np.full(snap.M, hi * 0.06), "|", ms=6, color="0.35",
            alpha=0.7, label="snapshot instants")
    ax.legend(loc="upper left", fontsize=7, framealpha=0.7)


def plot_profiles(ax, snap, powers, names, bins, ref_db, args, x_speed):
    """Every (bin, variant) curve on one axis: the point is to compare them.

    Colour is the bin and dash pattern the variant, because two conditioning variants
    often land within a dB of each other over most of the axis and differ only where
    it matters. Drawn in reverse so the first variant given, usually the reference one,
    stays on top of the ones meant to be judged against it.
    """
    import matplotlib.pyplot as plt  # the backend is already chosen by main()

    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    styles = ["-", "--", ":", "-."]
    for i, (b, label) in enumerate(bins):
        col = colours[i % len(colours)]
        for j in reversed(range(len(powers))):
            y = 10.0 * np.log10(powers[j][b] + 1e-30) - ref_db
            tag = f"{label}, {b * snap.m_per_bin:.1f} m"
            ax.plot(x_speed, y, lw=1.4, color=col, ls=styles[j % len(styles)],
                    alpha=1.0 if j == 0 else 0.85,
                    label=tag if names[j] is None else f"{tag}  [{names[j]}]")

    ax.axvline(0.0, color="k", alpha=0.3, lw=0.8)
    if args.replicas:
        # The spectrum repeats every 1/T_TDD, so a real target and its copies are
        # indistinguishable on one map. Drawn from 0 Hz: the copies of a peak sit at
        # this spacing from it, wherever the peak is.
        step = (1.0 / snap.t_tdd) * (snap.lambda_m / 2.0 if args.x == "speed" else 1.0)
        k = 1
        while k * step <= x_speed.max():
            for s in (-1, 1):
                ax.axvline(s * k * step, color="crimson", ls=":", lw=0.8, alpha=0.5,
                           label="TDD replica spacing" if (k == 1 and s == 1) else None)
            k += 1

    ax.set_xlabel("speed [m/s]" if args.x == "speed" else "doppler [Hz]")
    ax.set_ylabel("dB below the unconditioned peak")
    ax.set_ylim(args.floor, 3.0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2, framealpha=0.7)


def plot_map(ax, snap, power, ref_db, args, x_speed, bins):
    y = np.arange(snap.n_bins) * snap.m_per_bin
    keep = slice(None)
    if args.max_range is not None:
        keep = slice(0, min(int(np.searchsorted(y, args.max_range)) + 1, snap.n_bins))
    img = 10.0 * np.log10(power + 1e-30) - ref_db
    mesh = ax.pcolormesh(x_speed, y[keep], img[keep, :], vmin=args.floor, vmax=3.0,
                         shading="nearest", cmap="viridis")
    ax.figure.colorbar(mesh, ax=ax, label="dB below the unconditioned peak")
    for b, _ in bins:
        ax.axhline(b * snap.m_per_bin, color="white", lw=0.8, ls="--", alpha=0.7)
    ax.set_xlabel("speed [m/s]" if args.x == "speed" else "doppler [Hz]")
    ax.set_ylabel("path length [m]")
    ax.set_title("the map these profiles are rows of", fontsize=9)


# ---------------------------------------------------------------------------


def parse_list(value, cast):
    return [cast(p.strip()) for p in str(value).split(",") if p.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="file written by nr_ue_sensing_dump_snapshots()")

    c = ap.add_argument_group("conditioning, as nr_ue_sensing_range_doppler() does it")
    c.add_argument("--clutter", default="kernel", choices=("none", "mean", "kernel"),
                   help="how the static scene is removed (default: kernel, the C default)")
    c.add_argument("--compare", metavar="MODES",
                   help="draw several clutter modes over each other, e.g. "
                        "--compare none,mean,kernel. Overrides --clutter.")
    c.add_argument("--trend", default=str(SLOW_TREND_DEGREE), metavar="DEG",
                   help=f"degree of the slow trend removed from every bin "
                        f"(default {SLOW_TREND_DEGREE}, 0 is the plain mean). "
                        f"Comma-separated draws one curve per degree.")
    c.add_argument("--max-paths", type=int, default=CLUTTER_MAX_PATHS,
                   help=f"static paths the kernel mode fits at most, the direct path "
                        f"included (default {CLUTTER_MAX_PATHS}, 1 for the direct path alone)")
    c.add_argument("--no-gain", action="store_true",
                   help="skip the pilot-count normalisation, to see what a varying "
                        "grant does to the Doppler axis when nothing corrects it")

    d = ap.add_argument_group("the Doppler axis")
    d.add_argument("--max-speed", type=float, default=10.0, metavar="M_S",
                   help="fastest target the axis has to cover, as --sensing-max-speed "
                        "(default 10). One TDD period of margin is added, as in C.")
    d.add_argument("--oversample", type=int, default=1, metavar="K",
                   help="points per resolution cell, over the 2 the C uses. Buys a "
                        "smooth curve, not resolution: 1/t_span is still the width of "
                        "the narrowest feature.")
    d.add_argument("--f-max", type=float, default=None, metavar="HZ",
                   help="set the edge of the axis directly, bypassing --max-speed")
    d.add_argument("--x", choices=("speed", "doppler"), default="speed",
                   help="horizontal axis unit (default: speed)")

    b = ap.add_argument_group("what to profile")
    b.add_argument("--bin", action="append", default=[], metavar="SPEC",
                   help="'los', 'peak', a bin index, or a path length like '30m'. "
                        "Repeatable; default: los and peak.")
    b.add_argument("--max-bin", type=int, default=None,
                   help="bins to keep from the dump (default: all of them)")

    p = ap.add_argument_group("output")
    p.add_argument("--map", action="store_true",
                   help="also draw the 2-D map the profiles are rows of")
    p.add_argument("--replicas", action="store_true",
                   help="mark the TDD replica spacing, where the copies of a peak sit")
    p.add_argument("--no-grants", action="store_true",
                   help="drop the n_m / a_m panel")
    p.add_argument("--max-range", type=float, default=None,
                   help="clip the range axis of --map, in metres")
    p.add_argument("--floor", type=float, default=-60.0,
                   help="lowest dB shown (default: -60)")
    p.add_argument("--fc", type=float, default=None, metavar="HZ",
                   help="override the carrier of the dump, used for the speed axis")
    p.add_argument("--plot", metavar="PNG", help="write a figure instead of showing it")
    args = ap.parse_args()

    snap = Snapshots(args.csv, max_bin=args.max_bin, fc=args.fc)
    # nr_ue_sensing.h has 10 slots per TDD period; Snapshots derives t_tdd from its own
    # constant, so it is set here from the header value the map was built with.
    mu = np.log2(snap.scs / 15000.0)
    snap.t_tdd = TDD_PERIOD_SLOTS * 1e-3 / (2.0 ** mu)
    snap.describe()

    freqs, f_max, t_span = doppler_grid(snap, args.max_speed, args.oversample, args.f_max)
    print(f"doppler axis   : +-{f_max:.1f} Hz (+-{snap.speed_ms(f_max):.2f} m/s), "
          f"{len(freqs)} points, {freqs[1] - freqs[0]:.2f} Hz apart")
    print(f"resolution     : {1.0 / t_span:.2f} Hz ({snap.speed_ms(1.0 / t_span):.3f} m/s) "
          f"over {t_span * 1e3:.1f} ms")
    print(f"TDD replicas   : every {1.0 / snap.t_tdd:.1f} Hz "
          f"({snap.speed_ms(1.0 / snap.t_tdd):.2f} m/s)")

    # The reference every curve is measured against: the peak of the unconditioned map,
    # gain applied. One number per file, independent of the variant, so the variants can
    # be compared on one axis and the clutter stage is visibly what pushes things down.
    res0, info0 = condition(snap, clutter="none", gain=not args.no_gain)
    ref_db = 10.0 * np.log10(doppler_power(res0, snap.t, freqs).max() + 1e-30)
    print(f"LoS bin        : {info0['bin_los']:.2f} "
          f"({info0['bin_los'] * snap.m_per_bin:.1f} m)")

    # One variant per (clutter mode, trend degree). A single pair leaves the curves
    # unlabelled, which is the common case and reads better.
    modes = parse_list(args.compare, str) if args.compare else [args.clutter]
    degrees = parse_list(args.trend, int)
    for m in modes:
        if m not in ("none", "mean", "kernel"):
            sys.exit(f"unknown clutter mode '{m}' in --compare")
    variants = [(m, deg) for m in modes for deg in degrees]
    named = len(variants) > 1

    powers = []
    for mode, deg in variants:
        print(f"\nclutter '{mode}', trend degree {deg if mode != 'none' else -1}:")
        res, info = condition(snap, clutter=mode, trend=deg, max_paths=args.max_paths,
                              gain=not args.no_gain, verbose=True)
        if mode == "kernel":
            print(f"  {len(info['paths'])} static path(s) removed, "
                  f"{info['n_q']} trend vector(s) projected out")
        else:
            print(f"  {info['n_q']} trend vector(s) projected out")
        powers.append(doppler_power(res, snap.t, freqs))

    specs = [parse_bin(s) for s in (args.bin or ["los", "peak"])]
    bins = resolve_bins(specs, snap, powers[0], info0["bin_los"])
    print("\nprofiled bins  : " + ", ".join(f"{b} ({b * snap.m_per_bin:.1f} m)"
                                            for b, _ in bins))

    import matplotlib
    if args.plot:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_panels = 1 + int(args.map) + int(not args.no_grants)
    heights = [3.0] + ([3.2] if args.map else []) + ([1.6] if not args.no_grants else [])
    fig, axes = plt.subplots(n_panels, 1, figsize=(11, 3.1 * sum(heights) / 3.0),
                             gridspec_kw=dict(height_ratios=heights), squeeze=False)
    axes = axes[:, 0]

    x = snap.speed_ms(freqs) if args.x == "speed" else freqs
    names = [f"{m}, K={d}" if named else None for m, d in variants]
    plot_profiles(axes[0], snap, powers, names, bins, ref_db, args, x)
    axes[0].set_title(f"Doppler profiles, {snap.M} snapshots over {t_span * 1e3:.1f} ms, "
                      f"{len(freqs)} points on the axis", fontsize=9)

    k = 1
    if args.map:
        plot_map(axes[k], snap, powers[0], ref_db, args, x, bins)
        k += 1
    if not args.no_grants:
        plot_grants(axes[k], snap)

    fig.tight_layout()
    if args.plot:
        fig.savefig(args.plot, dpi=130)
        print(f"\nwrote {args.plot}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
