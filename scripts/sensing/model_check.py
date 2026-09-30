#!/usr/bin/env python3
"""
Check the delay-response model against the measured snapshots.

The model, equation (10) of doc/sensing_model.tex:

    h_m[b] = sum_l  a_l * exp(j2pi nu_l t_m) * A_{n_m}(b - u_l)
                        * exp(j2pi c_m (b - u_l) / N)

with

    u_l  = tau_l * N * k_step * SCS        path delay, in bins
    a_m  = (k_first - k_first % k_step) / k_step
    c_m  = a_m + (n_m - 1) / 2             centroid of the measured span
    A_L  = zero-phase transform of the length-L Hann taper

Everything grant-dependent is in A_{n_m} and c_m, and both are read per snapshot
from the CSV, so nothing here is assumed constant.

GROUND TRUTH FROM THE RELAY

--target takes the SAME syntax as zmq_channel_relay.py, so the arguments can be
copied across unchanged, along with --fc and --distance-is:

    relay: --target 0:0:0dB --target 30.5:8.78:-40dB --target 74.4:-6.58:0.2
    here:  --target 0:0:0dB --target 30.5:8.78:-40dB --target 74.4:-6.58:0.2

Amplitudes are linear, or a level relative to the unit-gain LOS with a dB suffix. An
optional fourth field (Doppler spread in Hz) is accepted and ignored here.

Two conversions have to agree with the relay, and both are done here the same way
it does them in to_baseband():

    delay      tau = 2R/c   with --distance-is range (default, R is one way)
               tau = L/c    with --distance-is path
    Doppler    f_d = 2v / lambda,  lambda = c / fc,  positive v = closing

A delay bin here is c/(N*k_step*SCS) = c/fs metres, the relay's sample grid. The relay
applies delays fractionally, so u_l = tau*fs keeps its sub-bin part. Only a relay run
with --integer-delays puts every path on an integer bin; pass --integer-delays here too.

What is left unknown is a global delay offset, because the UE's timing loop decides
where the FFT window sits and therefore which bin the LOS lands in. It also absorbs
the relay's common look-ahead of its delay filters. It is one number shared by every
path; --bin-offset searches for it.

Usage:
    ./model_check.py /tmp/map.csv.snap.csv --target 0:0:1.0 --target 30.5:8.78:0.3
    ./model_check.py /tmp/map.csv.snap.csv --target ... --fix-amp
    ./model_check.py /tmp/map.csv.snap.csv --paths 4          # no priors at all

The kernel A_L can also be looked at on its own, with no data involved:

    ./model_check.py --plot_kernel_A --L 1638
    ./model_check.py --plot_kernel_A --L 1638,819,204 --with-rect
    ./model_check.py /tmp/map.csv.snap.csv --plot_kernel_A --L 1638   # N and m/bin from the dump
"""

import argparse

import numpy as np

C0 = 299_792_458.0

# ---------------------------------------------------------------------------
# the window transform
# ---------------------------------------------------------------------------


def D_L(L, u, N):
    """Dirichlet kernel of a length-L rectangular window, zero-phase.

    D_L(u) = sin(pi L u / N) / sin(pi u / N), real and even, D_L(0) = L.
    Singular where u is a multiple of N; the limit there is L*(-1)^(m(L-1)).
    """
    u = np.asarray(u, dtype=float)
    s = np.sin(np.pi * u / N)
    out = np.empty(u.shape, dtype=float)

    sing = np.abs(s) < 1e-12
    ok = ~sing
    out[ok] = np.sin(np.pi * L * u[ok] / N) / s[ok]

    m = np.rint(u[sing] / N)
    out[sing] = L * np.where(np.mod(m * (L - 1), 2) == 0, 1.0, -1.0)
    return out


def A_L(L, u, N):
    """Zero-phase transform of the symmetric Hann taper of length L.

    A_L(u) = 0.5 D_L(u) + 0.25 D_L(u + dw) + 0.25 D_L(u - dw),  dw = N/(L-1)
    with A_L(0) = (L-1)/2 exactly.
    """
    dw = N / (L - 1.0)
    return 0.5 * D_L(L, u, N) + 0.25 * D_L(L, u + dw, N) + 0.25 * D_L(L, u - dw, N)


def kernel_stats(L, N, span=64.0, n=200_001):
    """Peak, -3 dB width, first sidelobe and -40 dB reach of A_L, in bins.

    Measured off a fine grid rather than derived in closed form, so it stays honest
    for small L where the three shifted Dirichlets are not well separated.
    """
    u = np.linspace(0.0, span, n)                 # even in u, so half the axis is enough
    a = np.abs(A_L(L, u, N))
    peak = a[0]

    # -3 dB: first crossing of the half-power level on the way down
    half = np.where(a < peak / np.sqrt(2.0))[0]
    w3 = 2.0 * u[half[0]] if len(half) else np.nan

    # first null, then the tallest lobe beyond it
    below = np.where(a < peak * 1e-3)[0]
    if len(below):
        first_null = u[below[0]]
        tail = a[below[0]:]
        sll = 20 * np.log10(tail.max() / peak + 1e-30)
    else:
        first_null, sll = np.nan, np.nan

    # last place the kernel is still above -40 dB: how far one path reaches
    loud = np.where(a > peak * 1e-2)[0]
    reach = u[loud[-1]] if len(loud) else np.nan

    return dict(peak=peak, w3=w3, first_null=first_null, sll_db=sll, reach40=reach)


# ---------------------------------------------------------------------------
# the dump
# ---------------------------------------------------------------------------

META = ["snap", "t_sample", "t_rel_s", "k_step", "k_first", "n_pilots",
        "idft_size", "n_bins", "scs_hz", "fs_hz", "carrier_hz"]

# nr_ue_sensing.h: slots in one TDD period, and a slot is 1 ms / 2^mu. Must track
# NR_SENSING_TDD_PERIOD_SLOTS: it sets where the Doppler spectrum repeats, 1/T_TDD.
TDD_PERIOD_SLOTS = 10


class Snapshots:
    """One nr_ue_sensing_dump_snapshots() file."""

    def __init__(self, path, max_bin=None, fc=None):
        raw = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
        n_meta = len(META)

        for i, name in enumerate(META):
            setattr(self, name, raw[:, i])

        self.N = int(self.idft_size[0])
        self.n_bins = int(self.n_bins[0])
        self.M = raw.shape[0]
        self.scs = float(self.scs_hz[0])
        self.fs = float(self.fs_hz[0])
        self.fc = float(fc) if fc is not None else float(self.carrier_hz[0])

        vals = raw[:, n_meta:n_meta + 2 * self.n_bins]
        self.h = vals[:, 0::2] + 1j * vals[:, 1::2]      # [M, n_bins]

        # reduces the range here in self.h
        if max_bin is not None and max_bin < self.n_bins:
            self.n_bins = max_bin
            self.h = self.h[:, :max_bin]

        self.b = np.arange(self.n_bins, dtype=float)
        self.t = self.t_rel_s

        # per snapshot: pilots, lattice start, centroid of the measured span
        self.n_m = self.n_pilots.astype(int)
        offset = np.mod(self.k_first, self.k_step)
        self.a_m = (self.k_first - offset) / self.k_step
        self.c_m = self.a_m + (self.n_m - 1) / 2.0

        # metres per bin, c / (N * k_step * SCS), which is also c / fs
        self.m_per_bin = C0 / (self.N * self.k_step[0] * self.scs)

        self.lambda_m = C0 / self.fc
        mu = np.log2(self.scs / 15000.0)
        self.t_tdd = TDD_PERIOD_SLOTS * 1e-3 / (2.0 ** mu)

    def doppler_hz(self, v_ms):
        """Radial speed to Doppler, the relay's convention: f_d = 2v/lambda."""
        return 2.0 * v_ms / self.lambda_m

    def speed_ms(self, f_d):
        return f_d * self.lambda_m / 2.0

    def describe(self):
        print(f"snapshots      : {self.M}")
        print(f"bins kept      : {self.n_bins}  ({self.n_bins * self.m_per_bin:.0f} m "
              f"at {self.m_per_bin:.2f} m/bin)")
        print(f"idft size N    : {self.N}   k_step {int(self.k_step[0])}")
        print(f"time span      : {self.t[-1] * 1e3:.1f} ms  "
              f"(Doppler resolution {self.lambda_m / (2 * self.t[-1]):.2f} m/s)")
        print(f"carrier        : {self.fc / 1e9:.5f} GHz")
        print(f"TDD period     : {self.t_tdd * 1e3:.2f} ms")
        lo, hi = self.n_m.min(), self.n_m.max()
        print(f"pilots n_m     : {lo} .. {hi}" + ("   (FIXED allocation)" if lo == hi
                                                  else "   (VARYING allocation)"))
        print(f"lattice a_m    : {self.a_m.min():.0f} .. {self.a_m.max():.0f}")


# ---------------------------------------------------------------------------
# ground truth
# ---------------------------------------------------------------------------


def parse_target(spec):
    """'DISTANCE_M:VELOCITY_MPS:AMPLITUDE[:SPREAD_HZ]' -> (d, v, a), as zmq_channel_relay.py takes it.

    AMPLITUDE is linear, or a level relative to the LOS with a dB suffix ('-40dB').
    """
    try:
        parts = spec.split(":")
        if len(parts) not in (3, 4):
            raise ValueError
        float(parts[3]) if len(parts) == 4 else None
        a = parts[2].strip()
        amp = 10.0 ** (float(a[:-2]) / 20.0) if a.lower().endswith("db") else float(a)
        return (float(parts[0]), float(parts[1]), amp)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"bad --target '{spec}', expected DISTANCE_M:VELOCITY_MPS:AMPLITUDE[:SPREAD_HZ], "
            f"e.g. 30.5:8.78:0.3 or 30.5:8.78:-40dB")


def targets_to_bins(snap, targets, distance_is_path, integer_delays=False):
    """Relay targets -> (bins, doppler Hz, amplitude), mirroring to_baseband().

    One relay sample is one bin here. The relay delays fractionally, so the bin keeps
    its sub-bin part, unless it was run with --integer-delays.
    """
    u, nu, amp = [], [], []
    for (dist, vel, a) in targets:
        tau = dist / C0 if distance_is_path else 2.0 * dist / C0
        u.append(float(round(tau * snap.fs)) if integer_delays else tau * snap.fs)
        nu.append(snap.doppler_hz(vel))
        amp.append(a)
    return np.array(u), np.array(nu), np.array(amp)


# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------


def find_paths(snap, n_paths, min_sep=3, exclude=()):
    """Path positions in bins, from the incoherent average profile.

    Bins near a position in exclude are masked, so a path already known from a
    prior is not found a second time.
    """
    if n_paths <= 0:
        return np.zeros(0)

    p = np.mean(np.abs(snap.h) ** 2, axis=0)
    taken = np.zeros(snap.n_bins, dtype=bool)
    for e in exclude:
        c = int(round(e))
        taken[max(0, c - min_sep):c + min_sep + 1] = True

    u0 = []
    for _ in range(n_paths):
        cand = np.where(~taken, p, -np.inf)
        b = int(np.argmax(cand))
        if not np.isfinite(cand[b]):
            break
        # parabolic refinement in dB, the sub-bin offset matters a lot to A_L
        if 0 < b < snap.n_bins - 1:
            y0, y1, y2 = np.log(p[b - 1:b + 2] + 1e-30)
            denom = y0 - 2 * y1 + y2
            delta = float(np.clip(0.5 * (y0 - y2) / denom, -0.5, 0.5)) if abs(denom) > 1e-12 else 0.0
        else:
            delta = 0.0
        u0.append(b + delta)
        taken[max(0, b - min_sep):b + min_sep + 1] = True

    return np.array(u0)


def design(snap, u, nu):
    """Model matrix G[(m,b), l] of equation (10), without the amplitudes."""
    M, B, L = snap.M, snap.n_bins, len(u)
    G = np.empty((M * B, L), dtype=complex)

    for l in range(L):
        du = np.broadcast_to(snap.b[None, :] - u[l], (M, B))
        K = np.empty((M, B), dtype=complex)
        # A_L depends on n_m, so it is built per distinct allocation
        for n_val in np.unique(snap.n_m):
            rows = snap.n_m == n_val
            K[rows] = (A_L(int(n_val), du[rows], snap.N)
                       * np.exp(1j * 2 * np.pi * snap.c_m[rows, None] * du[rows] / snap.N))
        K = K * np.exp(1j * 2 * np.pi * nu[l] * snap.t)[:, None]
        G[:, l] = K.ravel()

    return G


def fit(snap, u, nu, amp=None):
    """Least squares for the amplitudes; returns (alpha, model, nmse).

    With amp given, the relative amplitudes are held at the relay's values and only
    one global complex scale is free: the strictest test, since then every number in
    the model except that scale comes from the simulator. Without it each complex
    amplitude is free, which also absorbs path loss, the precoder and the AGC.
    """
    G = design(snap, u, nu)
    y = snap.h.ravel()

    if amp is None:
        alpha, *_ = np.linalg.lstsq(G, y, rcond=None)
    else:
        v = G @ amp
        alpha = amp * (np.vdot(v, y) / np.vdot(v, v))

    model = (G @ alpha).reshape(snap.M, snap.n_bins)
    nmse = np.sum(np.abs(snap.h - model) ** 2) / np.sum(np.abs(snap.h) ** 2)
    return alpha, model, nmse


def refine(snap, u, nu, amp=None, axis="u", span=0.5, steps=11, rounds=3):
    """Coordinate search on the path delays (axis 'u') or speeds (axis 'nu')."""
    x = (u if axis == "u" else nu).copy()
    best = fit(snap, u, nu, amp)[2]

    for _ in range(rounds):
        for l in range(len(x)):
            for g in x[l] + np.linspace(-span, span, steps):
                trial = x.copy()
                trial[l] = g
                e = fit(snap, trial if axis == "u" else u,
                        nu if axis == "u" else trial, amp)[2]
                if e < best:
                    best, x = e, trial
        span /= 4.0

    return x, best


def search_offset(snap, u, nu, amp=None, lo=-4, hi=32):
    """Global delay offset in bins, from the UE's FFT window placement.

    One number shared by every path, so it is searched rather than fitted: the
    timing loop decides which bin the LOS lands in, and the relay's delays are all
    measured from that.
    """
    best, best_off = np.inf, 0
    for off in range(lo, hi + 1):
        e = fit(snap, u + off, nu, amp)[2]
        if e < best:
            best, best_off = e, off
    return best_off, best


# ---------------------------------------------------------------------------
# the kernel on its own
# ---------------------------------------------------------------------------


def plot_kernel_A(Ls, N, span, out=None, with_rect=False, m_per_bin=None):
    """A_L alone, for a few L, with no data anywhere near it.

    L is the number of pilot subcarriers in the grant (n_m in the model), N the IDFT
    size. Only the ratio L/N sets the shape: the mainlobe is about N/L bins wide, so
    a smaller grant smears every path over more bins. The peak height (L-1)/2 is the
    only part that depends on L alone, and it is divided out here so the shapes can
    be compared.
    """
    import matplotlib
    if out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print(f"\nA_L with N = {N}\n")
    print(f"{'L':>7}  {'L/N':>7}  {'peak':>9}  {'-3dB bins':>10}  {'null bins':>10}  "
          f"{'sidelobe':>9}  {'-40dB bins':>11}")
    for L in Ls:
        s = kernel_stats(L, N, span=max(span, 64.0))
        print(f"{L:7d}  {L / N:7.3f}  {s['peak']:9.1f}  {s['w3']:10.3f}  "
              f"{s['first_null']:10.3f}  {s['sll_db']:8.1f}dB  {s['reach40']:11.3f}")

    u = np.linspace(-span, span, 4001)
    fig, ax = plt.subplots(figsize=(9, 5))
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, L in enumerate(Ls):
        col = colours[i % len(colours)]
        a = A_L(L, u, N)
        peak = A_L(L, np.zeros(1), N)[0]          # (L-1)/2

        # normalised magnitude in dB: the sidelobes, and how far they carry
        ax.plot(u, 20 * np.log10(np.abs(a) / peak + 1e-12), lw=1.3, color=col,
                label=f"L = {L}  (N/L = {N / L:.2f} bins)")

        if with_rect:
            d = D_L(L, u, N)
            ax.plot(u, 20 * np.log10(np.abs(d) / L + 1e-12), lw=0.8, ls=":",
                    color=col, alpha=0.6, label=f"L = {L}, rect (no taper)")

    ax.set_title(f"A_L, magnitude normalised to the peak (L-1)/2     N = {N}")
    ax.set_ylim(-80, 3)
    ax.set_xlabel("u = b - u_l [bins]")
    ax.set_ylabel("dB")
    ax.axvline(0.0, color="k", alpha=0.2, lw=0.8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    if m_per_bin:
        sec = ax.secondary_xaxis("top", functions=(lambda x: x * m_per_bin,
                                                   lambda x: x / m_per_bin))
        sec.set_xlabel("path length [m]")

    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=130)
        print(f"\nwrote {out}")
    else:
        plt.show()


# ---------------------------------------------------------------------------


def parse_L_list(values):
    """--L 1638 --L 400  or  --L 1638,819,400  -> [1638, 819, 400]."""
    out = []
    for v in values:
        for piece in str(v).split(","):
            piece = piece.strip()
            if piece:
                out.append(int(piece))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?",
                    help="file written by nr_ue_sensing_dump_snapshots(). Not needed with "
                         "--plot_kernel_A, which touches no data at all.")

    k = ap.add_argument_group("the kernel on its own (no data needed)")
    k.add_argument("--plot_kernel_A", action="store_true",
                   help="plot A_L by itself and stop, to see what the taper transform does "
                        "for a given grant size. Everything below tunes it.")
    k.add_argument("--L", action="append", default=[], metavar="L",
                   help="pilot count(s) to draw, repeatable and comma-separated, e.g. "
                        "--L 1638,819,204 (default 1638). This is n_m in the model.")
    k.add_argument("--kernel-N", type=int, default=None, metavar="N",
                   help="IDFT size (default 2048, or the one in the csv if given)")
    k.add_argument("--kernel-span", type=float, default=8.0, metavar="BINS",
                   help="half-width of the delay axis, in bins (default 8)")
    k.add_argument("--with-rect", action="store_true",
                   help="also draw the untapered Dirichlet D_L, to see what the Hann buys")

    g = ap.add_argument_group("ground truth (same syntax as zmq_channel_relay.py)")
    g.add_argument("--target", type=parse_target, action="append", default=[],
                   metavar="DISTANCE_M:VELOCITY_MPS:AMPLITUDE[:SPREAD_HZ]",
                   help="a simulated target, repeatable. Copy the relay's arguments "
                        "verbatim, including the LOS '0:0:1.0' or '0:0:0dB'.")
    g.add_argument("--integer-delays", action="store_true",
                   help="the relay ran with --integer-delays: round every target to a whole bin")
    g.add_argument("--fc", type=float, default=None, metavar="HZ",
                   help="carrier used to turn velocities into Doppler. Defaults to the "
                        "carrier_hz in the dump; override if the relay used another one.")
    g.add_argument("--distance-is", choices=("range", "path"), default="range",
                   help="how to read --target distances, exactly as the relay reads them "
                        "(default: range, i.e. one way, delay 2R/c)")
    g.add_argument("--fix-amp", action="store_true",
                   help="hold the relative amplitudes at the relay's values, leaving one "
                        "global complex scale free")

    f = ap.add_argument_group("fitting")
    f.add_argument("--paths", type=int, default=0,
                   help="extra paths found from the data, taken as static. Used on its own "
                        "when there are no priors; default 0, or 3 if no --target is given.")
    f.add_argument("--max-bin", type=int, default=64,
                   help="bins to fit over, from bin 0 (default 64)")
    f.add_argument("--bin-offset", default="auto",
                   help="global delay offset in bins from the UE's FFT window placement. "
                        "'auto' searches for it (default), or give an integer, or 'none'.")
    f.add_argument("--refine", action="store_true",
                   help="coordinate search on the path delays")
    f.add_argument("--refine-doppler", action="store_true",
                   help="coordinate search on the path speeds, to see where the data puts "
                        "them relative to the priors")
    f.add_argument("--flip", action="store_true",
                   help="mirror the delay axis, if the IDFT sign convention is the other one")

    p = ap.add_argument_group("output")
    p.add_argument("--plot", metavar="PNG", help="write a figure instead of showing it")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--kernel", choices=("scaled", "bare", "off"), default="scaled",
                   help="overlay the taper transform A_L of each path on the profiles, to "
                        "see how far it reaches. 'scaled' weights it by the fitted |alpha_l|, "
                        "so the curves add up to the model (default); 'bare' draws the shape "
                        "alone, its peak put on the measured curve; 'off' drops it.")
    args = ap.parse_args()

    # ---- the kernel on its own: no csv, no fit, nothing else runs
    if args.plot_kernel_A:
        Ls = parse_L_list(args.L) or [1638]
        if any(L < 2 for L in Ls):
            ap.error("--L must be at least 2")
        N, m_per_bin = args.kernel_N, None
        if args.csv:                              # borrow the geometry from a real dump
            s = Snapshots(args.csv, max_bin=args.max_bin, fc=args.fc)
            N = N or s.N
            m_per_bin = s.m_per_bin
        plot_kernel_A(Ls, N or 2048, args.kernel_span, out=args.plot,
                      with_rect=args.with_rect, m_per_bin=m_per_bin)
        return

    if not args.csv:
        ap.error("give a csv, or --plot_kernel_A to look at the kernel on its own")

    snap = Snapshots(args.csv, max_bin=args.max_bin, fc=args.fc)
    if args.flip:
        snap.h = np.conj(snap.h)
    snap.describe()

    # ---- paths: priors first, then whatever else the data shows, taken as static
    u, nu, amp = targets_to_bins(snap, args.target, args.distance_is == "path", args.integer_delays)
    n_known = len(u)

    n_extra = args.paths if (args.paths or args.target) else 3
    extra = find_paths(snap, n_extra, exclude=u)
    if len(extra):
        u = np.concatenate([u, extra])
        nu = np.concatenate([nu, np.zeros(len(extra))])
        amp = np.concatenate([amp, np.full(len(extra), np.nan)])

    if len(u) == 0:
        ap.error("no paths: give --target or --paths")

    if args.fix_amp and not np.all(np.isfinite(amp)):
        ap.error("--fix-amp needs every path to come from --target (drop --paths)")
    amp_fixed = amp.astype(complex) if args.fix_amp else None

    # ---- global delay offset from the UE's FFT window
    if args.bin_offset == "auto" and n_known:
        off, _ = search_offset(snap, u, nu, amp_fixed)
        print(f"\nbin offset     : {off} bins ({off * snap.m_per_bin:.2f} m), searched")
    elif args.bin_offset in ("auto", "none"):
        off = 0
    else:
        off = int(args.bin_offset)
        print(f"\nbin offset     : {off} bins ({off * snap.m_per_bin:.2f} m), given")
    u = u + off

    # ---- report the scene the model will be built from
    print(f"\n{'path':>4}  {'bin':>7}  {'path m':>8}  {'speed m/s':>10}  {'f_d Hz':>8}  "
          f"{'amp':>6}   source")
    for l in range(len(u)):
        v = snap.speed_ms(nu[l])
        a_s = f"{amp[l]:6.3f}" if np.isfinite(amp[l]) else "     -"
        print(f"{l:4d}  {u[l]:7.2f}  {u[l] * snap.m_per_bin:8.2f}  {v:10.2f}  {nu[l]:8.2f}  "
              f"{a_s}   {'relay' if l < n_known else 'from data'}")

    alpha, model, nmse = fit(snap, u, nu, amp_fixed)
    print(f"\nNMSE           : {nmse:.5f}   ({10 * np.log10(nmse):.1f} dB)"
          + ("   [amplitudes held at the relay's values]" if args.fix_amp else ""))

    if args.refine:
        u, _ = refine(snap, u, nu, amp_fixed, axis="u")
        alpha, model, nmse = fit(snap, u, nu, amp_fixed)
        print(f"refined bins   : {np.round(u, 3)}   NMSE {nmse:.5f} "
              f"({10 * np.log10(nmse):.1f} dB)")

    if args.refine_doppler:
        nu, _ = refine(snap, u, nu, amp_fixed, axis="nu",
                       span=snap.doppler_hz(1.0), steps=21)
        alpha, model, nmse = fit(snap, u, nu, amp_fixed)
        print(f"refined speeds : {np.round([snap.speed_ms(x) for x in nu], 3)} m/s   "
              f"NMSE {nmse:.5f} ({10 * np.log10(nmse):.1f} dB)")

    # ---- per-snapshot residual, the interesting one when the allocation varies
    per_snap = (np.sum(np.abs(snap.h - model) ** 2, axis=1)
                / np.sum(np.abs(snap.h) ** 2, axis=1))
    print(f"\nper snapshot   : best {per_snap.min():.5f}, worst {per_snap.max():.5f}")

    if args.no_plot:
        return

    import matplotlib
    if args.plot:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 1, figsize=(9, 10))

    def db(x):
        return 20 * np.log10(np.abs(x) + 1e-9)

    worst = int(np.argmax(per_snap))
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for k, (idx, name) in enumerate([(0, "first snapshot"),
                                     (worst, f"worst snapshot ({worst})")]):
        x = snap.b * snap.m_per_bin
        ax[k].plot(x, db(snap.h[idx]), label="measured", lw=1.4)
        ax[k].plot(x, db(model[idx]), "--", label="model", lw=1.2)

        # each path's taper transform on its own: the phase ramp and the Doppler
        # term are unit modulus, so the magnitude of a path's contribution is
        # exactly |alpha_l| * |A_L(b - u_l)|. Drawing it shows how far the
        # sidelobes of one path reach into the bins of the others.
        if args.kernel != "off":
            for l in range(len(u)):
                a_l = A_L(int(snap.n_m[idx]), snap.b - u[l], snap.N)
                if args.kernel == "scaled":
                    curve = db(alpha[l] * a_l)
                else:
                    peak = db(snap.h[idx][int(np.clip(round(u[l]), 0, snap.n_bins - 1))])
                    curve = db(a_l / A_L(int(snap.n_m[idx]), np.zeros(1), snap.N)[0]) + peak
                ax[k].plot(x, curve, lw=2.5, alpha=0.85, ls="-",
                           color=colours[(l + 2) % len(colours)],   # 0,1 are the two curves above
                           label=f"A_L path {l}" + ("" if l < n_known else " (data)"))

        for l in range(n_known):
            ax[k].axvline(u[l] * snap.m_per_bin, color="k", alpha=0.25, lw=0.8)
        ax[k].set_title(f"{name}   n_m = {snap.n_m[idx]}   NMSE {per_snap[idx]:.5f}")
        ax[k].set_xlabel("path length [m]")
        ax[k].set_ylabel("dB")
        # the kernel nulls dive far below anything measured; keep the useful range
        floor = db(snap.h[idx]).max() - 80.0
        ax[k].set_ylim(bottom=floor)
        ax[k].legend(fontsize=8, ncol=2)
        ax[k].grid(alpha=0.3)

    fig.tight_layout()
    if args.plot:
        fig.savefig(args.plot, dpi=130)
        print(f"\nwrote {args.plot}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
