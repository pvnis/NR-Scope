#!/usr/bin/env python3
# SPDX-License-Identifier: LicenseRef-CSSL-1.0
"""Plot range-Doppler maps dumped by the UE sensing pipeline.

The maps are produced by nr_ue_sensing_dump_map(), which appends one CSV line per
map to the file named by --sensing-dump:

    frame,slot,aarx,ports,layer,k_step,k_offset,n_snapshots,t_span_s,m_per_bin,
    f_max_hz,carrier_hz,n_bins,n_freq,n_marker,n_aoa,bin_los,power...,marker...,aoa...

The payload is row-major with the delay bin outer and the Doppler bin inner, i.e.
power[b * n_freq + f]. Everything needed to scale both axes is in the line itself:

    path length of bin b  = b * m_per_bin                        metres
    Doppler of bin f      = -f_max_hz + f * 2*f_max_hz/(n_freq-1) Hz
    speed                 = doppler * lambda / 2, lambda = c / carrier_hz

bin_los is where the direct path sits on that axis, to sub-bin resolution. The delay
axis is referenced to the UE's FFT window rather than to the gNB's transmit instant,
so the direct path lands wherever the sync loop parked the window and not at bin 0.
The map keeps the absolute axis; the AoA block below reports excess range over the
direct path, so bin_los * m_per_bin is what converts between the two.

Markers, when present, are what --sensing-tdd-detect decided, four numbers each:
range_m, speed_ms, snr_dB, verdict. Verdict 0 is a target the detector kept, 1 a peak
it threw out as a replica of another one, 2 residue of a peak already accepted. They
are drawn on the map because a plot of the survivors alone looks the same whether the
detector ran or not: the rejections are the visible part of its work.

The speed axis reaches --sensing-max-speed plus one TDD period, so a car sits on the
map at its true speed. Because of TDD the spectrum repeats every 1/T_TDD (200 Hz, or
8.8 m/s at 3.41 GHz), so every target also appears as evenly spaced copies. Those are
what the detector rejects as replicas; without --sensing-tdd-detect they are all drawn
as ordinary cells.

The AoA block, when present, is what --sensing-antenna-avg measured across the Rx
array at the target cells, five numbers each: range_m, speed_ms,
angle_deg, power_dB (SNR over the CFAR noise), quality_dB. range_m here is excess over the
direct path, not absolute path length, which is what bistatic localisation needs; it
is drawn at range_m + bin_los * m_per_bin to land on the map's own axis. The angle is
0 at the array normal. quality_dB is the spatial peak over the mean of the spatial spectrum, which
is 10*log10(Mr) for a single plane wave, so 6 dB on four antennas; cells well under
that held more than one arrival and their angle is drawn faded. With
--sensing-tdd-detect the AoA cells are the detector's targets (true speed, can be
outside the map's speed axis); without it the AoA finds its own peaks on the averaged
map with CFAR.

--carrier overrides carrier_hz, for dumps written before that column existed.

One line is one map for one Rx antenna and one measurement stream, identified by
(ports, layer, k_step, k_offset). ports is the DM-RS port bitmap of the allocation
and layer picks one of its columns, so at rank 2 the same slot produces two maps,
one per layer, alike in every other field. They are different beams, H w_0 and
H w_1, and are never combined here. A rank change also changes ports and therefore
starts a separate stream, because layer 0 of a rank 2 grant is not the same beam as
the single layer of a rank 1 grant.

Collect data:
    ./nr-uesoftmodem ... --sensing-symbols 384 --sensing-dump /tmp/map.csv

Plot it:
    ./plot_range_doppler.py /tmp/map.csv
    ./plot_range_doppler.py /tmp/map.csv --index -1 --max-range 150
    ./plot_range_doppler.py /tmp/map.csv --layer 1 --x speed -o layer1.png
    ./plot_range_doppler.py /tmp/map.csv --list

Browse every map in one window, nothing written to disk:
    ./plot_range_doppler.py /tmp/map.csv --browse
    ./plot_range_doppler.py /tmp/map.csv --browse --index 40 --trail 8 --max-range 50

    <- / ->      previous / next map
    up / down    10 maps forward / back
    Home / End   first / last map
    q            quit

Play every map back like a video, nothing written to disk:
    ./plot_range_doppler.py /tmp/map.csv --play
    ./plot_range_doppler.py /tmp/map.csv --play --fps 15 --trail 8 --max-range 50

    space        play / pause
    <- / ->      step one map while paused
    -  /  +      slower / faster
    drag slider  scrub to any map
    q            quit

Playback loops back to the first map after the last. Unlike --browse it redraws on fixed
axes (the colour bar and layout are built once), so it stays smooth frame to frame.

With --trail K the targets the detector accepted in the K maps before the current one are
drawn faded on top of it. Without localisation this is the quickest way to tell a person
from clutter residue: a real target leaves a trail that moves smoothly in range and speed
from one map to the next, residue appears somewhere else every time.

--list marks with "!!" every map on which the detector rejected at least one peak as
a TDD replica, and says so again, louder, when the rejected replica was the strongest
cell of the map: that is the case where reading the map by its argmax reports an
artefact as a target.
"""

import argparse
import sys

import numpy as np

C_M_PER_S = 299792458.0

# Header fields ahead of the payload, in dump order. Every layout that has existed
# still parses: the current one and, going back, the same without the noise reference,
# without the grant statistics, without the LoS bin, without the stream identity,
# without the AoA block, without detection markers, and an early one without
# carrier_hz. They are told apart by whether the payload length works out, so an old
# file still plots, with --carrier supplying the wavelength where the dump has none.
_META = ("frame", "slot", "aarx", "ports", "layer", "k_step", "k_offset",
         "n_snapshots", "t_span_s", "m_per_bin", "f_max_hz", "carrier_hz",
         "n_bins", "n_freq", "n_marker", "n_aoa", "bin_los",
         "n_pilots_min", "n_pilots_max", "n_positions", "noise_ref")
_META_INT = {"frame", "slot", "aarx", "ports", "layer", "k_step", "k_offset",
             "n_snapshots", "n_bins", "n_freq", "n_marker", "n_aoa",
             "n_pilots_min", "n_pilots_max", "n_positions"}

# Layouts predating the noise reference, when power was the transform's own numbers
# rather than multiples of the map's noise floor. Levels are then not comparable
# between maps of such a file, which is what noise_ref was added to fix.
_META_NO_NOISE = tuple(n for n in _META if n != "noise_ref")

# Layouts predating the grant statistics of the window.
_META_NO_GRANT = tuple(n for n in _META_NO_NOISE if n not in ("n_pilots_min", "n_pilots_max", "n_positions"))

# Layouts predating the LoS bin, when the AoA range axis was still absolute and the
# range a cell reported was its own path length rather than its excess over the
# direct path.
_META_NO_LOS = tuple(n for n in _META_NO_GRANT if n != "bin_los")

# Layouts predating the stream identity, when a map was told apart by k_step alone
# and the layers of a rank 2 grant would have been indistinguishable.
_META_NO_STREAM = tuple(n for n in _META_NO_LOS if n not in ("ports", "layer", "k_offset"))
_META_NO_AOA = tuple(n for n in _META_NO_STREAM if n != "n_aoa")
_META_NO_MARKERS = tuple(n for n in _META_NO_AOA if n != "n_marker")
_META_LEGACY = tuple(n for n in _META_NO_MARKERS if n != "carrier_hz")

# verdict codes, mirroring nr_tdd_verdict_t. UNCERTAIN (a target with an equally strong
# mirror at -v) is listed by --list but deliberately not drawn on the map.
VERDICT_TARGET, VERDICT_REPLICA, VERDICT_DUPLICATE, VERDICT_UNCERTAIN = 0, 1, 2, 3

# A single plane wave gives a spatial peak-to-mean of Mr, so 6 dB on four antennas.
# Below this the cell holds more than one arrival, or nothing but noise.
AOA_QUALITY_MIN_DB = 4.0


class Map:
    """One range-Doppler map and the parameters needed to scale its axes."""

    def __init__(self, fields):
        for meta in (_META, _META_NO_NOISE, _META_NO_GRANT, _META_NO_LOS, _META_NO_STREAM,
                     _META_NO_AOA, _META_NO_MARKERS, _META_LEGACY):
            self._read_meta(meta, fields)
            if self.payload.size == (self.n_bins * self.n_freq
                                     + 4 * self.n_marker + 5 * self.n_aoa):
                break
        else:
            raise ValueError("payload does not match n_bins * n_freq under any layout")

        n_power = self.n_bins * self.n_freq
        n_mark = 4 * self.n_marker
        # power[b * n_freq + f] -> [bin, doppler]
        self.power = self.payload[:n_power].reshape(self.n_bins, self.n_freq)
        # markers[i] = range_m, speed_ms, snr_dB, verdict
        self.markers = self.payload[n_power:n_power + n_mark].reshape(self.n_marker, 4)
        # aoa[i] = range_m, speed_ms, angle_deg, power_dB, quality_dB
        self.aoa = self.payload[n_power + n_mark:].reshape(self.n_aoa, 5)

    def _read_meta(self, meta, fields):
        self.carrier_hz = None
        self.n_marker = 0
        self.n_aoa = 0
        # Absent from dumps written before the direct path was located. Zero is what
        # those dumps assumed: the AoA range axis was the absolute one, so leaving it
        # at zero makes the overlay below land where it always did.
        self.bin_los = 0.0
        # absent from dumps written before the grant statistics were recorded
        self.n_pilots_min = None
        self.n_pilots_max = None
        self.n_positions = None
        # absent from dumps written before power was scaled to the noise floor. None
        # rather than 1.0 so a reader can tell "this map was never scaled" from "this
        # map was scaled by one", which a degenerate map reports as 0.0.
        self.noise_ref = None
        # absent from dumps written before the stream identity was recorded
        self.ports = None
        self.layer = None
        self.k_offset = None
        for name, raw in zip(meta, fields):
            # int(float(...)) rather than int(...): when a layout guess is wrong an
            # integer field lands on a power value in scientific notation, and that
            # has to fall through to the payload length check rather than raise here.
            setattr(self, name, int(float(raw)) if name in _META_INT else float(raw))
        self.payload = np.asarray(fields[len(meta):], dtype=np.float64)

    def markers_of(self, verdict):
        """(speeds, ranges, snr_dB) of every marker carrying this verdict."""
        if self.n_marker == 0:
            return np.empty(0), np.empty(0), np.empty(0)
        keep = self.markers[:, 3].astype(int) == verdict
        m = self.markers[keep]
        return m[:, 1], m[:, 0], m[:, 2]

    def speed_step(self, override=None):
        """Spacing of the speed grid, in m/s."""
        v = self.speeds(override)
        return float(abs(v[1] - v[0])) if v.size > 1 else 0.0

    @property
    def n_replica_rejected(self):
        return int(self.markers_of(VERDICT_REPLICA)[0].size)

    def detection_lines(self, override=None):
        """What --sensing-tdd-detect decided, as lines for --list.

        Empty when the dump carries no markers, which is how a run without the
        detector is told apart from one where it ran and rejected nothing."""
        if not self.n_marker:
            return []
        tv, tr, tsnr = self.markers_of(VERDICT_TARGET)
        rv, rr, _ = self.markers_of(VERDICT_REPLICA)
        dv = self.markers_of(VERDICT_DUPLICATE)[0]

        parts = [f"{tv.size} target(s)"]
        parts.append(f"{rv.size} TDD replica(s) rejected" if rv.size else "no replicas rejected")
        if dv.size:
            parts.append(f"{dv.size} residue")
        uv, ur, usnr = self.markers_of(VERDICT_UNCERTAIN)
        if uv.size:
            parts.append(f"{uv.size} uncertain (mirrored at -v, not drawn)")
        lines = ["tdd-detect: " + ", ".join(parts)]
        for v, r, snr in zip(tv, tr, tsnr):
            lines.append(f"target {r:8.1f} m {v:+7.2f} m/s {snr:5.1f} dB")
        for v, r, snr in zip(uv, ur, usnr):
            lines.append(f"uncert {r:8.1f} m {v:+7.2f} m/s {snr:5.1f} dB")

        # The case worth seeing: the strongest cell of the map is something the
        # detector threw out, so reading the map by argmax alone reports a replica as
        # a target. One bin in range and two Doppler cells is the same tolerance the
        # detector uses to call a candidate the same object.
        pr, pv = self.peak(override)
        tol_v = max(2.0 * self.speed_step(override), 1e-9)
        for v, r in zip(rv, rr):
            if abs(r - pr) <= self.m_per_bin and abs(v - pv) <= tol_v:
                lines.append(f"!! strongest cell is a rejected replica: argmax alone "
                             f"would report {pr:.1f} m / {pv:+.2f} m/s")
                break
        return lines

    def aoa_lines(self):
        """What the AoA stage measured, as lines for --list.

        Empty when the dump carries no angles, which is how a run without
        --sensing-antenna-avg is told apart from one where every cell was too weak."""
        if not self.n_aoa:
            return []
        lines = [f"aoa: {self.n_aoa} cell(s), direct path at bin {self.bin_los:.2f} "
                 f"({self.los_offset_m:.1f} m)"]
        for r, v, ang, pwr, q in self.aoa:
            # A low quality figure is the honest case to show rather than hide: it
            # means the four antennas did not agree on one plane wave, so the angle
            # next to it is not one direction but a blend of several.
            warn = "" if q >= AOA_QUALITY_MIN_DB else "  (!) low quality"
            # r is excess over the direct path; the absolute path length that puts
            # it on the map axis follows in brackets
            lines.append(f"cell  +{r:7.1f} m ({r + self.los_offset_m:7.1f} abs) "
                         f"{v:+7.2f} m/s {ang:+7.1f} deg "
                         f"SNR {pwr:5.1f} dB, q {q:4.1f} dB{warn}")
        return lines

    @property
    def signal(self):
        """Short name of the measurement stream, for plot titles and --list."""
        if self.ports is None:
            return f"k_step {self.k_step}"
        rank = bin(self.ports).count("1")
        
        # 0xff marks a map summed over every layer, mirroring aarx = -1 for antennas
        avg_lay = "avg_lay" if self.layer == 0xFF else "Layer" + str(self.layer)

        return f"DMRS port {self.ports:x}, {avg_lay}/{rank}, k-step {self.k_step}"

    @property
    def ranges(self):
        """Path length of each delay bin, in metres."""
        return np.arange(self.n_bins) * self.m_per_bin

    @property
    def los_offset_m(self):
        """Where the direct path sits on the map's own range axis, in metres.

        The delay axis is referenced to the UE's FFT window, not to the gNB's
        transmit instant, so the direct path lands wherever the sync loop parked
        the window. That offset is what the AoA stage subtracts to report an
        excess range, and what has to be added back to plot a cell on this axis."""
        return self.bin_los * self.m_per_bin

    @property
    def dopplers(self):
        """Doppler of each frequency bin, in Hz. Matches the C grid exactly."""
        if self.n_freq == 1:
            return np.zeros(1)
        return np.linspace(-self.f_max_hz, self.f_max_hz, self.n_freq)

    def carrier(self, override=None):
        """Carrier in Hz: the override if given, else the dumped value."""
        if override:
            return override
        if not self.carrier_hz:
            sys.exit("this dump carries no carrier_hz, pass --carrier to scale the speed axis")
        return self.carrier_hz

    def speeds(self, override=None):
        """Doppler grid converted to radial speed, in m/s."""
        return self.dopplers * (C_M_PER_S / self.carrier(override)) / 2.0

    @property
    def doppler_resolution_hz(self):
        return 1.0 / self.t_span_s if self.t_span_s > 0 else float("nan")

    """
    Normalization:
    
    In simulation, the scene is deterministic, so the peak power is stable and can be used as
    reference for normalization. In real scene, the peak is usually clutter residue, so making it
    the reference does not make sense since we want to remove it by clutter removal. Also
    its better to make a detection algorithm for "something above the reference of that much...", 
    and not the opposite saying we detect a target is we are "below the reference of that much...".
    """

    def power_dB(self, norm_value):
        # """Normalised to the strongest cell of this map."""
        # peak = self.power.max()
        # if peak <= 0.0:
        #     return np.full_like(self.power, -np.inf)
        # with np.errstate(divide="ignore"):
        #     return 10.0 * np.log10(self.power / peak)
        
        norm_factor = self.power.max() if norm_value == "peak" else np.median(self.power[self.power > 0])
        return 10 * np.log10(self.power / norm_factor)

    def peak(self, override=None):
        """(range_m, speed_m_s) of the strongest cell."""
        b, f = np.unravel_index(np.argmax(self.power), self.power.shape)
        return self.ranges[b], self.speeds(override)[f]

    def label(self, override=None):
        # aarx -1 marks a map averaged over every Rx antenna, see --sensing-antenna-avg,
        # and -2 one built from chain 0 minus weighted chain 1 (spatial_null)
        rx = "rx null" if self.aarx == -2 else "rx avg" if self.aarx < 0 else f"rx{self.aarx}"
        return (f"{rx} {self.signal}: "
                f"{self.n_snapshots} snapshots over {self.t_span_s * 1e3:.1f} ms, "
                f"{self.n_bins} bins x {self.n_freq} doppler, "
                f"{self.m_per_bin:.2f} m/bin, f_max +-{self.f_max_hz:.0f} Hz, "
                f"{self.carrier(override) / 1e9:.2f} GHz")


def load(path):
    maps = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                maps.append(Map(line.split(",")))
            except ValueError as e:
                print(f"{path}:{lineno}: skipped ({e})", file=sys.stderr)
    if not maps:
        sys.exit(f"{path}: no usable maps found")
    return maps


def draw_map(fig, ax, m, args, trail=(), cax=None):
    """Draw one map on ax, with its markers and AoA cells, and its colour bar on fig.

    trail: maps shown before this one, whose accepted targets are overlaid faded, oldest
    faintest (see --trail).
    cax: reuse this axes for the colour bar instead of carving a new one out of ax every
    call. --play clears and redraws it each frame, so one fixed axes keeps the layout
    steady and the redraw cheap; without it (the one-shot and --browse paths) a new colour
    bar is made next to ax as before."""
    vmax = 0.0 if args.norm == "peak" else 40.0

    x = m.speeds(args.carrier) if args.x == "speed" else m.dopplers
    # Markers and AoA cells are stored as speeds (m/s). In --x doppler mode they
    # have to be turned back into Hz, or they land at the wrong place on the axis.
    hz_per_ms = 2.0 / (C_M_PER_S / m.carrier(args.carrier))
    to_x = (lambda v: np.asarray(v)) if args.x == "speed" else (lambda v: np.asarray(v) * hz_per_ms)
    y = m.ranges
    img = m.power_dB(args.norm)

    keep = slice(None)
    if args.max_range is not None:
        n = int(np.searchsorted(y, args.max_range)) + 1
        keep = slice(0, min(n, m.n_bins))

    mesh = ax.pcolormesh(x, y[keep], img[keep, :], vmin=args.floor, vmax=vmax,
                         shading="nearest", cmap="viridis")
    if cax is not None:
        fig.colorbar(mesh, cax=cax, label="dB below peak")
    else:
        fig.colorbar(mesh, ax=ax, label="dB below peak")

    # zero Doppler is where every static return should sit
    ax.axvline(0.0, color="white", lw=0.6, ls=":", alpha=0.7)

    r, v = m.peak(args.carrier)
    ax.plot([to_x(v)], [r], marker="o", mfc="none", mec="white", ms=13, mew=1.0,
            alpha=0.8, label=f"strongest cell {r:.1f} m / {v:+.2f} m/s")

    if not args.no_markers and trail:
        # Targets of the previous maps, oldest faintest. Drawn first, so the current
        # map's own verdicts stay on top.
        labelled = False
        for age, old in enumerate(reversed(list(trail)), 1):
            tv, tr, _ = old.markers_of(VERDICT_TARGET)
            if tv.size:
                ax.plot(to_x(tv), tr, linestyle="none", marker="s", mfc="lime", mec="none",
                        ms=6, alpha=max(0.15, 0.6 * (1.0 - (age - 1) / len(trail))),
                        label=None if labelled else f"targets of the {len(trail)} previous map(s)")
                labelled = True

    if not args.no_markers and m.n_marker:
        # Rejections first, so an accepted target drawn at the same place stays
        # on top: the interesting case is a target sitting under its own replica.
        vx, vy, _ = m.markers_of(VERDICT_DUPLICATE)
        if vx.size:
            ax.plot(to_x(vx), vy, linestyle="none", marker="+", mec="orange", ms=9, mew=1.4,
                    label=f"residue of an accepted peak ({vx.size})")
        vx, vy, _ = m.markers_of(VERDICT_REPLICA)
        if vx.size:
            ax.plot(to_x(vx), vy, linestyle="none", marker="x", mec="red", ms=9, mew=1.8,
                    label=f"rejected: TDD replica ({vx.size})")
        vx, vy, snr = m.markers_of(VERDICT_TARGET)
        if vx.size:
            ax.plot(to_x(vx), vy, linestyle="none", marker="s", mfc="none", mec="lime",
                    ms=13, mew=2.0, label=f"target ({vx.size})")
            for x0, y0, s0 in zip(vx, vy, snr):
                ax.annotate(f"{y0:.1f} m\n{x0:+.2f} m/s\n{s0:.0f} dB",
                            (to_x(x0), y0), textcoords="offset points", xytext=(9, 9),
                            fontsize=7, color="lime")

    if not args.no_aoa and m.n_aoa:
        # Drawn below the detector's markers, so the detector's verdict stays
        # legible on top. With --sensing-tdd-detect these cells are the detector's
        # targets and can sit outside the map's speed axis (true speed of a fast
        # target); matplotlib widens the axis to show them.
        # aoa range is excess over the direct path, the map axis is absolute
        # path length, so the offset goes back on before they can be drawn
        # together. Zero on dumps predating bin_los, where both were absolute.
        ar, av = m.aoa[:, 0] + m.los_offset_m, to_x(m.aoa[:, 1])
        good = m.aoa[:, 4] >= AOA_QUALITY_MIN_DB
        for sel, tag, style in (
                (good, "", dict(mec="deepskyblue", mew=1.6)),
                (~good, ", low quality", dict(mec="steelblue", mew=0.9, alpha=0.5))):
            if sel.any():
                ax.plot(av[sel], ar[sel], linestyle="none", marker="D", mfc="none",
                        ms=8, label=f"aoa cell ({int(sel.sum())}){tag}", **style)
        for (r0, v0, ang0, _, q0), y0 in zip(m.aoa, ar):
            col = "deepskyblue" if q0 >= AOA_QUALITY_MIN_DB else "steelblue"
            ax.annotate(f"{ang0:+.0f}°", (to_x(v0), y0),
                        textcoords="offset points", xytext=(-26, -4),
                        fontsize=7, color=col)

    ax.legend(loc="upper right", fontsize=7, framealpha=0.7)

    ax.set_xlabel("speed [m/s]" if args.x == "speed" else "doppler [Hz]")
    ax.set_ylabel("path length [m]")
    ax.set_title(m.label(args.carrier), fontsize=8)


def _info_lines(m, i, n, args):
    """The lines shown under an interactive map: which map it is, its grants, and the
    detector's and AoA verdicts. Shared by --browse and --play."""
    lines = [f"map {i + 1}/{n}   frame.slot {m.frame}.{m.slot}"]
    if m.n_pilots_min is not None:
        lines.append(f"grants: {m.n_pilots_min} to {m.n_pilots_max} pilots, "
                     f"{m.n_positions} distinct position(s) a_m")
    lines += m.detection_lines(args.carrier) + m.aoa_lines()
    return lines


# keys the viewer uses, and the step each one moves by
_BROWSE_STEPS = {"right": 1, "left": -1, "up": 10, "down": -10, "pageup": 10, "pagedown": -10}


def browse(maps, args, start=0):
    """Show the maps one at a time in a single window and step through them with the keys.

    Nothing is written to disk: each key press redraws the figure with another map."""
    import matplotlib.pyplot as plt

    # These keys belong to the viewer here; by default matplotlib's toolbar also uses
    # left/right for its zoom history and Home to reset the view.
    for key, taken in (("keymap.back", ("left",)), ("keymap.forward", ("right",)),
                       ("keymap.home", ("home",))):
        plt.rcParams[key] = [k for k in plt.rcParams[key] if k not in taken]

    fig = plt.figure(figsize=(10, 7))
    state = {"i": start % len(maps)}

    def show():
        i = state["i"]
        m = maps[i]
        fig.clf()
        ax = fig.add_subplot(111)
        trail = maps[max(0, i - args.trail):i] if args.trail > 0 else []
        draw_map(fig, ax, m, args, trail)

        info = _info_lines(m, i, len(maps), args)
        fig.text(0.01, 0.01, "\n".join(info), fontsize=7, family="monospace", va="bottom")
        fig.suptitle("<- / -> previous / next map    up / down 10 maps    Home / End first / last    q quit",
                     fontsize=8)
        bottom = min(0.45, 0.03 + 0.017 * len(info))
        fig.tight_layout(rect=(0, bottom, 1, 0.97))
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key in _BROWSE_STEPS:
            state["i"] = min(max(state["i"] + _BROWSE_STEPS[event.key], 0), len(maps) - 1)
        elif event.key == "home":
            state["i"] = 0
        elif event.key == "end":
            state["i"] = len(maps) - 1
        else:
            return
        show()

    fig.canvas.mpl_connect("key_press_event", on_key)
    show()
    plt.show()
    return fig, on_key, state


# how far --play's -/+ keys can push the frame rate
_PLAY_FPS_MIN, _PLAY_FPS_MAX = 0.5, 60.0


def play(maps, args, start=0, fps=5.0):
    """Play the maps back like a video in one window, nothing written to disk.

    Everything that does not change between maps (figure, axes, colour-bar axes, slider)
    is built once; each frame only clears the plot and the colour bar and redraws them, so
    playback stays smooth where --browse, which rebuilds the whole figure per key press,
    does not. A timer advances the frame, space pauses, the arrows step while paused, -/+
    change the rate, and the slider scrubs."""
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    # The arrows step frames here; by default matplotlib's toolbar steals left/right for
    # its zoom history.
    for key, taken in (("keymap.back", ("left",)), ("keymap.forward", ("right",))):
        plt.rcParams[key] = [k for k in plt.rcParams[key] if k not in taken]

    fig = plt.figure(figsize=(10, 7.5))
    # Fixed geometry: the plot, its colour bar, and the scrubber never move, so no
    # per-frame tight_layout and no drift. The bottom strip is left for the info text.
    ax = fig.add_axes((0.08, 0.34, 0.80, 0.56))
    cax = fig.add_axes((0.90, 0.34, 0.015, 0.56))
    sax = fig.add_axes((0.08, 0.20, 0.80, 0.025))
    info = fig.text(0.01, 0.01, "", fontsize=7, family="monospace", va="bottom")

    state = {"i": start % len(maps), "playing": True, "fps": float(fps)}

    slider = Slider(sax, "map", 1, len(maps), valinit=state["i"] + 1, valstep=1)

    def title():
        rate = f"{state['fps']:g} fps"
        mode = "playing" if state["playing"] else "paused"
        fig.suptitle(f"space play / pause ({mode})    <- / -> step    - / + speed ({rate})"
                     f"    drag to scrub    q quit", fontsize=8)

    def render():
        i = state["i"]
        m = maps[i]
        ax.cla()
        cax.cla()
        trail = maps[max(0, i - args.trail):i] if args.trail > 0 else []
        draw_map(fig, ax, m, args, trail, cax=cax)
        info.set_text("\n".join(_info_lines(m, i, len(maps), args)))
        title()
        fig.canvas.draw_idle()

    def goto(i, from_slider=False):
        state["i"] = min(max(i, 0), len(maps) - 1)
        if not from_slider:
            # set_val would re-enter goto through on_changed; mute it for the round trip
            slider.eventson = False
            slider.set_val(state["i"] + 1)
            slider.eventson = True
        render()

    slider.on_changed(lambda v: goto(int(round(v)) - 1, from_slider=True))

    timer = fig.canvas.new_timer(interval=int(1000.0 / state["fps"]))

    def tick():
        if state["playing"]:
            goto(0 if state["i"] + 1 >= len(maps) else state["i"] + 1)

    timer.add_callback(tick)

    def set_fps(new):
        state["fps"] = min(max(new, _PLAY_FPS_MIN), _PLAY_FPS_MAX)
        timer.interval = int(1000.0 / state["fps"])
        title()
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key == " ":
            state["playing"] = not state["playing"]
            title()
            fig.canvas.draw_idle()
        elif event.key == "right":
            state["playing"] = False
            goto(state["i"] + 1)
        elif event.key == "left":
            state["playing"] = False
            goto(state["i"] - 1)
        elif event.key in ("+", "="):
            set_fps(state["fps"] * 1.5)
        elif event.key == "-":
            set_fps(state["fps"] / 1.5)

    fig.canvas.mpl_connect("key_press_event", on_key)
    render()
    timer.start()
    plt.show()
    return fig, state


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="file written via --sensing-dump")
    ap.add_argument("--index", type=int, action="append",
                    help="plot only this map, 0-based, negative counts from the end "
                         "(repeatable; default: the last one)")
    ap.add_argument("--k-step", type=int,
                    help="keep only maps on this pilot spacing (2 = rank 1 DMRS "
                         "type 1, 4 = despread rank 2 type 1, 6 = despread type 2)")
    ap.add_argument("--ports", type=lambda v: int(v, 0),
                    help="keep only maps of this DM-RS port bitmap, e.g. 1 for rank "
                         "1 on port 1000, 3 for rank 2 on ports 1000+1001")
    ap.add_argument("--layer", type=int,
                    help="keep only maps of this layer of the allocation")
    ap.add_argument("--aarx", type=int, help="keep only maps from this Rx antenna")
    ap.add_argument("--x", choices=("speed", "doppler"), default="speed",
                    help="horizontal axis unit (default: speed)")
    ap.add_argument("--carrier", type=float,
                    help="override the carrier frequency in Hz; by default the "
                         "carrier_hz column of the dump is used")
    ap.add_argument("--max-range", type=float,
                    help="clip the range axis at this many metres")
    ap.add_argument("--floor", type=float, default=-40.0,
                    help="lowest dB shown (default: -40)")
    ap.add_argument("--no-markers", action="store_true",
                    help="do not draw what --sensing-tdd-detect decided, even when "
                         "the dump carries it")
    ap.add_argument("--norm", choices=("peak", "median"), default="peak",
                help="dB reference: the strongest cell (default) or the map median")
    ap.add_argument("--no-aoa", action="store_true",
                    help="do not draw the angle of arrival cells, even when the dump "
                         "carries them")
    ap.add_argument("--list", action="store_true",
                    help="print a summary of every map and exit")
    ap.add_argument("--browse", action="store_true",
                    help="show the maps one at a time in one window and step through them "
                         "with the arrow keys; nothing is written to disk. Starts at the "
                         "first --index if given, else at the first map")
    ap.add_argument("--play", action="store_true",
                    help="play the maps back like a video in one window; nothing is written "
                         "to disk. space pauses, the arrows step while paused, -/+ change "
                         "the rate, and a slider scrubs. Starts at the first --index if "
                         "given, else at the first map")
    ap.add_argument("--fps", type=float, default=5.0,
                    help="starting frame rate for --play (default 5); adjust live with -/+")
    ap.add_argument("--trail", type=int, default=5,
                    help="with --browse or --play, also draw faded the targets of this many "
                         "previous maps, to see whether a detection moves like a real target "
                         "(default 5, 0 to turn off)")
    ap.add_argument("-o", "--out", help="write the figure here instead of showing it")
    args = ap.parse_args()

    maps = load(args.csv)
    if args.k_step is not None:
        maps = [m for m in maps if m.k_step == args.k_step]
    if args.ports is not None:
        maps = [m for m in maps if m.ports == args.ports]
    if args.layer is not None:
        maps = [m for m in maps if m.layer == args.layer]
    if args.aarx is not None:
        maps = [m for m in maps if m.aarx == args.aarx]
    if not maps:
        sys.exit("no maps left after filtering")

    if args.list:
        for i, m in enumerate(maps):
            r, v = m.peak(args.carrier)
            # flagged in a column at the left, so maps where the detector actually
            # threw something out can be picked out by eye in a long listing
            flag = "!!" if m.n_replica_rejected else "  "
            idx = f"[{i}]"
            print(f"{idx} {flag} {m.label(args.carrier)}, "
                  f"doppler resolution {m.doppler_resolution_hz:.1f} Hz, "
                  f"peak at {r:.1f} m / {v:+.2f} m/s")
            pad = " " * (len(idx) + 4)
            if m.n_pilots_min is not None:
                print(pad + f"grants: {m.n_pilots_min} to {m.n_pilots_max} pilots, "
                            f"{m.n_positions} distinct position(s) a_m")
            for line in m.detection_lines(args.carrier):
                print(pad + line)
            for line in m.aoa_lines():
                print(pad + line)
        return

    if args.browse or args.play:
        if args.browse and args.play:
            sys.exit("--browse and --play are two different viewers; pick one")
        mode = "--browse" if args.browse else "--play"
        if args.out:
            sys.exit(f"{mode} shows the maps in a window; it cannot be combined with --out")
        if args.play and args.fps <= 0:
            sys.exit("--fps must be positive")
        start = args.index[0] if args.index else 0
        if not -len(maps) <= start < len(maps):
            sys.exit(f"--index out of range, {len(maps)} map(s) available")
        if args.browse:
            browse(maps, args, start)
        else:
            play(maps, args, start % len(maps), args.fps)
        return

    wanted = args.index if args.index else [-1]
    try:
        chosen = [maps[i] for i in wanted]
    except IndexError:
        sys.exit(f"--index out of range, {len(maps)} map(s) available")

    import matplotlib
    if args.out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(chosen), squeeze=False,
                             figsize=(6.5 * len(chosen), 5.5))
    for ax, m in zip(axes[0], chosen):
        draw_map(fig, ax, m, args)

    fig.tight_layout()
    if args.out:
        fig.savefig(args.out, dpi=140)
        print(f"wrote {args.out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
