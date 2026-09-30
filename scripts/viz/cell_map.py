#!/usr/bin/env python3
"""Draw the static cell map of a cell NR-Scope has decoded off the air.

Reads a cell_<date>_<time>_pci<N>.json summary written by the recorder (see
RunRecorder::record_cell_summary) and draws where the SS/PBCH block, CORESET#0
and SIB1 PDSCH sit in the resource grid, in resource blocks relative to Point A,
with a legend row of descriptive blocks (sizes / frequencies) above it. Every
value comes from what the sniffer decoded from the commercial gNB; nothing is
read from a gNB-side log.

Usage:
    python3 cell_map.py cell/cell_2026-09-26_16-17-15_pci632.json   # -> .png beside it
    python3 cell_map.py <cell.json> -o out.png                      # choose the path
    python3 cell_map.py <cell.json> --show                          # GUI window (needs X)
"""

import argparse
import json
import os
import sys

import matplotlib

# The backend is selected in main() before any figure is built: Agg when saving
# (the default, and what works over a plain SSH session), an interactive one
# only when --show is asked for and one is actually available.
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# One colour per logical channel. Kept deliberately small and flat.
C_BWP = "#e9edf2"      # initial DL BWP / carrier background
C_SSB = "#f3d9a0"      # SS/PBCH block envelope
C_PSS = "#9fb8d6"
C_SSS = "#7fbf7f"
C_PBCH = "#f0e2b6"
C_CORESET0 = "#c9743a"  # CORESET#0 (carries SIB1's PDCCH / SI-RNTI DCI)
C_SIB1 = "#5b8a72"      # SIB1 PDSCH (the payload the DCI points to)
C_POINTA = "#c0392b"
C_TEXT = "#1a1a1a"

NSC_PER_RB = 12
SSB_NOF_SC = 240      # 20 RBs
SSB_NOF_SYMB = 4


def load_summary(path):
    with open(path) as f:
        return json.load(f)


def hz_to_crb(freq_hz, point_a_hz, common_scs_khz):
    """Common-SCS resource-block index of an absolute frequency, from Point A."""
    return (freq_hz - point_a_hz) / (NSC_PER_RB * common_scs_khz * 1e3)


def ssb_first_symbol(pattern, ssb_idx):
    """First OFDM symbol of an SS/PBCH block within the half radio frame.

    TS 38.213 sec 4.1. Pattern A/C share the {2,8}+14n layout, B uses
    {4,8,16,20}+28n. Returned in symbols from the start of the half frame; an
    SSB spans this symbol and the three that follow.
    """
    if pattern in ("A", "C"):
        base, step = [2, 8], 14
    elif pattern == "B":
        base, step = [4, 8, 16, 20], 28
    else:  # FR2 patterns D/E, not expected here
        base, step = [4, 8, 16, 20], 28
    k = len(base)
    return base[ssb_idx % k] + step * (ssb_idx // k)


SYMB_PER_SLOT = 14


def _draw_carrier(ax, carrier_lo_rb, carrier_bw):
    """Carrier band, spanning the whole slot. No labels: those go in the legend."""
    if carrier_bw:
        ax.add_patch(Rectangle((-0.5, carrier_lo_rb), SYMB_PER_SLOT, carrier_bw,
                               facecolor=C_BWP, edgecolor="#b8c0cc", lw=1, zorder=0))


def _draw_ssb(ax, s0, ssb_lo_rb, ssb_h_rb):
    """SS/PBCH block: PSS / PBCH / SSS+PBCH / PBCH over 4 symbols from s0."""
    # Ticks sit on symbol centres, so symbol s spans s-0.5 .. s+0.5.
    s0 = s0 - 0.5
    ax.add_patch(Rectangle((s0, ssb_lo_rb), SSB_NOF_SYMB, ssb_h_rb,
                           facecolor=C_SSB, edgecolor="#c9a24a", lw=1.2, zorder=2))
    mid_lo = ssb_lo_rb + ssb_h_rb * (56 / SSB_NOF_SC)   # centre 127 subcarriers
    mid_h = ssb_h_rb * (127 / SSB_NOF_SC)
    ax.add_patch(Rectangle((s0, mid_lo), 1, mid_h,
                           facecolor=C_PSS, edgecolor="#4a6fa5", lw=0.8, zorder=3))
    for k in (1, 2, 3):
        ax.add_patch(Rectangle((s0 + k, ssb_lo_rb), 1, ssb_h_rb,
                               facecolor=C_PBCH, edgecolor="#c9a24a", lw=0.6, zorder=3))
    ax.add_patch(Rectangle((s0 + 2, mid_lo), 1, mid_h,
                           facecolor=C_SSS, edgecolor="#3f8f3f", lw=0.8, zorder=4))


def _draw_coreset0(ax, s0, dur, cs0_lo_rb, cs0_h_rb):
    s0 = s0 - 0.5  # symbol s spans s-0.5 .. s+0.5
    ax.add_patch(Rectangle((s0, cs0_lo_rb), dur, cs0_h_rb,
                           facecolor=C_CORESET0, edgecolor="#8a4a24", lw=1.2, zorder=3))


def _draw_sib1(ax, s0, dur, lo_rb, h_rb):
    # Below the SSB (zorder 2+): the PDSCH is rate-matched around the SSB REs,
    # so where they overlap in time the SSB should show on top.
    s0 = s0 - 0.5  # symbol s spans s-0.5 .. s+0.5
    ax.add_patch(Rectangle((s0, lo_rb), dur, h_rb,
                           facecolor=C_SIB1, edgecolor="#3c5e4d", lw=1.2, zorder=1.6))


def _draw_legend(ax, entries):
    """A row of descriptive blocks (colour swatch + text) above the grid."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    n = len(entries)
    col_w = 1.0 / n
    for i, e in enumerate(entries):
        x0 = i * col_w
        ax.add_patch(Rectangle((x0 + 0.006, 0.04), col_w - 0.012, 0.92,
                               facecolor="#fafafa", edgecolor="#dcdcdc", lw=1))
        ax.add_patch(Rectangle((x0 + 0.03, 0.70), 0.028, 0.20,
                               facecolor=e["color"], edgecolor="#333", lw=0.8))
        ax.text(x0 + 0.072, 0.80, e["title"], fontsize=9.5, fontweight="bold",
                va="center", ha="left")
        ax.text(x0 + 0.03, 0.60, "\n".join(e["lines"]), fontsize=8.2,
                va="top", ha="left", linespacing=1.5, color="#333")


def build_figure(c):
    common_scs = c["common_scs_khz"]
    ssb_scs = c["ssb_scs_khz"]

    # Point A in Hz, recovered from CORESET#0's known lowest subcarrier and its
    # offset (in common-SCS RBs) from Point A: both are sniffer-computed.
    point_a_hz = (
        c["coreset0_lower_freq_hz"]
        - c["coreset0_offset_rb"] * NSC_PER_RB * common_scs * 1e3
    )

    # SSB envelope in common-SCS RBs.
    ssb_lo_hz = c["ssb_center_freq_hz"] - (SSB_NOF_SC / 2) * ssb_scs * 1e3
    ssb_lo_rb = hz_to_crb(ssb_lo_hz, point_a_hz, common_scs)
    ssb_h_rb = SSB_NOF_SC * ssb_scs / (NSC_PER_RB * common_scs)

    # CORESET#0 frequency extent.
    cs0_lo_rb = c["coreset0_offset_rb"]
    cs0_h_rb = c["coreset0_bw_rb"]
    cs0_dur = max(1, c["coreset0_duration_symbols"])

    # Carrier band. Its position relative to Point A is offsetToCarrier (in the
    # carrier's own SCS, i.e. common-SCS RBs here), NOT offsetToPointA -- that
    # one positions the SSB. With offsetToCarrier = 0 the carrier starts at
    # Point A. carrier_bw is in common-SCS RBs.
    carrier_lo_rb = c.get("carrier_offset_to_carrier", 0)
    carrier_bw = c.get("carrier_bw_rb", 0) or 0

    # --- exact time-domain placement ---
    l_ssb = ssb_first_symbol(c.get("ssb_pattern", "C"), c.get("ssb_idx", 0))
    ssb_slot = l_ssb // SYMB_PER_SLOT
    ssb_sym = l_ssb % SYMB_PER_SLOT
    cs0_slot = c.get("coreset0_slot_n0", 0)
    cs0_sym = c.get("coreset0_first_symbol", 0)
    sfn_txt = "odd SFN" if c.get("coreset0_sfn_c", 0) else "even SFN"

    # SIB1 PDSCH, if this capture decoded it (prb count 0 => not captured).
    sib1_nof_prb = c.get("sib1_nof_prb", 0) or 0
    have_sib1 = sib1_nof_prb > 0
    sib1_lo_rb = cs0_lo_rb + c.get("sib1_prb_start", 0)   # PRBs are on the CORESET#0 grid
    sib1_h_rb = sib1_nof_prb
    sib1_sym = c.get("sib1_symbol_start", 0)
    sib1_dur = max(1, c.get("sib1_nof_symbols", 1))
    sib1_slot = c.get("sib1_slot_idx", cs0_slot)

    y_top = max(carrier_lo_rb + carrier_bw, ssb_lo_rb + ssb_h_rb, cs0_lo_rb + cs0_h_rb,
                (sib1_lo_rb + sib1_h_rb) if have_sib1 else 0) + 3

    # --- descriptive blocks for the legend row (sizes / frequencies) ---
    rb_hz = NSC_PER_RB * common_scs * 1e3
    cs0_hi_hz = c["coreset0_lower_freq_hz"] + cs0_h_rb * rb_hz
    carrier_lo_hz = point_a_hz + carrier_lo_rb * rb_hz
    entries = [
        {"color": C_SSB, "title": "SS/PBCH block", "lines": [
            f"{int(round(ssb_h_rb))} RB x {SSB_NOF_SYMB} sym",
            f"centre {c['ssb_center_freq_hz']/1e6:.3f} MHz",
            f"pattern {c.get('ssb_pattern','?')} - SSB#{c.get('ssb_idx',0)}",
            f"kSSB {c['k_ssb']} - slot {ssb_slot}, sym {ssb_sym}-{ssb_sym+3}"]},
        {"color": C_CORESET0, "title": "CORESET#0", "lines": [
            f"{cs0_h_rb} RB x {cs0_dur} sym",
            f"{c['coreset0_lower_freq_hz']/1e6:.3f}-{cs0_hi_hz/1e6:.3f} MHz",
            f"idx {c['coreset0_idx']} - Type0-PDCCH",
            f"slot {cs0_slot} ({sfn_txt}), sym {cs0_sym}-{cs0_sym+cs0_dur-1}"]},
    ]
    if have_sib1:
        sib1_lo_hz = point_a_hz + sib1_lo_rb * rb_hz
        sib1_hi_hz = sib1_lo_hz + sib1_h_rb * rb_hz
        entries.append({"color": C_SIB1, "title": "SIB1 PDSCH", "lines": [
            f"{sib1_h_rb} RB x {sib1_dur} sym",
            f"{sib1_lo_hz/1e6:.3f}-{sib1_hi_hz/1e6:.3f} MHz",
            f"slot {sib1_slot}, sym {sib1_sym}-{sib1_sym+sib1_dur-1}",
            "scheduled by SI-RNTI DCI"]})
    else:
        entries.append({"color": C_SIB1, "title": "SIB1 PDSCH",
                        "lines": ["not captured", "(no SI-RNTI DCI decoded)"]})
    entries.append({"color": C_BWP, "title": "Carrier", "lines": [
        f"{carrier_bw} RB @ {common_scs} kHz",
        f"{carrier_lo_hz/1e6:.3f}-{(carrier_lo_hz + carrier_bw*rb_hz)/1e6:.3f} MHz",
        f"Point A {point_a_hz/1e6:.3f} MHz"]})

    # One grid panel per distinct slot, so blocks sit at their true symbol index.
    used_slots = {ssb_slot, cs0_slot}
    if have_sib1:
        used_slots.add(sib1_slot)
    slots = sorted(used_slots)

    fig = plt.figure(figsize=(max(9.0, 5.0 * len(slots)), 9))
    gs = fig.add_gridspec(2, len(slots), height_ratios=[1.05, 4.2],
                          left=0.09, right=0.98, top=0.94, bottom=0.07,
                          hspace=0.18, wspace=0.08)
    legend_ax = fig.add_subplot(gs[0, :])
    _draw_legend(legend_ax, entries)

    axes = []
    for i, slot in enumerate(slots):
        ax = fig.add_subplot(gs[1, i], sharey=axes[0] if axes else None)
        axes.append(ax)
        first = i == 0
        _draw_carrier(ax, carrier_lo_rb, carrier_bw)
        if slot == ssb_slot:
            _draw_ssb(ax, ssb_sym, ssb_lo_rb, ssb_h_rb)
        if slot == cs0_slot:
            _draw_coreset0(ax, cs0_sym, cs0_dur, cs0_lo_rb, cs0_h_rb)
        if have_sib1 and slot == sib1_slot:
            _draw_sib1(ax, sib1_sym, sib1_dur, sib1_lo_rb, sib1_h_rb)

        ax.axhline(0, color=C_POINTA, lw=1.6, zorder=5)
        if first:
            ax.text(-0.3, 0, "Point A", ha="left", va="bottom",
                    fontsize=8, color=C_POINTA)
        ax.set_xlim(-0.5, SYMB_PER_SLOT - 0.5)
        ax.set_ylim(-3, y_top)
        ax.set_xticks(range(0, SYMB_PER_SLOT, 2))
        ax.set_xlabel(f"Slot {slot}")
        for b in range(SYMB_PER_SLOT + 1):
            ax.axvline(b - 0.5, color="#dfe3e8", lw=0.5, zorder=0)
        if not first:
            ax.tick_params(labelleft=False)

    axes[0].set_ylabel(f"Resource blocks from Point A ({common_scs} kHz)")
    fig.suptitle(f"Cell map  -  PCI {c['pci']}  -  {carrier_bw} RB @ {common_scs} kHz",
                 fontsize=11, fontweight="bold")
    return fig


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cell_json", help="cell_*.json written by the recorder")
    ap.add_argument("-o", "--out",
                    help="PNG path (default: alongside the JSON)")
    ap.add_argument("--show", action="store_true",
                    help="open an interactive window instead of saving "
                         "(needs a working GUI backend; unreliable over SSH)")
    args = ap.parse_args(argv)

    # Pick the backend before building the figure. Switching after a figure
    # exists is unsupported and is what breaks Agg->Tk over `ssh -X`.
    interactive = False
    if args.show and not args.out:
        for backend in ("QtAgg", "TkAgg", "GTK3Agg"):
            try:
                matplotlib.use(backend, force=True)
                interactive = True
                break
            except Exception:
                continue
        if not interactive:
            print("no interactive backend available, saving to a file instead",
                  file=sys.stderr)
    if not interactive:
        matplotlib.use("Agg", force=True)

    c = load_summary(args.cell_json)
    fig = build_figure(c)

    if interactive:
        plt.show()
    else:
        out = args.out or os.path.splitext(args.cell_json)[0] + ".png"
        fig.savefig(out, dpi=150)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
