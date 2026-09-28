#!/usr/bin/env python3
"""Draw the static cell map of a cell NR-Scope has decoded off the air.

Reads a cell_<date>_<time>_pci<N>.json summary written by the recorder (see
RunRecorder::record_cell_summary) and draws where the SS/PBCH block, CORESET#0,
SIB1 region and initial DL BWP sit in the resource grid, in resource blocks
relative to Point A. Every value comes from what the sniffer decoded from the
commercial gNB; nothing is read from a gNB-side log.

Usage:
    python3 cell_map.py cell/cell_2026-09-26_16-17-15_pci632.json
    python3 cell_map.py <cell.json> -o out.png        # save instead of show
"""

import argparse
import json
import sys

import matplotlib

matplotlib.use("Agg")  # switched to an interactive backend below if showing
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# One colour per logical channel. Kept deliberately small and flat.
C_BWP = "#e9edf2"      # initial DL BWP / carrier background
C_SSB = "#f3d9a0"      # SS/PBCH block envelope
C_PSS = "#9fb8d6"
C_SSS = "#7fbf7f"
C_PBCH = "#f0e2b6"
C_CORESET0 = "#c9743a"  # CORESET#0 (also carries SIB1's PDCCH)
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


def riv_to_start_len(riv, n_bwp):
    """Decode a locationAndBandwidth RIV into (start_rb, len_rb)."""
    if n_bwp <= 0:
        return 0, 0
    len_minus_1 = riv // n_bwp
    start = riv % n_bwp
    if len_minus_1 <= n_bwp // 2:
        length = len_minus_1 + 1
    else:
        length = n_bwp - len_minus_1 + 1
        start = n_bwp - 1 - start
    return start, length


SYMB_PER_SLOT = 14


def _draw_background(ax, slot, carrier_lo_rb, carrier_bw, bwp_start, bwp_len,
                     common_scs, label_carrier):
    """Carrier and initial-DL-BWP bands, spanning the whole slot."""
    if carrier_bw:
        ax.add_patch(Rectangle((-0.5, carrier_lo_rb), SYMB_PER_SLOT, carrier_bw,
                               facecolor=C_BWP, edgecolor="#b8c0cc", lw=1, zorder=0))
        if label_carrier:
            ax.text(SYMB_PER_SLOT - 0.3, carrier_lo_rb + carrier_bw - 1,
                    f"carrier {carrier_bw} RB @ {common_scs} kHz",
                    ha="right", va="top", fontsize=8, color="#6b7280")
    if bwp_len:
        ax.add_patch(Rectangle((-0.5, carrier_lo_rb + bwp_start), SYMB_PER_SLOT, bwp_len,
                               facecolor="none", edgecolor="#8aa0bf", lw=1.4,
                               ls="--", zorder=1))
        if label_carrier:
            ax.text(-0.3, carrier_lo_rb + bwp_start + bwp_len, "initial DL BWP",
                    ha="left", va="bottom", fontsize=8, color="#5b6b82")


def _draw_ssb(ax, s0, ssb_lo_rb, ssb_h_rb):
    """SS/PBCH block: PSS / PBCH / SSS+PBCH / PBCH over 4 symbols from s0."""
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
    for k, lab in ((0, "PSS"), (1, "PBCH"), (2, "SSS"), (3, "PBCH")):
        ax.text(s0 + k + 0.5, ssb_lo_rb + ssb_h_rb / 2, lab,
                rotation=90, ha="center", va="center", fontsize=7.5, color=C_TEXT)


def _draw_coreset0(ax, s0, dur, cs0_lo_rb, cs0_h_rb):
    ax.add_patch(Rectangle((s0, cs0_lo_rb), dur, cs0_h_rb,
                           facecolor=C_CORESET0, edgecolor="#8a4a24", lw=1.2, zorder=3))
    ax.text(s0 + dur / 2, cs0_lo_rb + cs0_h_rb / 2,
            f"CORESET#0\n{cs0_h_rb} RB x {dur} sym",
            ha="center", va="center", fontsize=8, color="white", zorder=4)


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

    # Carrier / initial DL BWP. carrier_bw is in common-SCS RBs from Point A,
    # shifted by offsetToPointA (which counts 15 kHz RBs).
    carrier_lo_rb = c["offset_to_point_a_rb"] * 15.0 / common_scs
    carrier_bw = c.get("carrier_bw_rb", 0) or 0
    bwp_start, bwp_len = riv_to_start_len(c.get("init_dl_bwp_riv", 0), carrier_bw)

    # --- exact time-domain placement ---
    l_ssb = ssb_first_symbol(c.get("ssb_pattern", "C"), c.get("ssb_idx", 0))
    ssb_slot = l_ssb // SYMB_PER_SLOT
    ssb_sym = l_ssb % SYMB_PER_SLOT
    cs0_slot = c.get("coreset0_slot_n0", 0)
    cs0_sym = c.get("coreset0_first_symbol", 0)
    sfn_txt = "odd SFN" if c.get("coreset0_sfn_c", 0) else "even SFN"

    y_top = max(carrier_lo_rb + carrier_bw, ssb_lo_rb + ssb_h_rb, cs0_lo_rb + cs0_h_rb) + 3

    # One panel per distinct slot, so blocks sit at their true symbol index.
    same_slot = ssb_slot == cs0_slot
    slots = [ssb_slot] if same_slot else sorted({ssb_slot, cs0_slot})
    fig, axes = plt.subplots(1, len(slots), figsize=(5.0 * len(slots), 8),
                             sharey=True, squeeze=False)
    axes = axes[0]

    for ax, slot in zip(axes, slots):
        first = ax is axes[0]
        _draw_background(ax, slot, carrier_lo_rb, carrier_bw, bwp_start, bwp_len,
                         common_scs, label_carrier=first)
        if slot == ssb_slot:
            _draw_ssb(ax, ssb_sym, ssb_lo_rb, ssb_h_rb)
        if slot == cs0_slot:
            _draw_coreset0(ax, cs0_sym, cs0_dur, cs0_lo_rb, cs0_h_rb)

        ax.axhline(0, color=C_POINTA, lw=1.6, zorder=5)
        if first:
            ax.text(-0.3, 0, "Point A (CRB 0)", ha="left", va="bottom",
                    fontsize=8, color=C_POINTA)
        ax.set_xlim(-0.5, SYMB_PER_SLOT - 0.5)
        ax.set_ylim(-3, y_top)
        ax.set_xticks(range(0, SYMB_PER_SLOT, 2))
        ax.set_xlabel(f"OFDM symbol in slot {slot}")
        for b in range(SYMB_PER_SLOT + 1):
            ax.axvline(b - 0.5, color="#dfe3e8", lw=0.5, zorder=0)

    axes[0].set_ylabel(f"Resource blocks from Point A ({common_scs} kHz)")
    fig.suptitle(
        f"Cell map  -  PCI {c['pci']}  -  SSB {c['ssb_center_freq_hz']/1e6:.3f} MHz "
        f"(pattern {c.get('ssb_pattern','?')}), kSSB {c['k_ssb']}\n"
        f"CORESET0 idx {c['coreset0_idx']} -> slot {cs0_slot} ({sfn_txt}), "
        f"sym {cs0_sym}..{cs0_sym + cs0_dur - 1}   |   "
        f"SSB#{c.get('ssb_idx',0)} -> slot {ssb_slot}, sym {ssb_sym}..{ssb_sym + 3}",
        fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cell_json", help="cell_*.json written by the recorder")
    ap.add_argument("-o", "--out", help="save PNG here instead of showing")
    args = ap.parse_args(argv)

    c = load_summary(args.cell_json)
    fig = build_figure(c)

    if args.out:
        fig.savefig(args.out, dpi=150)
        print(f"wrote {args.out}")
    else:
        try:
            matplotlib.use("TkAgg", force=True)
        except Exception:
            pass
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
