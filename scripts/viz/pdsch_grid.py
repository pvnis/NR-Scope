#!/usr/bin/env python3
"""Draw the PDSCH resource grid NR-Scope decoded: which REs carry DM-RS, which carry data.

For a few consecutive slots, every resource element (RE) of the carrier is
coloured by what the decoded DCIs put there:

  green   PDSCH DM-RS
  white   PDSCH data
  grey    not scheduled to the UE (control region, other slots, unused PRBs)
  hatched DM-RS symbol subcarriers left empty (CDM group without data)

Two panels share the time axis (OFDM symbols, slot boundaries marked):
  - overview: the whole carrier, one cell per PRB and symbol, to see the size
    and position of each grant and the DM-RS symbols changing from slot to slot;
  - zoom: a few PRBs at RE resolution, to see the DM-RS comb.

Everything comes from NR-Scope's DCI CSV (prbs, time_start, time_length,
dmrs_symbols, nof_dmrs_cdm_groups); the DM-RS check in the same CSV confirms the
pilots are really there on air. DM-RS configuration type 1 is assumed (the
comb is every other subcarrier, CDM group 0 on the even ones), which is what
this cell uses.

Usage:
    python3 pdsch_grid.py                         # newest DCIs/ file, busiest second, 10 slots
    python3 pdsch_grid.py DCIs/dci_<run>_pci1.csv --at 512.4 --slots 6
    python3 pdsch_grid.py --zoom-prb 30 --zoom-nprb 6 --show
"""

import argparse
import csv
import glob
import os
from collections import Counter, defaultdict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SLOTS_PER_FRAME = 20  # 30 kHz SCS
NSYMB = 14
NRE = 12

# RE categories and colours (value in the grid -> colour)
UNUSED, DATA, DMRS, EMPTY = 0, 1, 2, 3
C_UNUSED = "#e6e6e3"
C_DATA = "#ffffff"
C_DMRS = "#008300"
C_EMPTY = "#e6e6e3"   # drawn with a hatch on top, so it reads apart from UNUSED
C_GRID = "#c9c9c4"
C_SLOT = "#4a4a46"
C_TEXT = "#2b2b28"
C_MUTED = "#6b6b66"


def parse_prbs(s):
    """'0-272' or '0-9,20-29' -> set of PRB indices."""
    out = set()
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        a, _, b = part.partition("-")
        out.update(range(int(a), int(b or a) + 1))
    return out


def load(path):
    rows = [r for r in csv.DictReader(open(path)) if r["direction"] == "DL" and r["prbs"] and r["time_length"]]
    if not rows:
        raise SystemExit(f"{path}: no downlink grants with a PDSCH allocation")
    rows.sort(key=lambda r: float(r["timestamp"]))
    # Continuous slot index: SFN wraps every 10.24 s, unwrap with the timestamps
    t0, s0 = float(rows[0]["timestamp"]), int(rows[0]["sfn"]) * SLOTS_PER_FRAME + int(rows[0]["slot"])
    period = 1024 * SLOTS_PER_FRAME
    grants = []
    for r in rows:
        s = int(r["sfn"]) * SLOTS_PER_FRAME + int(r["slot"]) + int(r["k"] or 0)
        k = round(((float(r["timestamp"]) - t0) - (s - s0) / (SLOTS_PER_FRAME * 100)) / 10.24)
        symbols = [int(x) for x in r.get("dmrs_symbols", "").split()] if r.get("dmrs_symbols") else None
        grants.append(dict(abs=k * period + s, sfn=int(r["sfn"]), slot=int(r["slot"]), rnti=int(r["rnti"]),
                           prbs=parse_prbs(r["prbs"]), S=int(r["time_start"]), L=int(r["time_length"]),
                           dmrs=symbols, cdm=int(r["nof_dmrs_cdm_groups"] or 1), t=float(r["timestamp"]),
                           add_pos=r["dmrs_add_pos"], typea=r["dmrs_typeA_pos"]))
    return grants


def dmrs_symbols_type_a(typea, add_pos, S, L):
    """TS 38.211 table 7.4.1.1.2-3, single symbol, for CSVs without dmrs_symbols."""
    ld, l0, add = S + L, int(typea), int(add_pos)
    if ld <= 7 or add == 0:
        extra = ()
    elif ld <= 9:
        extra = (7,)
    elif ld <= 11:
        extra = (9,) if add == 1 else (6, 9)
    elif ld == 12:
        extra = {1: (9,), 2: (6, 9)}.get(add, (5, 8, 11))
    else:
        extra = {1: (11,), 2: (7, 11)}.get(add, (5, 8, 11))
    return [l0, *extra]


def pick_start(grants, at):
    if at:
        sfn, slot = (int(x) for x in at.split("."))
        for g in grants:
            if g["sfn"] == sfn and g["slot"] == slot:
                return g["abs"]
        raise SystemExit(f"no downlink grant at SFN {at}")
    # Default: the first slot of the busiest second, where the carrier is loaded
    t0 = grants[0]["t"]
    busiest = Counter(int(g["t"] - t0) for g in grants).most_common(1)[0][0]
    return min(g["abs"] for g in grants if int(g["t"] - t0) == busiest)


def build_grid(grants, start, nslots, nof_prb):
    """RE grid [subcarrier, symbol] over nslots slots, and the grants drawn."""
    grid = np.full((nof_prb * NRE, nslots * NSYMB), UNUSED, dtype=np.int8)
    by_slot = defaultdict(list)
    for g in grants:
        if start <= g["abs"] < start + nslots:
            by_slot[g["abs"]].append(g)
    for abs_slot, gs in by_slot.items():
        base = (abs_slot - start) * NSYMB
        for g in gs:
            dmrs = g["dmrs"] if g["dmrs"] is not None else dmrs_symbols_type_a(g["typea"], g["add_pos"], g["S"], g["L"])
            for prb in g["prbs"]:
                if prb >= nof_prb:
                    continue
                rows = slice(prb * NRE, (prb + 1) * NRE)
                for l in range(g["S"], min(g["S"] + g["L"], NSYMB)):
                    if l in dmrs:
                        col = np.full(NRE, EMPTY if g["cdm"] >= 2 else DATA, dtype=np.int8)
                        col[0::2] = DMRS  # type 1, CDM group 0: even subcarriers
                        grid[rows, base + l] = col
                    else:
                        grid[rows, base + l] = DATA
    return grid, by_slot


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dci_csv", nargs="?", help="default: the newest file in DCIs/")
    ap.add_argument("--at", help="first slot to draw, as SFN.slot (default: start of the busiest second)")
    ap.add_argument("--slots", type=int, default=10, help="consecutive slots to draw (default 10, half a frame)")
    ap.add_argument("--carrier-prb", type=int, default=273)
    ap.add_argument("--zoom-prb", type=int, default=0, help="first PRB of the RE-resolution zoom (default 0)")
    ap.add_argument("--zoom-nprb", type=int, default=4, help="PRBs in the zoom (default 4)")
    ap.add_argument("-o", "--output", help="PNG path (default: beside the CSV)")
    ap.add_argument("--show", action="store_true")
    a = ap.parse_args()

    path = a.dci_csv or max(glob.glob(os.path.join(ROOT, "DCIs", "dci_*.csv")), key=os.path.getmtime)
    grants = load(path)
    start = pick_start(grants, a.at)
    grid, by_slot = build_grid(grants, start, a.slots, a.carrier_prb)

    if not a.show:
        matplotlib.use("Agg")
    cmap = ListedColormap([C_UNUSED, C_DATA, C_DMRS, C_EMPTY])
    nsym = a.slots * NSYMB
    fig, (ax_all, ax_zoom) = plt.subplots(
        2, 1, figsize=(max(9, 0.28 * nsym + 2), 10.5), gridspec_kw=dict(height_ratios=[1.35, 1]), sharex=True)
    fig.patch.set_facecolor("#fcfcfb")

    # Overview: one cell per PRB and symbol. A PRB on a DM-RS symbol is drawn
    # green (it holds DM-RS on half its subcarriers); the comb is in the zoom.
    prb_grid = grid.reshape(a.carrier_prb, NRE, nsym).max(axis=1)
    ax_all.imshow(prb_grid, cmap=cmap, vmin=0, vmax=3, aspect="auto", origin="lower", interpolation="nearest",
                  extent=(0, nsym, 0, a.carrier_prb))
    ax_all.set_ylabel("PRB (carrier)", color=C_TEXT)
    ax_all.set_title("Whole carrier: grant size and position, DM-RS symbols (one cell = 1 PRB x 1 symbol)",
                     fontsize=10, color=C_MUTED, loc="left", pad=46)

    # Zoom: RE resolution with a thin grid
    z0, zn = a.zoom_prb, min(a.zoom_nprb, a.carrier_prb - a.zoom_prb)
    sub = grid[z0 * NRE:(z0 + zn) * NRE, :]
    ax_zoom.imshow(sub, cmap=cmap, vmin=0, vmax=3, aspect="auto", origin="lower", interpolation="nearest",
                   extent=(0, nsym, z0 * NRE, (z0 + zn) * NRE))
    ys, xs = np.nonzero(sub == EMPTY)
    for y, x in zip(ys, xs):
        ax_zoom.add_patch(plt.Rectangle((x, z0 * NRE + y), 1, 1, fill=False, hatch="////", edgecolor=C_MUTED, lw=0))
    for x in range(nsym + 1):
        ax_zoom.axvline(x, color=C_GRID, lw=0.4)
    for y in range(z0 * NRE, (z0 + zn) * NRE + 1):
        ax_zoom.axhline(y, color=C_GRID, lw=1.0 if y % NRE == 0 else 0.4)
    ax_zoom.set_yticks([z0 * NRE + NRE * i + NRE / 2 for i in range(zn)])
    ax_zoom.set_yticklabels([f"PRB {z0 + i}" for i in range(zn)])
    ax_zoom.set_ylabel("subcarriers", color=C_TEXT)
    ax_zoom.set_title(f"Zoom on PRBs {z0}-{z0 + zn - 1}: one cell = 1 RE (subcarrier x symbol), the DM-RS comb",
                      fontsize=10, color=C_MUTED, loc="left")

    # Slot boundaries and labels on both panels
    for i in range(a.slots):
        abs_slot = start + i
        gs = by_slot.get(abs_slot, [])
        label = f"SFN {gs[0]['sfn']}.{gs[0]['slot']}" if gs else "no DL grant\n(uplink slot,\nor UE not scheduled)"
        if gs:
            g = gs[0]
            dm = g["dmrs"] if g["dmrs"] is not None else dmrs_symbols_type_a(g["typea"], g["add_pos"], g["S"], g["L"])
            label += f"\n{len(g['prbs'])} PRB, symb {g['S']}-{g['S'] + g['L'] - 1}\nDM-RS {' '.join(map(str, dm))}"
        ax_all.text(i * NSYMB + NSYMB / 2, a.carrier_prb + 3, label, ha="center", va="bottom", fontsize=8.5, color=C_TEXT)
        for ax in (ax_all, ax_zoom):
            ax.axvline(i * NSYMB, color=C_SLOT, lw=1.6)
    for ax in (ax_all, ax_zoom):
        ax.axvline(nsym, color=C_SLOT, lw=1.6)
        ax.tick_params(colors=C_MUTED, labelsize=8)
        for sp in ax.spines.values():
            sp.set_visible(False)
    ax_zoom.set_xticks([i + 0.5 for i in range(nsym)])
    ax_zoom.set_xticklabels([str(i % NSYMB) for i in range(nsym)], fontsize=7)
    ax_zoom.set_xlabel("OFDM symbol within the slot", color=C_TEXT)

    rnti = Counter(g["rnti"] for gs in by_slot.values() for g in gs).most_common(1)
    who = f"UE 0x{rnti[0][0]:04x}" if rnti else "no UE"
    fig.suptitle(f"PDSCH resource grid decoded by NR-Scope, {who}, {a.slots} slots from "
                 f"{os.path.basename(path)}", x=0.01, ha="left", fontsize=12, color=C_TEXT, y=0.995)
    fig.legend(handles=[Patch(facecolor=C_DMRS, label="PDSCH DM-RS"),
                        Patch(facecolor=C_DATA, edgecolor=C_GRID, label="PDSCH data"),
                        Patch(facecolor=C_EMPTY, edgecolor=C_MUTED, hatch="////", label="empty (CDM group without data)"),
                        Patch(facecolor=C_UNUSED, label="not scheduled to the UE")],
               loc="upper right", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.99, 0.975))
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    out = a.output or os.path.splitext(path)[0] + f"_grid_{start % (1024 * SLOTS_PER_FRAME)}.png"
    fig.savefig(out, dpi=140, facecolor=fig.get_facecolor())
    print(f"wrote {out}")
    if a.show:
        plt.show()


if __name__ == "__main__":
    main()
