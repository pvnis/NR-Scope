#!/usr/bin/env python3
"""Zoom into CORESET#0 and show where the SIB1 PDCCH sits on the REG grid.

The general cell map (cell_map.py) can't show the control channel: CORESET#0 is
tens of RBs inside a carrier of hundreds, so the PDCCH is a sliver. This view
zooms to CORESET#0 alone and draws its resource-element groups (REGs, 1 RB x 1
symbol), colouring the control channel elements (CCEs) that carry the SIB1
PDCCH. The CCE-to-REG-bundle mapping follows TS 38.211 7.3.2.2, reproducing the
interleaving, so the highlighted REGs land on their true positions -- which is
why an interleaved CORESET scatters one PDCCH across the band.

Reads the same cell_*.json as cell_map.py. Every value is sniffer-decoded.

Usage:
    python3 pdcch_map.py cell/cell_2026-09-29_15-06-41_pci1.json   # -> .png beside it
    python3 pdcch_map.py <cell.json> -o out.png
    python3 pdcch_map.py <cell.json> --show
"""

import argparse
import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

C_REG = "#eef1f4"       # a REG not used by the SIB1 PDCCH
C_REG_EDGE = "#c7ced6"
C_PDCCH = "#c9743a"     # a REG carrying the SIB1 PDCCH
C_TEXT = "#33393f"

REGS_PER_CCE = 6        # fixed by the standard


def load_summary(path):
    with open(path) as f:
        return json.load(f)


def cce_to_regs(n_rb, n_symb, reg_bundle, interleaved, interleaver, shift):
    """Map each CCE index to its REG indices, per TS 38.211 7.3.2.2.

    REGs are numbered time-first: REG(rb, sym) = rb*n_symb + sym. A REG bundle
    of L REGs is {iL .. iL+L-1}; a CCE is 6 REGs = 6/L bundles. For interleaved
    mapping the bundle picked for logical position x is
        f(x) = (r*C + c + shift) mod N_bundle,  x = c*R + r,  C = N_reg/(L*R).
    """
    n_reg = n_rb * n_symb
    n_bundle = n_reg // reg_bundle
    n_cce = n_reg // REGS_PER_CCE
    bundles_per_cce = REGS_PER_CCE // reg_bundle

    if interleaved and interleaver:
        R = interleaver
        C = n_reg // (reg_bundle * R)

        def f(x):
            r = x % R
            c = x // R
            return (r * C + c + shift) % n_bundle
    else:
        def f(x):
            return x

    mapping = {}
    for j in range(n_cce):
        regs = []
        for k in range(bundles_per_cce):
            b = f(bundles_per_cce * j + k)
            regs.extend(range(b * reg_bundle, (b + 1) * reg_bundle))
        mapping[j] = regs
    return mapping, n_cce


def build_figure(c):
    n_rb = c["coreset0_bw_rb"]
    n_symb = max(1, c["coreset0_duration_symbols"])
    reg_bundle = c.get("coreset0_reg_bundle_size", 6) or 6
    interleaved = bool(c.get("coreset0_interleaved", 0))
    interleaver = c.get("coreset0_interleaver_size", 0)
    shift = c.get("coreset0_shift_index", 0)
    crb0 = c["coreset0_offset_rb"]     # CORESET#0 lowest RB, from Point A

    mapping, n_cce = cce_to_regs(n_rb, n_symb, reg_bundle,
                                 interleaved, interleaver, shift)
    reg_to_cce = {reg: cce for cce, regs in mapping.items() for reg in regs}

    agg = c.get("sib1_pdcch_agg_level", 0) or 0
    ncce = c.get("sib1_pdcch_ncce", 0)
    sib1_cces = set(range(ncce, ncce + agg)) if agg else set()

    # frequency (RB) on x, time (symbol) on y -> a wide, short grid.
    fig_w = min(16.0, max(8.0, 0.26 * n_rb + 2))
    fig, ax = plt.subplots(figsize=(fig_w, 2.4 + 0.7 * n_symb))

    for rb in range(n_rb):
        for sym in range(n_symb):
            reg = rb * n_symb + sym
            cce = reg_to_cce.get(reg)
            in_sib1 = cce in sib1_cces
            ax.add_patch(Rectangle((rb, sym), 1, 1,
                                   facecolor=C_PDCCH if in_sib1 else C_REG,
                                   edgecolor=C_REG_EDGE, lw=0.5,
                                   zorder=2 if in_sib1 else 1))
            if cce is not None and n_rb <= 60:
                ax.text(rb + 0.5, sym + 0.5, str(cce),
                        ha="center", va="center", fontsize=6,
                        color="white" if in_sib1 else C_TEXT, zorder=3)

    ax.set_xlim(0, n_rb)
    ax.set_ylim(0, n_symb)
    ax.set_xlabel(f"RB within CORESET#0  (CRB {crb0}..{crb0 + n_rb - 1} from Point A)")
    ax.set_ylabel("symbol")
    ax.set_yticks([s + 0.5 for s in range(n_symb)])
    ax.set_yticklabels([str(c["coreset0_first_symbol"] + s) for s in range(n_symb)])
    ax.set_xticks(range(0, n_rb + 1, 6))

    interl_txt = (f"interleaved (L={reg_bundle}, R={interleaver}, shift={shift})"
                  if interleaved else f"non-interleaved (L={reg_bundle})")
    if agg:
        pdcch_txt = (f"SIB1 PDCCH: aggregation level {agg} CCE"
                     f"{'s' if agg > 1 else ''}, first CCE {ncce} "
                     f"(CCEs {sorted(sib1_cces)})")
    else:
        pdcch_txt = "SIB1 PDCCH: not captured (no SI-RNTI DCI location)"
    fig.suptitle(
        f"CORESET#0 detail  -  PCI {c['pci']}  -  idx {c['coreset0_idx']}  -  "
        f"{n_rb} RB x {n_symb} sym, {n_cce} CCE, {interl_txt}\n{pdcch_txt}",
        fontsize=9.5)
    fig.subplots_adjust(top=0.80, bottom=0.20, left=0.06, right=0.99)
    return fig


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cell_json", help="cell_*.json written by the recorder")
    ap.add_argument("-o", "--out", help="PNG path (default: alongside the JSON, _pdcch)")
    ap.add_argument("--show", action="store_true",
                    help="open an interactive window (needs a GUI backend)")
    args = ap.parse_args(argv)

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
        out = args.out or os.path.splitext(args.cell_json)[0] + "_pdcch.png"
        fig.savefig(out, dpi=150)
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
