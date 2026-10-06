#!/usr/bin/env python3
"""Check the time-domain slots the SIB1 decoder was handed (NRSCOPE_DUMP_SIB_SLOTS).

A SIB1 that never decodes looks the same from outside as a cell not transmitting one:
the PDCCH search reports energy with no correlation either way. This answers which half
of the pipeline is at fault, from the samples that actually reached the FFT.

Per slot it reports:
  cp_corr   cyclic-prefix correlation at the assumed symbol boundaries. A 30 kHz slot at
            122.88 Msps is 1x(4096+352) + 13x(4096+288) samples; a correctly aligned
            OFDM slot correlates its CP against the symbol tail at ~1.0. Low means the
            samples are not aligned OFDM at all.
  best_off  the sample offset whose CP correlation is highest, searched +-64. Non-zero
            means the slot handed over starts in the wrong place, which no amount of
            frequency mapping downstream can fix.
  ssb_bin   where the strongest narrow band sits, in subcarriers from the carrier
            centre. The SSB should be at (ssb_freq - rx_center) / scs: -1368 for the
            Benetel cell (3408.96 against 3450 MHz at 30 kHz).

  sib_slot_check.py /tmp/sib_slots.bin
  sib_slot_check.py /tmp/sib_slots_4ant.bin --expect-ssb-bin -1368
"""
import argparse
import sys

import numpy as np

NFFT = 4096
CP_LONG = 352   # symbol 0 of each half subframe
CP_NORM = 288


def symbol_starts(n_symbols=14):
    """Sample offset and CP length of each symbol of a 30 kHz slot at 122.88 Msps."""
    out, pos = [], 0
    for l in range(n_symbols):
        cp = CP_LONG if l == 0 else CP_NORM
        out.append((pos, cp))
        pos += cp + NFFT
    return out


def cp_corr(x, offset=0):
    """Mean |CP correlation| over the slot's symbols, at a sample offset."""
    vals = []
    for pos, cp in symbol_starts():
        a0 = pos + offset
        if a0 < 0 or a0 + cp + NFFT > len(x):
            continue
        a = x[a0:a0 + cp]
        b = x[a0 + NFFT:a0 + NFFT + cp]
        d = np.sqrt(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2))
        if d > 0:
            vals.append(abs(np.vdot(b, a)) / d)
    return float(np.mean(vals)) if vals else 0.0


def read_records(path):
    with open(path, "rb") as f:
        raw = f.read()
    off, recs = 0, []
    while off + 12 <= len(raw):
        sfn, slot_idx, n = np.frombuffer(raw, dtype=np.uint32, count=3, offset=off)
        off += 12
        need = int(n) * 8  # cf_t = 2 x float32
        if off + need > len(raw):
            break
        x = np.frombuffer(raw, dtype=np.float32, count=int(n) * 2, offset=off)
        off += need
        recs.append((int(sfn), int(slot_idx), (x[0::2] + 1j * x[1::2]).astype(np.complex64)))
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--expect-ssb-bin", type=int, default=-1368,
                    help="subcarriers from centre where the SSB should be")
    ap.add_argument("--search", type=int, default=64, help="sample offsets searched either side")
    a = ap.parse_args()

    recs = read_records(a.dump)
    if not recs:
        print(f"no records in {a.dump}")
        return 1
    print(f"{len(recs)} slot(s), {len(recs[0][2])} samples each\n")
    print(" sfn slot    rms_dBFS  cp_corr  best_off  cp@best  ssb_bin  ssb_dB_over_mean")

    aligned = 0
    for sfn, slot_idx, x in recs:
        rms = 10 * np.log10(max(np.mean(np.abs(x) ** 2), 1e-30))
        c0 = cp_corr(x)
        offs = range(-a.search, a.search + 1)
        best = max(offs, key=lambda o: cp_corr(x, o))
        cbest = cp_corr(x, best)

        # Average spectrum over the slot's symbols, at the best alignment
        acc = np.zeros(NFFT)
        n = 0
        for pos, cp in symbol_starts():
            s = pos + best + cp
            if s + NFFT <= len(x):
                acc += np.abs(np.fft.fftshift(np.fft.fft(x[s:s + NFFT]))) ** 2
                n += 1
        acc /= max(n, 1)
        # strongest 240-subcarrier band (the SSB is 20 RB wide)
        k = np.convolve(acc, np.ones(240) / 240, mode="same")
        peak = int(np.argmax(k)) - NFFT // 2
        ssb_db = 10 * np.log10(max(k.max(), 1e-30) / max(acc.mean(), 1e-30))

        ok = cbest > 0.7 and best == 0
        aligned += ok
        print(f"{sfn:4d} {slot_idx:4d}  {rms:9.1f}  {c0:7.3f}  {best:+8d}  {cbest:7.3f}  "
              f"{peak:+7d}  {ssb_db:14.1f}{'' if ok else '   <--'}")

    print(f"\n{aligned}/{len(recs)} slot(s) are aligned OFDM at offset 0")
    print(f"Expected SSB bin: {a.expect_ssb_bin:+d}")
    print("\nReading: cp_corr near 1.0 at best_off 0 means the samples handed to the FFT are")
    print("properly aligned OFDM, so the fault is downstream (mapping or descrambling).")
    print("A non-zero best_off, or a low cp_corr everywhere, means the slot itself is wrong,")
    print("so the fault is upstream in the chain-0 path through the slot queue.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
