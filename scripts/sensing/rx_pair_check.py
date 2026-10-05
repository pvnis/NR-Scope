#!/usr/bin/env python3
"""Compare two Rx chains from a raw two-channel capture (rx_samples_to_file, sc16).

Reports, per chain: carrier power, noise floor (quiet TDD symbols) and SNR, and
the power across the carrier. Between chains: power ratio, delay, coherence and
how steady the relative phase stays over time.

  rx_pair_check.py usrp_samples.0.dat usrp_samples.1.dat --rate 122.88e6 --bw 98.28e6
"""
import argparse
import numpy as np


def load(path, n_max):
    raw = np.fromfile(path, dtype=np.int16, count=2 * n_max if n_max else -1)
    return (raw[0::2] + 1j * raw[1::2]).astype(np.complex64) / 32768.0


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ch0")
    ap.add_argument("ch1")
    ap.add_argument("--rate", type=float, default=122.88e6)
    ap.add_argument("--bw", type=float, default=98.28e6, help="carrier width centred on the tuned frequency")
    ap.add_argument("--block", type=int, default=4384, help="block length, about one OFDM symbol at 30 kHz SCS, 122.88 Msps")
    ap.add_argument("--nfft", type=int, default=4096)
    ap.add_argument("--max-samples", type=int, default=0)
    a = ap.parse_args()

    x = [load(a.ch0, a.max_samples), load(a.ch1, a.max_samples)]
    n = min(len(x[0]), len(x[1]))
    x = [v[:n] for v in x]
    print(f"{n} samples per chain = {n / a.rate * 1e3:.1f} ms")

    for i, v in enumerate(x):
        pk = np.max(np.abs(v))
        clip = np.mean((np.abs(v.real) > 0.99) | (np.abs(v.imag) > 0.99))
        print(f"rx{i}: peak {db(pk ** 2):.1f} dBFS, rms {db(np.mean(np.abs(v) ** 2)):.1f} dBFS, clipped {clip * 100:.3f} %")

    # Block powers. Quiet TDD symbols (UL, guard) give the noise floor, busy ones the DL.
    nb = n // a.block
    blk = [v[: nb * a.block].reshape(nb, a.block) for v in x]
    pw = [np.mean(np.abs(b) ** 2, axis=1) for b in blk]
    ref = pw[0] + pw[1]
    quiet = ref <= np.percentile(ref, 5)
    busy = ref >= np.percentile(ref, 80)
    print("\nper chain, from block powers (busy = top 20 % of blocks, quiet = bottom 5 %)")
    for i in range(2):
        pb, pq = np.mean(pw[i][busy]), np.mean(pw[i][quiet])
        print(f"rx{i}: busy {db(pb):6.1f} dBFS  quiet {db(pq):6.1f} dBFS  SNR {db(pb / pq - 1):5.1f} dB")
    print(f"busy-block power rx1 vs rx0: {db(np.mean(pw[1][busy]) / np.mean(pw[0][busy])):+.1f} dB, "
          f"quiet-block: {db(np.mean(pw[1][quiet]) / np.mean(pw[0][quiet])):+.1f} dB")

    # Spectra over busy blocks, by sub-band.
    f = np.fft.fftshift(np.fft.fftfreq(a.nfft, 1 / a.rate))
    idx = np.flatnonzero(busy)
    seg = []
    for k in idx:
        s = k * a.block
        if s + a.nfft <= n:
            seg.append(s)
    w = np.hanning(a.nfft).astype(np.float32)
    S = [np.fft.fftshift(np.fft.fft(np.stack([v[s:s + a.nfft] * w for s in seg]), axis=1), axes=1) for v in x]
    P = [np.mean(np.abs(s) ** 2, axis=0) for s in S]
    C = np.mean(S[0] * np.conj(S[1]), axis=0)
    coh = np.abs(C) / np.sqrt(P[0] * P[1])
    in_band = np.abs(f) < a.bw / 2 * 0.98
    edges = np.linspace(-a.bw / 2, a.bw / 2, 11)
    print("\nsub-band     rx0 dB   rx1 dB   rx1-rx0   coherence")
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (f >= lo) & (f < hi) & in_band
        print(f"{lo / 1e6:+6.1f} MHz  {db(np.mean(P[0][m])):7.1f}  {db(np.mean(P[1][m])):7.1f}  "
              f"{db(np.mean(P[1][m]) / np.mean(P[0][m])):+7.1f}   {np.mean(coh[m]):.3f}")
    print(f"whole carrier: coherence {np.mean(coh[in_band]):.3f}")

    # Delay between chains: slope of the cross-spectrum phase over the carrier.
    ph = np.unwrap(np.angle(C[in_band]))
    slope = np.polyfit(f[in_band], ph, 1)[0]
    print(f"delay rx1 against rx0 (positive = rx1 later): {slope / (2 * np.pi) * 1e9:+.2f} ns  "
          f"({slope / (2 * np.pi) * a.rate:+.3f} samples)")

    # Relative phase and power ratio over time, busy blocks only, delay removed.
    rot = np.exp(-1j * slope * f[in_band])
    c = np.sum(S[0][:, in_band] * np.conj(S[1][:, in_band]) * rot, axis=1)
    e0 = np.sum(np.abs(S[0][:, in_band]) ** 2, axis=1)
    e1 = np.sum(np.abs(S[1][:, in_band]) ** 2, axis=1)
    phase = np.unwrap(np.angle(c))
    ratio = db(e1 / e0)
    t = np.array(seg) / a.rate * 1e3
    print(f"\nover {t[-1] - t[0]:.0f} ms of busy blocks:")
    print(f"relative phase: mean {np.degrees(np.angle(np.mean(c / np.abs(c)))):+.1f} deg, "
          f"spread (std) {np.degrees(np.std(phase - np.polyval(np.polyfit(t, phase, 1), t))):.1f} deg, "
          f"drift {np.degrees(np.polyfit(t, phase, 1)[0]) * 1e3:+.2f} deg/s")
    print(f"power ratio rx1-rx0: median {np.median(ratio):+.1f} dB, 5-95 % {np.percentile(ratio, 5):+.1f} .. {np.percentile(ratio, 95):+.1f} dB")
    corr = np.abs(np.mean(c)) / np.sqrt(np.mean(e0) * np.mean(e1))
    print(f"wideband correlation |rho| {corr:.3f}")


if __name__ == "__main__":
    main()
