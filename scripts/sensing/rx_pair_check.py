#!/usr/bin/env python3
"""Compare the Rx chains of a raw multi-channel capture (rx_samples_to_file, sc16).

Takes two to four chain files, chain 0 first, and reports everything against chain 0.

Per chain: peak and rms level, clipping, carrier power, noise floor (quiet TDD symbols)
and SNR. Against chain 0: power ratio, delay, coherence across the carrier, and how
steady the relative phase stays over time.

The last block is the one the four-chain AoA needs. The X410's chains share a sample
clock, but each ZBX daughterboard settles its LO at an arbitrary phase on every retune,
and chains 2 and 3 are on the second board: until those offsets are known and constant,
the steering vector in nr_ue_aoa.h describes an array that does not exist. Feed one
source into every input through equal-length cables, run this, and if the phase is steady
(small spread, little drift) paste the printed line into NR_AOA_CAL_DEG. If it drifts
within a run, the calibration has to become a per-map estimate instead of a constant.

  rx_pair_check.py usrp_samples.0.dat usrp_samples.1.dat --rate 122.88e6 --bw 98.28e6
  rx_pair_check.py usrp_samples.{0,1,2,3}.dat --rate 122.88e6
"""
import argparse
import numpy as np

# Chains the sensing library handles, NR_SENSING_MAX_RX in nr_ue_sensing.h
MAX_RX = 4


def load(path, n_max):
    raw = np.fromfile(path, dtype=np.int16, count=2 * n_max if n_max else -1)
    return (raw[0::2] + 1j * raw[1::2]).astype(np.complex64) / 32768.0


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chans", nargs="+", help=f"chain files, chain 0 first, 2 to {MAX_RX} of them")
    ap.add_argument("--rate", type=float, default=122.88e6)
    ap.add_argument("--bw", type=float, default=98.28e6, help="carrier width centred on the tuned frequency")
    ap.add_argument("--block", type=int, default=4384, help="block length, about one OFDM symbol at 30 kHz SCS, 122.88 Msps")
    ap.add_argument("--nfft", type=int, default=4096)
    ap.add_argument("--max-samples", type=int, default=0)
    ap.add_argument("--bands", type=int, default=10, help="sub-bands the carrier is split into")
    a = ap.parse_args()

    if len(a.chans) < 2:
        ap.error("at least two chain files are needed to compare anything")
    if len(a.chans) > MAX_RX:
        ap.error(f"{len(a.chans)} chains given, the sensing pipeline handles {MAX_RX}")

    x = [load(p, a.max_samples) for p in a.chans]
    n = min(len(v) for v in x)
    x = [v[:n] for v in x]
    nc = len(x)
    print(f"{nc} chains, {n} samples each = {n / a.rate * 1e3:.1f} ms")

    for i, v in enumerate(x):
        pk = np.max(np.abs(v))
        clip = np.mean((np.abs(v.real) > 0.99) | (np.abs(v.imag) > 0.99))
        print(f"rx{i}: peak {db(pk ** 2):.1f} dBFS, rms {db(np.mean(np.abs(v) ** 2)):.1f} dBFS, clipped {clip * 100:.3f} %")

    # Block powers. Quiet TDD symbols (UL, guard) give the noise floor, busy ones the DL.
    # Busy and quiet are decided on the sum over the chains, so every chain is measured on
    # the same blocks and the ratios between them stay meaningful.
    nb = n // a.block
    blk = [v[: nb * a.block].reshape(nb, a.block) for v in x]
    pw = [np.mean(np.abs(b) ** 2, axis=1) for b in blk]
    ref = sum(pw)
    quiet = ref <= np.percentile(ref, 5)
    busy = ref >= np.percentile(ref, 80)
    print("\nper chain, from block powers (busy = top 20 % of blocks, quiet = bottom 5 %)")
    for i in range(nc):
        pb, pq = np.mean(pw[i][busy]), np.mean(pw[i][quiet])
        rel = db(pb / np.mean(pw[0][busy]))
        print(f"rx{i}: busy {db(pb):6.1f} dBFS  quiet {db(pq):6.1f} dBFS  SNR {db(pb / pq - 1):5.1f} dB"
              f"  busy vs rx0 {rel:+5.1f} dB")

    # Spectra over busy blocks, by sub-band.
    f = np.fft.fftshift(np.fft.fftfreq(a.nfft, 1 / a.rate))
    seg = [k * a.block for k in np.flatnonzero(busy) if k * a.block + a.nfft <= n]
    if not seg:
        print("\nno busy block long enough for one FFT; nothing further to report")
        return
    w = np.hanning(a.nfft).astype(np.float32)
    S = [np.fft.fftshift(np.fft.fft(np.stack([v[s:s + a.nfft] * w for s in seg]), axis=1), axes=1) for v in x]
    P = [np.mean(np.abs(s) ** 2, axis=0) for s in S]
    # cross-spectrum of every chain against chain 0, and the coherence from it
    C = [np.mean(S[0] * np.conj(S[i]), axis=0) for i in range(nc)]
    coh = [np.abs(C[i]) / np.sqrt(P[0] * P[i]) for i in range(nc)]
    in_band = np.abs(f) < a.bw / 2 * 0.98
    edges = np.linspace(-a.bw / 2, a.bw / 2, a.bands + 1)

    head = "\nsub-band     rx0 dB" + "".join(f"   rx{i}-rx0  coh{i}" for i in range(1, nc))
    print(head)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (f >= lo) & (f < hi) & in_band
        row = f"{lo / 1e6:+6.1f} MHz  {db(np.mean(P[0][m])):7.1f}"
        for i in range(1, nc):
            row += f"  {db(np.mean(P[i][m]) / np.mean(P[0][m])):+8.1f}  {np.mean(coh[i][m]):.3f}"
        print(row)
    print("whole carrier: coherence " + "  ".join(f"rx{i} {np.mean(coh[i][in_band]):.3f}" for i in range(1, nc)))

    # Delay of each chain against chain 0: slope of the cross-spectrum phase over the carrier.
    print()
    delay = [0.0] * nc
    for i in range(1, nc):
        ph = np.unwrap(np.angle(C[i][in_band]))
        delay[i] = np.polyfit(f[in_band], ph, 1)[0]
        print(f"delay rx{i} against rx0 (positive = rx{i} later): {delay[i] / (2 * np.pi) * 1e9:+.2f} ns  "
              f"({delay[i] / (2 * np.pi) * a.rate:+.3f} samples)")

    # Relative phase and power ratio over time, busy blocks only, each chain's delay removed.
    t = np.array(seg) / a.rate * 1e3
    e = [np.sum(np.abs(S[i][:, in_band]) ** 2, axis=1) for i in range(nc)]
    print(f"\nover {t[-1] - t[0]:.0f} ms of busy blocks:")
    cal_deg = [0.0] * nc
    for i in range(1, nc):
        rot = np.exp(-1j * delay[i] * f[in_band])
        c = np.sum(S[0][:, in_band] * np.conj(S[i][:, in_band]) * rot, axis=1)
        phase = np.unwrap(np.angle(c))
        fit = np.polyfit(t, phase, 1)
        ratio = db(e[i] / e[0])
        # The offset rx{i} carries in excess of rx0, which is what NR_AOA_CAL_DEG removes:
        # nr_ue_aoa_cell_vector() multiplies chain a by exp(-j*cal[a]). The cross-spectrum
        # is S0 conj(Si), so its angle is phi0 - phi_i and the offset is the negative of it.
        cal_deg[i] = -np.degrees(np.angle(np.mean(c / np.abs(c))))
        rho = np.abs(np.mean(c)) / np.sqrt(np.mean(e[0]) * np.mean(e[i]))
        print(f"rx{i}: offset vs rx0 {cal_deg[i]:+7.1f} deg, spread (std) "
              f"{np.degrees(np.std(phase - np.polyval(fit, t))):5.1f} deg, "
              f"drift {np.degrees(fit[0]) * 1e3:+8.2f} deg/s, "
              f"power {np.median(ratio):+5.1f} dB "
              f"({np.percentile(ratio, 5):+.1f} .. {np.percentile(ratio, 95):+.1f}), |rho| {rho:.3f}")

    print("\nper-chain phase offsets for nr_ue_aoa.h, if the spread and drift above are small:")
    cal = [f"{v:.1f}" for v in cal_deg] + ["0.0"] * (MAX_RX - nc)
    tail = "" if nc == MAX_RX else f"  // rx{nc} and up were not captured, left at 0"
    print("#define NR_AOA_CAL_DEG {" + ", ".join(cal) + "}" + tail)


if __name__ == "__main__":
    main()
