#!/usr/bin/env python3
"""Find the 5G NR cells (SSB frequency and PCI) around: for the field test on commercial
networks, where the SSB frequencies are not known in advance.

For every centre frequency: capture --ms of IQ with UHD's rx_samples_to_file (X410, DPDK),
correlate the 3 PSS at every sync-raster frequency (GSCN) of the captured band, and for
the strong ones read the SSS two symbols later (equalized with the PSS) to get N_ID1,
hence the PCI. The radio's reference error (--ppm, the X410's is about -4) is corrected,
with the sign of NR-Scope's freq_offset. Prints one line per cell found, strongest first,
and the NR-Scope radio settings to use for it.

usage: sudo field_scan.py [--centers 3350e6,3450e6,3550e6,3650e6,3750e6] [--rate 122.88e6]
         [--scs 30] [--ms 40] [--gain 35] [--ppm -4.02] [--args UHD_ARGS] [--ant TX/RX0]
         [--file capture.sc16 --fc HZ]   (analyse an existing capture instead)
"""
import argparse
import os
import subprocess
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--centers", default="3350e6,3450e6,3550e6,3650e6,3750e6")
ap.add_argument("--rate", type=float, default=122.88e6)
ap.add_argument("--scs", type=float, default=30, help="SSB subcarrier spacing, kHz")
ap.add_argument("--ms", type=float, default=40)
ap.add_argument("--gain", type=float, default=35)
ap.add_argument("--ppm", type=float, default=-4.02, help="radio reference error (NR-Scope freq_offset_ppm)")
ap.add_argument("--args", default="type=x4xx,master_clock_rate=245760000,addr=192.168.20.2,mgmt_addr=192.168.1.190,"
                                  "use_dpdk=1234")
ap.add_argument("--ant", default="TX/RX0")
ap.add_argument("--rx", default="/usr/local/lib/uhd/examples/rx_samples_to_file")
ap.add_argument("--out", default="/tmp/field_scan")
ap.add_argument("--file", help="analyse this capture (sc16 or fc32) instead of capturing")
ap.add_argument("--fc", type=float, help="centre frequency of --file")
ap.add_argument("--min-ratio", type=float, default=8.0, help="PSS peak/median to report")
ap.add_argument("--all", action="store_true", help="also list PSS peaks whose SSS does not fit (side lobes)")
a = ap.parse_args()
fs, scs = a.rate, a.scs * 1e3
nfft = int(round(fs / scs))
cp = int(round(nfft * 144 / 2048))  # normal CP (the long one only at symbols 0/7 of a half subframe)


def pss_seq(nid2):
    x = [0, 1, 1, 0, 1, 1, 1]
    for i in range(127 - 7):
        x.append((x[i + 4] + x[i]) % 2)
    return np.array([1 - 2 * x[(m + 43 * nid2) % 127] for m in range(127)], dtype=np.float32)


def sss_seq(nid1, nid2):
    x0 = [1, 0, 0, 0, 0, 0, 0]
    x1 = [1, 0, 0, 0, 0, 0, 0]
    for i in range(127 - 7):
        x0.append((x0[i + 4] + x0[i]) % 2)
        x1.append((x1[i + 1] + x1[i]) % 2)
    m0 = 15 * (nid1 // 112) + 5 * nid2
    m1 = nid1 % 112
    return np.array([(1 - 2 * x0[(n + m0) % 127]) * (1 - 2 * x1[(n + m1) % 127]) for n in range(127)], np.float32)


k127 = np.arange(127) - 64
pss_t = []
for nid2 in range(3):
    X = np.zeros(nfft, complex)
    X[k127 % nfft] = pss_seq(nid2)
    pss_t.append(np.fft.ifft(X).astype(np.complex64))
SSS = {(n1, n2): sss_seq(n1, n2) for n2 in range(3) for n1 in range(336)}


def raster(fc):
    half = fs / 2 * 0.8
    if fc < 3e9:
        c = [N * 1.2e6 + M * 50e3 for N in range(int((fc - half) // 1.2e6), int((fc + half) // 1.2e6) + 2)
             for M in (1, 3, 5)]
    else:
        c = [3e9 + N * 1.44e6 for N in range(int((fc - half - 3e9) // 1.44e6), int((fc + half - 3e9) // 1.44e6) + 2)]
    return [f for f in c if abs(f - fc) < half - 127 * scs / 2]


def analyse(path, fc):
    """one FFT of the whole capture; per raster candidate, the +-3.84 MHz around it (the
    PSS/SSS are 127 subcarriers) taken from that spectrum and brought back to the time
    domain at 7.68 Msps (30 kHz: 256-point symbols), where the PSS is correlated and the
    SSS read"""
    n = int(fs * a.ms / 1e3)
    raw = (np.fromfile(path, dtype=np.int16, count=2 * n).astype(np.float32) if not path.endswith(".fc32") else
           np.fromfile(path, dtype=np.float32, count=2 * n))
    x = (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)
    n = len(x)
    if n < nfft * 8:
        return []
    peak_dbfs = 20 * np.log10(np.abs(x).max() / (32768 if not path.endswith(".fc32") else 1) + 1e-12)
    X = np.fft.fft(x)
    fd = 256 * scs  # decimated rate: 256 subcarriers
    W = int(round(n * fd / fs))  # bins kept per candidate
    nd = 256
    cpd = int(round(nd * 144 / 2048))
    tmpl = []
    for nid2 in range(3):
        T = np.zeros(nd, complex)
        T[k127 % nd] = pss_seq(nid2)
        tmpl.append(np.fft.ifft(T).astype(np.complex64))
    found = []
    for f in raster(fc):
        # a signal at f shows at f - fc + the radio's error; NR-Scope's freq_offset = ppm * f
        off = f - fc + a.ppm * 1e-6 * f
        b0 = int(round(off / fs * n))
        idx = (b0 + np.arange(-W // 2, W - W // 2)) % n
        Yd = np.fft.ifftshift(X[idx])  # DC (the candidate) at bin 0
        y = np.fft.ifft(Yd).astype(np.complex64) * (W / n)
        L = 1 << int(np.ceil(np.log2(len(y) + nd)))
        Yf = np.fft.fft(y, L)
        for nid2 in range(3):
            c = np.abs(np.fft.ifft(Yf * np.conj(np.fft.fft(tmpl[nid2], L))))[:len(y)]
            p = int(c.argmax())
            ratio = c[p] / (np.median(c) + 1e-9)
            if ratio < a.min_ratio:
                continue
            s0, s2 = p, p + 2 * (nd + cpd)  # PSS: SSB symbol 0, SSS: symbol 2
            if s2 + nd > len(y):
                continue
            P = np.fft.fft(y[s0:s0 + nd])[k127 % nd]
            pfit = np.abs(np.vdot(pss_seq(nid2), P)) / (np.abs(P).sum() + 1e-9)
            if pfit < 0.5:  # peak/median alone is fooled by TDD (the median falls in the UL gaps)
                continue
            S = np.fft.fft(y[s2:s2 + nd])[k127 % nd]
            H = P / pss_seq(nid2)
            Se = S * np.conj(H) / (np.abs(H) ** 2 + 1e-9)
            # magnitude: a residual CFO turns the SSS against the PSS two symbols earlier
            best = max(range(336), key=lambda n1: np.abs(np.vdot(SSS[(n1, nid2)], Se)))
            sc = np.abs(np.vdot(SSS[(best, nid2)], Se)) / (np.abs(Se).sum() + 1e-9)
            found.append((ratio, f, 3 * best + nid2, nid2, sc, peak_dbfs))
    found.sort(reverse=True)
    out, seen = [], set()
    for r in found:
        key = (round(r[1]), r[2])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


os.makedirs(a.out, exist_ok=True)
jobs = [(a.file, a.fc)] if a.file else []
if not a.file:
    for fc in [float(v) for v in a.centers.split(",")]:
        path = os.path.join(a.out, "scan_%.0fMHz.sc16" % (fc / 1e6))
        cmd = [a.rx, "--args", a.args, "--rate", str(fs), "--freq", str(fc), "--gain", str(a.gain), "--ant", a.ant,
               "--duration", str(a.ms / 1e3 + 0.02), "--type", "short", "--file", path, "--setup", "0.5"]
        print("capturing %.1f MHz ..." % (fc / 1e6), flush=True)
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired:
            print("  capture timed out")
            continue
        if r.returncode != 0 or not os.path.exists(path):
            print("  capture failed: %s" % (r.stderr.strip().splitlines()[-1:] or ["?"])[0])
            continue
        jobs.append((path, fc))
cells = []
for path, fc in jobs:
    for r in analyse(path, fc):
        cells.append(r + (fc,))
cells.sort(reverse=True)
if not a.all:
    cells = [c for c in cells if c[4] >= 0.5]
print("\n%-12s %-5s %-9s %-8s %-8s %s" % ("SSB MHz", "PCI", "PSS ratio", "SSS fit", "peak dBFS", "captured at"))
for ratio, f, pci, nid2, sc, pk, fc in cells:
    print("%-12.3f %-5d %-9.1f %-8.2f %-8.1f %.1f MHz%s" % (f / 1e6, pci, ratio, sc, pk, fc / 1e6,
                                                         "  (SSS weak: PCI unsure)" if sc < 0.5 else ""))
if cells:
    ratio, f, pci, nid2, sc, pk, fc = cells[0]
    print("\nNR-Scope radio settings for the strongest (rx_center_freq: refine with the carrier line "
          "NR-Scope prints after SIB1):")
    print("  ssb_freq: %d\n  rx_center_freq: %d\n  pci: %d\n  freq_offset_ppm: %.2f\n  scs_index: %d" %
          (round(f), round(f), pci, a.ppm, 0 if a.scs == 15 else 1))
