"""OFDM demodulation of the raw X410 IQ (sc16, 122.88 Msps, 30 kHz SCS, FFT 4096), independent of NR-Scope.
Slot = 61440 samples: symbol 0 has a 352-sample CP, symbols 1-13 288. The carrier's subcarrier k sits at FFT bin
k - K0 + bin0 (K0 = half the carrier's subcarriers, bin0 found by matching)."""
import numpy as np

FS, NFFT, SLOT = 122.88e6, 4096, 61440
CP = [352] + [288] * 13
START = np.cumsum([0] + [c + NFFT for c in CP[:-1]])   # symbol starts within the slot


def read(path, s0, n):
    """n complex samples from sample s0 of an sc16 file"""
    x = np.fromfile(path, np.int16, 2 * n, offset=4 * int(s0)).astype(np.float32)
    return x[0::2] + 1j * x[1::2]


def symbol(x, start, cp, backoff=0):
    """FFT of one symbol whose CP starts at `start` (window = CP end - backoff), DC-centred bins"""
    w = start + cp - backoff
    return np.fft.fftshift(np.fft.fft(x[w:w + NFFT])) * np.exp(-2j * np.pi * (np.arange(NFFT) - NFFT // 2) * backoff / NFFT)
