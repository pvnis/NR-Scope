#!/usr/bin/env python3
"""A time slice of an interleaved int16 IQ capture, as float32 for srsRAN file radio (FIFO-safe).
usage: sc16_slice.py in.sc16 out.fc32 <start_s> <dur_s> [rate=122.88e6]"""
import sys, numpy as np
src, dst, t0, dur = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
rate = float(sys.argv[5]) if len(sys.argv) > 5 else 122.88e6
n0, n = int(round(t0 * rate)), int(round(dur * rate)); chunk = 1 << 21   # complex samples
with open(src, "rb") as fi, open(dst, "wb") as fo:
    fi.seek(n0 * 4)
    while n > 0:
        a = np.fromfile(fi, dtype=np.int16, count=2 * min(chunk, n))
        if a.size == 0: break
        fo.write((a.astype(np.float32) / 32768.0).tobytes()); n -= a.size // 2
