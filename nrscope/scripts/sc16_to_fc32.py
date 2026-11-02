#!/usr/bin/env python3
"""Convert an interleaved int16 IQ capture (.sc16, as in the GoldenSniffer dataset) to
the interleaved float32 (.fc32) that srsRAN's file radio reads (device_name: file).
usage: sc16_to_fc32.py in.sc16 out.fc32 [scale=1/32768]"""
import sys
import numpy as np

src, dst = sys.argv[1], sys.argv[2]
scale = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0 / 32768.0
chunk = 1 << 22
with open(src, "rb") as fi, open(dst, "wb") as fo:
    while True:
        a = np.fromfile(fi, dtype=np.int16, count=chunk)
        if a.size == 0:
            break
        fo.write((a.astype(np.float32) * scale).tobytes())  # not tofile(): it needs a seekable file (FIFO)
