"""Reader for NR-Scope's sensing.record_estimates file (nrscope_sensing_record.h)."""
import struct

import numpy as np


def records(path, max_records=None):
    """Yield dict(sfn, slot, slot_abs, aarx, ports, layer, dmrs_mask, symbols=[(l, k, H)])."""
    with open(path, "rb") as f:
        hdr = f.read(48)
        magic, ver, nant, fc, fs, ofdm, scs, nsc, _ = struct.unpack("<8sIIQdIIII", hdr)
        nbytes = (nsc + 7) // 8
        n = 0
        while True:
            rh = f.read(32)
            if len(rh) < 32:
                return
            m, aarx, ports, layer, nl, tdd, _, mask, _, sfn, slot, slot_abs = struct.unpack("<IHHBBBBHHIIQ", rh)
            if m != 0x4E525352:
                raise ValueError("bad record magic")
            syms = []
            for l in range(14):
                if not mask & (1 << l):
                    continue
                nv = struct.unpack("<I", f.read(4))[0]
                bits = np.unpackbits(np.frombuffer(f.read(nbytes), np.uint8), bitorder="little")[:nsc]
                k = np.nonzero(bits)[0]
                H = np.frombuffer(f.read(8 * nv), np.complex64)
                syms.append((l, k, H))
            yield dict(sfn=sfn, slot=slot, slot_abs=slot_abs, aarx=aarx, ports=ports, layer=layer,
                       dmrs_mask=mask, symbols=syms, n_sc=nsc, scs=scs, ofdm=ofdm)
            n += 1
            if max_records and n >= max_records:
                return
