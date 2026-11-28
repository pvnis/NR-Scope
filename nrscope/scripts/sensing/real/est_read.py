"""Reader for so_light's channel-estimate files (NRSEST1): per slot and chain, the symbols in the mask, each with its
valid subcarriers (packed bits) and complex values.

  read(path) -> dict(nsc, slots: {slot_abs: {symbol: X (4, nsc) complex64, NaN where not measured}})
  snapshots(path) -> t_s (n,) chunk time of each symbol, X (n, 4, nsc)"""
import struct
import numpy as np

HEAD = struct.Struct("<8sIIQdIIII")
REC = struct.Struct("<IHHBBBBHHIIQ")


def read(path):
    b = open(path, "rb").read(); h = HEAD.unpack_from(b, 0); nsc = h[7]; o = HEAD.size; nb = (nsc + 7) // 8; slots = {}
    while o < len(b):
        r = REC.unpack_from(b, o); o += REC.size; c, mask, a = r[1], r[7], r[11]
        for l in range(14):
            if not mask >> l & 1:
                continue
            n = struct.unpack_from("<I", b, o)[0]; o += 4
            valid = np.unpackbits(np.frombuffer(b, np.uint8, nb, o), bitorder="little")[:nsc].astype(bool); o += nb
            v = np.frombuffer(b, np.complex64, n, o); o += 8 * n
            X = slots.setdefault(a, {}).setdefault(l, np.full((4, nsc), np.nan + 0j, np.complex64)); X[c, valid] = v
    return dict(nsc=nsc, slots=slots)


def snapshots(path):
    d = read(path); T, X = [], []
    for a in sorted(d["slots"]):
        for l in sorted(d["slots"][a]):
            T.append(a * 5e-4 + l / 28e3); X.append(d["slots"][a][l])
    return np.array(T), np.array(X)
