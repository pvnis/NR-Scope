"""TRS / CSI-RS helpers: GSC1 records (raw REs of one symbol over the carrier), the 3GPP Gold sequence,
CSI-RS sequence init (TS 38.211 7.4.1.5.2), TRS (density 3, REs c, c+4, c+8) and 8-port CSI-RS (FD-CDM2 pairs on
subcarriers 0-7, density 1) despreading."""
import struct, os, re, numpy as np
PCI = int(os.environ.get("PCI", "632"))   # set_pci() / env PCI: the cell whose TRS and CSI-RS are despread

def set_pci(chunk_dir=None):
    """the PCI the chunk was decoded with (replay_blind_ch0/config.yaml), unless env PCI is set"""
    global PCI
    if "PCI" not in os.environ and chunk_dir:
        for d in ("replay_blind_ch0", "replay_ch0"):
            p = os.path.join(chunk_dir, d, "config.yaml")
            if os.path.exists(p):
                m = re.search(r"^\s*pci:\s*(\d+)", open(p).read(), re.M)
                if m: PCI = int(m.group(1)); break
    _cache.clear(); return PCI

def gold(cinit, n):
    Nc = 1600; L = Nc + n
    x1 = np.zeros(L + 31, np.uint8); x2 = np.zeros(L + 31, np.uint8); x1[0] = 1
    x2[:31] = [(cinit >> i) & 1 for i in range(31)]
    for k in range(L):
        x1[k + 31] = (x1[k + 3] + x1[k]) & 1
        x2[k + 31] = (x2[k + 3] + x2[k + 2] + x2[k + 1] + x2[k]) & 1
    return (x1[Nc:Nc + n] + x2[Nc:Nc + n]) & 1

def qpsk_seq(cinit, m):
    c = gold(cinit, 2 * m).astype(float)
    return ((1 - 2 * c[0::2]) + 1j * (1 - 2 * c[1::2])) / np.sqrt(2)

def csi_cinit(slot, l, nid=None):
    nid = PCI if nid is None else nid
    return ((1 << 10) * (14 * slot + l + 1) * (2 * nid + 1) + nid) % (1 << 31)

_cache = {}
def seq(slot, l, m):
    k = (slot, l, m)
    if k not in _cache: _cache[k] = qpsk_seq(csi_cinit(slot, l), m)
    return _cache[k]

def load_gsc(path):
    """-> dict (a, l) -> (first_crb, raw REs 12*nprb complex); a = sfn*20 + slot unwrapped like gs_all.load"""
    raw = open(path, "rb").read(); i = 0; out = {}; base = None
    while i + 12 <= len(raw):
        tag = raw[i:i + 4]
        if tag in (b"GSD1", b"GSD2", b"GSD3"):
            if i + 20 > len(raw): break      # skipped here; GSD3 = DM-RS type 2, 4 values per PRB
            sfn, sl, ns, rn, npb, occ, coh = struct.unpack_from("<HBBHHff", raw, i + 4); j = i + 20; per = {b"GSD1": 3, b"GSD2": 6, b"GSD3": 4}[tag] * npb * 8
            if j + npb + ns * (1 + per) > len(raw): break
            a = sfn * 20 + sl
            if base is None: base = a
            i = j + npb + ns * (1 + per)
        elif tag == b"GSE1":           # per-PRB power (energy check), not used here
            i += 10 + 4 * struct.unpack_from("<H", raw, i + 8)[0]
        elif tag == b"GSC1":
            sfn, sl, l, fc, npb = struct.unpack_from("<HBBHH", raw, i + 4); j = i + 12
            if j + 12 * npb * 8 > len(raw): break
            a = sfn * 20 + sl
            if base is None: base = a
            a = a + 20480 if a < base - 1000 else a
            out[(a, l)] = (fc, np.frombuffer(raw, np.complex64, 12 * npb, j).astype(complex))
            i = j + 12 * npb * 8
        else: i += 1
    return out

def trs_despread(fc, z, slot, l, c):
    """TRS on REs c, c+4, c+8 of each PRB, m = 3n + j -> (nprb, 3) channel"""
    npb = len(z) // 12; r = seq(slot, l, 3 * (fc + npb) + 4)
    n = fc + np.arange(npb)
    out = np.empty((npb, 3), complex)
    for j in range(3): out[:, j] = z[12 * np.arange(npb) + c + 4 * j] * np.conj(r[3 * n + j]) * 2
    return out

def csi8_despread(fc, z, slot, l=13):
    """8-port CSI-RS: subcarrier k in 0..7 uses m = 2n + (k & 1); ports 2i, 2i+1 = (z_2i +- z_2i+1) / 2 -> (nprb, 8)"""
    npb = len(z) // 12; r = seq(slot, l, 2 * (fc + npb) + 4); n = fc + np.arange(npb)
    zz = np.empty((npb, 8), complex)
    for k in range(8): zz[:, k] = z[12 * np.arange(npb) + k] * np.conj(r[2 * n + (k & 1)]) * 2
    P = np.empty((npb, 8), complex)
    for i in range(4):
        P[:, 2 * i] = (zz[:, 2 * i] + zz[:, 2 * i + 1]) / 2; P[:, 2 * i + 1] = (zz[:, 2 * i] - zz[:, 2 * i + 1]) / 2
    return P

def coh(v):
    """phase smoothness of a sequence along frequency: |sum z_k z*_{k-1}| / sum |z|^2"""
    v = np.ravel(v); return float(np.abs(np.sum(v[1:] * np.conj(v[:-1]))) / (np.sum(np.abs(v) ** 2) + 1e-30))


def trs_comb(G, n=20):
    """the TRS subcarrier offset c (REs c, c+4, c+8 of each PRB) of this cell: the one whose despread TRS is coherent
    across subcarriers, over up to n TRS symbols of one antenna's load_gsc() dict -> (c, coherence)"""
    keys = sorted(k for k in G if k[1] in (4, 8))[:n]
    if not keys: return 2, 0.0
    best = (2, -1.0)
    for c in range(4):
        v = []
        for a, l in keys:
            h = trs_despread(*G[(a, l)], a % 20, l, c).ravel()
            v.append(abs(np.sum(h[1:] * np.conj(h[:-1]))) / (np.sum(abs(h) ** 2) + 1e-30))
        m = float(np.median(v))
        if m > best[1]: best = (c, m)
    return best
