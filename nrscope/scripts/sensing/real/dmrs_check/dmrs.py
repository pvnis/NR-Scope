"""PDSCH DM-RS of one symbol, types 1 and 2, despread with nrseq (TS 38.211 7.4.1.1): comb delta (CDM group 0: even
subcarriers, 1: odd), k = 4n + 2k' + delta, r index 2n + k' with n counted from CRB 0 (carrier starts at CRB 0 here).
ports(sym, slot, l, n_id, n_scid, crbs) -> (len(crbs), 3, 4) per PRB, pair n, port 1000..1003 (FD-OCC +- per comb)"""
import numpy as np
import nrseq


def despread(sym, slot, l, n_id, n_scid, crbs, delta):
    r = nrseq.qpsk(nrseq.dmrs_cinit(slot, l, n_id, n_scid), 6 * (max(crbs) + 1))
    n = 3 * np.asarray(crbs)[:, None] + np.arange(3)[None, :]               # (P, 3)
    a = sym[4 * n + delta] * np.conj(r[2 * n]); b = sym[4 * n + 2 + delta] * np.conj(r[2 * n + 1])
    return a, b                                                               # the pair's two REs, despread


def ports(sym, slot, l, n_id, n_scid, crbs):
    out = np.zeros((len(crbs), 3, 4), complex)
    for g in (0, 1):
        a, b = despread(sym, slot, l, n_id, n_scid, crbs, g)
        out[:, :, 2 * g], out[:, :, 2 * g + 1] = (a + b) / 2, (a - b) / 2
    return out


def ports2(sym, slot, l, n_id, n_scid, crbs):
    """type 2 (TS 38.211 7.4.1.1.2): CDM group g on k = 6n + k' + 2g, r index 2n + k', n = 2 CRB + block ->
    (P, 2 blocks, 6): ports 1000..1005 (group g: FD-OCC + / - over k' = 0, 1)"""
    r = nrseq.qpsk(nrseq.dmrs_cinit(slot, l, n_id, n_scid), 4 * (max(crbs) + 1))
    n = 2 * np.asarray(crbs)[:, None] + np.arange(2)[None, :]
    out = np.zeros((len(crbs), 2, 6), complex)
    for g in range(3):
        a = sym[6 * n + 2 * g] * np.conj(r[2 * n]); b = sym[6 * n + 1 + 2 * g] * np.conj(r[2 * n + 1])
        out[:, :, 2 * g], out[:, :, 2 * g + 1] = (a + b) / 2, (a - b) / 2
    return out


def coherence(v):
    """|sum v_i conj(v_i+1)| / sum |v|^2 along the flattened pairs within PRBs (1: smooth channel, 0: noise)"""
    s = (v[:, 1:] * np.conj(v[:, :-1])).sum(); e = 0.5 * (np.abs(v[:, 1:]) ** 2 + np.abs(v[:, :-1]) ** 2).sum()
    return float(abs(s) / e) if e > 0 else 0.0
