#!/usr/bin/env python3
"""Is the DM-RS excess noise thermal or a signal? Across the 4 chains: the spatial covariance (eigenvalues, dB over
the smallest) of (1) the tablet's empty port 1003 (type 2, nSCID 0, PCI), (2) the TRS line residual, (3) FFT bins
outside the carrier (receiver noise), plus their powers per RE. Thermal noise: flat eigenvalues; a co-scheduled UE or
a neighbour cell: one or two dominant.

  python3 check3.py CHUNK_DIR CHUNK_INDEX [n_grants]"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ofdm, dmrs, nrseq

IQ = "/home/sens/isac_data/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq/iq_ch%d.sc16"
PCI, D, SLOT0, AOFF = 630, 410, 897 + 4384, 5985


def full_slot(k, a):
    s = int(round(k * 2.5 * ofdm.FS)) + SLOT0 + (a - AOFF) * ofdm.SLOT; out = np.zeros((4, 14, 4096), complex)
    for c in range(4):
        x = ofdm.read(IQ % c, s, ofdm.SLOT + 512).astype(np.complex128)
        for l in range(14):
            out[c, l] = ofdm.symbol(x, ofdm.START[l], ofdm.CP[l])
    return out


def eig(V):
    """V (4, n) -> eigenvalues of the 4x4 covariance in dB over the smallest, and the power per element dB"""
    R = V @ V.conj().T / V.shape[1]; e = np.sort(np.linalg.eigvalsh(R))[::-1]
    return np.round(10 * np.log10(e / e[-1]), 1), 10 * np.log10(np.trace(R).real / 4)


def main(cd, k, ng):
    rows = list(csv.DictReader(open(cd + "/b%d/replay_blind_ch0/golden_sniffer_grants.csv" % k)))[:ng]
    E, T, Gd = [], [], []
    for r in rows:
        a = int(r["sfn"]) * 20 + int(r["slot"]); slot = a % 20; F = full_slot(k, a); S = F[:, :, D:D + 3276]
        crbs = sorted({c for run in r["prb_runs"].split(";") for c in range(int(run.split("-")[0]), int(run.split("-")[1]) + 1)})
        E.append(np.array([dmrs.ports2(S[c, 2], slot, 2, PCI, 0, crbs)[:, :, 3].ravel() for c in range(4)]))
        Gd.append(np.concatenate([F[:, 2, 40:D - 40], F[:, 2, D + 3276 + 40:4096 - 40]], 1))   # outside the carrier, away from DC and edges
        if slot in (6, 7):
            cc = np.arange(273); n = 3 * cc[:, None] + np.arange(3)[None, :]
            rr = nrseq.qpsk(((1 << 10) * (14 * slot + 4 + 1) * (2 * PCI + 1) + PCI) % (1 << 31), 3 * 273 + 4)
            v = np.array([S[c, 4][12 * cc[:, None] + 4 * np.arange(3)[None, :]] * np.conj(rr[n]) for c in range(4)])
            if dmrs.coherence(v[0]) > 0.8:
                T.append((v[:, :, 0] - 2 * v[:, :, 1] + v[:, :, 2]) / np.sqrt(6))   # per-RE noise scale
    e1, p1 = eig(np.concatenate(E, 1)); print("empty port 1003 (pair value; per RE +3 dB): eig %s, power %.1f dB/value -> %.1f dB/RE" % (e1, p1, p1 + 3))
    if T:
        e2, p2 = eig(np.concatenate(T, 1)); print("TRS line residual (per RE):                 eig %s, power %.1f dB/RE (n %d)" % (e2, p2, len(T)))
    e3, p3 = eig(np.concatenate(Gd, 1)); print("bins outside the carrier (receiver noise):  eig %s, power %.1f dB/bin" % (e3, p3))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 40)
