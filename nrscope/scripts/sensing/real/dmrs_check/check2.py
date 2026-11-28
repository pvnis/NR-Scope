#!/usr/bin/env python3
"""Independent type-2 check of a UE's PDSCH DM-RS from the raw IQ: per port (1000-1005), scrambling ID and nSCID, the
power and the coherence between the PRB's two blocks; the noise per value from what must be empty (the unused OCC
code and CDM groups); the same noise on the TRS (residual of a line through its 3 values per PRB) and on REs
outside every grant. All in the same FFT units, chain 0.

  python3 check2.py CHUNK_DIR CHUNK_INDEX [n_grants] [nids, comma]"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grid, dmrs, nrseq

DESLOPE = os.environ.get("DESLOPE", "0") == "1"
IQ = "/home/sens/isac_data/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq/iq_ch%d.sc16"
PCI = 630


def main(cd, k, ng, nids):
    C0 = cd + "/b%d/replay_blind_ch0" % k; G = grid.Grid([IQ % c for c in range(4)], k * 2.5, 897 + 4384, 5985)
    rows = list(csv.DictReader(open(C0 + "/golden_sniffer_grants.csv")))[:ng]; A = {}
    add = lambda key, v: A.setdefault(key, []).append(v)
    for r in rows:
        a = int(r["sfn"]) * 20 + int(r["slot"]); slot = a % 20; S = G.slot(a)[0]
        crbs = sorted({c for run in r["prb_runs"].split(";") for c in range(int(run.split("-")[0]), int(run.split("-")[1]) + 1)})
        if DESLOPE:   # one delay per slot from the TRS-free estimate: port 1000's two blocks (6 subcarriers apart)
            for l in range(14):
                P0 = dmrs.ports2(S[2], slot, 2, PCI, 0, crbs)[:, :, 0]; z = np.sum(P0[:, 1] * np.conj(P0[:, 0]))
                S[l] = S[l] * np.exp(-1j * np.angle(z) / 6 * np.arange(3276))
                if l == 0:
                    add(("slope", 0), np.angle(z) / 6)
        for l in (2, 11):
            for nid in nids:
                for sc in (0, 1):
                    P = dmrs.ports2(S[l], slot, l, nid, sc, crbs)
                    for p in range(6):
                        add(("pow", l, nid, sc, p), np.mean(np.abs(P[:, :, p]) ** 2)); add(("coh", l, nid, sc, p), dmrs.coherence(P[:, :, p]))
        if slot in (6, 7):   # TRS slots (only every other frame carries it): keep the coherent ones
            for l in (4, 8):
                c = np.arange(273); n = 3 * c[:, None] + np.arange(3)[None, :]
                rr = nrseq.qpsk(((1 << 10) * (14 * slot + l + 1) * (2 * PCI + 1) + PCI) % (1 << 31), 3 * 273 + 4)
                v = S[l][12 * c[:, None] + 4 * np.arange(3)[None, :]] * np.conj(rr[n])
                if dmrs.coherence(v) > 0.8:
                    d2 = v[:, 0] - 2 * v[:, 1] + v[:, 2]; add(("trs_pow", l), np.mean(np.abs(v) ** 2)); add(("trs_noise", l), np.mean(np.abs(d2) ** 2) / 6)
        out = np.setdiff1d(np.arange(273), crbs)
        if len(out) >= 10:
            for l in range(14):
                add(("out", l), np.mean(np.abs(S[l][(12 * out[:, None] + np.arange(12)[None, :]).ravel()]) ** 2))
    db = lambda x: 10 * np.log10(np.mean(x)) if len(x) else float("nan")
    print("chunk %d, %d grants of rnti %s (type 2)" % (k, len(rows), rows[0]["rnti"]))
    for l in (2, 11):
        for nid in nids:
            for sc in (0, 1):
                print("  sym %2d nid %5d scid %d | power dB %s | coherence %s" % (l, nid, sc, " ".join("%6.1f" % db(A[("pow", l, nid, sc, p)]) for p in range(6)),
                                                                            " ".join("%.2f" % np.mean(A[("coh", l, nid, sc, p)]) for p in range(6))))
    if ("slope", 0) in A:
        sl = np.median(A[("slope", 0)]); print("  removed phase slope: median %.4f rad per subcarrier (= %.0f ns, %.0f samples)" % (sl, sl / (2 * np.pi * 30e3) * 1e9, sl / (2 * np.pi * 30e3) * 122.88e6))
    for l in (4, 8):
        if ("trs_pow", l) in A:
            print("  TRS sym %d: power %.1f dB, noise per value %.1f dB (line residual), n %d" % (l, db(A[("trs_pow", l)]), db(A[("trs_noise", l)]), len(A[("trs_pow", l)])))
    if ("out", 0) in A:
        print("  REs outside the grant, per symbol dB:", " ".join("%d:%.1f" % (l, db(A[("out", l)])) for l in range(14)))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 40, [int(x) for x in (sys.argv[4] if len(sys.argv) > 4 else "630,27109").split(",")])
