#!/usr/bin/env python3
"""What sits on the tablet's empty DM-RS code? Type 2, symbol 2, CDM group 1 (k = 6n + 2 + k'): the tablet's own
port 1002 uses the + code, so e = (a r_a* - b r_b*) / 2 (port 1003, our sequence r) holds only noise and anything
else on those REs. For a candidate transmitter with sequence s on its port 1002 (+) or 1003 (-), e = h w with
w = (s_a r_a* -+ s_b r_b*) / 2 known: h = e / w where |w| > 0.5, and its coherence across neighbouring blocks is ~1
for the right s, ~0 otherwise. Candidates: N_ID in 0..1007 and extra ones, nSCID 0 / 1. Chain 0.

  python3 scan.py CHUNK_DIR CHUNK_INDEX [n_grants] [extra nids]"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grid, nrseq

IQ = "/home/sens/isac_data/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq/iq_ch%d.sc16"
PCI = 630


def coh(h, ok):
    p = h[1:] * np.conj(h[:-1]); m = ok[1:] & ok[:-1]; e = 0.5 * (np.abs(h[1:]) ** 2 + np.abs(h[:-1]) ** 2)
    return float(np.abs(p[m].sum()) / max(e[m].sum(), 1e-30))


def main(cd, k, ng, extra):
    G = grid.Grid([IQ % 0], k * 2.5, 897 + 4384, 5985)
    rows = list(csv.DictReader(open(cd + "/b%d/replay_blind_ch0/golden_sniffer_grants.csv" % k)))[:ng]
    nids = list(range(1008)) + extra; cand = [(n, s) for n in nids for s in (0, 1)]; score = np.zeros((len(cand), 2)); l = 2
    for r in rows:
        a = int(r["sfn"]) * 20 + int(r["slot"]); slot = a % 20; S = G.slot(a)[0, l]
        crbs = np.array(sorted({c for run in r["prb_runs"].split(";") for c in range(int(run.split("-")[0]), int(run.split("-")[1]) + 1)}))
        n = (2 * crbs[:, None] + np.arange(2)[None, :]).ravel(); A, B = S[6 * n + 2], S[6 * n + 3]
        R = nrseq.qpsk(nrseq.dmrs_cinit(slot, l, PCI, 0), 2 * (n.max() + 1)); ra, rb = R[2 * n], R[2 * n + 1]
        e = 0.5 * (A * np.conj(ra) - B * np.conj(rb))
        Sq = nrseq.qpsk_many([nrseq.dmrs_cinit(slot, l, nid, sc) for nid, sc in cand], 2 * (n.max() + 1)); sa, sb = Sq[:, 2 * n], Sq[:, 2 * n + 1]
        for j, sg in enumerate((-1, 1)):   # candidate on port 1002 (+ code: w = sa ra* - sb rb*) or 1003 (- code)
            w = 0.5 * (sa * np.conj(ra) + sg * sb * np.conj(rb)); ok = np.abs(w) > 0.5; h = np.where(ok, e[None, :] / np.where(ok, w, 1), 0)
            score[:, j] += [coh(h[i], ok[i]) for i in range(len(cand))]
    score /= len(rows); o = np.argsort(-score.max(1))
    print("chunk %d, %d grants: coherence of the empty code under candidate sequences (median over candidates %.3f, p99 %.3f)" % (k, len(rows), np.median(score), np.percentile(score, 99)))
    for i in o[:10]:
        print("  N_ID %5d nSCID %d: as port 1002 %.3f, as port 1003 %.3f" % (*cand[i], *score[i]))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20, [int(x) for x in (sys.argv[4] if len(sys.argv) > 4 else "27109").split(",")])
