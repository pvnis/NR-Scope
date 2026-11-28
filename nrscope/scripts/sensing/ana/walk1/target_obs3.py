"""Observation 3 for the TARGET (walk1). From jp_cache (per window w, group b, antenna c, hypothesis h: matched sums r and norms n of
the high-passed DM-RS (dominant tablet beam) and TRS against the lidar steering). Per stream z[w,b] = sum_c conj(a_c) r/n
(a from the TRS, all windows). q[w,b] = z_dmrs conj(z_trs): the target's DM-RS/TRS ratio, window phase cancelled.
 coherence over windows per group  C_b = |sum_w q| / sum_w |q|   (1 = the ratio persists; noise only ~ 1/sqrt(Nw))
 adjacent-window coherence          A_b = |sum q_w conj q_w+1| / sum |q_w||q_w+1|
 at the truth range and at +25 m (noise / clutter reference); also the per-(w,b) SNR of each stream.
usage: target_obs3.py <chunks>"""
import sys, numpy as np
for k in [int(x) for x in sys.argv[1].split(",")]:
    Z = np.load("/mnt/mydrive/field/ana_walk1/jp_cache_b%d.npz" % k); v = Z["valid"]
    rd, nd, rt, nt, sd2, st2 = Z["rd"][:, v], Z["nd"][:, v], Z["rt"][:, v], Z["nt"][:, v], float(Z["sd2"]), float(Z["st2"])
    gt = rt[0] / np.maximum(nt[0], 1e-9)[..., None]                      # [W, B, 4]
    a = np.linalg.eigh(np.einsum("wbc,wbd->cd", gt, gt.conj()))[1][:, -1]
    out = []
    for h in (0, 3):
        zd = np.einsum("wbc,c->wb", rd[h], a.conj()) / np.maximum(nd[h], 1e-9); zt = np.einsum("wbc,c->wb", rt[h], a.conj()) / np.maximum(nt[h], 1e-9)
        ok = (nd[h] > 0) & (nt[h] > 0); q = np.where(ok, zd * zt.conj(), 0)
        snr_d = np.abs(zd) ** 2 / (sd2 / np.maximum(nd[h], 1e-9)); snr_t = np.abs(zt) ** 2 / (st2 / np.maximum(nt[h], 1e-9))
        C = np.abs(q.sum(0)) / np.maximum(np.abs(q).sum(0), 1e-30)
        A = np.abs((q[:-1] * q[1:].conj()).sum(0)) / np.maximum((np.abs(q[:-1]) * np.abs(q[1:])).sum(0), 1e-30)
        # TRS against itself across groups (does the TRS target echo keep its per-group phase?): zt / (window mean over groups)
        kt = (zt * ok).sum(1, keepdims=True); Ct = np.abs((zt * kt.conj()).sum(0)) / np.maximum(np.abs(zt * kt.conj()).sum(0), 1e-30)
        kd = (zd * ok).sum(1, keepdims=True); Cd = np.abs((zd * kd.conj()).sum(0)) / np.maximum(np.abs(zd * kd.conj()).sum(0), 1e-30)
        out.append("%s: q coherence median %.2f, adjacent %.2f | per-group vs own band mean: TRS %.2f, DM-RS %.2f | SNR per (w,b) median DM-RS %.1f dB, TRS %.1f dB"
                   % ("truth" if h == 0 else "+25 m", np.median(C), np.median(A), np.median(Ct), np.median(Cd), 10 * np.log10(np.median(snr_d[ok])), 10 * np.log10(np.median(snr_t[ok]))))
    print("chunk %d (%d windows, %d groups; noise-only coherence ~%.2f):\n   %s\n   %s" % (k, len(v), rd.shape[2], 1 / np.sqrt(len(v)), out[0], out[1]))
