"""Alignment again from a sync_walk.py npz, with each (range, Doppler) cell normalised against itself over the whole
recording (median / MAD over maps): a band that is there all the time (a fixed scatterer that flickers, foliage, a
road along one ellipse) cancels, a walker passing through does not. Prints the best offsets and writes the score curves
and a range-time plot of the normalised maps.  usage: resync.py <npz> <out.png> [label]"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = np.load(sys.argv[1]); T, Z, mpb, t_rec, TR = D["T"], D["Z"].astype(np.float32), float(D["mpb"]), float(D["t_rec"]), D["tr"]
o = np.argsort(T); T, Z = T[o], Z[o]; nf = Z.shape[2]; v = np.linspace(-1, 1, nf) * 650.173071 * (3e8 / 3.75e9) / 2; dv = v[1] - v[0]
med = np.median(Z, 0); mad = 1.4826 * np.median(np.abs(Z - med), 0) + 1e-6; Zn = (Z - med) / mad
ids = sorted(set(TR[:, 0].astype(int)))
def score(off, sgn, ZZ):
    s = []
    for i, t in enumerate(T):
        tt = t_rec + off + t
        for k in ids:
            m = TR[:, 0] == k; tk = TR[m, 1]
            if tk[0] <= tt <= tk[-1]:
                d = np.interp(tt, tk, TR[m, 2]); vv = sgn * np.interp(tt, tk, TR[m, 3])
                if abs(vv) < 0.2: continue
                b = int(round(d / mpb)); j = int(round((vv - v[0]) / dv))
                if 1 <= b < Z.shape[1] - 1 and 1 <= j < nf - 1: s.append(ZZ[i, b - 1:b + 2, j - 1:j + 2].max())
    return np.mean(s) if len(s) > 50 else np.nan
offs = np.arange(-5, 40, 0.1); R = {}
for sgn in (1, -1):
    R[sgn] = np.array([score(o_, sgn, Zn) for o_ in offs]); s = R[sgn]; i = int(np.nanargmax(s))
    print("%s sign %+d: best IQ start = recorder start %+.1f s, score %.2f; median %.2f std %.2f -> %.1f sigma" % (sys.argv[3] if len(sys.argv) > 3 else "", sgn, offs[i], s[i], np.nanmedian(s), np.nanstd(s), (s[i] - np.nanmedian(s)) / np.nanstd(s)))
best = offs[int(np.nanargmax(R[1]))]
fig, ax = plt.subplots(2, 1, figsize=(13, 8), gridspec_kw=dict(height_ratios=[1, 2.4]))
ax[0].plot(offs, R[1], label="Doppler sign +1"); ax[0].plot(offs, R[-1], label="sign -1", alpha=.6); ax[0].axvline(best, c="C3", ls="--"); ax[0].legend(fontsize=8)
ax[0].set_xlabel("IQ start - recorder start (s)"); ax[0].set_ylabel("mean z at truth (cell-normalised)")
mv = (np.abs(v) > 0.25) & (np.abs(v) < 3); RTI = Zn[:, :, mv].max(2)
ax[1].pcolormesh(T, np.arange(Z.shape[1]) * mpb, RTI.T, cmap="magma", vmin=0, vmax=8, shading="nearest")
for k in ids:
    m = TR[:, 0] == k; ax[1].plot(TR[m, 1] - t_rec - best, TR[m, 2], c="c", lw=1.2, alpha=.8)
ax[1].set_ylim(0, 80); ax[1].set_xlim(T.min(), T.max()); ax[1].set_ylabel("excess path (m)"); ax[1].set_xlabel("s into the IQ")
ax[1].set_title("range-time, each cell normalised over the recording; cyan = lidar walkers at %+.1f s" % best)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=100); print("wrote", sys.argv[2])
