"""Plots for the radio <-> lidar alignment (sync_walk.py npz): score vs offset; range-time intensity (max over the
moving Doppler bins) with the lidar truth; Doppler-time at the truth's range bin with the truth's bistatic speed.
usage: plot_sync.py <npz> <offset s> <out.png>"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = np.load(sys.argv[1]); off = float(sys.argv[2]); T, Z, mpb, t_rec = D["T"], D["Z"], float(D["mpb"]), float(D["t_rec"])
nf = Z.shape[2]; v = np.linspace(-1, 1, nf) * 650.173071 * (3e8 / 3.75e9) / 2   # m/s, common axis
TR = D["tr"]; ids = sorted(set(TR[:, 0].astype(int)))
o = np.argsort(T); T, Z = T[o], Z[o]
mv = (np.abs(v) > 0.25) & (np.abs(v) < 3)
RTI = Z[:, :, mv].max(2)                                            # maps x range bins
fig = plt.figure(figsize=(13, 11)); g = fig.add_gridspec(3, 1, height_ratios=[1, 2.2, 2.2])
a = fig.add_subplot(g[0]); a.plot(D["offs"], D["sp"][:, 0], label="Doppler sign +1"); a.plot(D["offs"], D["sm"][:, 0], label="sign -1", alpha=.6)
a.axvline(off, c="C3", ls="--"); a.set_xlabel("IQ start - recorder start (s)"); a.set_ylabel("mean z at truth"); a.legend(fontsize=8)
a = fig.add_subplot(g[1]); a.pcolormesh(T, np.arange(Z.shape[1]) * mpb, RTI.T, cmap="magma", vmin=0, vmax=6, shading="nearest")
for i in ids:
    s = TR[:, 0] == i; a.plot(TR[s, 1] - t_rec - off, TR[s, 2], lw=1.2, c="c", alpha=.8)
a.set_ylim(0, 80); a.set_ylabel("excess path (m)"); a.set_title("range-time: max z over 0.25 < |v| < 3 m/s; cyan = lidar walkers at offset %+.1f s" % off)
a.set_xlim(T.min(), T.max())
# Doppler-time at the main walker's range
a = fig.add_subplot(g[2]); s = TR[:, 0] == ids[0]; tt = TR[s, 1] - t_rec - off
dr = np.interp(T, tt, TR[s, 2], left=np.nan, right=np.nan); vt = np.interp(T, tt, TR[s, 3], left=np.nan, right=np.nan)
DT = np.full((len(T), nf), np.nan)
for k in range(len(T)):
    if np.isfinite(dr[k]):
        b = int(round(dr[k] / mpb))
        if 1 <= b < Z.shape[1] - 1: DT[k] = Z[k, b - 1:b + 2].max(0)
a.pcolormesh(T, v, DT.T, cmap="magma", vmin=0, vmax=6, shading="nearest"); a.plot(T, vt, c="c", lw=1.2, alpha=.8)
for i in ids[1:]:
    s2 = TR[:, 0] == i; a.plot(TR[s2, 1] - t_rec - off, TR[s2, 3], c="w", lw=.8, alpha=.5)
a.set_ylim(-3, 3); a.set_ylabel("bistatic speed (m/s)"); a.set_xlabel("s into the IQ"); a.set_xlim(T.min(), T.max())
a.set_title("Doppler-time at walker %d's range (+-1 bin); cyan = its lidar bistatic speed" % ids[0])
fig.tight_layout(); fig.savefig(sys.argv[3], dpi=100); print("wrote", sys.argv[3])
