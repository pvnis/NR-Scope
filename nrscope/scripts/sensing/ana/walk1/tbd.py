"""Track-before-detect on the radar maps (sync_walk.py npz): per map, the evidence per range bin is the best cell-
normalised power over the moving Doppler bins (0.2 < |v| < 3 m/s); a Viterbi path over range bins with at most MAXSTEP
bins between consecutive maps (and a small penalty per bin moved) gives the single most supported moving-target
track. Scored against the lidar walkers' excess path at the radio<->lidar offset: share of the walk the track is within
3 m (and 5 m) of a walker; the same against the truth shifted by +-15..30 s gives the chance level.
usage: tbd.py <npz> <offset s> [maxstep] [png]"""
import sys, numpy as np
D = np.load(sys.argv[1]); OFF = float(sys.argv[2]); MAXSTEP = int(sys.argv[3]) if len(sys.argv) > 3 else 1
T, Z, mpb, t_rec, TR = D["T"], D["Z"].astype(np.float32), float(D["mpb"]), float(D["t_rec"]), D["tr"]
o = np.argsort(T); T, Z = T[o], Z[o]; nf = Z.shape[2]; v = np.linspace(-1, 1, nf) * 650.173071 * 0.04
med = np.median(Z, 0); mad = 1.4826 * np.median(np.abs(Z - med), 0) + 1e-6; Zn = (Z - med) / mad
mv = (np.abs(v) > 0.2) & (np.abs(v) < 3); EV = Zn[:, :, mv].max(2)           # maps x range bins
R0, R1 = 4, 30                                                               # 10-73 m: off the direct path
ev0 = EV[:, R0:R1].copy(); n, m = ev0.shape; pen = 0.5; K = int(__import__("os").environ.get("K", "1")); PATHS = []
for _k in range(K):
  ev = ev0
  score = ev[0].copy(); back = np.zeros((n, m), int)
  for i in range(1, n):
      best = np.full(m, -1e9); arg = np.zeros(m, int)
      for d in range(-MAXSTEP, MAXSTEP + 1):
          sh = np.roll(score, d) - pen * abs(d)
          if d > 0: sh[:d] = -1e9
          if d < 0: sh[d:] = -1e9
          upd = sh > best; best[upd] = sh[upd]; arg[upd] = (np.arange(m) - d)[upd]
      score = best + ev[i]; back[i] = arg
  path = np.zeros(n, int); path[-1] = int(np.argmax(score))
  for i in range(n - 1, 0, -1): path[i - 1] = back[i, path[i]]

  PATHS.append(path.copy())
  for i_ in range(n): ev0[i_, max(0, path[i_] - 2):path[i_] + 3] = -5.0   # suppress it for the next track
RNG = [(p_ + R0) * mpb for p_ in PATHS]; rng = RNG[0]
ids = sorted(set(TR[:, 0].astype(int)))
def near(shift):
    w3 = w5 = tot = 0
    for ii, t in enumerate(T):
        te = t_rec + OFF + shift + t; ds = []
        for k in ids:
            mm = TR[:, 0] == k; tk = TR[mm, 1]
            if tk[0] <= te <= tk[-1]: ds += [abs(np.interp(te, tk, TR[mm, 2]) - R_[ii]) for R_ in RNG]
        if ds: tot += 1; w3 += min(ds) < 3; w5 += min(ds) < 5
    return w3 / max(tot, 1), w5 / max(tot, 1), tot
a3, a5, tot = near(0.0); ch = [near(s) for s in (-30, -25, -20, -15, 15, 20, 25, 30)]
print("%-26s %d maps: track within 3 m of a walker %.0f %% (chance %.0f +- %.0f %%), within 5 m %.0f %% (chance %.0f %%)"
      % (sys.argv[1].split("/")[-1], tot, 100 * a3, 100 * np.mean([c[0] for c in ch]), 100 * np.std([c[0] for c in ch]), 100 * a5, 100 * np.mean([c[1] for c in ch])))
if len(sys.argv) > 4:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(13, 4.5)); ax.pcolormesh(T, np.arange(Z.shape[1]) * mpb, EV.T, cmap="magma", vmin=0, vmax=8, shading="nearest")
    for k in ids:
        mm = TR[:, 0] == k; ax.plot(TR[mm, 1] - t_rec - OFF, TR[mm, 2], c="c", lw=1.2, alpha=.8)
    for q_, R_ in enumerate(RNG): ax.plot(T, R_, c=["lime", "yellow", "w"][q_ % 3], lw=1.5, label="track %d" % (q_ + 1)); ax.set_ylim(0, 80); ax.set_xlim(T.min(), T.max())
    ax.set_xlabel("s into the IQ"); ax.set_ylabel("excess path (m)"); ax.legend(loc="upper right", fontsize=8)
    ax.set_title("%s: moving-Doppler evidence, lidar walkers (cyan), extracted track (green)" % sys.argv[1].split("/")[-1])
    fig.tight_layout(); fig.savefig(sys.argv[4], dpi=100); print("wrote", sys.argv[4])
