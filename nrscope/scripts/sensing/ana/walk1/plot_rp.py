import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Z = np.load("/mnt/mydrive/field/ana_walk1/range_prof.npz"); D = Z["DEL"]; i0 = int(np.argmin(np.abs(D)))
S = np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz"); TR = S["tr"]; t_rec = float(S["t_rec"])
def wr(k):
    v = []
    for s in np.arange(0.1, 2.5, 0.2):
        te = t_rec + 20.2 + k * 2.5 + s
        for tid in (1, 4):
            m = TR[:, 0] == tid; tk = TR[m, 1]
            if tk[0] <= te <= tk[-1]: v.append(np.interp(te, tk, TR[m, 2])); break
    return np.mean(v), (v[-1] - v[0]) / 2.4
fig, ax = plt.subplots(2, 4, figsize=(15, 6.6), sharex=True, sharey="row")
for j, k in enumerate((3, 12, 13, 15)):
    r, vr = wr(k)
    for row, (lab, col, nm) in enumerate((("trs", "#1f5fa8", "TRS"), ("dm_anch", "#d1495b", "DM-RS (tablet, anchored)"))):
        a_ = ax[row, j]; p = Z["%s_%d" % (lab, k)]; f = Z["%s_flip_%d" % (lab, k)]; ref = p[i0]
        a_.plot(D, 10 * np.log10(p / ref), color=col, lw=1.8, label="walker's Doppler")
        a_.plot(D, 10 * np.log10(f / ref), color=col, lw=1.1, ls="--", alpha=0.8, label="Doppler reversed (control)")
        a_.axvline(0, color="k", lw=0.8); a_.axvline(-r, color="#7a7a7a", lw=0.8, ls=":")
        a_.set_ylim(-20, 13); a_.grid(alpha=0.3)
        if row == 0: a_.set_title("chunk %d: walker at %.0f m, %+.1f m/s in range" % (k, r, vr), fontsize=9.5)
        if j == 0: a_.set_ylabel("%s\nenergy rel. to walker cell [dB]" % nm, fontsize=9)
        if row == 1: a_.set_xlabel("range offset from the walker (lidar) [m]", fontsize=9)
ax[0, 0].legend(fontsize=7.5, loc="lower left"); ax[0, 0].text(-r if False else -18, 11, "dotted: absolute range 0\n(direct path leak)", fontsize=7, color="#555")
fig.suptitle("walk1: walker range profile, 0.2 s coherent windows. TRS: walker peak at 0 drops when his Doppler is reversed (it is him). DM-RS: no change -> the floor is thermal noise, the walker is below it.", fontsize=10)
fig.tight_layout(); fig.savefig("/mnt/mydrive/field/ana_walk1/range_prof_walk1.png", dpi=130); print("ok")
