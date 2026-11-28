"""Figure 1 analysis: raw DMRS (one complex number per snapshot) against TRS-anchored DMRS (one gain per PRG and precoder
epoch), 18 Nov continuous capture (Sunrise PCI 632), same processing otherwise (so_light as um_20261118b.sh, sensing
cfg_trs.yaml, antennas averaged by the sensing: antenna_avg true). Maps from fig1_run.sh in /mnt/mydrive/field/ana_fig1.

Scale: every map is put back in absolute units (power x noise_ref) and divided by the energy of the same window without
clutter removal (the cfg_trs_noclut.yaml map, summed over all cells): dB relative to the received energy, which the
direct path dominates. Same reference definition for both panels.
Metrics per map (no ground truth: the "target" is the strongest moving cell of the anchored map):
  (a) residue: median over cells outside the zero-Doppler band (|v| > V0) and beyond the direct path (> G bins)
  (b) target peak: the strongest such cell of the anchored map; the raw map is read at the same cell
  (c) peak minus local residue: (b) against the median of the +-8 bins x +-8 Doppler cells around it (its +-1 excluded)
  (d) ground truth: not available (no LiDAR for this capture)
A "clear mover" is (c) >= 15 dB in the anchored map, and found again (within 2 bins, 1 m/s) in a neighbouring map.
usage: fig1_ana.py <map dir> <out dir>"""
import sys, os, glob, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
LAM = 299792458.0 / 3.75e9; V0, G, CLEAR = 0.6, 2, 15.0


def read(f):
    M = []
    for l in open(f):
        if l.startswith("#"): continue
        x = l.rstrip("\n").split(","); fr, sl = int(x[0]), int(x[1]); nb, nf = int(x[12]), int(x[13])
        mpb, fmax, nref = float(x[9]), float(x[10]), float(x[20])
        P = np.array(x[21: 21 + nb * nf], float).reshape(nb, nf) * nref
        M.append(dict(t=fr * 0.01 + sl * 0.0005, nb=nb, nf=nf, mpb=mpb, fmax=fmax, P=P))
    return M


def metrics(P, ref, mpb, fmax, cell=None):
    nb, nf = P.shape; v = np.linspace(-fmax, fmax, nf) * LAM / 2
    Z = 10 * np.log10(np.maximum(P, 1e-30) / ref)
    mov = (np.abs(v)[None, :] > V0) & (np.arange(nb)[:, None] > G)
    res = float(np.median(Z[mov]))
    if cell is None:
        Zm = np.where(mov, Z, -np.inf); cell = np.unravel_index(int(np.argmax(Zm)), Z.shape)
    b, f = cell; pk = float(Z[b, f])
    w = np.zeros_like(mov); w[max(0, b - 8): b + 9, max(0, f - 8): f + 9] = True; w[max(0, b - 1): b + 2, max(0, f - 1): f + 2] = False
    loc = float(np.median(Z[w & mov])) if (w & mov).any() else res
    return dict(residue_db=res, peak_db=pk, peak_minus_local_db=pk - loc, range_m=float(b * mpb), speed_ms=float(v[f]), cell=[int(b), int(f)]), Z, v


rows = []; keep = {}
for b in range(48):
    S = {}
    for var in ("anch", "raw"):
        fc, fn = os.path.join(D, "b%d_%s_clut.map2d.csv" % (b, var)), os.path.join(D, "b%d_%s_noclut.map2d.csv" % (b, var))
        if not (os.path.exists(fc) and os.path.exists(fn)): break
        S[var] = (read(fc), read(fn))
    if len(S) < 2: continue
    A, An = S["anch"]; R, Rn = S["raw"]
    t0 = A[0]["t"]
    tt_ = lambda L: np.array([(m["t"] - t0) % 10.24 for m in L])
    tA, tAn, tR, tRn = tt_(A), tt_(An), tt_(R), tt_(Rn)
    near = lambda tl, x: int(np.argmin(np.abs(tl - x))) if len(tl) else -1
    for i in range(len(A)):
        jn, jr = near(tAn, tA[i]), near(tR, tA[i]); jrn = near(tRn, tR[jr]) if jr >= 0 else -1
        if min(jn, jr, jrn) < 0 or abs(tAn[jn] - tA[i]) > 0.03 or abs(tR[jr] - tA[i]) > 0.03 or abs(tRn[jrn] - tR[jr]) > 0.03: continue
        if A[i]["P"].shape[0] != R[jr]["P"].shape[0]: continue
        ra, rr = An[jn]["P"].sum(), Rn[jrn]["P"].sum()
        va = np.linspace(-A[i]["fmax"], A[i]["fmax"], A[i]["nf"]); vr = np.linspace(-R[jr]["fmax"], R[jr]["fmax"], R[jr]["nf"])
        PR = np.array([np.interp(va, vr, row) for row in R[jr]["P"]])          # raw map on the anchored Doppler grid
        ma, Za, v = metrics(A[i]["P"], ra, A[i]["mpb"], A[i]["fmax"])
        mr, Zr, _ = metrics(PR, rr, A[i]["mpb"], A[i]["fmax"], cell=tuple(ma["cell"]))
        mr_own, _, _ = metrics(PR, rr, A[i]["mpb"], A[i]["fmax"])
        rows.append(dict(chunk=b, i=i, t_s=b * 2.5 + tA[i], anch=ma, raw_same_cell=mr, raw_own_peak=mr_own))
        keep[(b, i)] = (Za, Zr, v, A[i]["mpb"])
# clear movers: anchored (c) >= CLEAR and confirmed by a neighbouring map of the same chunk
by = {(r["chunk"], r["i"]): r for r in rows}
for r in rows:
    a = r["anch"]; conf = False
    for di in (-1, 1):
        q = by.get((r["chunk"], r["i"] + di))
        if q and abs(q["anch"]["range_m"] - a["range_m"]) <= 2 * 2.44 and abs(q["anch"]["speed_ms"] - a["speed_ms"]) <= 1.0 and q["anch"]["peak_minus_local_db"] >= CLEAR - 3:
            conf = True
    r["clear_mover"] = bool(a["peak_minus_local_db"] >= CLEAR and conf)
q = lambda x: dict(median=float(np.median(x)), p25=float(np.percentile(x, 25)), p75=float(np.percentile(x, 75)))
typ = dict(n_maps=len(rows), n_chunks=len({r["chunk"] for r in rows}),
           residue_db=dict(raw=q([r["raw_same_cell"]["residue_db"] for r in rows]), anch=q([r["anch"]["residue_db"] for r in rows])),
           peak_minus_local_db_at_anchored_target=dict(raw=q([r["raw_same_cell"]["peak_minus_local_db"] for r in rows]), anch=q([r["anch"]["peak_minus_local_db"] for r in rows])),
           share_clear_mover_anch=float(np.mean([r["clear_mover"] for r in rows])),
           share_raw_also_clear_at_that_cell=float(np.mean([r["raw_same_cell"]["peak_minus_local_db"] >= CLEAR for r in rows if r["clear_mover"]])) if any(r["clear_mover"] for r in rows) else None)
# chosen map: a confirmed clear mover with the largest anchored-vs-raw contrast at the target cell
cand = [r for r in rows if r["clear_mover"]]
best = max(cand, key=lambda r: r["anch"]["peak_minus_local_db"] - r["raw_same_cell"]["peak_minus_local_db"]) if cand else max(rows, key=lambda r: r["anch"]["peak_minus_local_db"])
Za, Zr, v, mpb = keep[(best["chunk"], best["i"])]
rng = np.arange(Za.shape[0]) * mpb
np.savez(os.path.join(OUT, "fig1_panel_raw.npz"), range_m=rng, speed_ms=v, map_db=Zr)
np.savez(os.path.join(OUT, "fig1_panel_anchored.npz"), range_m=rng, speed_ms=v, map_db=Za)
for nm, Z in (("raw", Zr), ("anchored", Za)):
    with open(os.path.join(OUT, "fig1_panel_%s.csv" % nm), "w") as fo:
        fo.write("# dB relative to the window's received energy; rows: range_m; columns: speed_ms\nrange_m," + ",".join("%.3f" % x for x in v) + "\n")
        for k in range(Z.shape[0]): fo.write("%.2f," % rng[k] + ",".join("%.2f" % x for x in Z[k]) + "\n")
# range-time alternative: ~10 s around the chosen map, per map the strongest moving cell of each range bin
tc = best["t_s"]; sel = sorted([r for r in rows if abs(r["t_s"] - tc) <= 5.0], key=lambda r: r["t_s"])
RT = {}
for nm, k in (("raw", 1), ("anchored", 0)):
    cols = []
    for r in sel:
        Z = keep[(r["chunk"], r["i"])][k]; vv = keep[(r["chunk"], r["i"])][2]
        cols.append(np.max(np.where(np.abs(vv)[None, :] > V0, Z, -np.inf), axis=1))
    nb = min(len(c) for c in cols); RT[nm] = np.array([c[:nb] for c in cols]).T
tt = np.array([r["t_s"] for r in sel])
np.savez(os.path.join(OUT, "fig1_rangetime.npz"), t_s=tt, range_m=np.arange(RT["raw"].shape[0]) * mpb, raw_db=RT["raw"], anchored_db=RT["anchored"])
out = dict(capture="20261118_154424_sunrise_632_4ch_cont_tabletdl_iq", ground_truth="not available (no LiDAR for this capture)",
           processing=dict(so_light="MODE=tab MINPRB=48 USE_DMRS=1 USE_CSI=1 (as um_20261118b.sh); raw: PREC_IGNORE=1",
                           sensing="cfg_trs.yaml (clutter removal on, antenna_avg: maps averaged over the 4 antennas)",
                           reference="each map divided by the energy of the same window without clutter removal (cfg_trs_noclut.yaml)",
                           zero_doppler_band_ms=V0, guard_bins=G, range_max_m=float(rng[-1])),
           chosen=dict(chunk=best["chunk"], map=best["i"], t_s=best["t_s"], anchored=best["anch"], raw_same_cell=best["raw_same_cell"], raw_own_peak=best["raw_own_peak"]),
           typicality=typ, per_map=rows)
json.dump(out, open(os.path.join(OUT, "fig1_metrics.json"), "w"), indent=1)
# quick look
vmax = max(np.max(Za), -10); vmin = vmax - 45
fig, ax = plt.subplots(2, 2, figsize=(10, 7))
for k, (nm, Z) in enumerate((("raw DMRS", Zr), ("TRS-anchored DMRS", Za))):
    im = ax[0, k].imshow(Z.T, origin="lower", aspect="auto", extent=[rng[0], rng[-1], v[0], v[-1]], vmin=vmin, vmax=vmax, cmap="magma")
    ax[0, k].set_title("%s, chunk %d map %d" % (nm, best["chunk"], best["i"])); ax[0, k].set_xlabel("bistatic range (m)"); ax[0, k].set_ylabel("speed (m/s)")
    ax[0, k].plot(best["anch"]["range_m"], best["anch"]["speed_ms"], "o", mfc="none", mec="c", ms=12)
    Rt = RT["raw" if k == 0 else "anchored"]
    ax[1, k].imshow(Rt, origin="lower", aspect="auto", extent=[tt[0], tt[-1], 0, Rt.shape[0] * mpb], vmin=vmin, vmax=vmax, cmap="magma")
    ax[1, k].set_xlabel("time (s)"); ax[1, k].set_ylabel("bistatic range (m)"); ax[1, k].set_title("%s, strongest moving cell per range" % nm)
fig.colorbar(im, ax=ax, label="dB rel. received energy", shrink=0.8)
fig.savefig(os.path.join(OUT, "fig1_quicklook.png"), dpi=110)
print(json.dumps(dict(chosen=out["chosen"], typicality=typ), indent=1))
