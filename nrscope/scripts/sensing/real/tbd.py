"""Track-before-detect on range x speed x angle cubes (cube.py), the TRS and DM-RS maps in one Bayesian recursion.

Per stream, each cell's whitened power z is first divided by its own slow background (exponential average of the
past maps of that stream, time constant bg_tau_s; >= 1), so cells that are always lit (static leakage, fans, trees)
do not count as evidence. The result is rank-whitened per map and range bin (over the moving speed x angle cells):
e = -ln(1 - rank / (N + 1)), exactly exponential under noise whatever the clutter's tails. Evidence: the likelihood
ratio of a Swerling-1 target of mean SNR s against exponential noise, LR = exp(e s / (1 + s)) / (1 + s) (E[LR] = 1 in
noise); cells too slow for the stream (|v| < its target_dv) are neutral (LR 1).
Recursion over the maps in time order (TRS and DM-RS interleaved, each stream with its own s and weight w):
  L_k(r, v, a) = LR_k^w * ((1 - q) + q * K[L_{k-1}](r + 2 v dt, v, a)),   q = survival^(dt / 0.16 s)
K: one-cell moves in range, speed and angle with probabilities p_r, p_v, p_a (stay otherwise); the range drifts as
the speed says (NR-Scope's convention R = R0 - 2 v t). In noise E[LR] = 1, so E[L] stays 1 (no drift); along a target
L grows geometrically. A gap longer than gap_s restarts L at 1. Detections of map k: local maxima over range x angle
of max_v log L_k above threshold (at most max_dets), each with its argmax speed. Written as a capture in NR-Scope's
map format (markers verdict 0 + AoA entries), so the paper tracker and scorer read it.

tbd.output "map_peaks": no recursion, each map's own (whitened) cube, its peaks over map_threshold (z) as detections:
a replacement for NR-Scope's chain-0 detector on a stream (tbd_track.py --as-rx2 fuses it with NR-Scope's TRS).
Tracks straight from the evidence (tbd.output "paths"): an offline Viterbi over the maps with an off state and an on
state per range x speed x angle cell (|v| <= path_v_max_ms), per map evidence the cell's log LR minus on_cost, one-cell
moves per map at move_cost (range, speed, angle) per cell, birth_cost / death_cost to switch on / off. The best path is taken, its neighbourhood (+-supp
cells) suppressed along its on-segments, and the next found, up to max_paths; segments shorter than min_maps or with
an evidence sum under min_score are dropped. Written as tbd_paths.json (per map: [id, excess, angle, speed]), which
score.py turns into the tracks the paper scorer scores (tracker.external_tracks).

  run(cfg, caps, cal, out_dir) -> dict(maps, dets_per_map, score quantiles)"""
import os, sys
import numpy as np
from scipy.ndimage import maximum_filter, convolve1d
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cube

META = 21


def evidence(Z, snr):
    """log LR per cell from the rank-whitened score"""
    nr, nv, na = Z.shape; mv = (Z > 0).any((0, 2)); X = Z[:, mv, :].reshape(nr, -1); n = X.shape[1]
    rk = np.argsort(np.argsort(X, 1), 1) + 1; e = -np.log1p(-rk / (n + 1.0))
    out = np.zeros(Z.shape, np.float32); out[:, mv, :] = (e * (snr / (1 + snr)) - np.log1p(snr)).reshape(nr, mv.sum(), na)
    return out


def normalise(t, Z, good, tau):
    """each cell over its slow background (exponential average of the stream's past maps, floor 1)"""
    bg = None; last = None; out = []
    for i in range(len(t)):
        if not good[i]:
            out.append(None); continue
        z = Z[i].astype(np.float32)
        if bg is None:
            bg = np.ones_like(z)
        out.append(z / np.maximum(bg, 1.0))
        a = 1 - np.exp(-(t[i] - last) / tau) if last is not None else 0.1
        bg = (1 - a) * bg + a * z; last = t[i]
    return out


def drift(S, v, dt, r_step):
    """S(r, v, a) moved along range by each speed's drift over dt (R = R0 - 2 v t)"""
    out = np.zeros_like(S); nr = S.shape[0]; idx = np.arange(nr)
    for j, vj in enumerate(v):
        src = idx + 2 * vj * dt / r_step; i0 = np.floor(src).astype(int); f = (src - i0)[:, None]
        ok0, ok1 = (i0 >= 0) & (i0 < nr), (i0 + 1 >= 0) & (i0 + 1 < nr)
        a = np.where(ok0[:, None], S[np.clip(i0, 0, nr - 1), j], 0); b = np.where(ok1[:, None], S[np.clip(i0 + 1, 0, nr - 1), j], 0)
        out[:, j] = (1 - f) * a + f * b
    return out


def move(L, T):
    """K: one-cell moves in range, speed, angle (axes 0, 1, 2), the rest stays"""
    for ax, p in ((0, T["p_r"]), (1, T["p_v"]), (2, T["p_a"])):
        L = convolve1d(L, np.array([p, 1 - 2 * p, p], np.float32), axis=ax, mode="constant")
    return L


def recursion(seq, T, r, v):
    """seq: [(t, log LR cube or None, weight)] in time order -> per map (max_v log L (nr, na), argmax speed index)"""
    G = None; last = None; out = []
    for t, l, w in seq:
        if l is None:
            out.append(None); continue
        if G is None or t - last > T["gap_s"]:
            prior = np.zeros_like(l)
        else:
            dt = t - last; q = T["survival"] ** (dt / 0.16); M = float(G.max())
            K = move(drift(np.exp(G - M), v, dt, T["r_step_m"]), T)
            prior = np.logaddexp(np.log1p(-q), np.log(q) + M + np.log(np.maximum(K, 1e-30))).astype(np.float32)
        G = w * l + prior; last = t
        out.append((G.max(1), G.argmax(1)))
    return out


def shifted(X, off, fill):
    """X moved by the integer offset vector off (out-of-range cells = fill)"""
    Y = np.full_like(X, fill); src, dst = [], []
    for o, n in zip(off, X.shape):
        src.append(slice(max(-o, 0), n - max(o, 0))); dst.append(slice(max(o, 0), n - max(-o, 0)))
    Y[tuple(dst)] = X[tuple(src)]
    return Y


def viterbi(E, alive, T):
    """E (n, *cells) per map evidence, alive (n,) False after a gap -> on cell per map (flat index) or -1"""
    n, shp = E.shape[0], E.shape[1:]; costs = T["move_cost"]; B, Dc = T["birth_cost"], T["death_cost"]
    offs = [np.array(o) for o in np.ndindex(*(3,) * len(shp))]; offs = [o - 1 for o in offs]
    flat = np.arange(int(np.prod(shp))).reshape(shp)
    on = np.full(shp, -1e18); off = 0.0; bk_on = np.zeros((n, flat.size), np.int32); bk_off = np.zeros(n, np.int32)
    for k in range(n):
        if not alive[k]:
            on[:] = -1e18
        best = np.full(shp, -1e18); src = np.full(shp, -1, np.int64)
        for o in offs:
            c = shifted(on, o, -1e18) - float(np.dot(np.abs(o), costs)); u = c > best
            best[u] = c[u]; src[u] = shifted(flat, o, -1)[u]
        born = off - B; u = born > best; best[u] = born; src[u] = -1
        leave = on.max() - Dc; bk_off[k] = int(on.argmax()) if leave > off else -1; off = max(off, leave)
        on = best + E[k] - T["on_cost"]; bk_on[k] = src.ravel()
    cells = np.full(n, -1); st = int(on.argmax()) if on.max() - Dc > off else -1
    for k in range(n - 1, -1, -1):
        cells[k] = st
        st = bk_on[k, st] if st >= 0 else bk_off[k]
    return cells


def paths(seq, T, r, v, s):
    """seq: [(t, log LR cube or None, weight)] -> per map [[id, excess, angle_deg, speed]]; state range x speed x angle"""
    n = len(seq); vk = np.abs(v) <= T["path_v_max_ms"]; vs = v[vk]
    E = np.full((n, len(r), vk.sum(), len(s)), 0.0, np.float32); alive = np.ones(n, bool)
    for k, (t, l, w) in enumerate(seq):
        if l is not None:
            E[k] = w * l[:, vk, :]
        if k and t - seq[k - 1][0] > T["gap_s"]:
            alive[k] = False
    out = [[] for _ in range(n)]; tid = 0; Ew = E.copy(); sp = T["supp"]; shp = E.shape[1:]
    for _ in range(T["max_paths"]):
        cells = viterbi(Ew, alive, T); on = cells >= 0; k = 0
        while k < n:
            if not on[k]:
                k += 1; continue
            j = k
            while j < n and on[j] and (j == k or alive[j]):
                j += 1
            seg = range(k, j); sc = sum(float(E[q].ravel()[cells[q]] - T["on_cost"]) for q in seg)
            if j - k >= T["min_maps"] and sc >= T["min_score"]:
                tid += 1
                for q in seg:
                    ir, iv, ia = np.unravel_index(int(cells[q]), shp); out[q].append([tid, float(r[ir]), float(np.degrees(np.arcsin(s[ia]))), float(vs[iv])])
            for q in seg:
                ir, iv, ia = np.unravel_index(int(cells[q]), shp); Ew[q, max(0, ir - sp):ir + sp + 1, :, max(0, ia - sp):ia + sp + 1] = -1e3
            k = j
        if not on.any():
            break
    return out


def detections(Sv, T, r, v, s):
    Sm, iv = Sv; peak = (Sm == maximum_filter(Sm, size=(T["nms_r_cells"], T["nms_a_cells"]), mode="constant")) & (Sm >= T["threshold"])
    ir, ia = np.nonzero(peak); o = np.argsort(-Sm[ir, ia])[:T["max_dets"]]
    return [(r[ir[k]], v[iv[ir[k], ia[k]]], np.degrees(np.arcsin(s[ia[k]])), float(Sm[ir[k], ia[k]])) for k in o]


def row_with(line, dets):
    f = line.rstrip("\n").split(","); nb, nf, nm, na = int(f[12]), int(f[13]), int(f[14]), int(f[15])
    off = float(f[16]) * float(f[9]); f[14] = str(len(dets)); f[15] = str(len(dets))
    mk = ["%.3f,%.4f,%.2f,0" % (R + off, V, 10 * np.log10(max(sc, 1e-3))) for R, V, A, sc in dets]
    ao = ["%.3f,%.4f,%.2f,%.2f,6" % (R + off, V, A, 10 * np.log10(max(sc, 1e-3))) for R, V, A, sc in dets]
    return ",".join(f[:META + nb * nf] + mk + ao) + "\n"


def run(cfg, caps, cal, out_dir):
    T = cfg["tbd"]; seq, rows = [], []
    for st in T["streams"]:
        d = cube.cubes(cfg, st, caps[st], cal); P = T["per_stream"][st]
        if T["output"] == "map_peaks":
            for i in range(len(d["t"])):
                rows.append((d["t"][i], d["rows"][i], st, d["Z"][i].astype(np.float32) if d["good"][i] else None, P["weight"]))
            r, v, s = d["r"], d["v"], d["s"]; continue
        zn = normalise(d["t"], d["Z"], d["good"], T["bg_tau_s"])
        for i in range(len(d["t"])):
            rows.append((d["t"][i], d["rows"][i], st, evidence(zn[i], P["snr"]) if zn[i] is not None else None, P["weight"]))
        r, v, s = d["r"], d["v"], d["s"]
    rows.sort(key=lambda x: x[0])
    seq = [(t, l, w) for t, _, _, l, w in rows]
    if T["output"] == "map_peaks":   # no recursion: each map's own whitened cube, peaks over map_threshold (z, noise mean 1)
        os.makedirs(out_dir, exist_ok=True); head = next(l for l in open(os.path.join(caps[T["streams"][0]], "map2d.csv")) if l.startswith("#")); n_det = []
        with open(os.path.join(out_dir, "map2d.csv"), "w") as f, open(os.path.join(out_dir, "map2d.csv.time.csv"), "w") as ft:
            f.write(head); ft.write("unix_s\n")
            for t, line, st, z, _ in rows:
                dets = []
                if z is not None:
                    Zm = z.max(1); iv = z.argmax(1)
                    dets = detections((Zm, iv), dict(T, threshold=T["map_threshold"]), r, v, s)
                f.write(row_with(line, dets)); ft.write("%.4f\n" % t); n_det.append(len(dets))
        return dict(maps=len(rows), dets_per_map=float(np.mean(n_det)))
    if T["output"] == "paths":
        os.makedirs(out_dir, exist_ok=True); P = paths(seq, T, r, v, s)
        import json
        json.dump(dict(unix=[float(x[0]) for x in rows], span=[float(x[1].split(",")[8]) for x in rows], tracks=P), open(os.path.join(out_dir, "tbd_paths.json"), "w"))
        Sv = [None] * len(rows)
    else:
        Sv = recursion(seq, T, r, v)
    os.makedirs(out_dir, exist_ok=True); head = next(l for l in open(os.path.join(caps[T["streams"][0]], "map2d.csv")) if l.startswith("#"))
    n_det, top = [], []
    with open(os.path.join(out_dir, "map2d.csv"), "w") as f, open(os.path.join(out_dir, "map2d.csv.time.csv"), "w") as ft:
        f.write(head); ft.write("unix_s\n")
        for (t, line, st, _, _), sv in zip(rows, Sv):
            dets = detections(sv, T, r, v, s) if sv is not None else ([(e, vv, a, 10.0) for _, e, a, vv in P[len(n_det)]] if T["output"] == "paths" else [])
            f.write(row_with(line, dets)); ft.write("%.4f\n" % t); n_det.append(len(dets))
            if sv is not None:
                top.append(float(sv[0].max()))
    if T["output"] == "paths":
        return dict(maps=len(rows), paths=len({q[0] for p in P for q in p}), on_per_map=float(np.mean([len(p) for p in P])))
    return dict(maps=len(rows), dets_per_map=float(np.mean(n_det)), score_p50_p90_p99=np.percentile(top, [50, 90, 99]).round(1).tolist())
