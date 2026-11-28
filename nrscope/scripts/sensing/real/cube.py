"""Range x speed x angle power cubes for track-before-detect, from the channel estimates, on the maps' own windows.

For every map row of a capture (NR-Scope map2d.csv + its UTC times) the snapshots of that window (or of the
tbd.window_s seconds ending with it: longer coherent integration) are taken from the
chunk's estimate file and transformed with NR-Scope's conventions (sensing_offline.cc, nr_ue_aoa.c):
  echo  H_a[k](t) = g exp(-j 2 pi k df tau) exp(-j 2 pi R(t) / lam) exp(j pi a sin(theta)),  R = R0 - 2 v t
  range C(tau) = sum_k w_k H[k] exp(+j 2 pi k df tau), at NR-Scope's bins (excess path = bin * m_per_bin)
  speed D(v) = sum_i w_i C_i exp(-j 2 pi (2 v / lam) t_i)    (v > 0 approaching, as NR-Scope)
  angle x_a = D_a exp(-j cal_a),  P(theta) = |sum_a x_a exp(-j pi a sin(theta))|^2
        or, with tbd.whiten "map" / "chunk", whitened against the noise's spatial covariance R (the moving cells of
        the map / of every map of its chunk; the
        floor is gNB-borne distortion arriving along the direct path, dmrs_check):
        P(theta) = |a^H R^-1 x|^2 / (a^H R^-1 a)     (R + 1 % diagonal loading)
Every snapshot is first aligned to the window's static scene (as NR-Scope's nr_ue_sensing_align.c): one delay,
phase and gain per snapshot from chain 0 against the window mean, applied to all chains (keeps the inter-chain
phases). The static scene is then removed as the window's mean per chain and subcarrier. Each cube is whitened per range bin by
the median over the moving cells (exponential noise: mean = median / ln 2), so a cell holds z = power / noise.

  cubes(cfg, stream, cap_dir) -> dict(t (n,), r (nr,), v (nv,), s (na,) sin(angle), Z (n, nr, nv, na) float16, rows)"""
import json, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfg as C, capture, est_read

LAM = 299792458.0 / 3.75e9


def grids(T):
    r = np.arange(T["r_min_m"], T["r_max_m"] + 1e-9, T["r_step_m"])
    v = np.r_[-np.arange(T["v_max_ms"], T["v_min_ms"] - 1e-9, -T["v_step_ms"]), np.arange(T["v_min_ms"], T["v_max_ms"] + 1e-9, T["v_step_ms"])]
    s = np.sin(np.radians(np.arange(-T["angle_max_deg"], T["angle_max_deg"] + 1e-9, T["angle_step_deg"])))
    return r, np.sort(v), s


def rows(cap_dir):
    L = [l for l in open(os.path.join(cap_dir, "map2d.csv")) if not l.startswith("#")]
    U = [float(x) for x in open(os.path.join(cap_dir, "map2d.csv.time.csv")).read().split()[1:]]
    return L, np.array(U)


def chunk_of(cfg, line, unix, a0s):
    """(chunk k, window end slot a_end) of a map row: the chunk whose slot clock gives this UTC end time"""
    v = line.split(",", 9); a = int(v[0]) * 20 + int(v[1]); rel = unix - (capture.t_rec(cfg) + cfg["recording"]["iq_offset_s"])
    for k in (int(rel // 2.5), int(rel // 2.5) - 1, int(rel // 2.5) + 1):
        if k in a0s and abs(k * 2.5 + ((a - a0s[k]) % 20480) * 5e-4 - rel) < 2e-3:
            return k, a, float(v[8])
    return None, a, float(v[8])


def transform(X, T, a_end, span, r, v, cal):
    """the window's snapshots -> calibrated range x speed per chain D (4, nr, nv), or None"""
    a = np.floor(T / 5e-4 + 1e-6).astype(int); d = (a - a_end + 10240) % 20480 - 10240
    m = (d <= 0) & (d > -span / 5e-4 - 0.5)
    if m.sum() < 4:
        return None
    H = X[m].astype(np.complex128); t = T[m] - T[m].mean(); n = len(t)
    ok = ~np.isnan(H); H = np.where(ok, H, 0)
    cnt = ok.sum(0)
    kk = np.arange(H.shape[2])
    for _ in range(2):                                                                           # align to the static scene
        R0 = H[:, 0].sum(0) / np.maximum(cnt[0], 1); x = H[:, 0] * np.conj(R0)                   # (n, K)
        xv = x[:, ok[:, 0].all(0)] if ok[:, 0].all(0).sum() > 16 else x
        kv = kk[ok[:, 0].all(0)] if ok[:, 0].all(0).sum() > 16 else kk
        dk = np.diff(kv).mean() if len(kv) > 1 else 1.0
        sl = np.angle((xv[:, 1:] * np.conj(xv[:, :-1])).sum(1)) / dk                             # phase per subcarrier
        g = (xv * np.exp(-1j * sl[:, None] * kv[None, :])).sum(1) / max(float((np.abs(R0[kv]) ** 2).sum()), 1e-30)
        g = np.where(np.abs(g) > 1e-12, g, 1.0)
        H = H * (np.exp(-1j * sl[:, None] * kk[None, :]) / g[:, None])[:, None, :]
    mu = H.sum(0) / np.maximum(cnt, 1); H = np.where(ok, H - mu, 0)                              # static scene out
    k = np.flatnonzero(ok.any((0, 1))); wk = np.hanning(len(k) + 2)[1:-1]
    A = wk[:, None] * np.exp(2j * np.pi * k[:, None] * 30e3 * (r[None, :] / 299792458.0))       # (K, nr)
    Cr = H[:, :, k] @ A                                                                          # (n, 4, nr)
    wt = np.hanning(n + 2)[1:-1]
    B = wt[:, None] * np.exp(-2j * np.pi * (2 * v[None, :] / LAM) * t[:, None])                  # (n, nv)
    return np.einsum("iar,iv->arv", Cr, B) * np.exp(-1j * np.radians(cal))[:, None, None]        # (4, nr, nv)


def noise_cov(D, mv):
    Y = D[:, :, mv].reshape(4, -1); return Y @ Y.conj().T / Y.shape[1]


def power(D, s, mv, R=None):
    """angle stage (Bartlett, or whitened against R) and per-range-bin whitening -> z (nr, nv, na)"""
    S = np.exp(-1j * np.pi * np.arange(4)[:, None] * s[None, :])                                 # conj steering (4, na)
    if R is not None:
        R = R / (np.trace(R).real / 4) + 0.01 * np.eye(4)
        Am = S.conj(); W = np.linalg.inv(R) @ Am; W = W / np.sqrt(np.real(np.sum(Am.conj() * W, 0)))[None, :]
        S = W.conj()                                                                             # w^H x, w = R^-1 a / sqrt(a^H R^-1 a)
    P = np.abs(np.einsum("arv,ag->rvg", D, S)) ** 2                                              # (nr, nv, na)
    noise = np.median(P[:, mv, :].reshape(P.shape[0], -1), 1) / np.log(2)
    Z = P / np.maximum(noise, 1e-30)[:, None, None]; Z[:, ~mv, :] = 0
    return Z


def cubes(cfg, stream, cap_dir, cal):
    T_ = cfg["tbd"]; r, v, s = grids(T_); L, U = rows(cap_dir); rec = cfg["recording"]
    a0s = {int(re.search(r"b(\d+)", f).group(1)): json.load(open(os.path.join(rec["lidar"], "grid", f)))["a0"]
           for f in os.listdir(os.path.join(rec["lidar"], "grid")) if re.match(r"b\d+\.json$", f)}
    by = {}
    for i, (l, u) in enumerate(zip(L, U)):
        k, a, span = chunk_of(cfg, l, u, a0s)
        if k is not None:
            by.setdefault(k, []).append((i, a, span))
    Z = np.zeros((len(L), len(r), len(v), len(s)), np.float16); good = np.zeros(len(L), bool)
    est = C.est_dir(cfg, stream)
    for k, items in sorted(by.items()):
        f = os.path.join(est, "b%d.bin" % k)
        if not os.path.exists(f):
            continue
        T, X = est_read.snapshots(f); mv = np.abs(v) >= cfg["streams"][stream]["target_dv"]
        Ds = [(i, transform(X, T, a, T_["window_s"] or span, r, v, cal)) for i, a, span in items]
        Ds = [(i, D_) for i, D_ in Ds if D_ is not None]
        Rc = sum(noise_cov(D_, mv) for _, D_ in Ds) / max(len(Ds), 1) if T_["whiten"] == "chunk" and Ds else None
        for i, D_ in Ds:
            R = noise_cov(D_, mv) if T_["whiten"] == "map" else Rc
            Z[i] = np.minimum(power(D_, s, mv, R), 6e4); good[i] = True
    return dict(t=U, r=r, v=v, s=s, Z=Z, good=good, rows=L)
