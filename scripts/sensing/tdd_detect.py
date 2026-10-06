# SPDX-License-Identifier: LicenseRef-CSSL-1.0
"""Python port of the TDD target detector, nrscope/src/libs/sensing/nr_ue_tdd_detect.c,
plus the +-v mirror rejection that nrscope_sensing.c applies after it.

Method: Henninger et al., "Target Detection for ISAC with TDD Transmission",
arXiv:2504.19260v2. The TDD slow-time sampling paints every target with a comb of
replicas at multiples of 1/T_TDD. A CFAR candidate is tested by coherently subtracting
its full modelled response (replicas included, from the actual sample times) and
checking that the replica positions lose power: a real target takes its comb down with
it, a replica does not. See the C header for the full discussion; every function here
mirrors the C one of the same name, with the same conventions:
  delay:   IDFT, X[q] = sum_p x[p] exp(+j 2 pi p q / N)
  Doppler: C[f] = sum_i h_i exp(-j 2 pi f t_i)

Used by replay_samples.py on the same conditioned rx0 samples the live pipeline feeds it.
"""
import math

import numpy as np

# nr_ue_tdd_detect.h
NR_TDD_MAX_TARGETS = 16
NR_TDD_MAX_CANDIDATES = 128
NR_TDD_MAX_REJECTED = 32
NR_TDD_MIRROR_DB = 6.0

VERDICT_TARGET, VERDICT_REPLICA, VERDICT_DUPLICATE, VERDICT_UNCERTAIN = 0, 1, 2, 3


def cfg_default():
    """nr_tdd_cfg_default()."""
    return dict(pfa=1e-6, cfar_guard_bin=2, cfar_guard_freq=1, cfar_train_bin=4,
                cfar_train_freq=8, gamma=0.2, n_sidelobes=1, ell_bin=2.0, ell_freq=3.0,
                focus_zoom=16, focus_span=2, max_targets=8, max_iter=64)


def _rot(turns):
    """cos/sin of -2 pi turns, phase reduced to one turn first (nr_tdd_rot)."""
    frac = turns - np.trunc(turns)
    return np.exp(-2j * np.pi * frac)


class Obs:
    """nr_tdd_obs_t. h: [n_snap, n_bins] complex conditioned samples (as the C reads
    them, float32), t: [n_snap] seconds, win: frequency taper of length win_len,
    trend_q: [n_q, n_snap] orthonormal removed subspace or None."""

    def __init__(self, h, t, idft_size, win_len, win_start, win, m_per_bin,
                 n_freq, f_max_hz, lambda_m, t_tdd_s, trend_q, slow_filter=None):
        self.h = np.asarray(h, dtype=np.complex64).astype(np.complex128)
        self.t = np.asarray(t, dtype=np.float64)
        self.n_snap, self.n_bins = self.h.shape
        self.idft_size, self.win_len, self.win_start = idft_size, win_len, win_start
        self.win = None if win is None else np.asarray(win, dtype=np.float32).astype(np.float64)
        self.m_per_bin, self.n_freq, self.f_max_hz = m_per_bin, n_freq, f_max_hz
        self.lambda_m, self.t_tdd_s = lambda_m, t_tdd_s
        self.trend_q = trend_q if trend_q is not None and len(trend_q) else None
        # any other linear slow-time filter the samples went through (e.g. the sliding ECA)
        self.slow_filter = slow_filter
        self.df = 2.0 * f_max_hz / (n_freq - 1) if n_freq > 1 else 0.0
        self.freqs = -f_max_hz + self.df * np.arange(n_freq)

    def freq_of(self, j):
        return -self.f_max_hz + j * self.df if self.n_freq > 1 else 0.0

    def index_of(self, f_hz):
        return (f_hz + self.f_max_hz) / self.df if self.df > 0 else 0.0


# ---------------------------------------------------------------- PSFs

def filtered_tone(obs, f0_hz):
    """exp(+j 2 pi f0 t) minus its projection on the slow-trend subspace."""
    d = np.conj(_rot(f0_hz * obs.t))
    if obs.trend_q is not None:
        for q in obs.trend_q:
            d = d - (q @ d) * q
    if obs.slow_filter is not None:
        d = obs.slow_filter(d)
    return d


def psf_doppler_filtered(obs, d, f_hz):
    """(1/n) sum_i d_i exp(-j 2 pi f t_i), f scalar or array."""
    f = np.atleast_1d(np.asarray(f_hz, dtype=np.float64))
    E = _rot(np.outer(f, obs.t))                       # (nf, n)
    out = (E @ d) / obs.n_snap if obs.n_snap > 0 else np.zeros(len(f), complex)
    return out if np.ndim(f_hz) else out[0]


def psf_range(obs, u_bins):
    """W_R(u) = exp(j 2 pi a u / N) sum_j win[j] exp(j 2 pi j u / N) / sum(win)."""
    u = np.atleast_1d(np.asarray(u_bins, dtype=np.float64))
    j = np.arange(obs.win_len)
    w = obs.win if obs.win is not None else np.ones(obs.win_len)
    S = np.conj(_rot(np.outer(u, j) / obs.idft_size)) @ w     # sum w exp(+j...)
    a = np.conj(_rot(obs.win_start * u / obs.idft_size))
    wsum = w.sum()
    out = S * a / wsum if wsum > 0 else np.zeros(len(u), complex)
    return out if np.ndim(u_bins) else out[0]


# ---------------------------------------------------------------- map stages

def periodogram(obs):
    """Unnormalised complex map C[b, f] = sum_i h[i, b] exp(-j 2 pi f t_i), stored as
    float32 complex like the C nr_tdd_cmap_t."""
    E = _rot(np.outer(obs.freqs, obs.t))               # (nf, n)
    return (obs.h.T @ E.T).astype(np.complex64)        # (n_bins, nf)


def cell_at(obs, b, f_hz):
    return obs.h[:, b] @ _rot(f_hz * obs.t)


def cfar_power(power, cfg, max_out=NR_TDD_MAX_CANDIDATES):
    """2D CA-CFAR on local maxima, strongest first (nr_tdd_cfar_power). Returns lists
    (bins, freqs, snr_dB)."""
    P = np.asarray(power, dtype=np.float64)
    nb, nf = P.shape
    gb, gf = cfg["cfar_guard_bin"], cfg["cfar_guard_freq"]
    hb, hf = gb + cfg["cfar_train_bin"], gf + cfg["cfar_train_freq"]

    # local maximum: no 8-neighbour strictly larger
    pad = np.pad(P, 1, constant_values=-np.inf)
    is_max = P > 0
    for db in (-1, 0, 1):
        for df in (-1, 0, 1):
            if db or df:
                is_max &= ~(pad[1 + db:1 + db + nb, 1 + df:1 + df + nf] > P)

    # clipped rectangle sums by summed-area table
    S = np.zeros((nb + 1, nf + 1))
    S[1:, 1:] = P.cumsum(0).cumsum(1)

    def rect(b, f, rb, rf):
        b0, b1 = np.maximum(b - rb, 0), np.minimum(b + rb, nb - 1) + 1
        f0, f1 = np.maximum(f - rf, 0), np.minimum(f + rf, nf - 1) + 1
        return S[b1, f1] - S[b0, f1] - S[b1, f0] + S[b0, f0], (b1 - b0) * (f1 - f0)

    bs, fs = np.nonzero(is_max)                        # scan order: b outer, f inner
    so, no = rect(bs, fs, hb, hf)
    sg, ng = rect(bs, fs, gb, gf)
    total, nt = so - sg, no - ng
    p = P[bs, fs]
    ok = nt >= 8
    noise = np.where(ok, total / np.maximum(nt, 1), 0.0)
    ok &= noise > 0
    alpha = np.where(ok, nt * (cfg["pfa"] ** (-1.0 / np.maximum(nt, 1)) - 1.0), np.inf)
    ok &= p > alpha * noise
    bs, fs, snr = bs[ok], fs[ok], 10.0 * np.log10(p[ok] / noise[ok])
    order = np.argsort(-snr, kind="stable")[:max_out]   # ties keep scan order, as the C insert
    return list(bs[order]), list(fs[order]), list(snr[order])


def focus(obs, cfg, bin0, f0):
    """Off-grid Doppler (fine direct transform), parabolic delay, model amplitude."""
    f_centre = obs.freq_of(f0)
    zoom = max(cfg["focus_zoom"], 1)
    span = max(cfg["focus_span"], 1)
    step = obs.df / zoom
    fgrid = f_centre + step * np.arange(-span * zoom, span * zoom + 1)
    fgrid = fgrid[(fgrid >= -obs.f_max_hz) & (fgrid <= obs.f_max_hz)]
    vals = obs.h[:, bin0] @ _rot(np.outer(obs.t, fgrid))
    best_f = float(fgrid[int(np.argmax(np.abs(vals) ** 2))]) if len(fgrid) else f_centre

    d_hat = float(bin0)
    if 1 <= bin0 and bin0 + 1 < obs.n_bins:
        y = [math.log(max(abs(cell_at(obs, bin0 - 1 + k, best_f)) ** 2, 1e-30)) for k in range(3)]
        den = y[0] - 2.0 * y[1] + y[2]
        if abs(den) > 1e-12:
            d_hat = bin0 + min(0.5, max(-0.5, 0.5 * (y[0] - y[2]) / den))

    c = cell_at(obs, bin0, best_f)
    w = psf_range(obs, bin0 - d_hat) * psf_doppler_filtered(obs, filtered_tone(obs, best_f), best_f)
    amp = c / w if abs(w) ** 2 > 1e-12 else c
    return dict(bin=d_hat, range_m=d_hat * obs.m_per_bin, f_hz=best_f,
                speed_ms=best_f * obs.lambda_m / 2.0, amp=amp, snr_dB=0.0, iteration=0,
                verdict=VERDICT_TARGET)


def psf_subtract(obs, tgt, cmap):
    """Coherent removal of amp * W_R(b - d) * W_D'(f; f0), replicas included."""
    rr = psf_range(obs, np.arange(obs.n_bins) - tgt["bin"])
    dd = psf_doppler_filtered(obs, filtered_tone(obs, tgt["f_hz"]), obs.freqs)
    cmap -= (tgt["amp"] * np.outer(rr, dd)).astype(np.complex64)


def _ellipse_mean(cmap, bin_c, freq_c, eb, ef):
    nb, nf = cmap.shape
    b = np.arange(math.floor(bin_c - eb), math.ceil(bin_c + eb) + 1)
    f = np.arange(math.floor(freq_c - ef), math.ceil(freq_c + ef) + 1)
    b, f = b[(b >= 0) & (b < nb)], f[(f >= 0) & (f < nf)]
    if not len(b) or not len(f):
        return None
    B, F = np.meshgrid(b, f, indexing="ij")
    inside = ((B - bin_c) / eb) ** 2 + ((F - freq_c) / ef) ** 2 <= 1.0
    if not inside.any():
        return None
    p = np.abs(cmap[B[inside], F[inside]].astype(np.complex128)) ** 2
    return float(p.mean())


def sidelobe_check(before, after, obs, cfg, tgt):
    """Eq. (23): every in-grid replica must drop by at least gamma."""
    if obs.t_tdd_s <= 0:
        return True
    df_rep = 1.0 / obs.t_tdd_s
    n_tested = 0
    for k in range(1, cfg["n_sidelobes"] + 1):
        for sign in (-1, 1):
            f_k = tgt["f_hz"] + sign * k * df_rep
            if f_k < -obs.f_max_hz or f_k > obs.f_max_hz:
                continue
            j_k = obs.index_of(f_k)
            pb = _ellipse_mean(before, tgt["bin"], j_k, cfg["ell_bin"], cfg["ell_freq"])
            pa = _ellipse_mean(after, tgt["bin"], j_k, cfg["ell_bin"], cfg["ell_freq"])
            if pb is None or pa is None:
                continue
            n_tested += 1
            if pa > (1.0 - cfg["gamma"]) * pb:
                return False
    return n_tested > 0


def detect(obs, cfg=None, max_out=NR_TDD_MAX_TARGETS):
    """Algorithm 1: periodogram, CFAR, CLEAN loop of focus / subtract / check.
    Returns (targets, rejected); rejected holds at most NR_TDD_MAX_REJECTED entries."""
    cfg = cfg or cfg_default()
    if obs.n_snap < 8 or obs.n_bins < 3 or obs.n_freq < 3:
        return [], []
    work = periodogram(obs)
    cb, cf, cs = cfar_power(np.abs(work.astype(np.complex128)) ** 2, cfg)
    out, rej = [], []

    def note(b, f, snr, it, verdict):
        if len(rej) < NR_TDD_MAX_REJECTED:
            fr = obs.freq_of(f)
            rej.append(dict(bin=float(b), range_m=b * obs.m_per_bin, f_hz=fr,
                            speed_ms=fr * obs.lambda_m / 2.0, snr_dB=snr, iteration=it,
                            verdict=verdict))

    p = it = 0
    while p < len(cb) and len(out) < max_out and len(out) < cfg["max_targets"] and it < cfg["max_iter"]:
        it += 1
        dup = any(abs(cb[p] - o["bin"]) <= cfg["ell_bin"]
                  and abs(obs.freq_of(cf[p]) - o["f_hz"]) <= cfg["ell_freq"] * obs.df for o in out)
        if dup:
            note(cb[p], cf[p], cs[p], it, VERDICT_DUPLICATE)
            p += 1
            continue
        cand = focus(obs, cfg, cb[p], cf[p])
        cand["snr_dB"], cand["iteration"] = cs[p], it
        saved = work.copy()
        psf_subtract(obs, cand, work)
        if sidelobe_check(saved, work, obs, cfg, cand):
            out.append(cand)
            cb, cf, cs = cfar_power(np.abs(work.astype(np.complex128)) ** 2, cfg)
            p = 0
        else:
            note(cb[p], cf[p], cs[p], it, VERDICT_REPLICA)
            work = saved
            p += 1
    return out, rej


def _lround(x):
    """C lround: halves away from zero (Python round() goes to even)."""
    return int(math.floor(x + 0.5)) if x >= 0 else -int(math.floor(-x + 0.5))


def mirror_reject(targets, rejected, power, f_max_hz):
    """Set aside a target whose mirror cell at -f on the same range is within
    NR_TDD_MIRROR_DB of it (nrscope_sensing.c). Each side is the strongest cell within
    one bin and one Doppler cell; a target within a couple of cells of 0 Hz is kept."""
    nb, nf = power.shape
    df = 2.0 * f_max_hz / (nf - 1) if nf > 1 else 0.0
    kept = []
    for t in targets:
        b0 = _lround(t["bin"])
        f0 = _lround((t["f_hz"] + f_max_hz) / df) if df > 0 else -1
        fm = nf - 1 - f0
        uncertain = False
        if 0 <= b0 < nb and 0 <= f0 < nf and abs(f0 - fm) > 2:
            bs = slice(max(b0 - 1, 0), min(b0 + 1, nb - 1) + 1)
            p_own = power[bs, max(f0 - 1, 0):min(f0 + 1, nf - 1) + 1].max()
            p_mir = power[bs, max(fm - 1, 0):min(fm + 1, nf - 1) + 1].max()
            uncertain = p_mir > 0 and p_own < p_mir * 10.0 ** (NR_TDD_MIRROR_DB / 10.0)
        if uncertain:
            if len(rejected) < NR_TDD_MAX_REJECTED:
                rejected.append(dict(t, verdict=VERDICT_UNCERTAIN))
        else:
            kept.append(t)
    return kept, rejected
