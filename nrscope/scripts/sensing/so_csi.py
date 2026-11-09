"""CSI-RS-based precoder correction ("CSI stitching"), shared by gs_to_nrsr.py (VAR=csi) and so_simprec.py.
For every DM-RS snapshot y (NA x 816, 3 CDM-pair values per PRB) and PRG b, with H_b the 8-port channel of the nearest
CSI-RS (NA x 4 PRBs x 8 ports, flat in a PRB): the PRG's precoder weights w_b = argmin |y_b - H_b w|, the burst's
reference beam w_ref = principal eigenvector of sum w_b w_b^H / |w_b|^2, and per (epoch, PRG) one complex number
d = sum_t <H_b(t) w_ref, y_b(t)> / sum_t |H_b(t) w_ref|^2 over the epoch's snapshots; y_b / d puts every PRG of every
epoch on the frame of the same beam w_ref, with no edge extrapolation and no anchor between epochs."""
import numpy as np

def csi_fix(A, EP, CS, NG, NA, max_dt=0.01, sub=4):
    """A: [(t, a, pm, Y)], EP: epoch per snapshot, CS: [(t, P)] with P: NA x N x 8.
    Returns (list of X like Y, kept mask per snapshot, stats)"""
    TC = np.array([c[0] for c in CS]); near = []
    for t, a, pm, Y in A:
        j = int(np.argmin(np.abs(TC - t))) if len(TC) else -1
        near.append(j if j >= 0 and abs(TC[j] - t) <= max_dt else -1)
    def Hb(j, b): return np.repeat(CS[j][1][:, 4 * b:4 * b + 4, :][:, :, None, :], 3, 2).reshape(-1, 8)
    Rw = np.zeros((8, 8), complex); res = []
    for i in range(0, len(A), sub):
        t, a, pm, Y = A[i]; j = near[i]
        if j < 0: continue
        for b in np.flatnonzero(pm):
            H = Hb(j, b); y = Y[:, 12 * b:12 * b + 12].reshape(-1)
            w, *_ = np.linalg.lstsq(H, y, rcond=None); n2 = np.vdot(w, w).real
            if n2 > 0:
                Rw += np.outer(w, w.conj()) / n2; res.append(np.sum(np.abs(y - H @ w) ** 2) / np.sum(np.abs(y) ** 2))
    ev, V = np.linalg.eigh(Rw); wref = V[:, -1]
    num = {}; den = {}
    for i, (t, a, pm, Y) in enumerate(A):
        j = near[i]
        if j < 0: continue
        for b in np.flatnonzero(pm):
            p = Hb(j, b) @ wref; y = Y[:, 12 * b:12 * b + 12].reshape(-1); k = (EP[i], b)
            num[k] = num.get(k, 0) + np.vdot(p, y); den[k] = den.get(k, 0) + np.vdot(p, p).real
    out = []; kept = []
    for i, (t, a, pm, Y) in enumerate(A):
        X = np.zeros_like(Y); ok = near[i] >= 0
        for b in np.flatnonzero(pm):
            k = (EP[i], b)
            if ok and k in num and den[k] > 0 and num[k] != 0:
                X[:, 12 * b:12 * b + 12] = Y[:, 12 * b:12 * b + 12] / (num[k] / den[k])
        out.append(X); kept.append(ok)
    stats = dict(span_residual=float(np.median(res)) if res else np.nan, wref_share=float(ev[-1] / max(ev.sum(), 1e-30)),
                 kept=float(np.mean(kept)) if kept else 0.0)
    return out, kept, stats


def csi_fix8(A, EP, CS, NG, NA, max_dt=0.01, ridge=1e-3, sub=4):
    """The full 8-weight version. Per (epoch, PRG) the weights w = argmin sum_t |y_b(t) - H_b(t) w|^2 (all the epoch's
    snapshots, H_b(t) from the CSI-RS nearest each one, small ridge); the burst's virtual beam w_ref (principal
    direction of the w's); each snapshot rebuilt as what w_ref would have given, the minimum-norm update of H that
    agrees with the measurement:  x_b(t) = H_b w_ref + (y_b(t) - H_b w) (w^H w_ref) / |w|^2."""
    TC = np.array([c[0] for c in CS]); near = []
    for t, a, pm, Y in A:
        j = int(np.argmin(np.abs(TC - t))) if len(TC) else -1
        near.append(j if j >= 0 and abs(TC[j] - t) <= max_dt else -1)
    def Hb(j, b): return np.repeat(CS[j][1][:, 4 * b:4 * b + 4, :][:, :, None, :], 3, 2).reshape(-1, 8)
    G = {}; r = {}
    for i, (t, a, pm, Y) in enumerate(A):
        j = near[i]
        if j < 0: continue
        for b in np.flatnonzero(pm):
            H = Hb(j, b); y = Y[:, 12 * b:12 * b + 12].reshape(-1); k = (EP[i], b)
            G[k] = G.get(k, 0) + H.conj().T @ H; r[k] = r.get(k, 0) + H.conj().T @ y
    W = {}
    for k in G:
        lam = ridge * np.trace(G[k]).real / 8
        W[k] = np.linalg.solve(G[k] + lam * np.eye(8), r[k])
    Rw = sum(np.outer(w, w.conj()) / max(np.vdot(w, w).real, 1e-30) for w in W.values())
    ev, V = np.linalg.eigh(Rw); wref = V[:, -1]
    out = []; kept = []; res = []
    for i, (t, a, pm, Y) in enumerate(A):
        X = np.zeros_like(Y); j = near[i]; ok = j >= 0
        for b in np.flatnonzero(pm):
            k = (EP[i], b)
            if not ok or k not in W: continue
            H = Hb(j, b); w = W[k]; y = Y[:, 12 * b:12 * b + 12].reshape(-1); n2 = np.vdot(w, w).real
            if n2 <= 0: continue
            e = y - H @ w; res.append(np.sum(np.abs(e) ** 2) / np.sum(np.abs(y) ** 2))
            X[:, 12 * b:12 * b + 12] = (H @ wref + e * (np.vdot(w, wref) / n2)).reshape(NA, 12)
        out.append(X); kept.append(ok)
    return out, kept, dict(span_residual=float(np.median(res)) if res else np.nan, wref_share=float(ev[-1] / max(ev.sum(), 1e-30)),
                           kept=float(np.mean(kept)) if kept else 0.0)


def dejump(Y, pm):
    """Remove NR-Scope's timing-window jumps from one DM-RS snapshot (NA x 816, 3 values per PRB, 4 subcarriers apart):
    the phase slope across the band, measured only between neighbouring values inside a PRG (a precoder jump at a PRG
    boundary does not enter), dominated by the direct path; the snapshot is deramped so the direct path sits at zero
    delay. The same for every snapshot, so window jumps (whole-sample delay steps) disappear."""
    m = np.repeat(pm, 12); idx = np.arange(Y.shape[1]); inside = (idx[1:] % 12) != 0          # pair (i-1, i) in one PRG
    ok = m[1:] & m[:-1] & inside
    z = np.sum((Y[:, 1:] * np.conj(Y[:, :-1]))[:, ok])
    s = np.angle(z) if z != 0 else 0.0
    return Y * np.exp(-1j * s * idx)[None, :], s

def dejump_csi(P):
    """The same for an 8-port CSI-RS snapshot (NA x N PRBs x 8, one value per PRB, 12 subcarriers apart): slope per PRB
    over all antennas and ports, deramped on the DM-RS index scale (3 values per PRB, the CSI value at the PRB centre)."""
    z = np.sum(P[:, 1:, :] * np.conj(P[:, :-1, :]))
    s = (np.angle(z) / 3) if z != 0 else 0.0
    q = 3 * np.arange(P.shape[1]) + 1
    return P * np.exp(-1j * s * q)[None, :, None], s


def glue_per_prg(A, EP, X0, max_slots=4):
    """Glue precoder epochs PRG by PRG. X0: each snapshot already stitched with its epoch's template (frequency
    coherent inside the epoch). At an epoch change, between the last snapshot of the old epoch (already glued) and the
    first of the new one (<= max_slots apart, the channel nearly unchanged), one complex number per shared PRG maps the
    new onto the old; it is applied to every snapshot of the new epoch. PRGs not shared take the common number over the
    shared ones (one number, as before); with nothing shared (or too far apart) the whole epoch keeps the previous
    epoch's numbers. Removes the per-PRG phase jumps a precoder change leaves after stitching."""
    NG_ = A[0][2].size; D = np.ones(NG_, complex); out = [X0[0]]
    for i in range(1, len(A)):
        if EP[i] != EP[i - 1]:
            pm = A[i][2] & A[i - 1][2]; Dn = D.copy()
            if pm.any() and A[i][1] - A[i - 1][1] <= max_slots:
                prev = out[i - 1]; num_all = 0; den_all = 0
                for b in np.flatnonzero(pm):
                    sl = slice(12 * b, 12 * b + 12); u = prev[:, sl]; v = X0[i][:, sl]; den = np.vdot(u, u).real
                    if den > 0:
                        Dn[b] = np.vdot(u, v) / den; num_all += np.vdot(u, v); den_all += den
                if den_all > 0:
                    g = num_all / den_all
                    for b in np.flatnonzero(A[i][2] & ~pm): Dn[b] = g
            D = Dn
        X = X0[i] / np.repeat(D, 12)[None, :]; X[:, ~np.repeat(A[i][2], 12)] = 0; out.append(X)
    return out


def ref_fix(A, EP, REF, max_dt=0.01, align_delay=False):
    """Per (epoch, PRG) one complex number onto an external full-band reference (the TRS): d = sum_t <R_b, y_b(t)> /
    sum_t |R_b|^2 over the epoch's snapshots, R = the reference snapshot nearest each one (<= max_dt); y_b / d. No
    chaining across PRGs (no stitching error random walk) and no anchor between epochs. REF: [(t, R)], R: NA x 816.
    align_delay: first line the reference's delay up with each snapshot (relative phase slope inside PRGs), so the
    PRG numbers carry no delay offset between the two signals."""
    TR = np.array([r[0] for r in REF]); num = {}; den = {}; near = []
    for t, a, pm, Y in A:
        j = int(np.argmin(np.abs(TR - t))) if len(TR) else -1
        near.append(j if j >= 0 and abs(TR[j] - t) <= max_dt else -1)
    idx = np.arange(REF[0][1].shape[1]) if REF else None
    inside = ((idx[1:] % 12) != 0) if REF else None
    for i, (t, a, pm, Y) in enumerate(A):
        j = near[i]
        if j < 0: continue
        R = REF[j][1]
        if align_delay:   # line the reference's delay up with this snapshot: relative slope inside PRGs
            m = np.repeat(pm, 12); ok = m[1:] & m[:-1] & inside; P_ = Y * np.conj(R)
            z = np.sum((P_[:, 1:] * np.conj(P_[:, :-1]))[:, ok])
            if z != 0: R = R * np.exp(1j * np.angle(z) * idx)[None, :]
        for b in np.flatnonzero(pm):
            sl = slice(12 * b, 12 * b + 12); k = (EP[i], b)
            num[k] = num.get(k, 0) + np.vdot(R[:, sl], Y[:, sl]); den[k] = den.get(k, 0) + np.vdot(R[:, sl], R[:, sl]).real
    out = []; kept = []
    for i, (t, a, pm, Y) in enumerate(A):
        X = np.zeros_like(Y); ok = False
        for b in np.flatnonzero(pm):
            k = (EP[i], b)
            if k in num and den[k] > 0 and num[k] != 0:
                X[:, 12 * b:12 * b + 12] = Y[:, 12 * b:12 * b + 12] / (num[k] / den[k]); ok = True
        out.append(X); kept.append(ok)
    return out, kept


def los_to_zero(S, NF=4096):
    """Put range zero on the direct path, PER ANTENNA: the strongest peak of each antenna's average delay power profile
    over the burst (circular, parabolic refinement), removed from that antenna by a phase ramp. Per antenna because
    each antenna of the 5 Nov recording was replayed by its own NR-Scope run, with its own window position: their
    direct paths sat up to ~40 m apart. S: [(t, a, pm, X, l)], X: NA x 816 (4 subcarriers apart).
    Returns (S shifted, offsets in bins of NF over the 816-sample grid, one per antenna)."""
    P = 0
    for t, a, pm, X, l in S[::3]:
        m = np.repeat(pm, 12)
        if m.sum() < 48: continue
        idx = np.flatnonzero(m); w = np.zeros(X.shape[1]); w[idx[0]:idx[-1] + 1] = np.hanning(idx[-1] - idx[0] + 1)
        P = P + np.abs(np.fft.ifft(X * (w * m)[None, :], NF, axis=1)) ** 2
    if np.isscalar(P): return S, np.zeros(1)
    b0 = []
    for Pa in P:
        k = int(np.argmax(Pa)); y0, y1, y2 = Pa[k - 1], Pa[k], Pa[(k + 1) % NF]; d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
        v = k + d; b0.append(v - NF if v > NF / 2 else v)
    b0 = np.array(b0); ramp = np.exp(2j * np.pi * b0[:, None] * np.arange(S[0][3].shape[1])[None, :] / NF)
    return [(t, a, pm, X * ramp, l) for t, a, pm, X, l in S], b0


def load_wsh(path):
    """GSW1 records of one dump file -> {a: window_shift_total (samples)}, a = sfn*20 + slot unwrapped like gs_all.load"""
    import struct
    raw = open(path, "rb").read(); out = {}; i = 0; base = None
    while True:
        i = raw.find(b"GSW1", i)
        if i < 0: break
        sfn, sl, z, ws = struct.unpack_from("<HBBq", raw, i + 4); a = sfn * 20 + sl
        if base is None: base = a
        a = a + 20480 if a < base - 1000 else a
        out[a] = ws; i += 16
    return out

PH_PER_SAMPLE = 2 * np.pi * 120e3 / 122.88e6   # phase step per DM-RS/TRS grid index (4 subcarriers) per sample of shift

def undo_wsh(Y, shifts, sign=1.0):
    """Y: NA x M on the 4-subcarrier grid; shifts: NA window shifts (samples) for this slot. Multiplies antenna c by
    exp(sign * j * PH_PER_SAMPLE * shift_c * q), undoing the window move (sign found empirically)."""
    q = np.arange(Y.shape[1])
    return Y * np.exp(1j * sign * PH_PER_SAMPLE * np.asarray(shifts, float)[:, None] * q[None, :])

def slope_of(Y, pm):
    m = np.repeat(pm, 12); idx = np.arange(Y.shape[1]); ok = m[1:] & m[:-1] & ((idx[1:] % 12) != 0)
    z = np.sum((Y[:, 1:] * np.conj(Y[:, :-1]))[:, ok]); return np.angle(z) if z != 0 else 0.0

def remove_drift(A):
    """After the window moves are undone, the direct path drifts slowly with the clock offset: fit a line to the
    per-snapshot phase slope over time (unwrapped, robust) and remove that line only. A: [(t, a, pm, Y, ...)]."""
    t = np.array([x[0] for x in A]); s = np.unwrap(np.array([slope_of(x[3], x[2]) for x in A]))
    p = np.polyfit(t, s, 1); r = s - np.polyval(p, t); keep = np.abs(r) < 3 * np.std(r) + 1e-9
    p = np.polyfit(t[keep], s[keep], 1); q = np.arange(A[0][3].shape[1])
    out = [tuple([x[0], x[1], x[2], x[3] * np.exp(-1j * np.polyval(p, x[0]) * q)[None, :]] + list(x[4:])) for x in A]
    return out, p, float(np.std(s[keep] - np.polyval(p, t[keep])))


def block_epochs(A, block_s, SL, q, NA, NG):
    """Precoder 'epochs' as fixed time blocks of block_s seconds instead of detected epochs: every snapshot of a block
    is aligned to the block's running template (delay slope search over SL, one phase) and added in, whatever its
    coherence, so one set of per-PRG stitching gains serves the whole block (no independent per-epoch estimation
    noise at every precoder change). A: [(t, a, pm, Y)]. Returns (EP labels, {label: (sum, count)})."""
    gc = 12 * np.arange(NG) + 5.5; ESL = np.exp(-1j * np.outer(SL, gc)); EP = []; TPL = {}; t0 = A[0][0]
    for t, a, pm, Y in A:
        e = int((t - t0) // block_s); m12 = np.repeat(pm, 12)[None, :]
        if e not in TPL:
            TPL[e] = (Y * m12, pm.astype(int))
        else:
            S_, CN = TPL[e]; sh = (CN > 0) & pm
            if sh.sum() >= 2:
                r = np.sum((np.conj(S_ / np.maximum(np.repeat(CN, 12), 1)[None, :]) * Y).reshape(NA, NG, 12), axis=(0, 2))[sh]
                c = np.abs(ESL[:, sh] @ r); b = int(np.argmax(c)); s_ = SL[b]; p0 = np.angle(ESL[b, sh] @ r)
                TPL[e] = (S_ + Y * np.exp(-1j * (s_ * q + p0))[None, :] * m12, CN + pm)
            else:
                TPL[e] = (S_ + Y * m12, CN + pm)
        EP.append(e)
    return EP, TPL


def inject_raw(A, targets, t0, fc=3.75e9, ula=0.5):
    """Inject targets into raw snapshots BEFORE any processing, with sensing_offline's model (an echo is the measured
    channel delayed, Doppler shifted, attenuated): y_k <- y_k (1 + g e^{-j2pi (k - K/2) df tau} e^{-j2pi R/lambda}
    e^{j2pi ula a sin(angle)}), on measured values only. A: [(t, a, pm, Y, ...)], Y: NA x 816 on the 4-subcarrier grid
    (subcarrier k = 12 PRB + 4 n + 1). targets: [(range_m, speed_ms, level_db, angle_deg)], R(t) = range - 2 speed (t-t0).
    Returns (A with the targets, truth rows [(a, target, range, speed)])."""
    lam = 3e8 / fc; kk = 12 * (np.arange(816) // 3) + 4 * (np.arange(816) % 3) + 1; out = []; truth = []
    for x in A:
        t, a, pm, Y = x[0], x[1], x[2], x[3]; m = 1.0 + 0j
        for i, (r0, v, lv, ang) in enumerate(targets):
            R = r0 - 2 * v * (t - t0); tau = R / 3e8
            g = 10 ** (lv / 20) * np.exp(-2j * np.pi * R / lam) * np.exp(-2j * np.pi * (kk - 1638) * 30e3 * tau)
            m = m + g[None, :] * np.exp(2j * np.pi * ula * np.arange(Y.shape[0]) * np.sin(np.radians(ang)))[:, None]
            truth.append((a, i, R, v))
        out.append(tuple([t, a, pm, Y * m] + list(x[4:])))
    return out, truth

def lowrank_fill(S, win=64, rank=4, iters=30):
    """Fill the missing PRGs of each snapshot from the others, per antenna and per window of win consecutive snapshots:
    the snapshots x 816 matrix of a window is completed with a rank-limited model (alternating least squares on the
    observed entries). S: [(t, a, pm, X, l)]; returns the same list with full-band X and full PRG masks (the filled
    values are an estimate; the measured ones are kept)."""
    out = list(S); n = len(S)
    for w0 in range(0, n, win):
        idx = list(range(w0, min(n, w0 + win)))
        if len(idx) < 2 * rank: continue
        Mk = np.array([np.repeat(S[i][2], 12) for i in idx]); NA_ = S[idx[0]][3].shape[0]; Xf = [S[i][3].copy() for i in idx]
        for c in range(NA_):
            M = np.array([S[i][3][c] for i in idx]); O = Mk
            U = np.linalg.svd(np.where(O, M, 0), full_matrices=False)[0][:, :rank]; V = np.zeros((rank, M.shape[1]), complex)
            for _ in range(iters):
                for k in range(M.shape[1]):
                    o = O[:, k]
                    if o.sum() >= rank: V[:, k] = np.linalg.lstsq(U[o], M[o, k], rcond=None)[0]
                for j in range(M.shape[0]):
                    o = O[j]
                    if o.sum() >= rank: U[j] = np.linalg.lstsq(V[:, o].T, M[j, o], rcond=None)[0]
            F = U @ V
            for jj, i in enumerate(idx): Xf[jj][c] = np.where(O[jj], M[jj], F[jj])
        for jj, i in enumerate(idx):
            t, a, pm, X, l = S[i]; out[i] = (t, a, np.ones_like(pm), Xf[jj], l)
    return out
