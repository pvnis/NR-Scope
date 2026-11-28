"""Synthetic test of the wideband beam anchoring (BA_PH=wide) algorithm, outside so_light: snapshots of one UE whose
grant covers runs of PRGs; per PRG a fixed beam shape S (4 antennas x 12 values, a few paths); per transmission one
random phase; an echo with its own delay and arrival angle, Doppler; white noise. The algorithm (as in so_light) gives
Y -> (Y e^{-j theta} - S_beam) / g; the echo's survival is measured against the ideal (echo with no anchoring loss)."""
import numpy as np, collections as C
rng = np.random.default_rng(5); NG = 68; L = 12 * NG; NA = 4
kk = (12 * (np.arange(L) // 12)).astype(float) + np.tile(np.r_[0:12] // 3 * 4 + 1, NG)[:L]   # approx subcarrier of each value
DF = 30e3; LAM = 0.08
def chan(paths):   # paths: (gain, delay s, angle rad) -> NA x L
    H = np.zeros((NA, L), complex)
    for g, tau, ang in paths:
        H += g * np.exp(-2j * np.pi * kk * DF * tau)[None, :] * np.exp(1j * np.pi * np.sin(ang) * np.arange(NA))[:, None]
    return H
static_trs = chan([(1, 0, 0.3), (0.3, 0.4e-6, -0.6), (0.2, 0.9e-6, 0.9)])
static_beam = chan([(1.6, 0, 0.3), (0.15, 0.4e-6, -0.6), (0.35, 0.9e-6, 0.9)])   # the beam weights the paths differently
def echo(t, eps=10 ** (-30 / 20), tau=100e-9, ang=-0.8, v=1.0):   # physical echo: own delay, angle, Doppler
    R = 30 - 2 * v * t
    return eps * np.exp(-2j * np.pi * R / LAM) * chan([(1, tau, ang)])
T = np.arange(0, 2.4, 0.0005)[::2]; E = []
runs_ = [(2, 30), (34, 60)]
for t in T:
    pm = np.zeros(NG, bool)
    for a_, b_ in runs_: pm[a_:b_] = True
    th = rng.uniform(-np.pi, np.pi)
    Y = np.exp(1j * th) * (static_beam + echo(t) * 1.6) + 10 ** (-20 / 20) * (rng.standard_normal((NA, L)) + 1j * rng.standard_normal((NA, L))) / np.sqrt(2) * np.abs(static_beam).mean()
    E.append((t, pm, Y))
def runs(m):
    d = np.diff(np.r_[0, m.astype(int), 0]); s = np.flatnonzero(d == 1); e = np.flatnonzero(d == -1); return list(zip(s, e))
# --- the algorithm (copy of the BA_PH=wide logic)
cth = 0.97; ttl = 0.5; TP = C.defaultdict(list); TPL = {}; own = {}; THETA = {}; nb = 0
for i, (t, pm, Y) in enumerate(E):
    for s0_, e0_ in runs(pm):
        cand = {}; z = 0
        for b in range(s0_, e0_):
            y = Y[:, 12 * b:12 * b + 12].ravel(); ny = np.linalg.norm(y); TP[b] = [Tt for Tt in TP[b] if t - Tt[2] <= ttl]
            best = None; bc = cth
            for Tt in TP[b]:
                c_ = abs(np.vdot(Tt[0], y)) / (np.linalg.norm(Tt[0]) * ny + 1e-30)
                if c_ >= bc: bc = c_; best = Tt
            cand[b] = best
            if best is not None: z += np.vdot(best[0] / best[1], y)
        th = np.angle(z) if z != 0 else 0.0
        for b in range(s0_, e0_):
            y = Y[:, 12 * b:12 * b + 12].ravel() * np.exp(-1j * th); Tt = cand[b]
            if Tt is None: Tt = [np.zeros_like(y), 0, t, nb]; TPL[nb] = Tt; nb += 1; TP[b].append(Tt)
            Tt[0] = Tt[0] + y; Tt[1] += 1; Tt[2] = t; own[(i, b)] = Tt[3]; THETA[(i, b)] = th
S = {k: Tt[0].reshape(NA, 12) / Tt[1] for k, Tt in TPL.items()}
print("beams:", nb)
# output and echo survival per snapshot (on measured PRGs), against the ideal echo in TRS units: echo scaled by 1/g
gam = []
for i, (t, pm, Y) in enumerate(E):
    d = []; ideal = []
    for b in np.flatnonzero(pm):
        k = own[(i, b)]; sl = slice(12 * b, 12 * b + 12); Sb = S[k]; Rb = static_trs[:, sl]; g = np.vdot(Rb, Sb) / np.vdot(Rb, Rb).real
        y = Y[:, sl] * np.exp(-1j * THETA[(i, b)])
        d.append(((y - Sb) / g).ravel()); ideal.append((echo(t)[:, sl] * 1.6 / g).ravel())
    d = np.concatenate(d); ideal = np.concatenate(ideal)
    gam.append(np.vdot(ideal, d) / np.vdot(ideal, ideal).real)
gam = np.array(gam)
print("echo survival |gamma| median %.2f; phase spread MAD %.2f rad; coherent mean |<gamma>| / <|gamma|> = %.2f"
      % (np.median(abs(gam)), 1.4826 * np.median(np.abs(np.angle(gam / np.exp(1j * np.angle(np.mean(gam)))))), abs(np.mean(gam)) / np.mean(abs(gam))))
