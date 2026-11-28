"""TRS cleaning test (walk1): (a) today's static removal (+-0.25 s moving mean) vs (b) per-snapshot fitted static removed first
(R_i - c_i Rbar, c_i = <Rbar, R_i>/|Rbar|^2; or with the top-r subspace of the snapshots), then the moving mean.
Per variant: residual/thermal, and the walker range profile (walker cell / floor at +30..+60 m, with the reversed-Doppler control).
usage: trs_clean.py <chunks> <offset>"""
import sys, os, io, json, contextlib, numpy as np
KS = [int(x) for x in sys.argv[1].split(",")]; OFF = float(sys.argv[2]); HP = 0.25; WIN = 0.2
CHD = "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks"; GD = "/home/sens/GIOVANNI/real/cubotron_24_11/lidar_24_11_one_target_back_and_forth/grid"
Dn = np.load("/mnt/mydrive/field/ana_walk1/sync_walk1.npz"); TR = Dn["tr"]; t_rec = float(Dn["t_rec"]); LAM = 3e8 / 3.75e9; DF = 30e3
def walker_at(te):
    for tid in (1, 4):
        m = TR[:, 0] == tid; tk = TR[m, 1]
        if m.any() and tk[0] <= te <= tk[-1]: return np.interp(te, tk, TR[m, 2])
def hp(T, V):
    C = np.concatenate([np.zeros((1,) + V.shape[1:], complex), np.cumsum(V, 0)]); lo = np.searchsorted(T, T - HP); hi = np.searchsorted(T, T + HP, "right")
    return V - (C[hi] - C[lo]) / (hi - lo)[:, None, None]
DL = np.r_[0.0, np.arange(30, 61, 3.0)]; RES = {}
for k in KS:
    os.environ.update(USE_DMRS="0", USE_CSI="0", MINPRB="48")
    sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/e2e/tc.bin"]; g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(open(g["__file__"]).read(), "so_light", "exec"), g)
    TRS, KKT = g["TRS"], g["KK_T"]; a0 = json.load(open("%s/b%d.json" % (GD, k)))["a0"]
    t_ = np.array([x[0] for x in TRS]); R0 = TRS[0][3]; po = np.polyfit(t_, np.unwrap([np.angle(np.vdot(R0, x[3])) for x in TRS]), 3)
    V = np.array([x[3] for x in TRS]) * np.exp(-1j * np.polyval(po, t_))[:, None, None]
    T = np.array([k * 2.5 + ((x[1] - a0) % 20480) * 5e-4 + x[2] / 28e3 for x in TRS]); A = np.array([x[1] for x in TRS])
    th = np.mean([np.mean(np.abs(V[j] - V[j + 1]) ** 2) / 2 for j in range(len(V) - 1) if A[j] == A[j + 1]])
    Rb = V.mean(0); X = V.reshape(len(V), -1)
    variants = {"moving mean": V}
    c = (X @ Rb.ravel().conj()) / np.vdot(Rb.ravel(), Rb.ravel()).real; variants["fit static per snapshot"] = V - c[:, None, None] * Rb[None]
    for r in (2, 4):
        U, s, Vh = np.linalg.svd(X - X.mean(0), full_matrices=False); P = np.concatenate([Rb.ravel()[None] / np.linalg.norm(Rb), Vh[:r]])
        Q, _ = np.linalg.qr(P.T); Xc = X - (X @ Q.conj()) @ Q.T; variants["fit static + top-%d subspace per snapshot" % r] = Xc.reshape(V.shape)
    kt = KKT - 1638; M = np.exp(-2j * np.pi * kt[None, :] * DF * DL[:, None] / 3e8)
    for nm, W in variants.items():
        Dt = hp(T, W); E = {}; nW = int(2.5 // WIN) + 1; acc = np.zeros((nW, len(DL) + 1, 4), complex)
        for j, t in enumerate(T):
            dr = walker_at(t_rec + OFF + t); w = int((t - k * 2.5) // WIN)
            if dr is None or not 0 <= w < nW: continue
            base = Dt[j] * np.exp(2j * np.pi * kt * DF * dr / 3e8)[None, :]
            acc[w, :len(DL)] += ((base * np.exp(2j * np.pi * dr / LAM)) @ M.conj().T).T; acc[w, len(DL)] += (base * np.exp(-2j * np.pi * dr / LAM)) @ M[0].conj()
        En = np.sum(np.abs(acc) ** 2, (0, 2)); fl = En[1:len(DL)].mean()
        r_ = RES.setdefault(nm, [0, 0, 0, 0]); r_[0] += En[0]; r_[1] += En[len(DL)]; r_[2] += fl; r_[3] += 10 * np.log10(np.mean(np.abs(Dt) ** 2) / th) / len(KS)
for nm, (tr_, fl_, fl, rt) in RES.items():
    print("%-42s residual/thermal %+.1f dB | walker cell / floor %.1f, reversed Doppler / floor %.1f" % (nm, rt, tr_ / fl, fl_ / fl))
