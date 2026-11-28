"""Delay-sparse precoder with a free complex scalar per PRG (template frames are arbitrary per PRG):
y_b ~ s_b H_b U phi(f_b), phi_q = exp(j2pi f_b q dtau), q=-K..K. Scalar-free shape residual 1-|<y,Hw>|^2/(|y|^2|Hw|^2).
CV: 4 folds over PRGs. Controls: neighbour's own v (b-1), random w, same model with PRG frequencies shuffled.  usage: dsfit2.py <chunks>"""
import sys, numpy as np
sf = lambda y, z: 1 - abs(np.vdot(z, y)) ** 2 / (np.vdot(y, y).real * np.vdot(z, z).real + 1e-30)
def fit(HB, YB, fb, K, idx, dtau, it=25, lam=1e-3):
    q = np.arange(-K, K + 1); PHI = np.exp(2j * np.pi * fb[:, None] * q[None, :] * dtau); s = np.ones(len(fb), complex)
    for _ in range(it):
        A = np.concatenate([s[j] * np.kron(PHI[j][None, :], HB[j]) for j in idx]); y = np.concatenate([YB[j] for j in idx])
        AhA = A.conj().T @ A; u = np.linalg.solve(AhA + lam * np.trace(AhA).real / len(AhA) * np.eye(len(AhA)), A.conj().T @ y); U = u.reshape(len(q), 8).T
        for j in idx: z = HB[j] @ (U @ PHI[j]); s[j] = np.vdot(z, YB[j]) / max(np.vdot(z, z).real, 1e-30)
    return U, PHI
rng = np.random.default_rng(0)
for k in [int(x) for x in sys.argv[1].split(",")]:
    Z = np.load("/mnt/mydrive/field/ana_walk1/hy_b%d.npz" % k); HB, YB, B, Vb = Z["HB"], Z["YB"], Z["B"], Z["Vb"]; n = len(B)
    fb = B * 48 * 30e3; dtau = 1 / ((B.max() - B.min() + 1) * 48 * 30e3)
    nb = [sf(YB[j], HB[j] @ Vb[j - 1]) for j in range(1, n) if B[j] - B[j - 1] == 1]
    rnd = [sf(YB[j], HB[j] @ (rng.standard_normal(8) + 1j * rng.standard_normal(8))) for j in range(n)]
    own = [sf(YB[j], HB[j] @ Vb[j]) for j in range(n)]
    line = "chunk %d: own free v %.1f dB | neighbour's v %.1f dB | random w %.1f dB |" % (k, 10 * np.log10(np.median(own)), 10 * np.log10(np.median(nb)), 10 * np.log10(np.median(rnd)))
    for K in (0, 2, 4, 8, 16):
        out = []
        for lab, f_ in (("true f", fb), ("shuffled f", rng.permutation(fb))):
            e = []
            for f in range(4):
                tr = [j for j in range(n) if j % 4 != f]; U, PHI = fit(HB, YB, f_, K, tr, dtau)
                e += [sf(YB[j], HB[j] @ (U @ PHI[j])) for j in range(n) if j % 4 == f]
            out.append("%s %.1f" % (lab, 10 * np.log10(np.median(e))))
        line += "\n     K=%2d held-out PRGs: %s dB" % (K, ", ".join(out))
    print(line, flush=True)
