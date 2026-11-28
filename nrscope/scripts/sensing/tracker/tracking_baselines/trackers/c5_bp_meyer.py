"""C5: belief-propagation multitarget tracking for an unknown number of targets and several
sensors (Meyer, Braca, Willett & Hlawatsch, IEEE TSP 2017; tutorial: Meyer et al., Proc. IEEE
2018, Sec. VIII-IX). Potential targets (PTs) with an existence probability r. Per scan (one
receiver's map): legacy PTs predicted (survival p_s per second); one new PT per measurement;
association by the SPA iterations of Williams & Lau 2014 (n_it), with
  beta_i(m) = r pd g_im / (mu_c f_c),  beta_i(0) = 1 - r + r (1 - pd),
  xi_m(0) = 1 + mu_n f_n / (mu_c f_c),  xi_m(i) = 1;
then the legacy beliefs (missed + each measurement, moment matched: the Gaussian form of the
paper's particle beliefs) and existences, and the new PTs' r = (xi_m(0) - 1) / (xi_m(0) + sum_i phi_im).
Sensors in time order, each its own scan (the sequential schedule, eq. 47-49: new PTs of one
receiver are legacy PTs of the next). Declared when r > r_th; pruned when r < r_prune.
mu_n / (mu_c f_c) is taken constant: new targets per map mu_n, spread like the clutter."""
import numpy as np

from . import common as C


def spa(beta, beta0, xi0, n_it):
    """phi[i, m] (PT -> measurement) and nu[m, i] (measurement -> PT) messages."""
    n, M = beta.shape
    nu = np.ones((M, n))
    phi = np.zeros((n, M))
    for _ in range(n_it):
        bn = beta * nu.T
        phi = beta / (beta0[:, None] + bn.sum(1, keepdims=True) - bn)
        s = phi.sum(0)
        nu = (1.0 / (xi0[:, None] + s[:, None] - phi.T))
    return phi, nu


def run(frames, q=0.35, p_s=0.9, p_d=0.9, clutter=2.0, mu_n=0.05, n_it=20, r_th=0.5, r_prune=1e-3, gate_p=0.999):
    pts, out, nid = [], [], 1
    for t, meas, rb, sc, geo, rx in frames:
        R = C.Rz(rb, sc)
        for p in pts:
            dt = t - p["t"]
            p["r"] *= p_s ** dt
            p["x"], p["P"] = C.predict(p["x"], p["P"], dt, q)
            p["t"] = t
        zs = [C.zvec(m) for m in meas]
        n, M = len(pts), len(meas)
        pd = np.array([p_d if C.in_fov(p["x"], geo) else 0.0 for p in pts])
        r = np.array([p["r"] for p in pts])
        beta = np.zeros((n, M))
        cache = {}
        for i, p in enumerate(pts):
            if pd[i] == 0:
                continue
            for j, z in enumerate(zs):
                nu, S, H = C.innov(p["x"], p["P"], z, geo, R)
                if float(nu @ np.linalg.solve(S, nu)) < C.gate(len(z), gate_p):
                    cache[i, j] = (nu, S, H)
                    beta[i, j] = r[i] * pd[i] * C.gauss(nu, S) / C.kappa(clutter, len(z))
        beta0 = 1 - r + r * (1 - pd)
        ratio = np.array([mu_n / clutter if m["xy"] is not None else 0.0 for m in meas])
        xi0 = 1 + ratio
        if n and M:
            phi, nuv = spa(beta, beta0, xi0, n_it)
        else:
            phi, nuv = np.zeros((n, M)), np.ones((M, n))
        for i, p in enumerate(pts):
            miss = r[i] * (1 - pd[i])
            ws, xs, Ps = [miss], [p["x"]], [p["P"]]
            for j in range(M):
                if (i, j) in cache:
                    xu, Pu = C.update(p["x"], p["P"], *cache[i, j])
                    ws.append(beta[i, j] * nuv[j, i]), xs.append(xu), Ps.append(Pu)
            num = sum(ws)
            p["r"] = num / (num + 1 - r[i])
            if num > 0:
                p["x"], p["P"] = C.moment_match(ws, xs, Ps)
            p["decl"] = p["r"] > r_th
        new = []
        for j, m in enumerate(meas):
            if ratio[j] > 0:
                rn = ratio[j] / (ratio[j] + 1 + (phi[:, j].sum() if n else 0.0))
                x0, P0 = C.birth(m, geo)
                new.append(dict(id=nid, x=x0, P=P0, t=t, r=rn, decl=rn > r_th))
                nid += 1
        pts = [p for p in pts + new if p["r"] >= r_prune]
        out.append([(p["id"], *map(float, p["x"])) for p in pts if p["decl"]])
    return out
