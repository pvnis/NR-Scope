"""C2: Gaussian-mixture PHD filter (Vo & Ma 2006), EKF form (their Sec. IV), in the plane.
Survival p_s per second, detection p_d inside the receiver's field of view (0 outside), Poisson
clutter of `clutter` false detections per map, uniform. Measurement-driven birth (Ristic, Clark,
Vo & Vo 2012): every localised measurement of the scan adds a birth component of weight w_birth,
updated from the next scan on. Pruning (w < 1e-5), merging (Mahalanobis U = 4), at most 100
components. Estimates: components of weight > 0.5. The PHD has no track identities: labels are
propagated by component tagging (Panta, Clark & Vo 2009), a merge keeping the heaviest label.
Fused: the receivers' scans in time order, each its own update (iterated corrector, Mahler 2009)."""
import numpy as np

from . import common as C


def run(frames, q=0.35, p_s=0.9, p_d=0.9, clutter=2.0, w_birth=0.05, gate_p=0.999, prune=1e-5, U=4.0, jmax=100):
    comps, out, nid, t_prev = [], [], 1, None   # comps: [w, x, P, label]
    for t, meas, rb, sc, geo, rx in frames:
        R = C.Rz(rb, sc)
        dt = 0.0 if t_prev is None else t - t_prev
        t_prev = t
        for c in comps:
            c[0] *= p_s ** dt
            c[1], c[2] = C.predict(c[1], c[2], dt, q)
        pd = [p_d if C.in_fov(c[1], geo) else 0.0 for c in comps]
        new = [[w * (1 - d), x, P, l] for (w, x, P, l), d in zip(comps, pd)]
        for m in meas:
            z = C.zvec(m)
            upd = []
            for (w, x, P, l), d in zip(comps, pd):
                if d == 0:
                    continue
                nu, S, H = C.innov(x, P, z, geo, R)
                if float(nu @ np.linalg.solve(S, nu)) > C.gate(len(z), gate_p):
                    continue
                xu, Pu = C.update(x, P, nu, S, H)
                upd.append([w * d * C.gauss(nu, S), xu, Pu, l])
            tot = C.kappa(clutter, len(z)) + sum(u[0] for u in upd)
            for u in upd:
                u[0] /= tot
            new += upd
        for m in meas:   # births for the next scans
            if m["xy"] is not None:
                x0, P0 = C.birth(m, geo)
                new.append([w_birth, x0, P0, nid])
                nid += 1
        comps = merge([c for c in new if c[0] > prune], U, jmax)
        est, seen = [], set()
        for c in sorted(comps, key=lambda c: -c[0]):
            if c[0] > 0.5:
                if c[3] in seen:   # a label shared by two targets (one parent, two measurements): a new one
                    c[3] = nid
                    nid += 1
                seen.add(c[3])
                est.append((c[3], *map(float, c[1])))
        out.append(est)
    return out


def merge(comps, U, jmax):
    comps = sorted(comps, key=lambda c: -c[0])
    out = []
    while comps:
        w0, x0, P0, l0 = comps[0]
        Pi = np.linalg.inv(P0)
        near = [c for c in comps if float((c[1] - x0) @ Pi @ (c[1] - x0)) <= U]
        ids = {id(c) for c in near}
        comps = [c for c in comps if id(c) not in ids]
        w = sum(c[0] for c in near)
        x, P = C.moment_match([c[0] for c in near], [c[1] for c in near], [c[2] for c in near])
        out.append([w, x, P, l0])
    return out[:jmax]
