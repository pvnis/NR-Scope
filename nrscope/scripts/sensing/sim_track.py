"""Track-level evaluation of the ways of combining signals, on the REAL occupancy of the 5 Nov recording.
A continuous geometric scene per 2.45 s burst, a 150 ms range-Doppler map every STEP (100 ms), detections with their angle of
arrival (4 receive antennas), localised and tracked with the tracker of ~/GIOVANNI/sim/track.py (run_ekf_tracker,
default settings, imported unchanged), scored against the true car positions.
Scene: gNB T = (0, 0, 15 m), 8-port ULA along x; receiver R = (100, 60, 2 m), 4-antenna ULA along x, facing -y (the
street). 20 static reflectors at random building positions (-45..-25 dB; calibrated so that PRG stitching errs ~2 deg
per edge, measured 1.6 deg on the cell). 3 cars per burst on the street (y = 38 or 42 m, z = 1 m), 4-12 m/s either way,
-22..-28 dB. Every path: delay, carrier phase (Doppler, range migration), departure and arrival angle from the geometry.
Signals and precoders as in sim_sense.py (DM-RS where the recording had them, per-UE two-beam precoders with per-PRG and
per-epoch changes, TRS fixed beam, CSI-RS 8 unprecoded ports). Variants: ideal, DMRS oracle, TRS, TRS+DMRS@TRS,
TRS+incoh (DM-RS groups of >= 4 coherent pieces), CSI-RS, and with the tablet DMRS today, DMRS today+TRS.
Per map: detections as in sim_sense (>= 12 dB over the moving floor, within 15 dB of the strongest, |nu| >= 10 Hz,
Doppler grid -500..500 Hz, 0-390 m), excess range e, bistatic speed v = f_D lambda / 2, angle by a beam scan over the 4 antennas, position by
track.localize. Scores per variant: GOSPA (c = 5 m, mean over maps), confirmed tracks never within 6 m of a car (false
tracks per burst), share of car-maps in view with a confirmed track within 6 m (coverage), and the per-map
detection numbers (cars found, false detections).
usage: sim_track.py <usemap json dir> <out.json> <bursts...>   env SC (notab | tab), MIX (0.2), SNR (25), SEED"""
import sys, os, json, numpy as np, collections as C
UM, OUTF = sys.argv[1], sys.argv[2]; BUR = [int(x) for x in sys.argv[3:]]; saved = list(sys.argv)
os.environ["CH"] = "0,1,2,3"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
sys.path.insert(0, "/home/sens/GIOVANNI/sim"); import track as TK
SC = os.environ.get("SC", "notab"); SNR = float(os.environ.get("SNR", "25")); MIX = float(os.environ.get("MIX", "0.2"))
rng = np.random.default_rng(int(os.environ.get("SEED", "5")))
NT = 8; FD = np.arange(-500, 501, 2.5); MOV = np.abs(FD) >= 10; KMAX = 160; FC0 = 3.75e9; c0 = 3e8; LAM = c0 / FC0
FK = (12 * np.repeat(np.arange(N), 3) + np.tile([1, 5, 9], N)) * 30e3
FT = (12 * np.repeat(np.arange(N), 3) + np.tile([2, 6, 10], N)) * 30e3
FCS = (12 * np.arange(N) + 3.5) * 30e3
T = np.array([0.0, 0.0, 15.0]); R = np.array([100.0, 60.0, 2.0]); BORE_DEG = -90.0
AXR = np.array([-np.sin(np.deg2rad(BORE_DEG)), np.cos(np.deg2rad(BORE_DEG)), 0.0]); AXT = np.array([1.0, 0.0, 0.0])
LTR = float(np.linalg.norm(R - T)); WIN = 0.15; STEP = float(os.environ.get("STEP", "0.1")); TH = np.deg2rad(np.arange(-80, 80.5, 0.5))
def atx_s(s): return np.exp(1j * np.pi * np.arange(NT) * s) / np.sqrt(NT)
def arx_s(s): return np.exp(1j * np.pi * np.arange(NA) * s)
AR_SCAN = np.exp(1j * np.pi * np.outer(np.sin(TH), np.arange(NA)))          # angles x antennas

def geom(p):
    """excess range, sin(aod), sin(aoa) of a point"""
    dT, dR = p - T, p - R
    return np.linalg.norm(dT) + np.linalg.norm(dR) - LTR, dT @ AXT / np.linalg.norm(dT), dR @ AXR / np.linalg.norm(dR)

class Scene:
    def __init__(s):
        s.st = []   # static paths: amp, e, s_aod, s_aoa
        uT = (R - T) / LTR; s.st.append((1.0, 0.0, uT @ AXT, (-uT) @ AXR))
        while len(s.st) < 21:
            p = np.array([rng.uniform(-150, 300), rng.uniform(-60, 160), rng.uniform(2, 20)])
            if abs(p[1] - 40) < 8: continue
            e, sd, sa = geom(p)
            if e > 300 or e < 5: continue
            s.st.append((10 ** (rng.uniform(-45, -25) / 20) * np.exp(1j * rng.uniform(0, 2 * np.pi)), e, sd, sa))
        s.cars = []
        for _ in range(3):
            v = rng.uniform(4, 12) * rng.choice([-1, 1]); x0 = rng.uniform(-60, 260)
            s.cars.append(dict(x0=x0, v=v, y=rng.choice([38.0, 42.0]), amp=10 ** (rng.uniform(-28, -22) / 20) * np.exp(1j * rng.uniform(0, 2 * np.pi))))
        s.sig = 10 ** (-SNR / 20) / np.sqrt(2)
    def car_pos(s, c, t): return np.array([c["x0"] + c["v"] * t, c["y"], 1.0])
    def paths(s, t):
        P = list(s.st)
        for c in s.cars:
            e, sd, sa = geom(s.car_pos(c, t)); P.append((c["amp"], e, sd, sa))
        return P
    def obs(s, f, t, Wk):
        P = s.paths(t); amp = np.array([p[0] for p in P]); e = np.array([p[1] for p in P])
        At = np.array([atx_s(p[2]) for p in P]); Ar = np.array([arx_s(p[3]) for p in P]).T
        tau = e / c0
        M = np.exp(-2j * np.pi * np.outer(tau, FC0 + f)) * (np.conj(At) @ Wk) * amp[:, None]
        y = Ar @ M
        return (y + s.sig * (rng.standard_normal(y.shape) + 1j * rng.standard_normal(y.shape))).astype(np.complex64)
    def truth(s, t):
        out = []
        for c in s.cars:
            p = s.car_pos(c, t); e, _, sa = geom(p)
            if e < KMAX * BIN - 5 and (p - R) @ np.array([np.cos(np.deg2rad(BORE_DEG)), np.sin(np.deg2rad(BORE_DEG)), 0]) > 0:
                out.append((p[0], p[1]))
        return out

def rd(ts, Xs, masks, t0):
    """clutter removal (mean over snapshots where present), range IFFT (Hann on own span), NUDFT over time (Hann);
    returns complex map (FD x NA x KMAX)"""
    if len(ts) < 2: return None
    Ms = np.array(masks, float); Xa = np.array(Xs, complex)
    mean = (Xa * Ms[:, None, :]).sum(0) / np.maximum(Ms.sum(0), 1)[None, :]
    Rr = []
    for X, m in zip(Xa, masks):
        idx = np.flatnonzero(m); w = np.zeros(L)
        if len(idx) < 24: Rr.append(np.zeros((NA, KMAX), complex)); continue
        w[idx[0]:idx[-1] + 1] = np.hanning(idx[-1] - idx[0] + 1)
        Rr.append(np.fft.ifft((X - mean) * (w * m)[None, :], NF, axis=1)[:, :KMAX])
    Rr = np.array(Rr); tt = np.array(ts) - t0; wt = np.hanning(len(tt) + 2)[1:-1]
    E = np.exp(-2j * np.pi * np.outer(FD, tt)) * wt[None, :]
    return np.einsum("fn,nak->fak", E, Rr)

def detect(Mc):
    """Mc: list of complex maps (several for the CSI-RS ports) -> detections [(f, k, p_db, theta)]"""
    if Mc is None or (isinstance(Mc, list) and not Mc): return []
    Ms = Mc if isinstance(Mc, list) else [Mc]
    P = sum(np.sum(np.abs(M) ** 2, 1) for M in Ms)
    pk = np.zeros_like(P, bool)
    pk[1:-1, 1:-1] = P[1:-1, 1:-1] >= np.maximum.reduce([P[:-2, 1:-1], P[2:, 1:-1], P[1:-1, :-2], P[1:-1, 2:], P[:-2, :-2], P[2:, 2:]])
    pk &= MOV[:, None]; fl = np.median(P[MOV]) + 1e-30
    ii = np.argwhere(pk); v = P[pk]
    if not len(v): return []
    keep = (v >= fl * 10 ** 1.2) & (v >= v.max() / 10 ** 1.5)
    out = []
    for (f, k), pv in zip(ii[keep], v[keep]):
        sc = sum(np.abs(AR_SCAN.conj() @ M[f, :, k]) ** 2 for M in Ms)
        out.append((int(f), int(k), float(10 * np.log10(pv / fl)), float(TH[int(np.argmax(sc))])))
    return out

def fit1(ref, x, pm):
    m = np.repeat(pm, 12); r = ref[:, m]; v = x[:, m]; den = np.vdot(r, r).real
    if den <= 0 or m.sum() < 48: return np.nan, np.nan
    c = np.vdot(r, v) / den; return c, float(np.sum(np.abs(v - c * r) ** 2) / np.sum(np.abs(v) ** 2))

VAR = ["ideal", "DMRS oracle", "TRS", "TRS+DMRS@TRS", "TRS+incoh", "CSI-RS"] + (["DMRS today", "DMRS today+TRS"] if SC == "tab" else [])
Rz = C.defaultdict(list)
import time
def mem_avail_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable"): return int(line.split()[1]) / 1e6
    return 99
for bi in BUR:
    while mem_avail_gb() < 3.0: time.sleep(20)   # a burst needs up to ~2.2 GB: never push the machine into swapping
    J = json.load(open(os.path.join(UM, "um_b%d.json" % bi))); seg = np.array(J["seg"][SC]).reshape(-1, 5)
    kind = {i: r[1] for i, r in enumerate(J["rnti"])}; slots = {s[0]: s for s in J["slots"]}
    scn = Scene(); ues = {}
    def W_of(u, a, typ):
        if u not in ues:
            pu = np.array([rng.uniform(-200, 300), rng.uniform(-100, 250), 1.5]); d = (pu - T) / np.linalg.norm(pu - T)
            s1 = d @ AXT; ues[u] = dict(b1=atx_s(s1), b2=atx_s(np.clip(s1 + rng.uniform(0.15, 0.4), -1, 1)), next=-1, typ=typ)
        U = ues[u]
        if a >= U["next"]:
            U["next"] = a + int(rng.exponential(38)) + 1
            U["g"] = (rng.standard_normal(NG) + 1j * rng.standard_normal(NG)) * (MIX if U["typ"] != "wideband" else 0) + MIX * (rng.standard_normal() + 1j * rng.standard_normal())
            U["ph"] = np.exp(1j * rng.choice([np.pi / 2, -np.pi / 2], NG)) if U["typ"] != "wideband" else np.ones(NG)
            U["sc"] = np.exp(1j * rng.uniform(0, 2 * np.pi)) * rng.uniform(0.5, 2)
        Wb = (U["b1"][:, None] + U["b2"][:, None] * U["g"][None, :]) * U["ph"][None, :]
        return np.repeat(Wb / np.linalg.norm(Wb, axis=0, keepdims=True) * U["sc"], 12, 1)
    # ---- every snapshot of the burst
    m = np.isin(seg[:, 3], [0, 1, 3, 4, 5, 6, 9]); byslot = C.defaultdict(list)
    for a, lo, hi, cat, u in seg[m]: byslot[int(a)].append((int(lo), int(hi), int(u)))
    w1 = np.repeat(np.eye(NT)[:, :1], L, 1)
    pieces = []   # [t, u, pm, Y, Yoracle]
    for a in sorted(byslot):
        byue = C.defaultdict(lambda: np.zeros(NG, bool))
        for lo, hi, u in byslot[a]:
            if kind.get(u) == "blind": u = 1000 + min(lo // 46, 5)
            byue[u][lo // 4:hi // 4 + 1] = True
        for u, pm in byue.items():
            typ = "per-PRG" if (u >= 1000 or kind.get(u) in ("tablet", "per-PRG")) else "wideband"
            Wk = W_of(u, a, typ); mk = np.repeat(pm, 12)
            for l in (2, 11):
                t = a * 0.5e-3 + l / 28e3; Y = scn.obs(FK, t, Wk); Y[:, ~mk] = 0; Yo = scn.obs(FK, t, w1); Yo[:, ~mk] = 0
                g = joint_gains(Y.astype(complex), pm); X = (Y / np.repeat(g, 12)[None, :]).astype(np.complex64); X[:, ~mk] = 0
                # raw Y kept only for the tablet (today's pipeline needs its epochs); memory
                pieces.append([t, u, pm, Y if kind.get(u) == "tablet" else None, Yo, X])
    pieces.sort(key=lambda x: x[0])
    TRS = [(a * 0.5e-3 + l / 28e3, scn.obs(FT, a * 0.5e-3 + l / 28e3, w1)) for a in sorted(slots) if slots[a][1] for l in (4, 8)]
    CSI = [(a * 0.5e-3 + 13 / 28e3, [scn.obs(FCS, a * 0.5e-3 + 13 / 28e3, np.repeat(np.eye(NT)[:, p:p + 1], N, 1)) for p in range(NT)]) for a in sorted(slots) if slots[a][2]]
    dl = sorted(slots)
    tP = np.array([p[0] for p in pieces]); tT = np.array([x[0] for x in TRS]); tC = np.array([x[0] for x in CSI])
    full = np.ones(L, bool)
    # self-stitched pieces, TRS anchor (once per burst)
    XP = [p[5] for p in pieces]
    viaT = []
    for i, p in enumerate(pieces):
        if not len(tT): break
        j = int(np.argmin(np.abs(tT - p[0])))
        if abs(tT[j] - p[0]) <= 0.01:
            c, e = fit1(TRS[j][1], XP[i], p[2])
            if np.isfinite(e) and e < 0.3: viaT.append((i, c))
    Rz[(SC, "joined via TRS share")].append(len(viaT) / max(len(pieces), 1))
    EPall, _ = epochs([(p[0], 0, p[2], XP[i]) for i, p in enumerate(pieces)]) if pieces else ([], None)
    # today's pipeline (reference = UE with most PRG-snapshots over the burst), once per burst
    today = None
    if SC == "tab" and pieces:
        cnt = C.Counter()
        for p in pieces:
            if p[3] is not None: cnt[p[1]] += int(p[2].sum())
        Ru = cnt.most_common(1)[0][0] if cnt else None; ia = [i for i, p in enumerate(pieces) if p[1] == Ru]
    if SC == "tab" and pieces and ia:
        A = [(pieces[i][0], 0, pieces[i][2], pieces[i][3].astype(complex)) for i in ia]
        EP, TPL = epochs(A); EG = {k: joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0) for k, (S2, C2) in TPL.items()}
        X0 = []
        for n_, (t, a, pm, Y) in enumerate(A):
            Z = Y / np.repeat(EG[EP[n_]], 12)[None, :]; Z[:, ~np.repeat(pm, 12)] = 0; X0.append(Z)
        XA = [X0[0]]; gg = 1 + 0j
        for n_ in range(1, len(A)):
            if EP[n_] != EP[n_ - 1]:
                f, nn = frame(XA[n_ - 1], A[n_ - 1][2], X0[n_], A[n_][2]); gg = f if (np.isfinite(f) and nn >= 4 and f != 0) else 1 + 0j
            XA.append(X0[n_] / gg)
        TA = np.array([pieces[i][0] for i in ia]); tod = [(pieces[i][0], XA[n_], np.repeat(pieces[i][2], 12)) for n_, i in enumerate(ia)]
        for i, p in enumerate(pieces):
            if p[1] == Ru: continue
            j = int(np.argmin(np.abs(TA - p[0])))
            if abs(TA[j] - p[0]) <= 0.01:
                c, e = fit1(XA[j], XP[i], p[2] & A[j][2])
                if np.isfinite(e) and e < 0.3: tod.append((p[0], XP[i] / c, np.repeat(p[2], 12)))
        todT = []
        for t, H in TRS:
            j = int(np.argmin(np.abs(TA - t)))
            if abs(TA[j] - t) <= 0.01:
                c, e = fit1(H, XA[j], A[j][2])
                if np.isfinite(e) and e < 0.3: todT.append((t, H * c, full))
        today = (sorted(tod, key=lambda x: x[0]), todT)
    # ---- maps every STEP
    frames = {v: [] for v in VAR}; det_stats = C.defaultdict(list); tr_truth = []
    t_end = max([p[0] for p in pieces] + [dl[-1] * 0.5e-3])
    for w0 in np.arange(dl[0] * 0.5e-3, t_end - WIN, STEP):
        w1_ = w0 + WIN; tc = w0 + WIN / 2; maps = {}
        sel = lambda arr: np.flatnonzero((arr >= w0) & (arr < w1_))
        ts = [a * 0.5e-3 + l / 28e3 for a in dl if w0 <= a * 0.5e-3 < w1_ for l in (2, 11)]
        maps["ideal"] = rd(ts, [scn.obs(FK, t, w1) for t in ts], [full] * len(ts), w0)
        ip = sel(tP); it = sel(tT); ic = sel(tC)
        maps["DMRS oracle"] = rd([pieces[i][0] for i in ip], [pieces[i][4] for i in ip], [np.repeat(pieces[i][2], 12) for i in ip], w0)
        maps["TRS"] = rd([TRS[i][0] for i in it], [TRS[i][1] for i in it], [full] * len(it), w0)
        vs = [(i, c) for i, c in viaT if w0 <= pieces[i][0] < w1_]
        maps["TRS+DMRS@TRS"] = rd([TRS[i][0] for i in it] + [pieces[i][0] for i, c in vs], [TRS[i][1] for i in it] + [XP[i] / c for i, c in vs],
                                  [full] * len(it) + [np.repeat(pieces[i][2], 12) for i, c in vs], w0)
        inc = [maps["TRS"]] if maps["TRS"] is not None else []
        grp = C.defaultdict(list)
        for i in ip: grp[EPall[i]].append(i)
        for e, idx in grp.items():
            if len(idx) >= 4:
                Mg = rd([pieces[i][0] for i in idx], [XP[i] for i in idx], [np.repeat(pieces[i][2], 12) for i in idx], w0)
                if Mg is not None: inc.append(Mg)
        maps["TRS+incoh"] = inc   # list: powers added in detect
        maps["CSI-RS"] = [rd([CSI[i][0] for i in ic], [np.repeat(CSI[i][1][p], 3, 1)[:, :L] for i in ic], [full] * len(ic), w0) for p in range(NT)] if len(ic) >= 2 else None
        if today:
            tod, todT = today
            sl_ = [x for x in tod if w0 <= x[0] < w1_]; slT = [x for x in todT if w0 <= x[0] < w1_]
            maps["DMRS today"] = rd([x[0] for x in sl_], [x[1] for x in sl_], [x[2] for x in sl_], w0)
            maps["DMRS today+TRS"] = rd([x[0] for x in sl_ + slT], [x[1] for x in sl_ + slT], [x[2] for x in sl_ + slT], w0)
        truth = scn.truth(tc); tr_truth.append(truth)
        for v in VAR:
            dets = detect(maps.get(v))
            meas = []
            for f, k, pdb, th in dets:
                e = k * BIN; vel = FD[f] * LAM / 2
                p = TK.localize(e, np.rad2deg(th), T, R, BORE_DEG)
                meas.append(dict(e=e, v=vel, theta=th if p is not None else None, xy=None if p is None else np.array(p[:2]), p_db=pdb))
            frames[v].append((tc, meas, BIN, max((1 / WIN) * LAM / 2, TK.SPEED_FLOOR_MS)))
            # per-map detection score: a car is found if a detection is within 2 bins / 4 Hz of its (e, Doppler)
            tv = []
            for c in scn.cars:
                p = scn.car_pos(c, tc); e, _, _ = geom(p)
                pp = scn.car_pos(c, tc + 1e-3); e2, _, _ = geom(pp); fd = -(e2 - e) / 1e-3 / LAM
                if e < KMAX * BIN - 5 and abs(fd) < 495: tv.append((e / BIN, (fd + 500) / 2.5))
            hit = 0; used = set()
            for kb, fb in tv:
                ok = [n_ for n_, d in enumerate(dets) if abs(d[1] - kb) <= 2 and abs(d[0] - fb) <= 2]
                if ok: hit += 1; used |= set(ok)
            if tv: det_stats[v].append((hit / len(tv), len(dets) - len(used)))
    # ---- tracking and scores
    geo = (T, R, LTR, np.deg2rad(BORE_DEG))
    for v in VAR:
        out = TK.run_ekf_tracker(frames[v], geo)
        gos = []; near_any = C.defaultdict(bool); seen = set(); cov_n = cov_d = 0
        for trks, truth in zip(out, tr_truth):
            est = [(tq[1], tq[2]) for tq in trks]
            g = TK.xy_gospa(est, truth, 5.0); gos.append(g[0])
            for tq in trks:
                seen.add(tq[0])
                if any(np.hypot(tq[1] - x, tq[2] - y) <= 6.0 for x, y in truth): near_any[tq[0]] = True
            for x, y in truth:
                cov_d += 1; cov_n += any(np.hypot(tq[1] - x, tq[2] - y) <= 6.0 for tq in trks)
        Rz[(SC, v, "gospa")].append(float(np.mean(gos)))
        Rz[(SC, v, "false tracks")].append(sum(1 for i in seen if not near_any[i]))
        Rz[(SC, v, "coverage")].append(cov_n / max(cov_d, 1))
        if det_stats[v]:
            Rz[(SC, v, "found")].append(float(np.mean([d[0] for d in det_stats[v]]))); Rz[(SC, v, "false det")].append(float(np.mean([d[1] for d in det_stats[v]])))
    print("burst %d done (%s), %d pieces, %d maps" % (bi, SC, len(pieces), len(tr_truth)), flush=True)
json.dump({"|".join(k): v for k, v in Rz.items()}, open(OUTF, "w"))
print("%-16s %8s %10s %12s %10s %8s" % ("variant (" + SC + ")", "found", "false det", "false trk/b", "coverage", "GOSPA"))
for v in VAR:
    if Rz.get((SC, v, "gospa")):
        print("%-16s %7.0f%% %10.1f %12.1f %9.0f%% %8.1f" % (v, 100 * np.mean(Rz.get((SC, v, "found"), [0])), np.mean(Rz.get((SC, v, "false det"), [0])),
              np.mean(Rz[(SC, v, "false tracks")]), 100 * np.mean(Rz[(SC, v, "coverage")]), np.mean(Rz[(SC, v, "gospa")])))
