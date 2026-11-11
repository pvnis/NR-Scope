"""Light network: one sensing stream from the TRS and the DM-RS of every UE except (optionally) the tablet, written as a
record_estimates file for NR-Scope's sensing_offline.
- TRS symbols 4/8 (slots 5/6, every 20 ms, full band, fixed beam): each TRS value (REs 2, 6, 10 of a PRB) is written on
  the DM-RS pair grid (REs 4n and 4n+2), so TRS and DM-RS share one pilot lattice and NR-Scope makes ONE map of both.
- DM-RS: the blind records (rnti 0, all 4 antennas, every PCI-scrambled DM-RS) minus the tablet's PRGs (MODE=notab), or
  plus the tablet (MODE=tab); each snapshot's PRGs put on the nearest TRS's reference (one number per snapshot and
  PRG; env ALIGN=1 delay-aligns the TRS to each snapshot first, off by default), or left out (USE_DMRS=0).
- Every antenna: NR-Scope's logged window moves undone (GSW1), straight-line drift removed, range zero on that
  antenna's direct path (measured on the TRS, applied to everything of that antenna).
- CSI-RS (USE_CSI=1): per occasion, the whole-band 8-port weights reproducing the nearest TRS, as an extra
  TRS-beam sample (CSI_DT 6 ms, kept when the fit residual < CSI_MAXRES 0.3); CSI_REF=1 also uses them as references
  for the DM-RS (nearest of TRS or CSI-RS-derived sample).
- INJ='range,speed,dB,angle;...': targets injected into every raw signal first (truth to <out>.truth.csv).
usage: so_light.py <burst dir> <out.bin>   env MODE (notab | tab), USE_DMRS (1), USE_CSI (0), CSI_REF (0), INJ, DIR"""
import sys, os, struct, numpy as np, collections as C
BD, OUTP = sys.argv[1].rstrip("/") + "/", sys.argv[2]; saved = list(sys.argv); os.environ["CH"] = "0,1,2,3"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gs_trs_lib import load_gsc, trs_despread
from so_csi import load_wsh, undo_wsh, slope_of, PH_PER_SAMPLE as PH
MODE = os.environ.get("MODE", "notab"); USE = os.environ.get("USE_DMRS", "1") == "1"; DIR = os.environ.get("DIR", "replay_blind_ch")
NSC = 3276; NFq = 4096; files = [BD + DIR + "%d/ue_dmrs.bin" % c for c in range(4)]
WS = [load_wsh(f) for f in files]
# targets injected into every RAW signal before any processing (env INJ='range,speed,dB,angle;...'), each at its own
# subcarriers, same echo model as sensing_offline; truth (one row per slot and target) to <out>.truth.csv
TG = [tuple(float(v) for v in x.split(",")) for x in os.environ["INJ"].split(";")] if os.environ.get("INJ", "") else []
T0 = [None]; TRUTH = {}
def inj(Y, t, a, kk):
    if not TG: return Y
    if T0[0] is None: T0[0] = t
    m = 1.0 + 0j
    for i, (r0, v, lv, ang) in enumerate(TG):
        R = r0 - 2 * v * (t - T0[0]); tau = R / 3e8; TRUTH[(a, i)] = (R, v)
        g = 10 ** (lv / 20) * np.exp(-2j * np.pi * R / (3e8 / 3.75e9)) * np.exp(-2j * np.pi * (kk - 1638) * 30e3 * tau)
        m = m + g[None, :] * np.exp(2j * np.pi * 0.5 * np.arange(Y.shape[0]) * np.sin(np.radians(ang)))[:, None]
    return Y * m
KK_D = 12 * (np.arange(816) // 3) + 4 * (np.arange(816) % 3) + 1; KK_T = KK_D + 1
def shifts(a): return [W.get(a, W[min(W, key=lambda b: abs(b - a))]) if W else 0 for W in WS]
# TRS
G = [load_gsc(f) for f in files]; gk = sorted(set(G[0]).intersection(*[set(g) for g in G[1:]]))
TRS = [(a * 0.5e-3 + l / 28e3, a, l, undo_wsh(inj(np.array([trs_despread(*G[c][(a, l)], a % 20, l, 2)[:N].ravel() for c in range(4)]),
                                               a * 0.5e-3 + l / 28e3, a, KK_T), shifts(a), -1.0))
       for a, l in gk if l in (4, 8)]
del G
# drift: straight line through the TRS slope (full band, fixed beam: the cleanest), per antenna
t_ = np.array([x[0] for x in TRS]); q = np.arange(L); full = np.ones(NG, bool)
qd = q - float(os.environ.get("WSH_Q0", "409.25"))   # phase ramps (delay, drift) are linear in frequency about DC
drift = []
for c in range(4):
    s = np.unwrap([slope_of(x[3][c:c + 1], full) for x in TRS]); p = np.polyfit(t_, s, 1); r = s - np.polyval(p, t_)
    k = np.abs(r) < 3 * np.std(r) + 1e-9; drift.append(np.polyfit(t_[k], s[k], 1))
def dedrift(t, Y): return Y * np.exp(-1j * np.array([np.polyval(drift[c], t) for c in range(4)])[:, None] * qd[None, :])
TRS = [(t, a, l, dedrift(t, R)) for t, a, l, R in TRS]
# range zero per antenna on the TRS direct path
P = 0
for t, a, l, R in TRS: P = P + np.abs(np.fft.ifft(R * np.hanning(L)[None, :], NFq, axis=1)) ** 2
b0 = []
for Pa in P:
    k = int(np.argmax(Pa)); y0, y1, y2 = Pa[k - 1], Pa[k], Pa[(k + 1) % NFq]; d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
    v = k + d; b0.append(v - NFq if v > NFq / 2 else v)
ramp = np.exp(2j * np.pi * np.array(b0)[:, None] * qd[None, :] / NFq)
TRS = [(t, a, l, R * ramp) for t, a, l, R in TRS]
print("TRS: %d symbols; direct path per antenna was at %s m" % (len(TRS), " ".join("%+.1f" % (b / NFq / 120e3 * 3e8) for b in b0)))
snaps = [(t, a, l, full, R, "trs") for t, a, l, R in TRS]
SIDE = [(a, l, "trs", -1, full, "ok") for t, a, l, R in TRS]   # every candidate symbol and its fate (env SIDE=<file>)
# CSI-RS as extra TRS-beam samples (env USE_CSI=1): per CSI-RS occasion (slots 2/10, symbol 13), the one set of 8 port
# weights for the whole band that best reproduces the nearest TRS (<= CSI_DT), i.e. the TRS beam through that occasion's
# ports; the CSI-RS times then carry TRS-like snapshots (more samples in time: less Doppler folding).
if os.environ.get("USE_CSI", "0") == "1":
    from gs_trs_lib import csi8_despread
    G = [load_gsc(f) for f in files]; gk = sorted(set(G[0]).intersection(*[set(g) for g in G[1:]]))
    TT0 = np.array([x[0] for x in TRS]); qc = 3 * np.arange(N) + 1; qd0 = float(os.environ.get("WSH_Q0", "409.25")); nc = 0; resid = []
    for a, l in gk:
        if l != 13 or a % 20 not in (2, 10): continue
        t = a * 0.5e-3 + 13 / 28e3; j = int(np.argmin(np.abs(TT0 - t)))
        if abs(TT0[j] - t) > float(os.environ.get("CSI_DT", "0.006")): SIDE.append((a, 13, "csi", -1, full, "far")); continue
        Pm = np.array([csi8_despread(*G[c][(a, l)], a % 20)[:N] for c in range(4)])            # 4 x N x 8
        if TG: Pm = np.stack([inj(Pm[:, :, p_], t, a, 12 * np.arange(N) + 3.5) for p_ in range(8)], 2)
        sh = np.array(shifts(a), float); Pm = Pm * np.exp(-1j * PH * sh[:, None, None] * (qc - float(os.environ.get("WSH_Q0", "409.25")))[None, :, None])  # window moves, about the carrier as undo_wsh
        Pm = Pm * np.exp(-1j * np.array([np.polyval(drift[c], t) for c in range(4)])[:, None, None] * (qc - qd0)[None, :, None])
        Pm = Pm * np.exp(2j * np.pi * np.array(b0)[:, None, None] * (qc - qd0)[None, :, None] / NFq)
        Rt = TRS[j][3].reshape(4, N, 3).mean(2)                                                # TRS per PRB
        Am = Pm.reshape(-1, 8); v, *_ = np.linalg.lstsq(Am, Rt.reshape(-1), rcond=None)
        S = (Am @ v).reshape(4, N); rr = np.sum(np.abs(Rt - S) ** 2) / np.sum(np.abs(Rt) ** 2); resid.append(rr)
        if rr > float(os.environ.get("CSI_MAXRES", "0.3")): SIDE.append((a, 13, "csi", -1, full, "fit")); continue
        snaps.append((t, a, 13, full, np.repeat(S, 3, 1), "csi")); nc += 1; SIDE.append((a, 13, "csi", -1, full, "ok"))
    del G
    print("CSI-RS as TRS-beam samples: %d of %d occasions kept, fit residual median %.3f" % (nc, len(resid), np.median(resid) if resid else np.nan))
def track_keys(E, dt=0.02, cmin=0.9):
    """Per-PRG precoder-epoch keys for the blind record, which mixes several UEs per symbol: each contiguous run of
    PRGs joins the track (seen within dt) it overlaps on >= 4 PRGs with the highest precoder coherence >= cmin (0.95
    below 6 shared PRGs), the continuity test of epochs(); otherwise it opens a new track. E: [(t, pm, Y)] in time."""
    TK = []; keys = []; nid = 0   # track: [template NA x NG x 12, count per PRG, last time, id]
    for t, pm, Y in E:
        Y3 = Y.reshape(NA, NG, 12); k = {}; TK = [T for T in TK if t - T[2] <= dt]
        for s0, e0 in runs(pm):
            rm = np.zeros(NG, bool); rm[s0:e0] = True; best = None; bc = 0.0
            for T in TK:
                sh = (T[1] > 0) & rm; ns = int(sh.sum())
                if ns < 4: continue
                r = np.sum(np.conj(T[0][:, sh] / T[1][sh][None, :, None]) * Y3[:, sh], axis=(0, 2))
                cc = np.abs(ESL[:, sh] @ r) / (np.abs(r).sum() + 1e-30); jj = int(np.argmax(cc))
                if cc[jj] >= (cmin if ns >= 6 else 0.95) and cc[jj] > bc: bc = cc[jj]; best = (T, SL[jj], np.angle(ESL[jj, sh] @ r))
            if best:
                T, s_, p0 = best; ph = np.exp(-1j * (s_ * q + p0)).reshape(NG, 12)
                T[0][:, rm] += Y3[:, rm] * ph[None, rm]; T[1][rm] += 1; T[2] = t
            else:
                T = [np.where(rm[None, :, None], Y3, 0), rm.astype(int), t, nid]; nid += 1; TK.append(T)
            for b in range(s0, e0): k[b] = ("trk", T[3])
        keys.append(k)
    return keys
# DM-RS
if USE:
    Dd = [load(f) for f in files]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
    cnt_ = C.Counter(k[1] for k in keys if k[1] != 0)
    tab = cnt_.most_common(1)[0][0] if cnt_ else None   # no UE decoded on all antennas: no tablet PRGs to remove
    tabp = C.defaultdict(lambda: np.zeros(NG, bool))
    for k in keys:
        if k[1] == tab: tabp[k[0]] |= np.logical_and.reduce([d[k][0] for d in Dd])
    # references for the DM-RS: the TRS, plus (CSI_REF=1) the CSI-RS-derived TRS-beam samples (closer in time)
    REFS = [(x[0], x[3]) for x in TRS]
    if os.environ.get("CSI_REF", "0") == "1": REFS += [(x[0], x[4]) for x in snaps if x[5] == "csi"]
    REFS.sort(key=lambda r: r[0])
    TT = np.array([r[0] for r in REFS]); inside = (q[1:] % 12) != 0; nd = 0
    PEND = []   # decoded-UE snapshots waiting for the epoch-averaged reference (TAB_EP=epoch)
    for k in sorted(keys):
        a, u, l = k
        if u != 0 and not (MODE == "tab" and u == tab): continue
        pm = np.logical_and.reduce([d[k][0] for d in Dd])
        if u == 0: pm = pm & ~tabp[a] if MODE == "notab" else pm & ~tabp[a]   # blind: the tablet's PRGs come from its own record
        if pm.sum() < 1: continue
        t = a * 0.5e-3 + l / 28e3; Y = inj(np.array([d[k][1] for d in Dd]).astype(complex), t, a, KK_D)
        Y = dedrift(t, undo_wsh(Y, shifts(a), -1.0)) * ramp
        j = int(np.argmin(np.abs(TT - t)))
        kind = "blind" if u == 0 else "dec"
        if abs(TT[j] - t) > 0.01: SIDE.append((a, l, kind, u, pm, "far")); continue
        R = REFS[j][1]; m = np.repeat(pm, 12); ok = m[1:] & m[:-1] & inside
        if os.environ.get("ALIGN", "0") == "1":   # off by default: the window moves are undone, so TRS and DM-RS share
            Pp = Y * np.conj(R)                       # one timing per antenna; aligning would move each snapshot's
            z = np.sum((Pp[:, 1:] * np.conj(Pp[:, :-1]))[:, ok])   # origin to its own (beam-dependent) mean delay
            if z != 0: R = R * np.exp(1j * np.angle(z) * q)[None, :]
        if kind == "dec" and os.environ.get("TAB_EP", "snap") in ("epoch", "hyb", "slot"): PEND.append((t, a, l, pm, Y, u)); continue
        if kind == "blind" and os.environ.get("BLIND_EP", "snap") in ("epoch", "hyb", "slot"): PEND.append((t, a, l, pm, Y, u)); continue
        X = np.zeros_like(Y)
        for b in np.flatnonzero(pm):
            sl = slice(12 * b, 12 * b + 12); dn = np.vdot(R[:, sl], Y[:, sl]); den = np.vdot(R[:, sl], R[:, sl]).real
            if den > 0 and dn != 0: X[:, sl] = Y[:, sl] / (dn / den)
        snaps.append((t, a, l, pm, X, "dmrs")); nd += 1; SIDE.append((a, l, kind, u, pm, "ok"))
    if PEND:   # reference numbers shared over a group of snapshots, per stream (TAB_EP: decoded UE, BLIND_EP: blind
        from so_csi import ref_fix, ref_fix_hyb   # record): epoch / hyb = detected precoder epochs, slot = per slot
        PEND.sort(key=lambda x: (x[5] == 0, x[0]))   # decoded UE first, then blind, each in time: epochs per kind
        md = lambda u: os.environ.get("BLIND_EP" if u == 0 else "TAB_EP", "snap")
        EP, _ = epochs([(t, a, pm, Y) for t, a, l, pm, Y, u in PEND])
        EP = [(u == 0, e) for e, (t, a, l, pm, Y, u) in zip(EP, PEND)]
        if os.environ.get("BLIND_TRACK", "0") == "1":   # blind record split into UE tracks: one epoch key per PRG
            ib = [k for k, x in enumerate(PEND) if x[5] == 0 and md(0) != "slot"]
            for k, kk in zip(ib, track_keys([(PEND[k][0], PEND[k][3], PEND[k][4]) for k in ib])): EP[k] = kk
        EP = [("slot", u == 0, a) if md(u) == "slot" else e for e, (t, a, l, pm, Y, u) in zip(EP, PEND)]
        if os.environ.get("HYB_MAXEP"):   # also cut epochs every HYB_MAXEP s
            me = float(os.environ["HYB_MAXEP"])
            EP = [{b: (v, int(x[0] // me)) for b, v in e.items()} if isinstance(e, dict) else (e, int(x[0] // me)) for e, x in zip(EP, PEND)]
        fit = ref_fix_hyb if "hyb" in (md(0), md(1)) else ref_fix   # hyb: d = c(t) * D(epoch, PRG)
        XR, kept = fit([(t, a, pm, Y) for t, a, l, pm, Y, u in PEND], EP, REFS, max_dt=0.01, align_delay=os.environ.get("TAB_ALIGN", "1") == "1")
        for (t, a, l, pm, Y, u), X, k in zip(PEND, XR, kept):
            if k: snaps.append((t, a, l, pm, X, "dmrs")); nd += 1; SIDE.append((a, l, "dec" if u else "blind", u, pm, "ok"))
            else: SIDE.append((a, l, "dec" if u else "blind", u, pm, "fit"))
        for nm, sel in (("tablet", 1), ("blind", 0)):   # how many snapshots share each per-PRG reference number
            ek = [((e[b] if isinstance(e, dict) else e), b) for e, x in zip(EP, PEND) if (x[5] != 0) == sel for b in np.flatnonzero(x[3])]
            if ek: print("shared reference, %s: %d PRG-snapshots, %.2f per number" % (nm, len(ek), len(ek) / len(set(ek))))
    del Dd
    print("DM-RS snapshots (%s): %d" % (MODE, nd))
# write: one record per slot and chain, every symbol present
byslot = C.defaultdict(dict)
if os.environ.get("GCAL", "0") == "1":   # per-4-PRB-group inter-antenna phase from the TRS within +-GCAL_HALF s
    from so_csi import gcal_prep, gcal_apply
    GP = gcal_prep([(x[0], x[3]) for x in TRS]); hw = float(os.environ.get("GCAL_HALF", "0.04"))
    snaps = [(t, a, l, pm, gcal_apply(t, X, GP, hw), kind) for t, a, l, pm, X, kind in snaps]
    print("inter-antenna phase calibrated per 4-PRB group on the TRS (+-%.0f ms)" % (1e3 * hw))
if os.environ.get("DUMPA"):   # diagnostic: the final snapshots (t, slot, symbol, mask, X, kind)
    import pickle; pickle.dump(snaps, open(os.environ["DUMPA"], "wb"))
for t, a, l, pm, X, kind in snaps:
    if l in byslot[a]:   # the tablet's record and the blind record of the other UEs share the symbol: disjoint PRGs, merge
        pm0, X0 = byslot[a][l]; m3 = np.repeat(pm & ~pm0, 12); X1 = X0.copy(); X1[:, m3] = X[:, m3]; byslot[a][l] = (pm0 | pm, X1)
    else: byslot[a][l] = (pm, X)
# MINPRB=<n> (or a list "0,24,48", one output <out>.m<n> each): leave out DM-RS symbols with fewer than n PRBs actually
# measured (the symbol's union of grants; TRS and CSI-RS are full band and never affected)
MINS = [int(x) for x in os.environ.get("MINPRB", "0").split(",")]
for MIN in MINS:
  OUTM = OUTP if len(MINS) == 1 else OUTP + ".m%d" % MIN
  with open(OUTM, "wb") as f:
    f.write(struct.pack("<8sIIQdIIII", b"NRSEST1", 1, 4, 3750000000, 122.88e6, 4096, 30000, NSC, 0))
    for a in sorted(byslot):
        syms = [l for l in sorted(byslot[a]) if l in (4, 8, 13) or 4 * int(byslot[a][l][0].sum()) >= MIN]
        if not syms: continue
        mask = sum(1 << l for l in syms)
        for c in range(4):
            f.write(struct.pack("<IHHBBBBHHIIQ", 0x4E525352, c, 1, 0, 1, 5, 0, mask, 0, (a // 20) % 1024, a % 20, a))
            for l in syms:
                pm, X = byslot[a][l]; valid = np.zeros(NSC, bool); val = np.zeros(NSC, complex)
                for qq in np.flatnonzero(np.repeat(pm, 4)):
                    for n in range(3):
                        for dd in (0, 2): kk = 12 * qq + 4 * n + dd; valid[kk] = True; val[kk] = X[c, 3 * qq + n]
                f.write(struct.pack("<I", int(valid.sum()))); f.write(np.packbits(valid, bitorder="little").tobytes()); f.write(val[valid].astype(np.complex64).tobytes())
if TG:
    with open(OUTP + ".truth.csv", "w") as ft:
        ft.write("t_s,frame,slot,slot_abs,target,active,range_m,speed_ms,level_db,angle_deg\n")
        for (a_, i_), (r_, v_) in sorted(TRUTH.items()): ft.write("0,%d,%d,%d,%d,1,%.3f,%.4f,%.1f,%.1f\n" % ((a_ // 20) % 1024, a_ % 20, a_, i_, r_, v_, TG[i_][2], TG[i_][3]))
print("wrote %s: %d slots" % (OUTP, len(byslot)))
if os.environ.get("SIDE"):
    import json
    with open(os.environ["SIDE"], "w") as fs_:
        for a_, l_, k_, u_, pm_, r_ in SIDE: fs_.write(json.dumps([int(a_), int(l_), k_, int(u_), np.flatnonzero(pm_).tolist(), r_]) + "\n")
