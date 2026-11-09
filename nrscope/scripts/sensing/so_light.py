"""Light network: one sensing stream from the TRS and the DM-RS of every UE except (optionally) the tablet, written as a
record_estimates file for NR-Scope's sensing_offline.
- TRS symbols 4/8 (slots 5/6, every 20 ms, full band, fixed beam): each TRS value (REs 2, 6, 10 of a PRB) is written on
  the DM-RS pair grid (REs 4n and 4n+2), so TRS and DM-RS share one pilot lattice and NR-Scope makes ONE map of both.
- DM-RS: the blind records (rnti 0, all 4 antennas, every PCI-scrambled DM-RS) minus the tablet's PRGs (MODE=notab), or
  plus the tablet (MODE=tab); each snapshot's PRGs put on the nearest TRS's reference (one number per snapshot and
  PRG, TRS delay-aligned first), or left out (USE_DMRS=0: TRS only).
- Every antenna: NR-Scope's logged window moves undone (GSW1), straight-line drift removed, range zero on that
  antenna's direct path (measured on the TRS, applied to everything of that antenna).
usage: so_light.py <burst dir> <out.bin>   env MODE (notab | tab), USE_DMRS (1), DIR (replay_blind_ch)"""
import sys, os, struct, numpy as np, collections as C
BD, OUTP = sys.argv[1].rstrip("/") + "/", sys.argv[2]; saved = list(sys.argv); os.environ["CH"] = "0,1,2,3"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gs_trs_lib import load_gsc, trs_despread
from so_csi import load_wsh, undo_wsh, slope_of
MODE = os.environ.get("MODE", "notab"); USE = os.environ.get("USE_DMRS", "1") == "1"; DIR = os.environ.get("DIR", "replay_blind_ch")
NSC = 3276; NFq = 4096; files = [BD + DIR + "%d/ue_dmrs.bin" % c for c in range(4)]
WS = [load_wsh(f) for f in files]
def shifts(a): return [W.get(a, W[min(W, key=lambda b: abs(b - a))]) if W else 0 for W in WS]
# TRS
G = [load_gsc(f) for f in files]; gk = sorted(set(G[0]).intersection(*[set(g) for g in G[1:]]))
TRS = [(a * 0.5e-3 + l / 28e3, a, l, undo_wsh(np.array([trs_despread(*G[c][(a, l)], a % 20, l, 2)[:N].ravel() for c in range(4)]), shifts(a), -1.0))
       for a, l in gk if l in (4, 8)]
del G
# drift: straight line through the TRS slope (full band, fixed beam: the cleanest), per antenna
t_ = np.array([x[0] for x in TRS]); q = np.arange(L); full = np.ones(NG, bool)
drift = []
for c in range(4):
    s = np.unwrap([slope_of(x[3][c:c + 1], full) for x in TRS]); p = np.polyfit(t_, s, 1); r = s - np.polyval(p, t_)
    k = np.abs(r) < 3 * np.std(r) + 1e-9; drift.append(np.polyfit(t_[k], s[k], 1))
def dedrift(t, Y): return Y * np.exp(-1j * np.array([np.polyval(drift[c], t) for c in range(4)])[:, None] * q[None, :])
TRS = [(t, a, l, dedrift(t, R)) for t, a, l, R in TRS]
# range zero per antenna on the TRS direct path
P = 0
for t, a, l, R in TRS: P = P + np.abs(np.fft.ifft(R * np.hanning(L)[None, :], NFq, axis=1)) ** 2
b0 = []
for Pa in P:
    k = int(np.argmax(Pa)); y0, y1, y2 = Pa[k - 1], Pa[k], Pa[(k + 1) % NFq]; d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
    v = k + d; b0.append(v - NFq if v > NFq / 2 else v)
ramp = np.exp(2j * np.pi * np.array(b0)[:, None] * q[None, :] / NFq)
TRS = [(t, a, l, R * ramp) for t, a, l, R in TRS]
print("TRS: %d symbols; direct path per antenna was at %s m" % (len(TRS), " ".join("%+.1f" % (b / NFq / 120e3 * 3e8) for b in b0)))
snaps = [(t, a, l, full, R, "trs") for t, a, l, R in TRS]
# DM-RS
if USE:
    Dd = [load(f) for f in files]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
    tab = C.Counter(k[1] for k in keys if k[1] != 0).most_common(1)[0][0]
    tabp = C.defaultdict(lambda: np.zeros(NG, bool))
    for k in keys:
        if k[1] == tab: tabp[k[0]] |= np.logical_and.reduce([d[k][0] for d in Dd])
    TT = np.array([x[0] for x in TRS]); inside = (q[1:] % 12) != 0; nd = 0
    for k in sorted(keys):
        a, u, l = k
        if u != 0 and not (MODE == "tab" and u == tab): continue
        pm = np.logical_and.reduce([d[k][0] for d in Dd])
        if u == 0: pm = pm & ~tabp[a] if MODE == "notab" else pm & ~tabp[a]   # blind: the tablet's PRGs come from its own record
        if pm.sum() < 1: continue
        t = a * 0.5e-3 + l / 28e3; Y = np.array([d[k][1] for d in Dd]).astype(complex)
        Y = dedrift(t, undo_wsh(Y, shifts(a), -1.0)) * ramp
        j = int(np.argmin(np.abs(TT - t)))
        if abs(TT[j] - t) > 0.01: continue
        R = TRS[j][3]; m = np.repeat(pm, 12); ok = m[1:] & m[:-1] & inside; Pp = Y * np.conj(R)
        z = np.sum((Pp[:, 1:] * np.conj(Pp[:, :-1]))[:, ok])
        if z != 0: R = R * np.exp(1j * np.angle(z) * q)[None, :]
        X = np.zeros_like(Y)
        for b in np.flatnonzero(pm):
            sl = slice(12 * b, 12 * b + 12); dn = np.vdot(R[:, sl], Y[:, sl]); den = np.vdot(R[:, sl], R[:, sl]).real
            if den > 0 and dn != 0: X[:, sl] = Y[:, sl] / (dn / den)
        snaps.append((t, a, l, pm, X, "dmrs")); nd += 1
    del Dd
    print("DM-RS snapshots (%s): %d" % (MODE, nd))
# write: one record per slot and chain, every symbol present
byslot = C.defaultdict(dict)
for t, a, l, pm, X, kind in snaps: byslot[a][l] = (pm, X)
with open(OUTP, "wb") as f:
    f.write(struct.pack("<8sIIQdIIII", b"NRSEST1", 1, 4, 3750000000, 122.88e6, 4096, 30000, NSC, 0))
    for a in sorted(byslot):
        syms = sorted(byslot[a]); mask = sum(1 << l for l in syms)
        for c in range(4):
            f.write(struct.pack("<IHHBBBBHHIIQ", 0x4E525352, c, 1, 0, 1, 5, 0, mask, 0, (a // 20) % 1024, a % 20, a))
            for l in syms:
                pm, X = byslot[a][l]; valid = np.zeros(NSC, bool); val = np.zeros(NSC, complex)
                for qq in np.flatnonzero(np.repeat(pm, 4)):
                    for n in range(3):
                        for dd in (0, 2): kk = 12 * qq + 4 * n + dd; valid[kk] = True; val[kk] = X[c, 3 * qq + n]
                f.write(struct.pack("<I", int(valid.sum()))); f.write(np.packbits(valid, bitorder="little").tobytes()); f.write(val[valid].astype(np.complex64).tobytes())
print("wrote %s: %d slots" % (OUTP, len(byslot)))
