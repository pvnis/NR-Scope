"""Synthetic DM-RS with simulated precoders on the REAL occupancy of the reference UE (tablet) of one burst, written as
record_estimates files for NR-Scope's sensing_offline (targets are injected there).
Scene: static only (direct path + 20 scatterers, -45..-25 dB, 10..120 m excess, calibrated like sim_sense), 8-port gNB
ULA, 4 receive antennas, noise at SNR (dB per pair value, direct path). The tablet's slots, DM-RS symbols and PRGs are
taken from the replay; gaps between its PRGs are kept (NR-Scope sensing.allow_gaps).
Precoder per snapshot: w_b = b1 + g_b b2 (MIX), times a +-90 deg phase per PRG, times one complex number per epoch;
a new epoch every ~19 ms (exponential), as measured on the cell.
Variants (one file each):
  oracle  the same snapshots through one fixed beam b1: no precoder change at all
  raw     the precoded snapshots as a receiver would see them
  csi     raw after the CSI-RS-based correction (so_csi.py), CSI-RS simulated (8 unprecoded ports, slots 2/10)
  proc2   as proc, but epochs glued PRG by PRG (so_csi.glue_per_prg)
  proc    raw after our processing (subframe step not simulated; epochs stitched + glued, as gs_to_nrsr.py VAR=proc)
  csi8    raw after the full 8-weight CSI-RS conversion to one virtual beam (so_csi.csi_fix8)
Timing-window jumps: env JUMPS=1 (steps of 1..16 samples, ~every 40 ms, on every signal); env DEJUMP=1 removes them
(so_csi.dejump: per-snapshot slope inside PRGs, on oracle/raw and the CSI-RS) before every processing.
Precoder direction diversity: env K extra beams with random per-PRG weights of scale MIX (K=0: two-beam model).
usage: so_simprec.py <burst dir> <out prefix>   env MIX (0.2), K (0), SNR (25), SEED"""
import sys, os, struct, numpy as np, collections as C
BD, OUTP = sys.argv[1].rstrip("/") + "/", sys.argv[2]; saved = list(sys.argv)
os.environ["CH"] = "0,1,2,3"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
MIX = float(os.environ.get("MIX", "0.2")); SNR = float(os.environ.get("SNR", "25")); rng = np.random.default_rng(int(os.environ.get("SEED", "1")))
NT = 8; NSC = 3276
FK = (12 * np.repeat(np.arange(N), 3) + np.tile([1, 5, 9], N)) * 30e3
def atx(s): return np.exp(1j * np.pi * np.arange(NT) * s) / np.sqrt(NT)
def arx(s): return np.exp(1j * np.pi * np.arange(NA) * s)
# scene
P = [(1.0, 0.0, 0.3, 0.2)] + [(10 ** (rng.uniform(-45, -25) / 20) * np.exp(1j * rng.uniform(0, 2 * np.pi)), rng.uniform(10, 120) / 3e8,
                              rng.uniform(-0.9, 0.9), rng.uniform(-0.9, 0.9)) for _ in range(20)]
E = np.array([np.exp(-2j * np.pi * FK * p[1]) for p in P]); At = np.array([atx(p[2]) for p in P]); Ar = np.array([arx(p[3]) for p in P]).T
amp = np.array([p[0] for p in P]); sig = 10 ** (-SNR / 20) / np.sqrt(2)
# NR-Scope timing-window jumps (env JUMPS=1): the FFT window moves by whole samples (8.14 ns at 122.88 Msps), measured on
# the cell from the TRS: about every other 20 ms TRS pair shows a step of 1..16 samples. Applied to every signal.
JUMPS = os.environ.get("JUMPS", "0") == "1"; TS = 1 / 122.88e6; _win = {"off": 0, "next": 0}
def win_delay(a):
    if JUMPS and a >= _win["next"]:
        if _win["next"]: _win["off"] = int(np.clip(_win["off"] + rng.choice([-1, 1]) * rng.integers(1, 17), -24, 24))
        _win["next"] = a + int(rng.exponential(80)) + 1
    return _win["off"] * TS
def obs(Wk, a=0, f=None):   # Wk: NT x 816
    y = Ar @ (E * (np.conj(At) @ Wk) * amp[:, None])
    if JUMPS: y = y * np.exp(-2j * np.pi * (FK if f is None else f) * win_delay(a))[None, :]
    return y + sig * (rng.standard_normal(y.shape) + 1j * rng.standard_normal(y.shape))
# the tablet's occupancy
D0 = load(BD + "replay_ch0/ue_dmrs.bin"); cnt = C.Counter(k[1] for k in D0 if k[1] != 0); U = cnt.most_common(1)[0][0]
occ = sorted((k[0], k[2], D0[k][0]) for k in D0 if k[1] == U); del D0
print("UE %#x: %d snapshots" % (U, len(occ)))
b1 = atx(0.25); b2 = atx(0.5); nxt = -1; snaps = {"oracle": [], "raw": []}
KB = int(os.environ.get("K", "0")); BX = [atx(s_) for s_ in rng.uniform(-0.6, 0.9, KB)]   # extra beams (direction diversity)
shares = []
for a, l, pm in occ:
    if a >= nxt:
        nxt = a + int(rng.exponential(38)) + 1
        g = (rng.standard_normal(NG) + 1j * rng.standard_normal(NG)) * MIX; ph = np.exp(1j * rng.choice([np.pi / 2, -np.pi / 2], NG))
        sc = np.exp(1j * rng.uniform(0, 2 * np.pi)) * rng.uniform(0.5, 2)
        GX = (rng.standard_normal((KB, NG)) + 1j * rng.standard_normal((KB, NG))) * MIX
    Wb = (b1[:, None] + b2[:, None] * g[None, :] + sum(BX[i][:, None] * GX[i][None, :] for i in range(KB))) * ph[None, :]
    Wb = Wb / np.linalg.norm(Wb, axis=0, keepdims=True) * sc
    if len(shares) < 3000 and pm.sum() >= 4:
        Wn = Wb[:, pm] / np.linalg.norm(Wb[:, pm], axis=0, keepdims=True); ev_, V_ = np.linalg.eigh(Wn @ Wn.conj().T)
        shares += list(np.abs(V_[:, -1].conj() @ Wn) ** 2)
    m = np.repeat(pm, 12); t = a * 0.5e-3 + l / 28e3
    for nm, W in (("raw", np.repeat(Wb, 12, 1)), ("oracle", np.repeat(b1[:, None], L, 1))):
        Y = obs(W, a); Y[:, ~m] = 0; snaps[nm].append((t, a, pm, Y, l))
# timing-window jumps removed first (env DEJUMP=1, so_csi.dejump), on every variant derived from the measurement
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from so_csi import dejump, dejump_csi
DJ = os.environ.get("DEJUMP", "0") == "1"
if DJ:
    for nm in ("oracle", "raw"):
        snaps[nm] = [(t, a, pm, dejump(Y, pm)[0], l) for t, a, pm, Y, l in snaps[nm]]
# our processing on raw (as gs_to_nrsr.py VAR=proc, without the subframe step, which is not simulated)
A = snaps["raw"]; EP, TPL = epochs([(t, a, pm, Y) for t, a, pm, Y, l in A])
EG = {k: joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0) for k, (S2, C2) in TPL.items()}
X0 = []
for i, (t, a, pm, Y, l) in enumerate(A):
    Z = Y / np.repeat(EG[EP[i]], 12)[None, :]; Z[:, ~np.repeat(pm, 12)] = 0; X0.append(Z)
XA = [X0[0]]; gg = 1 + 0j
for i in range(1, len(A)):
    if EP[i] != EP[i - 1]:
        f, n = frame(XA[i - 1], A[i - 1][2], X0[i], A[i][2]); gg = f if (np.isfinite(f) and n >= 4 and f != 0) else 1 + 0j
    XA.append(X0[i] / gg)
snaps["proc"] = [(t, a, pm, X, l) for (t, a, pm, Y, l), X in zip(A, XA)]
def _spread(S, EPs, across):
    out = []
    for i in range(1, len(S)):
        if (EPs[i] != EPs[i - 1]) != across or S[i][1] - S[i - 1][1] > 4: continue
        pm = S[i][2] & S[i - 1][2]
        if pm.sum() < 8: continue
        r = np.array([np.vdot(S[i - 1][3][:, 12 * b:12 * b + 12], S[i][3][:, 12 * b:12 * b + 12]) for b in np.flatnonzero(pm)]); c = np.mean(r / np.abs(r))
        out.append(np.degrees(np.std(np.angle(r / c))))
    return np.median(out) if out else np.nan
from so_csi import glue_per_prg
snaps["proc2"] = [(t, a, pm, X, l) for (t, a, pm, Y, l), X in zip(A, glue_per_prg([(t, a, pm, Y) for t, a, pm, Y, l in A], EP, X0))]
# the TRS check of the real data (trsref): every snapshot against a TRS-like beam (port 0) observed 0.5 ms later
if os.environ.get("TRSCHECK", "0") == "1":
    w0 = np.repeat(np.eye(NT)[:, :1], L, 1); full = np.ones(NG, bool)
    def cmp_(R, X, pm):
        b = np.flatnonzero(pm); r = np.array([np.vdot(R[:, 12 * k:12 * k + 12], X[:, 12 * k:12 * k + 12]) for k in b]); u = r / np.abs(r); c = np.mean(u)
        return abs(c), np.degrees(np.std(np.angle(u / c)))
    for nm in ("oracle", "raw", "proc"):
        v = np.array([cmp_(dejump(obs(w0, a), full)[0] if DJ else obs(w0, a), X, pm) for t, a, pm, X, l in snaps[nm][::10] if pm.sum() >= 16])
        print("TRS check, %-6s vs TRS-like beam: coherence %.2f, per-PRG phase spread %.0f deg" % (nm, np.median(v[:, 0]), np.median(v[:, 1])))
print("per-PRG gluing: spread across an epoch change %.0f deg" % _spread(snaps["proc2"], EP, True))
print("per-PRG phase spread after stitching+glue: across an epoch change %.0f deg, inside an epoch %.0f deg (real cell 50 / 12)" % (_spread(snaps["proc"], EP, True), _spread(snaps["proc"], EP, False)))
print("%d epochs; share of a PRG's precoder along its snapshot's main beam: median %.2f (real cell 0.38)" % (len(TPL), np.median(shares)))
# CSI-RS: 8 unprecoded ports, one value per PRB, slots 2 and 10 of every frame
FCS = (12 * np.arange(N) + 3.5) * 30e3; ECS = np.array([np.exp(-2j * np.pi * FCS * p[1]) for p in P])
a_lo, a_hi = occ[0][0], occ[-1][0]; CS = []
for a in range(a_lo, a_hi + 1):
    if a % 20 in (2, 10):
        Pm = np.zeros((NA, N, NT), complex)
        for pp in range(NT):
            y = Ar @ (ECS * (np.conj(At) @ np.repeat(np.eye(NT)[:, pp:pp + 1], N, 1)) * amp[:, None])
            if JUMPS: y = y * np.exp(-2j * np.pi * FCS * win_delay(a))[None, :]
            Pm[:, :, pp] = y + sig * (rng.standard_normal(y.shape) + 1j * rng.standard_normal(y.shape))
        CS.append((a * 0.5e-3 + 13 / 28e3, dejump_csi(Pm)[0] if DJ else Pm))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from so_csi import csi_fix
XC, kept, stc = csi_fix([(t, a, pm, Y) for t, a, pm, Y, l in A], EP, CS, NG, NA)
snaps["csi"] = [(t, a, pm, X, l) for (t, a, pm, Y, l), X, k in zip(A, XC, kept) if k]
from so_csi import csi_fix8
XC8, kept8, st8 = csi_fix8([(t, a, pm, Y) for t, a, pm, Y, l in A], EP, CS, NG, NA)
snaps["csi8"] = [(t, a, pm, X, l) for (t, a, pm, Y, l), X, k in zip(A, XC8, kept8) if k]
print("CSI-RS 8 weights: residual %.3f, reference beam %.0f%%" % (st8["span_residual"], 100 * st8["wref_share"]))
print("CSI-RS: residual %.3f, reference beam %.0f%%, kept %.0f%%" % (stc["span_residual"], 100 * stc["wref_share"], 100 * stc["kept"]))
def write(path, S):
    byslot = C.defaultdict(dict)
    for t, a, pm, X, l in S: byslot[a][l] = (pm, X)
    for a in byslot:
        if len(byslot[a]) > 1:
            cm = np.logical_and.reduce([v[0] for v in byslot[a].values()])
            byslot[a] = {l: (cm, X) for l, (pm, X) in byslot[a].items()} if cm.any() else {min(byslot[a]): byslot[a][min(byslot[a])]}
    with open(path, "wb") as f:
        f.write(struct.pack("<8sIIQdIIII", b"NRSEST1", 1, NA, 3750000000, 122.88e6, 4096, 30000, NSC, 0))
        for a in sorted(byslot):
            syms = sorted(byslot[a]); mask = sum(1 << l for l in syms)
            for c in range(NA):
                f.write(struct.pack("<IHHBBBBHHIIQ", 0x4E525352, c, 1, 0, 1, 5, 0, mask, 0, (a // 20) % 1024, a % 20, a))
                for l in syms:
                    pm, X = byslot[a][l]; valid = np.zeros(NSC, bool); val = np.zeros(NSC, complex)
                    for q in np.flatnonzero(np.repeat(pm, 4)):
                        for n in range(3):
                            for d in (0, 2): k = 12 * q + 4 * n + d; valid[k] = True; val[k] = X[c, 3 * q + n]
                    f.write(struct.pack("<I", int(valid.sum()))); f.write(np.packbits(valid, bitorder="little").tobytes()); f.write(val[valid].astype(np.complex64).tobytes())
for nm, S in snaps.items():
    write("%s_%s.bin" % (OUTP, nm), S); print("wrote %s_%s.bin" % (OUTP, nm))
