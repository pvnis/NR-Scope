"""Write GoldenSniffer DM-RS dumps as an NR-Scope sensing.record_estimates file (nrscope_sensing_record.h), so that
NR-Scope's own sensing pipeline (sensing_offline: alignment, clutter removal, Doppler, detector, AoA, map dump, target
injection) runs on them. One stream: the reference UE (the UE with the most DM-RS snapshots, the tablet on 5 Nov).
Variants (env VAR):
  raw     the DM-RS as dumped (port 1000, 3 CDM-pair values per PRB), precoder jumps included
  proc    our processing: subframe phase step removed, precoder epochs stitched across PRGs (4 antennas), epochs glued
          with the nearest-snapshot anchor
NR-Scope needs one contiguous block per snapshot: the longest contiguous PRG run is kept (env RUN=all keeps all).
Each CDM-pair value is written on both REs of its pair (subcarriers 12q + 4n and 12q + 4n + 2), only on the UE's PRGs.
  csi     CSI-RS-based correction (so_csi.py): per (epoch, PRG) one number onto the burst's reference beam
          through the nearest 8-port CSI-RS (slots 2/10), subframe step removed first; no stitching, no anchor
          env CSI8=1: the full 8-weight conversion instead (so_csi.csi_fix8)
env GLUE=prg (with VAR=proc): glue epochs PRG by PRG (so_csi.glue_per_prg) instead of one number
env LOS0 (1): put range zero on the direct path (burst median of the strongest delay bin) before writing
env BLOCK=<s> (with VAR=proc): one stitching template per time block of <s> seconds (so_csi.block_epochs), not per epoch
env WSH=1: undo NR-Scope's window moves exactly from the GSW1 records (per antenna), then a straight-line clock drift
env DEJUMP=1: remove NR-Scope's timing-window jumps first (so_csi.dejump / dejump_csi), for every variant
  trs     per (epoch, PRG) one number onto the nearest TRS (so_csi.ref_fix): no stitching chain, no anchor;
          env REF_ALIGN (1: delay-align the TRS to each snapshot first), REF_DT (0.01 s: farthest TRS)
usage: gs_to_nrsr.py <burst dir> <out estimates.bin>   env VAR (raw | proc | csi | trs), DIR (replay_blind_ch), CH (0,1,2,3)"""
import sys, os, struct, numpy as np, collections as C
BD, OUTF = sys.argv[1].rstrip("/") + "/", sys.argv[2]; saved = list(sys.argv)
os.environ.setdefault("CH", "0,1,2,3")
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
VAR = os.environ.get("VAR", "proc"); DIR = os.environ.get("DIR", "replay_blind_ch")
NSC = 3276

Dd = [load(BD + DIR + "%d/ue_dmrs.bin" % c) for c in CH]
keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
S = C.defaultdict(list)
for k in sorted(keys):
    if k[1] == 0: continue
    pm = np.logical_and.reduce([Dd[c][k][0] for c in range(NA)]); Y = np.array([Dd[c][k][1] for c in range(NA)]).astype(complex)
    if pm.sum() and np.all(np.isfinite(Y[:, np.repeat(pm, 12)])): S[k[1]].append((k[0] * 0.5e-3 + k[2] / 28e3, k[0], pm, Y, k[2]))
del Dd
U = max(S, key=lambda u: len(S[u])); A = sorted(S[U], key=lambda x: x[0])
print("reference UE %#x, %d snapshots" % (U, len(A)))
if os.environ.get("WSH", "0") == "1":   # undo NR-Scope's logged window moves exactly (GSW1), per antenna, then the clock drift
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from so_csi import load_wsh, undo_wsh, slope_of, remove_drift
    WS = [load_wsh(BD + DIR + "%d/ue_dmrs.bin" % c) for c in CH]
    def shifts(a): return [W.get(a, W[min(W, key=lambda b: abs(b - a))]) if W else 0 for W in WS]
    jit = {}
    for sg in (1.0, -1.0):   # the sign that makes consecutive snapshots' slopes steadiest
        sl = np.unwrap([slope_of(undo_wsh(Y, shifts(a), sg), pm) for t, a, pm, Y, l in A[::3]]); jit[sg] = float(np.median(np.abs(np.diff(sl))))
    SG = min(jit, key=jit.get)
    A = [(t, a, pm, undo_wsh(Y, shifts(a), SG), l) for t, a, pm, Y, l in A]
    A, pdr, rs = remove_drift(A)
    print("window moves undone (sign %+d; slope jitter %.4f vs %.4f rad), drift %.4f rad/s on the grid, residual %.4f rad" % (SG, jit[SG], jit[-SG], pdr[0], rs))
if os.environ.get("DEJUMP", "0") == "1":   # remove NR-Scope's timing-window jumps first (so_csi.dejump)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from so_csi import dejump
    A = [(t, a, pm, dejump(Y, pm)[0], l) for t, a, pm, Y, l in A]
if VAR == "proc":
    st = []
    for (t1, a1, p1, Y1, _), (t2, a2, p2, Y2, _) in zip(A[:-1], A[1:]):
        m = np.repeat(p1 & p2, 12)
        if m.sum() >= 120 and a2 - a1 <= 2 and a2 // 2 - a1 // 2 == 1:
            v = np.vdot(Y1[:, m], Y2[:, m]); st.append(v / abs(v))
    phs = float(np.angle(np.mean(st))) if st else 0.0
    A = [(t, a, pm, Y * np.exp(-1j * phs * (a // 2)), l) for t, a, pm, Y, l in A]
    if os.environ.get("BLOCK", ""):   # one stitching template per time block instead of per detected epoch
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from so_csi import block_epochs
        EP, TPL = block_epochs([(t, a, pm, Y) for t, a, pm, Y, l in A], float(os.environ["BLOCK"]), SL, q, NA, NG)
    else:
        EP, TPL = epochs([(t, a, pm, Y) for t, a, pm, Y, l in A])
    EG = {k: joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0) for k, (S2, C2) in TPL.items()}
    X0 = []
    for i, (t, a, pm, Y, l) in enumerate(A):
        Z = Y / np.repeat(EG[EP[i]], 12)[None, :]; Z[:, ~np.repeat(pm, 12)] = 0; X0.append(Z)
    XA = [X0[0]]; g = 1 + 0j
    for i in range(1, len(A)):
        if EP[i] != EP[i - 1]:
            f, n = frame(XA[i - 1], A[i - 1][2], X0[i], A[i][2]); g = f if (np.isfinite(f) and n >= 4 and f != 0) else 1 + 0j
        XA.append(X0[i] / g)
    if os.environ.get("GLUE", "") == "prg":
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from so_csi import glue_per_prg
        XA = glue_per_prg([(t, a, pm, Y) for t, a, pm, Y, l in A], EP, X0)
    A = [(t, a, pm, X, l) for (t, a, pm, Y, l), X in zip(A, XA)]
    print("subframe step %.1f deg, %d epochs" % (np.degrees(phs), len(TPL)))
# NR-Scope's sensing pipeline needs a uniform pilot lattice (one contiguous block of subcarriers per snapshot):
# keep each snapshot's longest contiguous PRG run (env RUN=all keeps everything)
def longest(pm):
    rr = runs(pm)
    if not rr: return pm
    s_, e_ = max(rr, key=lambda r: r[1] - r[0]); m = np.zeros_like(pm); m[s_:e_] = True; return m
KEEP_ALL = os.environ.get("RUN", "longest") == "all"
if VAR == "csi":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gs_trs_lib import load_gsc, csi8_despread
    from so_csi import csi_fix
    st = []
    for (t1, a1, p1, Y1, _), (t2, a2, p2, Y2, _) in zip(A[:-1], A[1:]):
        m = np.repeat(p1 & p2, 12)
        if m.sum() >= 120 and a2 - a1 <= 2 and a2 // 2 - a1 // 2 == 1:
            v = np.vdot(Y1[:, m], Y2[:, m]); st.append(v / abs(v))
    phs = float(np.angle(np.mean(st))) if st else 0.0
    A = [(t, a, pm, Y * np.exp(-1j * phs * (a // 2)), l) for t, a, pm, Y, l in A]
    Gs = [load_gsc(BD + DIR + "%d/ue_dmrs.bin" % c) for c in CH]; gk = sorted(set(Gs[0]).intersection(*[set(g) for g in Gs[1:]]))
    CS = [(a * 0.5e-3 + 13 / 28e3, np.array([csi8_despread(*Gs[c][(a, l)], a % 20)[:N] for c in range(NA)]) * np.exp(-1j * phs * (a // 2)))
          for a, l in gk if l == 13 and a % 20 in (2, 10)]
    if os.environ.get("DEJUMP", "0") == "1":
        from so_csi import dejump_csi
        CS = [(t, dejump_csi(P)[0]) for t, P in CS]
    del Gs
    EP, _ = epochs([(t, a, pm, Y) for t, a, pm, Y, l in A])
    if os.environ.get("CSI8", "0") == "1":
        from so_csi import csi_fix8
        XC, kept, stc = csi_fix8([(t, a, pm, Y) for t, a, pm, Y, l in A], EP, CS, NG, NA)
    else:
        XC, kept, stc = csi_fix([(t, a, pm, Y) for t, a, pm, Y, l in A], EP, CS, NG, NA)
    A = [(t, a, pm, X, l) for (t, a, pm, Y, l), X, k in zip(A, XC, kept) if k]
    print("CSI-RS: %d snapshots, residual in the 8-port span %.3f, reference beam %.0f%% of the weight energy, kept %.0f%%" %
          (len(CS), stc["span_residual"], 100 * stc["wref_share"], 100 * stc["kept"]))
if VAR == "trs":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gs_trs_lib import load_gsc, trs_despread
    from so_csi import ref_fix, dejump
    st = []
    for (t1, a1, p1, Y1, _), (t2, a2, p2, Y2, _) in zip(A[:-1], A[1:]):
        m = np.repeat(p1 & p2, 12)
        if m.sum() >= 120 and a2 - a1 <= 2 and a2 // 2 - a1 // 2 == 1:
            v = np.vdot(Y1[:, m], Y2[:, m]); st.append(v / abs(v))
    phs = float(np.angle(np.mean(st))) if st else 0.0
    A = [(t, a, pm, Y * np.exp(-1j * phs * (a // 2)), l) for t, a, pm, Y, l in A]
    Gs = [load_gsc(BD + DIR + "%d/ue_dmrs.bin" % c) for c in CH]; gk = sorted(set(Gs[0]).intersection(*[set(g) for g in Gs[1:]]))
    REF = []
    for a, l in gk:
        if l not in (4, 8): continue
        R = np.array([trs_despread(*Gs[c][(a, l)], a % 20, l, 2)[:N].ravel() for c in range(NA)]) * np.exp(-1j * phs * (a // 2))
        if os.environ.get("DEJUMP", "0") == "1": R = dejump(R, np.ones(NG, bool))[0]
        if os.environ.get("WSH", "0") == "1": R = undo_wsh(R, shifts(a), SG) * np.exp(-1j * np.polyval(pdr, a * 0.5e-3 + l / 28e3) * np.arange(R.shape[1]))[None, :]
        REF.append((a * 0.5e-3 + l / 28e3, R))
    del Gs
    EP, _ = epochs([(t, a, pm, Y) for t, a, pm, Y, l in A])
    XR, kept = ref_fix([(t, a, pm, Y) for t, a, pm, Y, l in A], EP, REF, max_dt=float(os.environ.get("REF_DT", "0.01")),
                       align_delay=os.environ.get("REF_ALIGN", "1") == "1")
    print("TRS reference: %d of %d snapshots have a TRS within %s s" % (sum(kept), len(kept), os.environ.get("REF_DT", "0.01")))
    A = [(t, a, pm, X, l) for (t, a, pm, Y, l), X, k in zip(A, XR, kept) if k]
    print("TRS reference: %d TRS symbols, kept %d snapshots" % (len(REF), len(A)))
if os.environ.get("LOS0", "1") == "1" and A:   # range zero on the direct path, not on the mean delay (so_csi.los_to_zero)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from so_csi import los_to_zero
    A, b0 = los_to_zero(A)
    print("direct path moved to range zero per antenna: it was at %s m" % " ".join("%+.1f" % (b / 4096 / 120e3 * 3e8) for b in b0))
byslot = C.defaultdict(dict)
for t, a, pm, X, l in A: byslot[a][l] = (pm if KEEP_ALL else longest(pm), X)
# both DM-RS symbols of a slot share one record: use the PRGs present on both
for a in byslot:
    if len(byslot[a]) > 1:
        cm = np.logical_and.reduce([v[0] for v in byslot[a].values()])
        byslot[a] = {l: (cm, X) for l, (pm, X) in byslot[a].items()} if cm.any() else {min(byslot[a]): byslot[a][min(byslot[a])]}
with open(OUTF, "wb") as f:
    f.write(struct.pack("<8sIIQdIIII", b"NRSEST1", 1, NA, 3750000000, 122.88e6, 4096, 30000, NSC, 0))
    for a in sorted(byslot):
        syms = sorted(byslot[a]); mask = sum(1 << l for l in syms)
        for c in range(NA):
            f.write(struct.pack("<IHHBBBBHHIIQ", 0x4E525352, c, 1, 0, 1, 5, 0, mask, 0, (a // 20) % 1024, a % 20, a))
            for l in syms:
                pm, X = byslot[a][l]
                valid = np.zeros(NSC, bool); val = np.zeros(NSC, complex)
                for q in np.flatnonzero(np.repeat(pm, 4)):          # PRBs of the UE's PRGs
                    for n in range(3):
                        for d in (0, 2):
                            k = 12 * q + 4 * n + d; valid[k] = True; val[k] = X[c, 3 * q + n]
                f.write(struct.pack("<I", int(valid.sum())))
                f.write(np.packbits(valid, bitorder="little").tobytes())
                f.write(val[valid].astype(np.complex64).tobytes())
print("wrote %s: %d slots x %d chains, %s" % (OUTF, len(byslot), NA, VAR))
