"""What the sensing pipeline would actually use, slot by slot and PRG by PRG, on a 4-antenna burst recording.
Same processing as gs_all.py, run per 150 ms window (one range-Doppler map):
  reference UE  = UE with the most PRG-snapshots in the window (scenario "tab": the tablet; scenario "notab": tablet removed)
  reference     : epoch templates + joint PRG stitching (4 antennas); PRGs cut off from the template's main run by a hole
                  have an unknown phase ("split"); PRGs missing inside a snapshot's span but present in its epoch
                  template are filled from other slots ("epoch fill")
  other UEs     : self-stitched per snapshot; anchor = reference snapshot within BRIDGE s sharing >= 4 PRGs;
                  merged if anchor fit residual < GATE, rejected otherwise; no anchor but touching the reference in the
                  same slot -> edge matching (fallback); else not joinable
  on air        : PCI DM-RS seen on antenna 0 (survey, symbol 2, strict) that no placed UE covers
  periodic      : TRS (symbols 4+8) and CSI-RS on symbol 13 detected by the survey
  blind (env DIR=replay_blind_ch, item 6e): the survey's PCI-scrambled DM-RS of every slot (rnti 0) minus the PRGs of
                  placed UEs in that slot, handled as one more UE ("blind"): merged into the reference like any other UE;
                  reference itself only when no placed UE is in the window
  TRS (6f): the TRS symbols themselves ("TRS, full band", category 8), and DM-RS that the pipeline could not join
                  coherently (split by a hole, rejected, not joinable) but that fits the nearest TRS within BRIDGE with one
                  number and a relative delay slope under GATE ("joined via TRS anchor", category 9)
usage: gs_usemap.py <rec dir> <out dir> [burst ...]   env CH (default 0,3), BRIDGE (0.01), GATE (0.3), DIR (replay_ch)"""
import sys, os, glob, json, numpy as np, collections as C
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from gs_trs_lib import load_gsc, trs_despread
sys.argv_saved = list(sys.argv)
REC, OUT = sys.argv[1], sys.argv[2]; ONLY = [int(x) for x in sys.argv[3:]]
os.environ.setdefault("CH", "0,3")
# reuse gs_all.py's helpers without running its main loop
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", REC, "/dev/null"]; exec(src); sys.argv = sys.argv_saved; OUT = sys.argv[2]
DIR = os.environ.get("DIR", "replay_ch"); BRIDGE = float(os.environ.get("BRIDGE", "0.01")); GATE = float(os.environ.get("GATE", "0.3")); WIN = 300; P = 20480
os.makedirs(OUT, exist_ok=True)
CAT = ["reference, stitched", "reference, split by a hole", "reference, epoch fill", "merged (anchor, gate passed)",
       "merged (edge only)", "rejected by gate", "not joinable", "on air, UE not placed", "TRS, full band", "joined via TRS anchor"]

NFS = 16384
def fit_eps_slope(ref, x, pm):
    """one number and a relative delay slope (NR-Scope timing-window jumps): x ~ c e^{-j s q} ref -> eps"""
    m = np.repeat(pm, 12)
    if m.sum() < 48: return np.nan
    z = np.sum(np.conj(ref) * x * m[None, :], 0)
    F = np.fft.fft(z, NFS); k = int(np.argmax(np.abs(F)))
    a_, b_, c_ = np.abs(F[k - 1]), np.abs(F[k]), np.abs(F[(k + 1) % NFS]); d = 0.5 * (a_ - c_) / (a_ - 2 * b_ + c_ + 1e-30)
    s_ = -2 * np.pi * (k + d) / NFS
    r = ref[:, m] * np.exp(-1j * s_ * q[m])[None, :]; v = x[:, m]; den = np.vdot(r, r).real
    if den <= 0: return np.nan
    cc = np.vdot(r, v) / den; return float(np.sum(np.abs(v - cc * r) ** 2) / np.sum(np.abs(v) ** 2))

def unwrap_to(a, ref):
    return a + round((ref - a) / P) * P

def hexbits(h):
    return np.unpackbits(np.frombuffer(bytes.fromhex(h + "0" * (len(h) % 2)), np.uint8))[:273].astype(bool)

def segs_of(mask_prg):
    return [(s, e) for s, e in runs(mask_prg)]

def do_window(S, wslots, R_ue, out, epstarts, stats, typ):
    """S: rnti -> list of (t,a,pm,Y) in the window. Fills out[(a,prg)] = (cat, ue)"""
    A = sorted(S[R_ue], key=lambda x: x[0])
    EP, TPL = epochs(A)
    EG = {}; MAIN = {}
    for k, (S2, C2) in TPL.items():
        EG[k] = joint_gains(S2 / np.maximum(np.repeat(C2, 12), 1)[None, :], C2 > 0)
        rr = runs(C2 > 0); m = np.zeros(NG, bool)
        if rr:
            s, e = max(rr, key=lambda x: x[1] - x[0]); m[s:e] = True
        MAIN[k] = m
    prev = None
    for i, (t, a, pm, Y) in enumerate(A):
        if EP[i] != prev: epstarts.append(a); prev = EP[i]
    XA = []
    for i, (t, a, pm, Y) in enumerate(A):
        X = Y / np.repeat(EG[EP[i]], 12)[None, :]; X[:, ~np.repeat(pm, 12)] = 0; XA.append((t, a, pm, X))
        m = MAIN[EP[i]]
        for b in np.flatnonzero(pm): out.setdefault((a, b), (0 if m[b] else 1, R_ue))
        idx = np.flatnonzero(pm)
        if len(idx):
            for b in range(idx[0], idx[-1] + 1):
                if not pm[b] and m[b]: out.setdefault((a, b), (2, R_ue))
    TA = np.array([x[0] for x in XA]); bya = {}
    for i, x in enumerate(XA): bya.setdefault(x[1], []).append(i)
    for u, lst in S.items():
        if u == R_ue: continue
        for t, a, pb, Yb in sorted(lst, key=lambda x: x[0]):
            if any((a, b) in out and out[(a, b)][0] in (3, 4) for b in np.flatnonzero(pb)): continue  # other DM-RS symbol already done
            g = joint_gains(Yb, pb); Z = Yb / np.repeat(g, 12)[None, :]; Z[:, ~np.repeat(pb, 12)] = 0
            nr = [k for k in np.flatnonzero(np.abs(TA - t) <= BRIDGE) if XA[k][1] != a or True]
            cs = [(k,) + frame(XA[k][3], XA[k][2], Z, pb) for k in nr]; cs = [x for x in cs if x[2] >= 4 and np.isfinite(x[1])]
            cat = None
            if cs:
                k2, f, n = min(cs, key=lambda x: (abs(TA[x[0]] - t) > 0.002, -x[2]))
                mm = np.repeat(XA[k2][2] & pb, 12); fr = nmse(Z[:, mm] / f, XA[k2][3][:, mm])
                cat = 3 if fr < GATE else 5
                stats["fitres"].append(round(float(fr), 3)); stats["dta"].append(round(float(abs(TA[k2] - t)) * 1e3, 2))
            elif a in bya:
                pa = np.zeros(NG, bool)
                for i in bya[a]: pa |= XA[i][2]
                touch = any(r1[1] == r2[0] or r2[1] == r1[0] for r1 in runs(pa) for r2 in runs(pb))
                cat = 4 if touch else 6
            else:
                cat = 6
            for b in np.flatnonzero(pb):
                if (a, b) not in out or out[(a, b)][0] == 6: out[(a, b)] = (cat, u)
            stats["b_" + typ.get(u, "?")][cat] += int(pb.sum())
    return len(TPL), [round((A[i][0] - A[0][0]) * 1e3, 1) for i in range(1)], EP, A

BURSTS = sorted(glob.glob(REC + "/b*/"), key=lambda s: int(s.rstrip("/").split("/b")[-1]))
SUMMARY = {"cats": CAT, "bridge_ms": BRIDGE * 1e3, "gate": GATE, "bursts": []}
for bi, bd in enumerate(BURSTS):
    if ONLY and bi not in ONLY: continue
    fs = [bd + DIR + "%d/ue_dmrs.bin" % c for c in CH]
    if not all(os.path.exists(f) and os.path.getsize(f) > 20e6 for f in fs): print("skip", bd, flush=True); continue
    Dd = [load(f) for f in fs]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
    S = C.defaultdict(list)
    for k in sorted(keys):
        pm = np.logical_and.reduce([Dd[c][k][0] for c in range(NA)]); Y = np.array([Dd[c][k][1] for c in range(NA)]).astype(complex)
        if pm.sum() and np.all(np.isfinite(Y[:, np.repeat(pm, 12)])): S[k[1]].append((k[0] * 0.5e-3 + k[2] / 28e3, k[0], pm, Y))
    del Dd
    if 0 in S:   # blind: drop the PRGs that a placed UE has in the same slot
        cov0 = C.defaultdict(lambda: np.zeros(NG, bool))
        for u in S:
            if u != 0:
                for x in S[u]: cov0[x[1]] |= x[2]
        S[0] = [(t, a, pm & ~cov0[a], Y) for t, a, pm, Y in S[0] if (pm & ~cov0[a]).sum() >= 1]
        if not S[0]: del S[0]
    cnt = sorted(((len(v), u) for u, v in S.items()), reverse=True)
    if not cnt: continue
    real = [u for n, u in cnt if u != 0]
    if not real: print("skip burst %d: no placed UE" % bi, flush=True); continue
    TAB = real[0]
    # precoder type per UE over the burst
    typ = {}
    for n, u in cnt:
        st = []
        for t, a, pm, Y in S[u][::max(1, len(S[u]) // 40)]: joint_gains(Y, pm, st)
        typ[u] = "per-PRG" if st and np.median(np.abs(np.degrees(np.angle(st)))) > 20 else ("wideband" if st else "?")
    typ[TAB] = "tablet"
    if 0 in S: typ[0] = "blind"
    allslots = sorted({x[1] for u in S for x in S[u]}); base = allslots[0]
    # survey (antenna 0): strict DM-RS bitmap of symbol 2 + periodic signals
    surv = {}; per = {}; hist = [base]
    for line in list(open(bd + DIR + "0/survey.csv"))[1:]:
        r = line.rstrip("\n").split(",")
        if len(r) < 40: continue
        a = unwrap_to(int(r[0]) * 20 + int(r[1]), float(np.median(hist[-64:]))); hist.append(a)
        if len(r[33]) == 69: surv[a] = hexbits(r[33])
        trs = csi = 0
        for e in r[39].split(";"):
            p = e.split(":")
            if len(p) == 4 and float(p[3]) >= 0.6:
                if int(p[0]) in (4, 8): trs = 1
                elif int(p[0]) == 13: csi += 1
        per[a] = (trs, csi)
    a0 = min(min(surv), base); a1 = max(max(surv), allslots[-1])
    Gs = [load_gsc(f) for f in fs]; gk = sorted(set(Gs[0]).intersection(*[set(g) for g in Gs[1:]]))
    TR = []
    for (a, l) in gk:
        if l not in (4, 8): continue
        H = np.array([trs_despread(*Gs[c][(a, l)], a % 20, l, 2)[:N].ravel() for c in range(NA)])
        pw = np.mean(np.abs(H) ** 2, 0).reshape(N, 3).mean(1); TR.append((a * 0.5e-3 + l / 28e3, a, H, pw > np.median(pw) / 10))
    del Gs
    TTR = np.array([x[0] for x in TR])
    J = {"rnti": [], "span": a1 - a0 + 1, "slots": [], "seg": {}, "win": {}, "ep": {}}
    UEI = {}
    def ui(u):
        if u not in UEI: UEI[u] = len(J["rnti"]); J["rnti"].append(["%#x" % u, typ.get(u, "?")])
        return UEI[u]
    ui(TAB)
    for a in sorted(set(surv) | set(per)):
        tr, cs_ = per.get(a, (0, 0)); J["slots"].append([a - a0, tr, cs_])
    BST = {"burst": bi, "tablet": "%#x" % TAB, "n_ues": len(S), "types": C.Counter(typ.values())}
    for sc in ("tab", "notab"):
        out = {}; eps = []; W = []
        stats = {"fitres": [], "dta": [], "b_wideband": [0] * 10, "b_per-PRG": [0] * 10, "b_?": [0] * 10, "b_tablet": [0] * 10, "b_blind": [0] * 10}
        for w0 in range(a0, a1 + 1, WIN):
            Sw = C.defaultdict(list)
            for u in S:
                if sc == "notab" and u == TAB: continue
                for x in S[u]:
                    if w0 <= x[1] < w0 + WIN: Sw[u].append(x)
            dl = [a for a in surv if w0 <= a < w0 + WIN]
            if not Sw:
                W.append([w0 - a0, -1, 0, len(dl), 0, 0, 0]); continue
            cand = [u for u in Sw if u != 0] or list(Sw)
            R = max(cand, key=lambda u: sum(int(x[2].sum()) for x in Sw[u]))
            ne, _, EPw, Aw = do_window(Sw, None, R, out, eps, stats, typ)
            # longest epoch in ms
            el = C.defaultdict(list)
            for e, x in zip(EPw, Aw): el[e].append(x[0])
            longest = max((max(v) - min(v)) * 1e3 for v in el.values())
            used = sum(4 for (a, b), (c, u) in out.items() if w0 <= a < w0 + WIN and c in (0, 2, 3, 4))
            ref_only = sum(4 for (a, b), (c, u) in out.items() if w0 <= a < w0 + WIN and c == 0)
            W.append([w0 - a0, ui(R), used, len(dl), ne, round(longest, 1), ref_only])
        # TRS anchor for what could not be joined coherently
        for u in S:
            if sc == "notab" and u == TAB: continue
            for t, a, pm, Y in S[u]:
                bs = [b for b in np.flatnonzero(pm) if (a, b) in out and out[(a, b)][1] == u and out[(a, b)][0] in (1, 5, 6)]
                if not bs or not len(TTR): continue
                j = int(np.argmin(np.abs(TTR - t)))
                if abs(TTR[j] - t) > BRIDGE: continue
                g = joint_gains(Y, pm); X = Y / np.repeat(g, 12)[None, :]; X[:, ~np.repeat(pm, 12)] = 0
                e = fit_eps_slope(TR[j][2], X, pm)
                if np.isfinite(e) and e < GATE:
                    for b in bs: out[(a, b)] = (9, u)
        # per-window band with the TRS: unique (slot, PRB) used by DM-RS (0, 2, 3, 4, 9) or a TRS symbol
        trs_prb = C.defaultdict(lambda: np.zeros(273, bool))
        for t, a, H, tcov in TR: trs_prb[a][:N] |= tcov
        usedp = C.defaultdict(lambda: np.zeros(273, bool))
        for (a, b), (c, u) in out.items():
            if c in (0, 2, 3, 4, 9): usedp[a][4 * b:4 * b + 4] = True
        for w in W:
            w0 = w[0] + a0
            w.append(int(sum((usedp[a] | trs_prb[a]).sum() for a in set(usedp) | set(trs_prb) if w0 <= a < w0 + WIN)))
        # on-air DM-RS not covered by any placed UE (in notab: the tablet's PRBs are treated as absent)
        cov = C.defaultdict(lambda: np.zeros(273, bool))
        for u in S:
            for x in S[u]: cov[x[1]][:272] |= np.repeat(x[2], 4)
        seg = []
        bycat = C.Counter()
        for t, a, H, tcov in TR:
            if a < a0 or a > a1: continue
            x = np.flatnonzero(tcov)
            if len(x): seg += [a - a0, int(x[0]), int(x[-1]), 8, -1]; bycat[8] += int(len(x)) // 2
        for a in sorted({k[0] for k in out}):
            row = C.OrderedDict()
            for b in range(NG):
                if (a, b) in out: row[b] = out[(a, b)]
            b = 0; keys_ = sorted(row)
            i = 0
            while i < len(keys_):
                b0 = keys_[i]; c, u = row[b0]; j = i
                while j + 1 < len(keys_) and keys_[j + 1] == keys_[j] + 1 and row[keys_[j + 1]] == (c, u): j += 1
                seg += [a - a0, 4 * b0, 4 * keys_[j] + 3, c, ui(u)]; bycat[c] += 4 * (j - i + 1); i = j + 1
        for a, bm in surv.items():
            m = bm & ~cov[a]
            x = np.flatnonzero(m)
            if not len(x): continue
            s = x[0]
            for k in range(1, len(x) + 1):
                if k == len(x) or x[k] != x[k - 1] + 1:
                    if x[k - 1] - s + 1 >= 2: seg += [a - a0, int(s), int(x[k - 1]), 7, -1]; bycat[7] += int(x[k - 1] - s + 1)
                    if k < len(x): s = x[k]
        J["seg"][sc] = seg; J["win"][sc] = W; J["ep"][sc] = sorted(set(e - a0 for e in eps))
        BST[sc] = {"bycat": dict(bycat), "fitres": stats["fitres"], "dta": stats["dta"],
                   "b": {k[2:]: v for k, v in stats.items() if k.startswith("b_")}, "win": W, "dl_slots": len(surv)}
    BST["rnti"] = J["rnti"]
    json.dump(J, open(OUT + "/um_b%d.json" % bi, "w"), separators=(",", ":"))
    json.dump(BST, open(OUT + "/sum_b%d.json" % bi, "w"), separators=(",", ":"), default=lambda o: dict(o) if isinstance(o, C.Counter) else int(o))
    print("burst %d: tablet %#x, %d UEs, tab cats %s, notab cats %s" % (bi, TAB, len(S), dict(BST["tab"]["bycat"]), dict(BST["notab"]["bycat"])), flush=True)
