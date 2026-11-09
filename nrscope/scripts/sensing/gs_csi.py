"""6f step 3: is every UE's PDSCH DM-RS "8 CSI-RS ports x 8 weights" per PRG? (the precoder w_b estimated from the CSI-RS)
Values: the blind DM-RS records (rnti 0, UE-independent, all 4 antennas); owner of each PRG from the decoded grants
on antennas 0 and 3: tablet, other placed UE, or unplaced (blind only). TRS symbols tested the same way.
Per PRG: y (4 antennas x 4 PRBs x 3 pairs) ~ sum_p H[a, p, PRB] w_p, H from the nearest CSI-RS (flat in a PRB),
8 unknown weights, LS residual. Compared with: the CSI-RS 500 ms away (stale), a random 8-dim basis level
(1 - 8/16 = 0.5), and one weight only (best single port). By time distance to the CSI-RS.
usage: gs_csi.py <rec dir> <out.json> <bursts>"""
import sys, os, json, numpy as np, collections as C
REC, OUTF = sys.argv[1], sys.argv[2]; ONLY = [int(x) for x in sys.argv[3:]]; saved = list(sys.argv)
os.environ["CH"] = "0,1,2,3"
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
sys.argv = ["gs_all.py", REC, "/dev/null"]; exec(src); sys.argv = saved
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from gs_trs_lib import *
DIR = "replay_blind_ch"

def fit8(y, P):
    """y: NA x 4 x 3, P: NA x 4 x 8 (CSI-RS per PRB) -> residual with 8 weights, and with the best single port"""
    A_ = np.repeat(P[:, :, None, :], 3, 2).reshape(-1, 8); yv = y.reshape(-1); e0 = np.sum(np.abs(yv) ** 2)
    if e0 <= 0: return np.nan, np.nan
    w, *_ = np.linalg.lstsq(A_, yv, rcond=None); r8 = np.sum(np.abs(yv - A_ @ w) ** 2) / e0
    r1 = min(np.sum(np.abs(yv - A_[:, p] * (np.vdot(A_[:, p], yv) / np.vdot(A_[:, p], A_[:, p]))) ** 2) / e0 for p in range(8))
    return float(r8), float(r1)

Rz = C.defaultdict(list)
BURSTS = sorted(glob.glob(REC + "/b*/"), key=lambda s: int(s.rstrip("/").split("/b")[-1]))
for bi, bd in enumerate(BURSTS):
    if ONLY and bi not in ONLY: continue
    fs = [bd + DIR + "%d/ue_dmrs.bin" % c for c in range(4)]
    if not all(os.path.exists(f) for f in fs): print("skip", bd); continue
    Dd = [load(f) for f in fs]
    # owners from antennas 0 and 3 (both decode every placed UE)
    own = {}
    cnt = C.Counter(k[1] for k in Dd[0] if k[1] != 0); tab = cnt.most_common(1)[0][0]
    for k, (pm, Y) in Dd[0].items():
        if k[1] == 0 or k not in Dd[3]: continue
        pm2 = pm & Dd[3][k][0]
        for b in np.flatnonzero(pm2): own[(k[0], b)] = "tablet" if k[1] == tab else "other"
    keys = [k for k in Dd[0] if k[1] == 0 and all(k in d for d in Dd[1:])]
    Gs = [load_gsc(f) for f in fs]; gk = sorted(set(Gs[0]).intersection(*[set(g) for g in Gs[1:]]))
    CS = []; TRl = []
    for (a, l) in gk:
        if l == 13: CS.append((a * 0.5e-3 + l / 28e3, np.array([csi8_despread(*Gs[c][(a, l)], a % 20)[:N] for c in range(4)])))
        elif l in (4, 8): TRl.append((a * 0.5e-3 + l / 28e3, np.array([trs_despread(*Gs[c][(a, l)], a % 20, l, 2)[:N] for c in range(4)])))
    del Gs
    TC = np.array([x[0] for x in CS])
    def test(t, y_of_b, bs, kind):
        j = int(np.argmin(np.abs(TC - t))); dt = abs(TC[j] - t)
        if dt > 0.01: return
        jn = int(np.argmin(np.abs(TC - t - 0.5))); stale = CS[jn][1] if abs(TC[jn] - t - 0.5) < 0.05 else None
        bin_ = "<1 ms" if dt <= 0.001 else "1-3 ms" if dt <= 0.003 else "3-10 ms"
        for b in bs:
            y = y_of_b(b); r8, r1 = fit8(y, CS[j][1][:, 4 * b:4 * b + 4, :])
            if not np.isfinite(r8): continue
            Rz[(kind, bin_, "8 ports")].append(r8); Rz[(kind, bin_, "1 port")].append(r1)
            if stale is not None and bin_ == "<1 ms":
                Rz[(kind, "<1 ms", "stale 8 ports")].append(fit8(y, stale[:, 4 * b:4 * b + 4, :])[0])
    for k in keys:
        a, _, l = k; t = a * 0.5e-3 + l / 28e3
        pm = np.logical_and.reduce([d[k][0] for d in Dd]); Y = np.array([d[k][1] for d in Dd]).astype(complex)
        bykind = C.defaultdict(list)
        for b in np.flatnonzero(pm): bykind[own.get((a, b), "unplaced")].append(b)
        for kind, bs in bykind.items():
            test(t, lambda b: Y[:, 12 * b:12 * b + 12].reshape(4, 4, 3), bs[::2], kind)
    for t, H in TRl[::2]:
        test(t, lambda b: H[:, 4 * b:4 * b + 4, :], list(range(0, NG, 3)), "TRS")
    print("burst %d: %d blind DM-RS records, %d CSI-RS, %d TRS" % (bi, len(keys), len(CS), len(TRl)), flush=True)
json.dump({"|".join(k): v for k, v in Rz.items()}, open(OUTF, "w"))
def q(k):
    v = np.array(Rz.get(k, []), float); return "%.2f (n %d)" % (np.median(v), len(v)) if len(v) else "-"
print("residual per PRG (median); random 8-dim basis = 0.50")
for kind in ("tablet", "other", "unplaced", "TRS"):
    for bin_ in ("<1 ms", "1-3 ms", "3-10 ms"):
        print("  %-8s %-7s 8 ports %-16s 1 port %-16s%s" % (kind, bin_, q((kind, bin_, "8 ports")), q((kind, bin_, "1 port")),
              ("  stale CSI-RS (500 ms) %s" % q((kind, "<1 ms", "stale 8 ports"))) if bin_ == "<1 ms" else ""))
