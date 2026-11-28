"""Per kind of DM-RS (tablet, other decoded UEs, blind), against the nearest TRS (<= 10 ms, antenna 0), per PRG:
  power   DM-RS power over the TRS power on the same REs (dB): how strongly that UE's beams reach the receiver
  fit     inside one PRG, the residual of one complex number mapping the TRS onto the DM-RS (noise + model)
  jump    between adjacent PRGs of one snapshot, the change of that number (phase in deg, |log ratio| in dB), after
          removing the snapshot's median step (a timing slope): how much the precoder changes from PRG to PRG
  same    share of adjacent PRG pairs whose numbers agree within 20 deg and 2 dB (the same precoder, probably)
Real or simulated burst dirs (same layout). usage: prg_stats.py <burst dir> [label] [tablet rnti]"""
import sys, os, numpy as np, collections as C
sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
os.environ.setdefault("CH", "0")
src = open(os.path.expanduser("/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/gs_all.py")).read().split("Rz=C.defaultdict(list)")[0]
saved = list(sys.argv); sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
from gs_trs_lib import load_gsc, trs_despread
bd = sys.argv[1].rstrip("/"); lab = sys.argv[2] if len(sys.argv) > 2 else bd
TABR = int(sys.argv[3], 0) if len(sys.argv) > 3 else 0x7304
f = bd + "/replay_blind_ch0/ue_dmrs.bin"
D = load(f); G = load_gsc(f)
TR = []
for (a, l), (fc, z) in G.items():
    if l in (4, 8) and a % 20 in (5, 6) and fc == 0 and len(z) >= 12 * 272:
        TR.append((a * 0.5e-3 + l / 28e3, trs_despread(fc, z[:12 * 272], a % 20, l, 2)[:272].ravel()))
TR.sort(key=lambda x: x[0]); TT = np.array([x[0] for x in TR])
S = C.defaultdict(lambda: dict(power=[], fit=[], jph=[], jdb=[], same=[], nprg=[]))
for (a, rn, l), (pm, Y) in D.items():
    t = a * 0.5e-3 + l / 28e3; j = int(np.argmin(np.abs(TT - t)))
    if abs(TT[j] - t) > 0.01 or pm.sum() < 1: continue
    T = TR[j][1]; kind = "tablet" if rn == TABR else ("blind" if rn == 0 else "other")
    k = S[kind]; gs = []
    for b in np.flatnonzero(pm):
        sl = slice(12 * b, 12 * b + 12); y, tr = Y[sl], T[sl]
        den = np.vdot(tr, tr).real
        if den <= 0 or not np.all(np.isfinite(y)): gs.append(np.nan); continue
        gb = np.vdot(tr, y) / den
        k["power"].append(10 * np.log10(np.vdot(y, y).real / den + 1e-30))
        k["fit"].append(np.sum(np.abs(y - gb * tr) ** 2) / max(np.vdot(y, y).real, 1e-30))
        gs.append(gb)
    k["nprg"].append(int(pm.sum()))
    idx = np.flatnonzero(pm); gs = np.array(gs)
    adj = [(gs[i + 1] / gs[i]) for i in range(len(idx) - 1) if idx[i + 1] == idx[i] + 1 and np.isfinite(gs[i]) and np.isfinite(gs[i + 1]) and gs[i] != 0]
    if adj:
        adj = np.array(adj); adj = adj / np.exp(1j * np.angle(np.median(adj.real) + 1j * np.median(adj.imag)))   # remove the slope
        ph = np.degrees(np.abs(np.angle(adj))); db = np.abs(20 * np.log10(np.abs(adj) + 1e-30))
        k["jph"] += ph.tolist(); k["jdb"] += db.tolist(); k["same"] += ((ph < 20) & (db < 2)).tolist()
print("== %s (%d TRS)" % (lab, len(TR)))
for kind in ("tablet", "other", "blind"):
    k = S[kind]
    if not k["power"]: print("  %-6s none" % kind); continue
    print("  %-6s PRGs/snapshot %4.1f | power vs TRS median %+5.1f dB | fit inside a PRG %.3f | adjacent PRG jump %5.1f deg %4.1f dB | same precoder %3.0f%%"
          % (kind, np.median(k["nprg"]), np.median(k["power"]), np.median(k["fit"]), np.median(k["jph"]) if k["jph"] else np.nan,
             np.median(k["jdb"]) if k["jdb"] else np.nan, 100 * np.mean(k["same"]) if k["same"] else np.nan))
