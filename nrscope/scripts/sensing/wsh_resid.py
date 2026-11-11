"""Residual phase slope between TRS occasions 20 ms apart vs the logged window move between them (after the undo)."""
import sys, os, numpy as np
here = "/home/sens/field/"; BD = sys.argv[1]
src = open(here + "so_light.py").read().split("# CSI-RS as extra TRS-beam samples")[0]
sys.argv = ["so_light.py", BD, "/dev/null"]; os.environ.pop("INJ", None); __file__ = here + "so_light.py"; exec(src)
T = {(a, l): R for t, a, l, R in TRS}; rows = []
for (a, l), R in T.items():
    if (a + 40, l) not in T: continue
    R2 = T[(a + 40, l)]; P = R2 * np.conj(R)
    for c in range(4):
        s = np.angle(np.sum(P[c, 1:] * np.conj(P[c, :-1])))            # residual phase per pilot step
        ds = shifts(a + 40)[c] - shifts(a)[c]
        rows.append((c, ds, s))
r = np.array(rows, float)
for c in range(4):
    v = r[r[:, 0] == c]; mv = v[:, 1] != 0
    if mv.sum() < 3: continue
    k = np.polyfit(v[mv, 1], v[mv, 2], 1)
    print("ant %d: moves %d | residual slope vs move: %.5f rad/pilot per sample (undo uses %.5f) | corr %.2f | residual std no-move %.4f, move %.4f" % (
        c, mv.sum(), k[0], PH, np.corrcoef(v[mv, 1], v[mv, 2])[0, 1], np.std(v[~mv, 2]), np.std(v[mv, 2])))
    hist = {}
    for ds, s in v[mv][:, 1:]: hist.setdefault(int(ds), []).append(s)
    print("   per move size (samples): " + " ".join("%+d:%+.3f(n%d)" % (d, np.median(x), len(x)) for d, x in sorted(hist.items()) if len(x) >= 3))
# per-group residual phase pattern of move pairs: where across the band does it break?
print("--- per-PRG residual phase (deg) vs PRG, median over pairs, by move size (ant 0)")
pat = {}
for (a, l), R in T.items():
    if (a + 40, l) not in T: continue
    R2 = T[(a + 40, l)]; ds = shifts(a + 40)[0] - shifts(a)[0]
    z = np.array([np.vdot(R[0, 12 * b:12 * b + 12], R2[0, 12 * b:12 * b + 12]) for b in range(NG)])
    z = z * np.exp(-1j * np.angle(z.sum())); pat.setdefault(ds, []).append(np.degrees(np.angle(z)))
for ds in sorted(pat):
    if len(pat[ds]) < 4: continue
    m = np.median(np.array(pat[ds]), 0)
    print("%+3d n%-3d " % (ds, len(pat[ds])) + " ".join("%+4.0f" % x for x in m[::4]))
print("--- per-pair coherence, ant 0 and all antennas jointly")
rows = []
for (a, l), R in sorted(T.items()):
    if (a + 40, l) not in T: continue
    R2 = T[(a + 40, l)]
    z0 = np.array([np.vdot(R[0, 12 * b:12 * b + 12], R2[0, 12 * b:12 * b + 12]) for b in range(NG)])
    zA = np.array([np.vdot(R[:, 12 * b:12 * b + 12], R2[:, 12 * b:12 * b + 12]) for b in range(NG)])
    zc = [abs(np.array([np.vdot(R[c, 12 * b:12 * b + 12], R2[c, 12 * b:12 * b + 12]) for b in range(NG)]).sum()) for c in range(4)]
    sa, sb = np.array(shifts(a)), np.array(shifts(a + 40))
    rows.append((a, l, abs(z0.sum()) / np.abs(z0).sum(), abs(zA.sum()) / np.abs(zA).sum(), int((sa != sb).any()), int(((sb - sa) != (sb - sa)[0]).any())))
r = np.array(rows)
for f, name in ((r[:, 4] == 0, "no move on any antenna"), ((r[:, 4] == 1) & (r[:, 5] == 0), "same move on all antennas"), (r[:, 5] == 1, "different moves per antenna")):
    if f.any(): print("%-28s n %3d | coh ant0 %.3f | coh 4 antennas jointly %.3f" % (name, f.sum(), np.median(r[f, 2]), np.median(r[f, 3])))
