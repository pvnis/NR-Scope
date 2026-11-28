import sys, os, numpy as np
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", MINPRB="48")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/stage.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
src = open("/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py").read().split("# DM-RS\nif USE:")[0]
exec(compile(src, "so_light_head", "exec"))
Dd = [load(f) for f in files]; keys = set(Dd[0]).intersection(*[set(d) for d in Dd[1:]])
import collections as CC
cnt = CC.Counter(k[1] for k in keys if k[1] != 0); tab = cnt.most_common(1)[0][0]; print("tablet rnti", hex(tab), cnt[tab], "records on all 4 antennas")
TK = sorted(k for k in keys if k[1] == tab)
def fr(X, R, b):
    sl = slice(12 * b, 12 * b + 12); x = X[:, sl].ravel(); r = R[:, sl].ravel(); g = np.vdot(r, x) / np.vdot(r, r)
    return np.sum(np.abs(x - g * r) ** 2) / (np.sum(np.abs(g * r) ** 2) + 1e-30)
def per_ant(X, R, b):   # per antenna separately
    sl = slice(12 * b, 12 * b + 12); out = []
    for c in range(4):
        x = X[c, sl]; r = R[c, sl]; g = np.vdot(r, x) / np.vdot(r, r); out.append(np.sum(np.abs(x - g * r) ** 2) / (np.sum(np.abs(g * r) ** 2) + 1e-30))
    return np.median(out)
stages = CC.defaultdict(list); pa = []
for k0, k1 in zip(TK, TK[1:]):
    if k1[0] - k0[0] != 1 or k1[2] != k0[2]: continue
    pm = np.logical_and.reduce([d[k0][0] for d in Dd]) & np.logical_and.reduce([d[k1][0] for d in Dd])
    if pm.sum() < 4: continue
    Y = [np.array([d[k][1] for d in Dd]).astype(complex) for k in (k0, k1)]; T_ = [k[0] * 0.5e-3 + k[2] / 28e3 for k in (k0, k1)]
    S = {"raw": Y}
    S["wsh"] = [undo_wsh(y, shifts(k[0]), -1.0) for y, k in zip(Y, (k0, k1))]
    S["drift"] = [dedrift(t, y) for t, y in zip(T_, S["wsh"])]
    S["ramp"] = [y * ramp for y in S["drift"]]
    for nm, (a_, b_) in S.items(): stages[nm] += [fr(b_, a_, b) for b in np.flatnonzero(pm)[:20]]
    pa += [per_ant(S["raw"][1], S["raw"][0], b) for b in np.flatnonzero(pm)[:20]]
    if len(stages["raw"]) > 3000: break
for nm, v in stages.items(): print("tablet, next slot, %-6s: residual median %.1f dB" % (nm, 10 * np.log10(np.median(v))))
print("tablet, next slot, raw, per antenna separately: %.1f dB" % (10 * np.log10(np.median(pa))))
print("first keys", TK[:6])
for k in TK[:4]: print(k, [int(d[k][0].sum()) for d in Dd], [np.round(np.abs(d[k][1][:6]), 3) for d in Dd[:1]])
