"""Page data of the usage map (same format as the 5 Nov page): per chunk um_b<k>.json (span, slots [a, trs, csi],
before [a, lo, hi] = DM-RS on air of every UE, sc{tab,notab}{seg, rnti, maps} from gen_um7.py re-based to the first
slot of the chunk) and data.json (agg, rows, bursts) for the page header and tables.
usage: build_um_page.py <chunks dir> <usemap dir (with json/f7_b*.json)> <out dir>"""
import sys, os, json, glob, numpy as np, collections as C
sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
CHD, UMD, ODIR = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(ODIR, exist_ok=True)
src = open("/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/gs_all.py").read().split("Rz=C.defaultdict(list)")[0]
saved = list(sys.argv); sys.argv = ["gs_all.py", "/dev/null", "/dev/null"]; exec(src); sys.argv = saved
from gs_trs_lib import load_gsc
agg = {"tab": [0] * 9, "notab": [0] * 9}; rows = []; bursts = []
for f in sorted(glob.glob(UMD + "/json/f7_b*.json"), key=lambda p: int(p.split("_b")[-1][:-5])):
    b = int(f.split("_b")[-1][:-5]); F = json.load(open(f))
    dump = CHD + "/b%d/replay_blind_ch0/ue_dmrs.bin" % b
    D = load(dump); G = load_gsc(dump)
    slots_all = sorted({k[0] for k in D} | {a for a, l in G})
    if not slots_all: continue
    A0 = slots_all[0]; span = slots_all[-1] - A0 + 1
    trs = {a for a, l in G if l in (4, 8)}; csi = {a for a, l in G if l == 13}
    slots = [[a - A0, int(a in trs), int(a in csi)] for a in slots_all]
    on = C.defaultdict(lambda: np.zeros(68, bool))
    for (a, u, l), v in D.items(): on[a] |= v[0]
    before = []
    for a in sorted(on):
        p = np.flatnonzero(on[a]); i = 0
        while i < len(p):
            j = i
            while j + 1 < len(p) and p[j + 1] == p[j] + 1: j += 1
            before += [int(a - A0), int(4 * p[i]), int(min(4 * p[j] + 3, 272))]; i = j + 1
    sc = {}; row = {"b": b, "hit": 1.0}
    for m in ("tab", "notab"):
        S = F[m]; seg = list(S["seg"])
        for k in range(0, len(seg), 6): seg[k] -= A0
        maps = [[t0 - A0, t1 - A0, n] for t0, t1, n in S["maps"]]
        sc[m] = {"seg": seg, "rnti": S["rnti"], "maps": maps}
        c = [int(S["cnt"].get(str(i), S["cnt"].get(i, 0))) for i in range(9)]
        for i in range(9): agg[m][i] += c[i]
        row[m] = {"c": c, "maps": len(maps), "snap": float(np.median([x[2] for x in maps])) if maps else 0.0}
    json.dump({"span": int(span), "slots": slots, "before": before, "sc": sc}, open(ODIR + "/um_b%d.json" % b, "w"), separators=(",", ":"))
    rows.append(row); bursts.append({"b": b, "span": int(span), "ues": int(any(r not in ("blind", "cell") for r in F["tab"]["rnti"]))})
    print("b%d span %d, %d DL slots, TRS %d, CSI %d, maps %d / %d" % (b, span, len(slots), len(trs), len(csi), row["tab"]["maps"], row["notab"]["maps"]), flush=True)
json.dump({"agg": agg, "rows": rows, "bursts": bursts}, open(ODIR + "/data.json", "w"), separators=(",", ":"))
print("wrote", ODIR)
