"""Usage map of the current pipeline (light stream, per-snapshot TRS reference) from so_light.py's SIDE file (every
candidate symbol and the converter's decision) and sensing_offline's NR_SENSING_TRACE (alignment drop A, random drop R,
pushed H, snapshots of each map M). Writes <out>/f7_b<b>.json per burst.
usage: gen_um7.py <dir with b*_{tab,notab}.side/.trace> <out dir>"""
import sys, os, json, glob, collections as C
D, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
SPS = 61440
# categories: 0 decoded-UE DM-RS in a map, 1 blind DM-RS in a map, 2 TRS in a map, 3 CSI-RS sample in a map,
# 4 dropped by the converter (no TRS within 10 ms / CSI-RS fit), 5 dropped by NR-Scope alignment, 6 NR-Scope random drop,
# 7 kept but in no map, 8 too narrow for NR-Scope (dropped before alignment)
KIND = {"dec": 0, "blind": 1, "trs": 2, "csi": 3}
for sf in sorted(glob.glob(D + "/b*_tab.side")):
    b = int(os.path.basename(sf)[1:].split("_")[0]); R = {"burst": b}
    for sc in ("tab", "notab"):
        side = [json.loads(x) for x in open(D + "/b%d_%s.side" % (b, sc))]
        fate = C.defaultdict(set); per_slot = C.defaultdict(set); maps = []; inmap = set()
        import re
        txt = open(D + "/b%d_%s.trace" % (b, sc)).read()   # two threads write: lines can interleave, parse by token
        for m_ in re.finditer(r"([AHR]) (\d+) (\d+)(?=\D|$)", txt):
            a, t = int(m_.group(2)), int(m_.group(3)); fate[t].add(m_.group(1)); per_slot[a].add(t)
        seen = set()
        for m_ in re.finditer(r"M (\d+)((?: \d+)+)", txt):
            ts = [int(x) for x in m_.group(2).split()]
            if len(ts) != int(m_.group(1)): continue
            inmap.update(ts)
            if (ts[0], len(ts)) not in seen: seen.add((ts[0], len(ts))); maps.append((None, ts[0], ts[-1], len(ts)))
        # symbol time in the record: slot start + 352 + 4384 l (sensing_offline's 4096-point symbols at 122.88 MHz)
        tl = {}; bad = 0
        for a, l, k, u, prgs, r in side:
            if r == "ok": tl[(a, l)] = a * SPS + 352 + 4384 * l
        rn = []; RI = {}
        def ui(u):
            if u not in RI: RI[u] = len(rn); rn.append("%#x" % u if u > 0 else ("blind" if u == 0 else "cell"))
            return RI[u]
        seg = []; cnt = C.Counter()
        for a, l, k, u, prgs, r in side:
            if r != "ok": c = 4
            else:
                t = tl.get((a, l))
                if t not in fate: c = 8        # never reached the alignment: too few pilots for NR-Scope (narrow grant)
                elif "H" in fate[t]: c = KIND[k] if t in inmap else 7
                elif "R" in fate[t]: c = 6
                elif "A" in fate[t]: c = 5
                else: c = 7
            # PRG list -> PRB runs
            i = 0
            while i < len(prgs):
                j = i
                while j + 1 < len(prgs) and prgs[j + 1] == prgs[j] + 1: j += 1
                lo, hi = 4 * prgs[i], min(4 * prgs[j] + 3, 272)
                seg += [a, l, lo, hi, c, ui(u)]; cnt[c] += hi - lo + 1; i = j + 1
        R[sc] = {"seg": seg, "rnti": rn, "cnt": dict(cnt), "unmatched_slots": bad,
                 "maps": [[t0 // SPS, t1 // SPS, n] for _, t0, t1, n in maps]}
    json.dump(R, open(OUT + "/f7_b%d.json" % b, "w"), separators=(",", ":"))
    print("b%d" % b, {sc: (R[sc]["cnt"], len(R[sc]["maps"]), R[sc]["unmatched_slots"]) for sc in ("tab", "notab")}, flush=True)
