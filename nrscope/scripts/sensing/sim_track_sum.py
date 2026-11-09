import json, glob, numpy as np, collections as C
for sc in ("notab", "tab"):
    R = C.defaultdict(list)
    for f in glob.glob("/home/sens/field/sim/trk_%s_[0-9].json" % sc):
        for k, v in json.load(open(f)).items(): R[k] += v
    if not R: continue
    nb = len(R.get("%s|ideal|gospa" % sc, []))
    print("== %s (%d bursts) ==" % ("without the tablet" if sc == "notab" else "with the tablet", nb))
    print("%-16s %8s %12s %14s %10s %8s" % ("variant", "found", "false det/map", "false trk/burst", "coverage", "GOSPA"))
    for v in ("ideal", "DMRS oracle", "TRS", "TRS+DMRS@TRS", "TRS+incoh", "CSI-RS", "DMRS today", "DMRS today+TRS"):
        if R.get("%s|%s|gospa" % (sc, v)):
            g = lambda k: R.get("%s|%s|%s" % (sc, v, k), [0])
            print("%-16s %7.0f%% %12.1f %14.1f %9.0f%% %8.1f" % (v, 100 * np.mean(g("found")), np.mean(g("false det")), np.mean(g("false tracks")), 100 * np.mean(g("coverage")), np.mean(g("gospa"))))
    print("DM-RS pieces joined via TRS: %.0f%%" % (100 * np.mean(R["%s|joined via TRS share" % sc])))
