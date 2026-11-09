import sys, json, glob, numpy as np, collections as C
for mix in ("0", "0.2"):
    for sc in ("notab", "tab"):
        R = C.defaultdict(list)
        for f in glob.glob("/home/sens/field/sim/%s_m%s_[0-9].json" % (sc, mix)):
            for k, v in json.load(open(f)).items(): R[k] += v
        if not R: continue
        print("== %s, precoder beam mix %s ==" % ("without the tablet" if sc == "notab" else "with the tablet", mix))
        for k, v in R.items():
            if "|diag|" in k: print("   diag: %s median %.3f" % (k.split("|")[2], np.median(v)))
        print("%-16s %8s %14s %10s %6s" % ("variant", "found", "false/window", "SNR dB", "n"))
        for nm in ("ideal", "DMRS oracle", "TRS", "TRS+DMRS@TRS", "TRS+incoh", "CSI-RS", "DMRS today", "DMRS today+TRS"):
            v = R.get("%s|%s|hit" % (sc, nm))
            if v: print("%-16s %7.0f%% %14.1f %10.1f %6d" % (nm, 100 * np.mean(v), np.mean(R["%s|%s|false" % (sc, nm)]), np.median(R["%s|%s|snr" % (sc, nm)]), len(v)))
        print("DM-RS pieces joined via TRS: %.0f%%" % (100 * np.mean(R["%s|joined via TRS share" % sc])))
