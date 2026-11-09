"""aggregate gs_trs.py outputs: gs_trs_sum.py <json> ..."""
import sys, json, numpy as np, collections as C
GATE = 0.3; NA = 4
LAGS = [("4 symbols", 0), ("1 slot", 1), ("20 ms", 40), ("40 ms", 80), ("100 ms", 200), ("500 ms", 1000)]
Rz = C.defaultdict(list)
for f in sys.argv[1:]:
    for k, v in json.load(open(f)).items():
        kk = tuple(k.split("|")); Rz[kk] += v
def q(k, f="%.2f"):
    v = np.array(Rz.get(k, []), float); v = v[np.isfinite(v)]
    return (f % np.median(v) + " (n %d)" % len(v)) if len(v) else "-"
print("A  1-|rho|^2 vs lag      TRS                 CSI-RS port 0          port 3               port 7")
for name, _ in LAGS: print("   %-10s %20s %20s %20s %20s" % (name, q(("A", "trs", name), "%.3f"), q(("A", "csi port 0", name), "%.3f"), q(("A", "csi port 3", name), "%.3f"), q(("A", "csi port 7", name), "%.3f")))
print("B  eps against the TRS (median, gate %.1f), and pass rate" % GATE)
for kd in ("tablet", "other", "blind"):
    v = np.array(Rz.get(("B", "trs", kd, "all"), [])); vt = np.array(Rz.get(("B", "tablet", kd, "all"), []))
    print("   %-7s TRS %s pass %s | " % (kd, q(("B", "trs", kd, "all")), "%.0f%%" % (100 * np.mean(v < GATE)) if len(v) else "-") +
          " ".join("%s %s" % (b, q(("B", "trs", kd, b))) for b in ("0-2 ms", "2-5 ms", "5-10 ms")) +
          (" | vs tablet anchor %s pass %.0f%%" % (q(("B", "tablet", kd, "all")), 100 * np.mean(vt < GATE)) if len(vt) else ""))
for kd in ("other", "blind"):
    a_ = sum(Rz.get(("C", kd, "PRG-snapshots"), [])); b_ = sum(Rz.get(("C", kd, "joined via TRS"), []))
    if a_: print("C  %-6s PRG-snapshots joined via TRS: %.0f%% of %d" % (kd, 100 * b_ / a_, a_))
for k in ("DM-RS via TRS", "TRS itself", "both"):
    v = np.array(Rz.get(("C", "band", k), []))
    if len(v): print("C  band per 150 ms window, %-14s median %.1f%%  p90 %.1f%%  zero %.0f%%" % (k, 100 * np.median(v), 100 * np.percentile(v, 90), 100 * np.mean(v < 1e-3)))
print("D  tablet DM-RS per PRG in the span of the 8 CSI-RS ports: residual near %s, CSI-RS 500 ms away %s; random-subspace level %.2f (4 antennas only)" %
      (q(("D", "near")), q(("D", "500 ms away")), 1 - 8 / (NA * 4)))
