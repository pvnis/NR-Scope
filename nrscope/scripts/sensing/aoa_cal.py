"""Chain phase calibration (aoa_cal_deg) from the TRS direct path: per chunk, the TRS static (so_light's corrections:
delay per antenna puts the direct path at range 0, drift, oscillator) range-gated at 0 (sum over the band) per chain;
expected phases of a lambda/2 ULA (axis = boresight + 90 deg, element 0 reference, NR-Scope's sign: exp(+j 2 pi d m s))
for the gNB direction in 3D; cal_m = arg(h_m / h_0) - 2 pi d m s. Checked on walk1 against the lidar walker: bias 29 deg -> 2.5 deg.
usage: aoa_cal.py <chunks dir> <session.json> [chunks, default: 5 spread over the recording] [cart dx,dy,dz]
Prints the per-chunk values and, last, the line "aoa_cal_deg: [...]" (exit 2 if the kept chunks disagree by > 15 deg at p90, or fewer than 3 are kept; chunks with an abnormal per-chain direct-path power pattern, > 3 dB, are dropped)."""
import sys, os, io, json, glob, re, contextlib, numpy as np
CHD = sys.argv[1]; S = json.load(open(sys.argv[2])); Li = S["lidar"]
ALL = sorted(int(re.match(r"b(\d+)$", os.path.basename(d)).group(1)) for d in glob.glob(CHD + "/b*") if re.match(r"b(\d+)$", os.path.basename(d)))
KS = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 and sys.argv[3] else [ALL[i] for i in np.unique(np.linspace(0, len(ALL) - 1, 5).astype(int))]
CO = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [-0.445, -0.018, -0.07]
KY = 6378137 * np.pi / 180; KX = KY * np.cos(np.radians(Li["lat"])); hd = np.radians(Li["heading"])
L2E = lambda x, y: x * np.array([np.cos(hd), -np.sin(hd)]) + y * np.array([np.sin(hd), np.cos(hd)])
RX = np.r_[L2E(CO[0], CO[1]), Li["h"] + CO[2]]; TX = np.array([(S["gnb"]["lon"] - Li["lon"]) * KX, (S["gnb"]["lat"] - Li["lat"]) * KY, S["gnb"]["h"]])
b = np.radians(90 - Li["heading"]); ax = np.array([-np.sin(b), np.cos(b), 0.0]); u = (TX - RX) / np.linalg.norm(TX - RX); s = float(u @ ax)
print("gNB from the receiver: sin(angle) on the array axis = %.3f (%.1f deg off boresight, elevation %.1f deg)" % (s, np.degrees(np.arcsin(s)), np.degrees(np.arcsin(u[2]))))
CAL = []; PWR = []
for k in KS:
    os.environ.update(USE_DMRS="0", USE_CSI="0", MINPRB="48")
    sys.argv = ["so_light.py", "%s/b%d" % (CHD, k), "/tmp/aoa_cal_%d.bin" % os.getpid()]; g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(open(g["__file__"]).read(), "so_light", "exec"), g)
    TRS = g["TRS"]; t_ = np.array([x[0] for x in TRS]); R0 = TRS[0][3]
    po = np.polyfit(t_, np.unwrap([np.angle(np.vdot(R0, x[3])) for x in TRS]), 3)
    Rb = np.mean([x[3] * np.exp(-1j * np.polyval(po, x[0])) for x in TRS], 0)          # 4 x L, direct path at range 0
    h = Rb.sum(1); m = np.arange(4)
    cal = np.degrees(np.angle(h / h[0]) - 2 * np.pi * 0.5 * m * s); cal = (cal + 180) % 360 - 180; CAL.append(cal); PWR.append(20 * np.log10(np.abs(h) / np.abs(h).max()))
    print("chunk %2d: direct-path power per chain %s dB; cal_deg %s" % (k, " ".join("%.1f" % v for v in 20 * np.log10(np.abs(h) / np.abs(h).max())), " ".join("%7.1f" % v for v in cal)), flush=True)
# robust: a chunk whose per-chain direct-path powers depart from the recording's median pattern by > 3 dB is dropped
# (walk1 chunk 4: one chain 22 dB down, something blocked that antenna); circular median of the rest; spread = p90
C = np.array(CAL); Pw = np.array(PWR); dev_p = np.abs(Pw - np.median(Pw, 0)).max(1); PDEV = float(os.environ.get("AOA_MAX_POWER_DEV_DB", "3.0")); SMAX = float(os.environ.get("AOA_MAX_SPREAD_DEG", "15")); keep = dev_p <= PDEV
for k_, d_ in zip(KS, dev_p):
    if d_ > PDEV: print("chunk %d dropped: direct-path power pattern off by %.1f dB" % (k_, d_))
Z = np.exp(1j * np.radians(C[keep])); med = np.degrees(np.angle(np.median(Z.real, 0) + 1j * np.median(Z.imag, 0)))
dev = np.degrees(np.abs(np.angle(np.exp(1j * np.radians(C[keep] - med)))))
spr = np.percentile(dev, 90, axis=0)
print("kept %d of %d chunks; deviation per chain p90 %s deg, max %s deg" % (keep.sum(), len(KS), " ".join("%.1f" % v for v in spr), " ".join("%.1f" % v for v in dev.max(0))))
print("aoa_cal_deg: [%s]" % ", ".join("%.2f" % v for v in med))
try: os.remove("/tmp/aoa_cal_%d.bin" % os.getpid())
except OSError: pass
sys.exit(2 if (spr.max() > SMAX or keep.sum() < 3) else 0)
