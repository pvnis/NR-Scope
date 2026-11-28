"""SNR budget of the tablet DM-RS after anchoring (one chunk): per PRG-symbol, thermal noise (same-slot symbol pair,
same beam: difference / 2 after the per-symbol phase), anchoring residual (symbol minus its beam's local mean over
+-0.25 s, same phase removal), static power. Runs on sens8 with so_light.py in-process (tablet only, BA_PH=prg)."""
import sys, os, io, contextlib, numpy as np, collections as CC
CH = sys.argv[1]
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="prg", BA_W="0", MINPRB="48", TAB_ONLY="1", NO_TRS_OUT="1")
sys.argv = ["so_light.py", CH, "/tmp/snr_tmp.bin"]
g = {"__file__": "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py", "__name__": "x"}; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(g["__file__"]).read(), "so_light", "exec"), g)
if os.path.exists("/tmp/snr_tmp.bin"): os.remove("/tmp/snr_tmp.bin")
E, own, Scnt = g["E"], g["own"], g["Scnt"]
# per-symbol common phase against the beams (as the wide anchoring estimates it)
Ssum = g["Ssum"]; psi = {}
for i, (t, a, l, pm, Y, u) in enumerate(E):
    z = 0
    for b in np.flatnonzero(pm):
        k = own.get((i, b))
        if k is not None and Scnt[k] >= 10: z += np.vdot(Ssum[k] / Scnt[k], Y[:, 12 * b:12 * b + 12])
    psi[i] = np.angle(z) if z != 0 else 0.0
# de-phased PRG values per (beam, PRG), in time order
ser = CC.defaultdict(list)
for (i, b), k in own.items():
    if Scnt[k] >= 10: ser[(k, b)].append(i)
stat, resid, therm, res1, res4 = [], [], [], [], []
idx = {(E[i][1], E[i][2], b): i for (i, b) in own}
for (k, b), ii in ser.items():
    ii.sort(key=lambda q: E[q][0]); T = np.array([E[q][0] for q in ii])
    V = np.array([E[q][4][:, 12 * b:12 * b + 12] * np.exp(-1j * psi[q]) for q in ii])
    for j, q in enumerate(ii):
        m = (abs(T - T[j]) <= 0.25); m[j] = False
        if m.sum() < 3: continue
        S = V[m].mean(0); stat.append(np.mean(abs(S) ** 2)); resid.append(np.mean(abs(V[j] - S) ** 2))
        al = np.vdot(S, V[j]) / np.vdot(S, S).real; res1.append(np.mean(abs(V[j] - al * S) ** 2))          # one complex gain per PRG-symbol
        A4 = (np.sum(np.conj(S) * V[j], 1) / np.sum(abs(S) ** 2, 1))[:, None]; res4.append(np.mean(abs(V[j] - A4 * S) ** 2))   # one per antenna
    # thermal: symbol 2 vs symbol 11 of the same slot, same beam and PRG (same precoder; a target barely moves in 0.3 ms)
    for q in ii:
        a, l = E[q][1], E[q][2]
        if l != 2: continue
        q2 = idx.get((a, 11, b))
        if q2 is None or own.get((q2, b)) != k: continue
        d = E[q][4][:, 12 * b:12 * b + 12] * np.exp(-1j * psi[q]) - E[q2][4][:, 12 * b:12 * b + 12] * np.exp(-1j * psi[q2])
        therm.append(np.mean(abs(d) ** 2) / 2)
stat, resid, therm, res1, res4 = map(np.array, (stat, resid, therm, res1, res4))
db = lambda x: 10 * np.log10(x)
print("chunk %s: %d PRG-symbols, %d same-slot pairs" % (os.path.basename(CH), len(resid), len(therm)))
print("static per RE           median %.1f dB" % db(np.median(stat)))
print("anchoring residual      median %.1f dB  (%.1f dB below static)" % (db(np.median(resid)), db(np.median(stat)) - db(np.median(resid))))
print("thermal (same-slot)     median %.1f dB  (%.1f dB below static)" % (db(np.median(therm)), db(np.median(stat)) - db(np.median(therm))))
print("residual over thermal: template %.1f dB | + 1 gain per PRG-symbol %.1f dB | + 1 gain per antenna %.1f dB (fit keeps 47/48 and 44/48 of the noise: expect +0.1 / -0.4 dB)" % (db(np.median(resid)) - db(np.median(therm)), db(np.median(res1)) - db(np.median(therm)), db(np.median(res4)) - db(np.median(therm))))
