"""Where does the DM-RS's random per-snapshot factor come from? Tablet DM-RS of one PRG and beam: complex factor alpha
against the beam's static shape (after the TRS oscillator model), compared between symbol 2 and symbol 11 of the
SAME slot (same transmission, same precoder) and between symbol 2 of consecutive slots (two transmissions).
Stable within a slot but random across slots = the gNB's per-transmission precoder factor; random within a slot too
= introduced by our chain."""
import sys, os, io, contextlib, numpy as np, collections as CC
os.environ.update(MODE="tab", USE_DMRS="1", USE_CSI="0", BA="1", BA_PH="cplx", BA_W="0", MINPRB="48", TAB_ONLY="1")
sys.argv = ["so_light.py", "/mnt/mydrive/field/20261124_140050_sunrise_632_4ch_cont_tabletdl_iq_chunks/b3", "/tmp/e2e/sp.bin"]
__file__ = "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing/so_light.py"; sys.path.insert(0, "/home/sens/GIOVANNI/NR-Scope/nrscope/scripts/sensing")
with contextlib.redirect_stdout(io.StringIO()): exec(compile(open(__file__).read(), "so_light", "exec"))
idx = {(E[i][1], E[i][2]): i for i in range(len(E))}          # (slot, symbol) -> snapshot index
same, nxt, mag_s, mag_n = [], [], [], []
for (a, l), i in idx.items():
    if l != 2: continue
    for (a2, l2), lab in (((a, 11), "same"), ((a + 1, 2), "next")):
        j = idx.get((a2, l2))
        if j is None: continue
        for b in np.flatnonzero(E[i][3] & E[j][3]):
            if (i, b) in ALPHA and (j, b) in ALPHA and own[(i, b)] == own[(j, b)]:
                r = ALPHA[(j, b)] / ALPHA[(i, b)]
                (same if lab == "same" else nxt).append(np.angle(r)); (mag_s if lab == "same" else mag_n).append(20 * np.log10(abs(r)))
for lab, v, m in (("symbol 2 -> 11, same slot", same, mag_s), ("symbol 2 -> 2, next slot", nxt, mag_n)):
    v = np.array(v); m = np.array(m)
    print("%-28s %6d PRG pairs: phase change MAD %.2f rad (|mean phasor| %.2f), amplitude change MAD %.2f dB"
          % (lab, len(v), 1.4826 * np.median(np.abs(np.angle(np.exp(1j * (v - np.angle(np.mean(np.exp(1j * v)))))))), abs(np.mean(np.exp(1j * v))), 1.4826 * np.median(np.abs(m - np.median(m)))))
