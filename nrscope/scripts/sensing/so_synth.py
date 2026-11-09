"""Sanity check for the NR-Scope chain: a synthetic static channel (direct path + 10 weak scatterers, no precoding,
4 antennas, full band) on the DL slots of a real burst, written as a record_estimates file like gs_to_nrsr.py
(pair values on REs 12q + 4n and + 2). Injected targets are then added by sensing_offline.
usage: so_synth.py <burst dir> <out.bin>"""
import sys, os, struct, numpy as np, csv
BD, OUTF = sys.argv[1].rstrip("/") + "/", sys.argv[2]
rng = np.random.default_rng(1); NSC = 3276; NA = 4; N = 273
slots = sorted({int(r[0]) * 20 + int(r[1]) for r in csv.reader(list(open(BD + "replay_ch0/survey.csv"))[1:]) if len(r) > 2})
base = slots[0]; slots = [a + 20480 if a < base - 1000 else a for a in slots]
k = np.arange(NSC); paths = [(1.0, 0.0, 0.2)] + [(10 ** (rng.uniform(-35, -15) / 20) * np.exp(1j * rng.uniform(0, 6.28)), rng.uniform(10, 120) / 3e8, rng.uniform(-1, 1)) for _ in range(10)]
H = sum(a * np.exp(1j * np.pi * np.arange(NA) * np.sin(th))[:, None] * np.exp(-2j * np.pi * k * 30e3 * tau)[None, :] for a, tau, th in paths)
valid = np.zeros(NSC, bool)
for q in range(N):
    for n in range(3): valid[12 * q + 4 * n] = valid[12 * q + 4 * n + 2] = True
with open(OUTF, "wb") as f:
    f.write(struct.pack("<8sIIQdIIII", b"NRSEST1", 1, NA, 3750000000, 122.88e6, 4096, 30000, NSC, 0))
    for a in slots:
        if a % 10 in (4, 9): continue   # DDDSU: no DL DM-RS in the U slots
        for c in range(NA):
            f.write(struct.pack("<IHHBBBBHHIIQ", 0x4E525352, c, 1, 0, 1, 5, 0, (1 << 2) | (1 << 11), 0, (a // 20) % 1024, a % 20, a))
            for l in (2, 11):
                y = H[c] + 10 ** (-25 / 20) / np.sqrt(2) * (rng.standard_normal(NSC) + 1j * rng.standard_normal(NSC))
                f.write(struct.pack("<I", int(valid.sum()))); f.write(np.packbits(valid, bitorder="little").tobytes()); f.write(y[valid].astype(np.complex64).tobytes())
print("wrote", OUTF, len(slots), "slots")
