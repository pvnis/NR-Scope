"""Echo transfer of the anchoring: two so_light outputs of the same chunk, with and without an injected target
(<inj.bin>.truth.csv). Their difference is the echo as processed; per snapshot-chain it is compared with the ideal
echo eps * static(k) * exp(-2j pi (R/lam + k df R/c)) (static = the TRS static channel of that chain, de-rotated as the
output): gain |gamma| (1 = kept), and gamma's phase over time (should follow the output's own rotation, i.e. be the
same for TRS and DM-RS snapshots at the same time).  usage: echo_tf.py <inj.bin> <noinj.bin> <level dB>"""
import sys, collections, csv, numpy as np
exec(open("/tmp/inj/mf_check.py").read().split("for path in sys.argv[1:]:")[0])
P1, P0, LV = sys.argv[1], sys.argv[2], float(sys.argv[3]); eps = 10 ** (LV / 20)
S1 = {(a, l, c): (vb, x) for a, l, c, vb, x in read(P1)}; S0 = {(a, l, c): (vb, x) for a, l, c, vb, x in read(P0)}
tr = {int(r["slot_abs"]): float(r["range_m"]) for r in csv.DictReader(open(P1 + ".truth.csv"))}; A = sorted(tr); rA = np.array([tr[q] for q in A])
kk = np.arange(NSC) - NSC // 2; G = collections.defaultdict(list)
for key in sorted(set(S1) & set(S0)):
    a, l, c = key; vb1, x1 = S1[key]; vb0, x0 = S0[key]; vb = vb1 & vb0
    if vb.sum() < 100: continue
    d = (x1 - x0)[vb]; rt = np.interp(a, A, rA)
    ideal = eps * x0[vb] * np.exp(-2j * np.pi * rt / LAM) * np.exp(-2j * np.pi * kk[vb] * DF * rt / 3e8)
    g = np.vdot(ideal, d) / np.vdot(ideal, ideal).real; coh = abs(np.vdot(ideal, d)) / (np.linalg.norm(ideal) * np.linalg.norm(d) + 1e-30)
    G["trs" if l in (4, 8) else "dmrs"].append((a, l, c, g, coh))
for k, v in G.items():
    g = np.array([x[3] for x in v]); co = np.array([x[4] for x in v])
    # phase consistency in time: gamma's phase against the TRS gamma of the same chain nearest in time
    print("%-4s %5d snapshot-chains: |gamma| median %.2f (p25 %.2f, p75 %.2f); shape coherence median %.2f; gamma phase MAD %.2f rad"
          % (k, len(v), np.median(abs(g)), np.percentile(abs(g), 25), np.percentile(abs(g), 75), np.median(co), 1.4826 * np.median(np.abs(np.angle(g / np.exp(1j * np.median(np.angle(g))))))))
# DM-RS gamma phase relative to the TRS gamma at the nearest time (same chain): coherence between the two streams
T = collections.defaultdict(list)
for a, l, c, g, co in G["trs"]: T[c].append((a, g))
rel = []
for a, l, c, g, co in G["dmrs"]:
    L = T[c]; j = int(np.argmin([abs(x[0] - a) for x in L])); rel.append(np.angle(g / L[j][1]))
rel = np.array(rel); print("DM-RS gamma phase against the nearest TRS gamma: median %.2f rad, MAD %.2f rad, |mean phasor| %.2f" % (np.median(rel), 1.4826 * np.median(np.abs(rel - np.median(rel))), abs(np.mean(np.exp(1j * rel)))))
