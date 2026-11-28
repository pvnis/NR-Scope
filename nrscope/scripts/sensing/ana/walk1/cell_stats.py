"""Statistics of a commercial cell from a decoded recording (chunks dir: b<k>/replay_blind_ch0, decoder with the joint
DM-RS test), for the paper's "cell at a glance" table. Channel 0 only. Chunks are decoded independently (2.45 s each),
so per-UE numbers are lower bounds (a UE is only seen once its PDCCH is found in that chunk).
usage: cell_stats.py <chunks dir> <label> <out.json>"""
import sys, os, glob, re, json, struct, csv, collections, numpy as np
CH, LAB, OUT = sys.argv[1:4]
NPRB = 273
def records(path):
    raw = open(path, "rb").read(); i = 0; R = []; C = []; P = {}
    while i + 12 <= len(raw):
        t = raw[i:i + 4]
        if t in (b"GSD1", b"GSD2", b"GSD3"):
            if i + 20 > len(raw): break
            sfn, sl, ns, rn, npb, o, c = struct.unpack_from("<HBBHHff", raw, i + 4); per = {b"GSD1": 3, b"GSD2": 6, b"GSD3": 4}[t] * npb * 8
            if i + 20 + npb + ns * (1 + per) > len(raw): break
            m = np.frombuffer(raw, np.uint8, npb, i + 20).astype(bool); j = i + 20 + npb; ls = []
            for _ in range(ns): ls.append(raw[j]); j += 1 + per
            R.append((sfn * 20 + sl, t == b"GSD3", rn, int(m.sum()), ls)); i = j
        elif t == b"GSC1":
            sfn, sl, l, fc, npb = struct.unpack_from("<HBBHH", raw, i + 4); C.append((sfn * 20 + sl, l, npb)); i += 12 + 12 * npb * 8
        elif t == b"GSE1":
            sfn, sl, _, npb = struct.unpack_from("<HBBH", raw, i + 4); P[sfn * 20 + sl] = np.frombuffer(raw, np.float32, npb, i + 10).copy(); i += 10 + 4 * npb
        else: i += 1
    return R, C, P
S = collections.defaultdict(list); rnti_all = collections.Counter(); per100 = []; sizes = collections.Counter(); fdra = collections.Counter()
dl_slots = 0; dm_slots = 0; prb_slot = []; full_slots = 0; t1_prb = t2_prb = 0; extra_sym = 0; nrec = 0; tab_prb = collections.Counter()
re_dmrs = 0.0; trs_sym = 0; csi_sym = 0; secs = 0.0; rate = []; dcis = 0
trs_slots = collections.Counter(); csi_slots = collections.Counter()
for d in sorted(glob.glob(os.path.join(CH, "b*")), key=lambda p: int(re.sub(r"\D", "", os.path.basename(p)))):
    o = os.path.join(d, "replay_blind_ch0")
    if not os.path.exists(os.path.join(o, "ue_dmrs.bin")): continue
    R, C, P = records(os.path.join(o, "ue_dmrs.bin"))
    if not R: continue
    aa = sorted(set(r[0] for r in R)); wr = aa[-1] - aa[0] > 10240   # a 2.45 s chunk spanning > half the SFN cycle wrapped
    uw = [x + 20480 if (wr and x < 10240) else x for x in aa]
    a0 = aa[0]; span = (max(uw) - min(uw) + 1) * 5e-4; secs += span
    # DCIs and UEs
    try:
        rows = list(csv.DictReader(open(os.path.join(o, "golden_sniffer.csv"))))
    except Exception: rows = []
    dcis += len(rows)
    win = collections.defaultdict(set)
    for r_ in rows:
        rn = r_.get("rnti", ""); rnti_all[rn] += 1
        if r_.get("size"): sizes[r_["size"]] += 1
        try: a = int(r_["sfn"]) * 20 + int(r_["slot"]); win[(a - a0) % 20480 // 200].add(rn)
        except Exception: pass
    per100 += [len(v) for v in win.values()]
    for line in open(os.path.join(o, "nrscope.log"), errors="ignore"):
        m = re.search(r"FDRA: (type \d|dynamic)[^;]*rbg config (\d)", line)
        if m: fdra["%s, rbg config %s" % (m.group(1), m.group(2))] += 1
        m = re.search(r"\[status\].* (\d+) slots/s", line)
        if m: rate.append(int(m.group(1)))
    # occupancy and DM-RS
    by = collections.defaultdict(list)
    for a, t2, rn, n, ls in R: by[a].append((t2, rn, n, ls))
    DLc = collections.Counter(a % 20 for a in P); DMc = collections.Counter(a % 20 for a in by)
    dl = {c for c in DLc if DMc[c] >= 0.05 * DLc[c]}
    for a, pw in P.items():
        if a % 20 not in dl: continue
        dl_slots += 1
    for a, L in by.items():
        if a % 20 not in dl: continue
        dm_slots += 1; tot = 0
        for t2, rn, n, ls in L:
            nrec += 1; tot += n
            if t2: t2_prb += n
            else: t1_prb += n
            extra_sym += len(ls) > 1; tab_prb[rn] += n
            re_dmrs += n * (4 if t2 else 6) * len(ls)          # REs of the CDM group carrying port 1000, per DM-RS symbol
        prb_slot.append(min(tot, NPRB)); full_slots += tot >= 0.9 * NPRB
    for a, l, npb in C:
        if l in (4, 8): trs_sym += 1; trs_slots[a % 40] += 1
        elif l == 13: csi_sym += 1; csi_slots[a % 20] += 1
del_tab = tab_prb.most_common(3)
tot_prb = t1_prb + t2_prb
out = {
    "label": LAB, "seconds": round(secs, 1), "dci": dcis, "dci_per_s": round(dcis / max(secs, 1e-9)),
    "rntis_seen": len([r for r in rnti_all if r]), "rntis_per_100ms_median": float(np.median(per100)) if per100 else None,
    "rntis_per_100ms_p90": float(np.percentile(per100, 90)) if per100 else None,
    "dl_slots_with_dmrs_pct": round(100 * dm_slots / max(dl_slots, 1), 1),
    "dmrs_prb_per_dl_slot_median": float(np.median(prb_slot)) if prb_slot else None, "dmrs_prb_per_dl_slot_mean": round(float(np.mean(prb_slot)), 1) if prb_slot else None,
    "full_band_slots_pct": round(100 * full_slots / max(dm_slots, 1), 1),
    "dmrs_type2_prb_share_pct": round(100 * t2_prb / max(tot_prb, 1), 1),
    "records_with_extra_dmrs_symbols_pct": round(100 * extra_sym / max(nrec, 1), 1),
    "top_rnti_prb_share_pct": [(r if r else "blind", round(100 * n / max(tot_prb, 1), 1)) for r, n in del_tab],
    "dci_sizes": dict(sizes.most_common(8)), "fdra_learned": dict(fdra.most_common()),
    "trs_symbols_per_s": round(trs_sym / max(secs, 1e-9), 1), "trs_slot_positions_mod40": {k: v for k, v in sorted(trs_slots.items()) if v >= 0.1 * max(trs_slots.values())},
    "csirs_symbols_per_s": round(csi_sym / max(secs, 1e-9), 1), "csirs_slot_positions_mod20": {k: v for k, v in sorted(csi_slots.items()) if v >= 0.1 * max(csi_slots.values())},
    "pilot_res_per_s_trs": round(trs_sym / max(secs, 1e-9) * NPRB * 3), "pilot_res_per_s_dmrs_port1000": round(re_dmrs / max(secs, 1e-9)),
    "decode_slots_per_s_median": float(np.median(rate)) if rate else None,
}
out["dmrs_over_trs_pilot_ratio"] = round(out["pilot_res_per_s_dmrs_port1000"] / max(out["pilot_res_per_s_trs"], 1), 1)
json.dump(out, open(OUT, "w"), indent=1); print(json.dumps(out, indent=1))
