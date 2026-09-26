#!/usr/bin/env python3
"""Rebuild the recording_mode CSVs (msg4/ and DCIs/) from a run that printed to the terminal.

Runs made before recording_mode existed, or with it off, only left their terminal
log. This parses that log, plus msg4_bytes.log and the RACH CSV the run wrote
next to the binary, into files with the same name pattern and columns as
nrscope/src/libs/run_recorder.cc writes, so statistics scripts can be built and
tested on them.

Columns the terminal never showed stay empty: coreset_id, and, for DCIs, the
timestamp, SFN and slot unless the run's logs_nrscope_*.csv rows line up with
the printed DCIs one for one (checked field by field, all or nothing).

usage: replay_log_to_recordings.py TERMINAL_LOG [--build-dir DIR] [--out-root DIR] [--pci N]
"""
import argparse
import csv
import datetime
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

DCI_COLUMNS = ("timestamp,pci,sfn,slot,direction,rnti,rnti_type,dci_format,ss_type,coreset_id,aggregation_level,cce,"
               "freq_alloc,time_alloc,dci_mcs,dci_ndi,dci_rv,harq_id,tpc,ports,dmrs_id,srs_request,"
               "k,mapping,time_start,time_length,prbs,nof_prb,nof_layers,"
               "dmrs_type,dmrs_add_pos,dmrs_len,dmrs_typeA_pos,nof_dmrs_cdm_groups,n_scid,beta_dmrs,"
               "modulation,mcs,tbs,code_rate,rv,ndi,nof_re,nof_bits,mcs_table,xoverhead,ca_variant,dci_bits,dci").split(",")
MSG4_COLUMNS = ("timestamp,pci,sfn,slot,tc_rnti,c_rnti,rrc_transaction_id,rrc_offset,nof_bytes,dci,"
                "msg4_bytes,master_cell_group").split(",")

ANSI = re.compile(r"\x1b\[[0-9;]*m")


def kv(text):
    """key=value pairs of a DCI line or a grant dump line."""
    return dict(re.findall(r"([A-Za-z0-9_\-]+)=(\S+)", text))


def as_int(v):
    return str(int(v, 16)) if v.startswith("0x") else v


def parse_dcis(lines):
    """Each "DCIDecoder -- Found DCI" line, with the PDSCH/PUSCH dump that follows it."""
    out = []
    i = 0
    while i < len(lines):
        m = re.match(r"DCIDecoder -- Found DCI: (.*)$", lines[i])
        if not m:
            i += 1
            continue
        dci_str = m.group(1).strip()
        d = kv(dci_str)
        cfg = {}
        j = i + 1
        # Older builds printed stray debug lines (carrier->nof_prb, new_node_timestamp) before the dump
        while j < len(lines) and j <= i + 3 and re.match(r"(carrier->nof_prb|new_node_timestamp):", lines[j]):
            j += 1
        if j < len(lines) and re.match(r"DCIDecoder -- P[DU]SCH_cfg:", lines[j]):
            j += 1
            while j < len(lines) and lines[j].startswith(" "):
                cfg.update(kv(lines[j]))
                j += 1
        out.append((dci_str, d, cfg))
        i = j
    return out


def dci_row(dci_str, d, cfg, pci):
    fmt = d.get("dci", "")
    row = dict.fromkeys(DCI_COLUMNS, "")
    rnti_key = next((k for k in d if k.endswith("-rnti")), None)
    row.update({
        "pci": pci,
        "direction": "DL" if fmt.startswith("1") else "UL",
        "rnti": as_int(d[rnti_key]) if rnti_key else "",
        "rnti_type": rnti_key.split("-")[0] if rnti_key else "",
        "dci_format": fmt,
        "ss_type": d.get("ss", ""),
        "aggregation_level": str(1 << int(d["L"])) if "L" in d else "",
        "cce": d.get("cce", ""),
        "freq_alloc": as_int(d.get("f_alloc", "")) if "f_alloc" in d else "",
        "time_alloc": as_int(d.get("t_alloc", "")) if "t_alloc" in d else "",
        "dci_mcs": d.get("mcs", ""),
        "dci_ndi": d.get("ndi", ""),
        "dci_rv": d.get("rv", ""),
        "harq_id": d.get("harq_id", ""),
        "tpc": d.get("tpc", ""),
        "ports": d.get("ports", ""),
        "dmrs_id": d.get("dmrs_id", ""),
        "srs_request": d.get("srs_request", d.get("srs_req", "")),
        "dci": dci_str,
        # ca_variant and dci_bits stay empty: the printed text cannot tell them. Every
        # DCI is printed through the CA-configured decoder, so even DCIs matched with
        # the normal sizes show cc=0.
    })
    if cfg:
        s, l = cfg.get("t_alloc", ":").split(":")
        f0, fn = cfg.get("f_alloc", ":").split(":")
        row.update({
            "k": cfg.get("k", ""), "mapping": cfg.get("mapping", ""), "time_start": s, "time_length": l,
            "prbs": f"{f0}-{int(f0) + int(fn) - 1}" if f0 and fn and int(fn) > 0 else "",
            "nof_prb": fn, "nof_layers": cfg.get("nof_layers", ""),
            "dmrs_type": cfg.get("type", ""), "dmrs_add_pos": cfg.get("add_pos", ""), "dmrs_len": cfg.get("len", ""),
            "dmrs_typeA_pos": cfg.get("typeA_pos", ""), "nof_dmrs_cdm_groups": cfg.get("nof_dmrs_cdm_grps", ""),
            "n_scid": cfg.get("n_scid", ""), "beta_dmrs": cfg.get("beta_dmrs", ""),
            "modulation": cfg.get("mod", ""), "mcs": cfg.get("mcs", ""), "tbs": cfg.get("tbs", ""),
            "code_rate": cfg.get("R", ""), "rv": cfg.get("rv", ""), "ndi": cfg.get("ndi", ""),
            "nof_re": cfg.get("nof_re", ""), "nof_bits": cfg.get("nof_bits", ""),
            "mcs_table": cfg.get("mcs_table", ""), "xoverhead": cfg.get("xoverhead", ""),
        })
    return row


def fill_time_from_csv(rows, csv_path, t0, t1):
    """Take timestamp/SFN/slot from the run's own DCI CSV, if its rows match the printed DCIs one for one."""
    if not os.path.exists(csv_path):
        return "no DCI CSV"
    logged = [r for r in csv.DictReader(open(csv_path))
              if r.get("timestamp", "timestamp") != "timestamp" and t0 <= float(r["timestamp"]) <= t1]
    if len(logged) != len(rows):
        return f"not filled: {len(logged)} CSV rows vs {len(rows)} printed DCIs"
    for r, c in zip(rows, logged):
        if (r["rnti"], r["dci_format"], r["dci_mcs"]) != (c["rnti"], c["dci_format"], c["mcs_index"]):
            return "not filled: CSV rows do not line up with the printed DCIs"
    for r, c in zip(rows, logged):
        r["timestamp"], r["sfn"], r["slot"] = c["timestamp"], c["system_frame_index"], c["slot_index"]
    return f"filled from {os.path.basename(csv_path)}"


def parse_master_cell_groups(lines):
    """The JSON after each "masterCellGroup:", compacted onto one line as the recorder writes it."""
    groups = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("masterCellGroup: {"):
            body = ["{"]
            i += 1
            while i < len(lines) and lines[i] != "}":
                body.append(lines[i].strip())
                i += 1
            body.append("}")
            groups.append("".join(body))
        i += 1
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("terminal_log")
    ap.add_argument("--build-dir", default=os.path.join(ROOT, "build", "nrscope", "src"))
    ap.add_argument("--out-root", default=ROOT)
    ap.add_argument("--pci", default="632")
    ap.add_argument("--log-name", default="logs_nrscope_x410", help="log_name of the run's config, without .csv")
    a = ap.parse_args()

    lines = [ANSI.sub("", l.rstrip("\n")) for l in open(a.terminal_log, errors="replace")]
    start = next((re.match(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)", l) for l in lines if "Cell search: Setting" in l), None)
    if not start:
        sys.exit("no 'Cell search: Setting' line: cannot tell when the run started")
    t_start = datetime.datetime.fromisoformat(start.group(1))
    stamp = t_start.strftime("%Y-%m-%d_%H-%M-%S")
    t0 = t_start.timestamp()
    t1 = os.path.getmtime(a.terminal_log) + 5

    # DCIs
    printed = parse_dcis(lines)
    if printed:
        rows = [dci_row(s, d, c, a.pci) for s, d, c in printed]
        how = fill_time_from_csv(rows, os.path.join(a.build_dir, a.log_name + ".csv"), t0, t1)
        os.makedirs(os.path.join(a.out_root, "DCIs"), exist_ok=True)
        path = os.path.join(a.out_root, "DCIs", f"dci_{stamp}_pci{a.pci}.csv")
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, DCI_COLUMNS)
            w.writeheader()
            w.writerows(rows)
        print(f"{path}: {len(rows)} DCIs (timestamp/sfn/slot {how})")

    # RRCSetups: bytes from msg4_bytes.log, JSON from the terminal, SFN/slot from the RACH CSV
    groups = parse_master_cell_groups(lines)
    tids = [l.split(":")[1].strip() for l in lines if l.startswith("rrc-TransactionIdentifier:")]
    setups = []
    mb = os.path.join(a.build_dir, "msg4_bytes.log")
    if os.path.exists(mb):
        for l in open(mb):
            m = re.match(r"(\S+) slot=(\d+) rnti=0x([0-9a-f]+) outcome=(\S+) rrc_offset=(\d+) nof_bytes=(\d+) "
                         r"dci=\{(.*)\} bytes=([0-9a-f]+)", l)
            if m and "rrc_setup" in m.group(4) and t0 <= float(m.group(1)) <= t1:
                setups.append(m)
    rach = {}
    rp = os.path.join(a.build_dir, a.log_name + "_rach.csv")
    if os.path.exists(rp):
        for r in csv.reader(open(rp)):
            if r and r[0] != "timestamp" and t0 <= float(r[0]) <= t1:
                rach.setdefault(r[3], r)
    if setups:
        if len(groups) != len(setups):
            print(f"warning: {len(setups)} RRCSetups in msg4_bytes.log but {len(groups)} printed masterCellGroups")
        os.makedirs(os.path.join(a.out_root, "msg4"), exist_ok=True)
        path = os.path.join(a.out_root, "msg4", f"msg4_{stamp}_pci{a.pci}.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(MSG4_COLUMNS)
            for n, m in enumerate(setups):
                rnti = str(int(m.group(3), 16))
                r = rach.get(rnti, ["", "", "", ""])
                w.writerow([m.group(1), a.pci, r[1], r[2], rnti, rnti, tids[n] if n < len(tids) else "", m.group(5), m.group(6),
                            m.group(7).strip(), m.group(8), groups[n] if n < len(groups) else ""])
        print(f"{path}: {len(setups)} RRCSetups")


if __name__ == "__main__":
    main()
