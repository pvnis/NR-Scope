# GoldenSniffer sensing: offline processing of commercial-cell DM-RS, TRS and CSI-RS

Offline analysis of what NR-Scope + GoldenSniffer dump from a commercial cell (first used on Sunrise SA n78,
PCI 632, 273 PRB, 30 kHz, TDD DDDSU, 5 Nov 2026, 4-antenna X410 recording). Python 3 + numpy (+ scipy for the
tracker). Most scripts take `<recording dir> <out> [burst indices]` and read the replays in `<rec>/b<n>/<DIR><c>/`.

## Input: what NR-Scope writes (`golden_sniffer.dump_file`, one file per replayed antenna)

| Record | Content | Enabled by |
|---|---|---|
| `GSD1` | one grant's PDSCH DM-RS, port 1000, 3 values per PRB (CDM pairs) per DM-RS symbol, PRB mask | `golden_sniffer.dump_rnti` |
| `GSD2` | the same with port 1001 after port 1000 | env `NRSCOPE_GS_DUMP_P1001` |
| `GSD*` with rnti 0 | blind DM-RS: the PCI-scrambled DM-RS of the whole slot (N_ID = PCI, n_SCID = 0), mask = strict survey map | `golden_sniffer.dump_blind` |
| `GSC1` | raw REs of one symbol over the carrier: CSI-RS (symbol 13) and TRS (symbols 4/8) | dump_file; TRS with `dump_blind` |
| `GSR1` | raw grid of symbols 2, 3, 11, 12 for every DL DCI of one size | env `NRSCOPE_GS_RAW_SIZE` |

## Algorithms (library: `gs_all.py`, imported by the others with `exec` of its top part)

- `deslope`, `joint_gains`: **PRG stitching**. Per 4-PRB boundary, line extrapolation of each side to the boundary on the
  delay-deramped estimate, one complex number over the 4 antennas, chained across the band.
- `epochs`: **precoder epochs**. A snapshot joins the running epoch template if it fits it (delay slope search, coherence
  >= 0.9); per-epoch templates are stitched once (`joint_gains` on the template).
- `edge_number`: **edge matching** across a hole or between two pieces (two UEs side by side).
- `frame`: one complex number between two snapshots on shared PRGs: the **anchor** (merging UEs, gluing epochs).
- T1-T6 in `gs_all.py`: multi-UE tests (same-PRG similarity, bridge distance, edge vs anchor, coverage).

## Scripts

| Script | What |
|---|---|
| `gs_all.py`, `gs_all_report.py` | library + multi-UE tests on a 4-antenna recording; report tables |
| `gs_edge_ue.py` | edge matching between UEs |
| `diag_ph.py` | finds the per-subframe phase step (~60 deg per 1 ms, NR-Scope) |
| `gs_omp.py` | OMP / sparse delay model in frequency (holes, two pieces): no gain over edge matching |
| `gs_omp_t.py` | **epoch gluing** in time: anchor vs Doppler-tone model (simultaneous OMP over Doppler, phase scan), subframe step correction |
| `gs_inc.py` | coherent merge vs incoherent sum of range-Doppler maps (controlled + real UEs) |
| `gs_usemap.py`, `gs_usemap_v6e.py`, `gs_usemap_v6f.py` | **usage map**: which DM-RS each step of the pipeline uses, per slot and PRG (v6e: + blind DM-RS, v6f: + TRS and the TRS anchor) |
| `gs_trs_lib.py` | Gold sequence, CSI-RS sequence init, TRS and 8-port CSI-RS despreading, `GSC1` reader |
| `gs_trs.py`, `gs_trs_sum.py` | TRS stability, TRS as anchor (one number + relative delay slope), light-network band, CSI-RS port span |
| `gs_csi.py` | every DM-RS kind (tablet, other UEs, blind) and the TRS in the span of the 8 CSI-RS ports (4 antennas) |
| `sim_sense.py`, `sim_sum.py` | simulation with known targets on the real occupancy, per-map scores |
| `sim_track.py`, `sim_track_sum.py` | the same as a continuous geometric scene, tracked with `~/GIOVANNI/sim/track.py` (env `TRACKDIR`), GOSPA / false tracks |
| `pages/` | usage-map page template + builder (`mk_usemap_page.py <dir>`), the math page, the 2-beam text |

## Recording and replay (sens6 paths are hard-coded)

| Script | What |
|---|---|
| `rec4c.sh` | 4-channel IQ bursts into RAM (2 GB/s does not fit any disk live), moved to disk; RAM wait, watchdog, `usrp-hwd` restart on the X410 |
| `replay.sh`, `replay4ch.sh`, `replay4ch_blind.sh` | one channel through NR-Scope (file radio, lossless); `_blind` adds `dump_blind` |
| `run_all2.sh`, `run_blind_all.sh`, `run_trs_all.sh`, `run_ch12.sh` | replay every burst x antenna, then the analyses |

Results and open items: `sens6:~/field/OPEN_ISSUES.md`. Heavy data: `sens6:/mnt/mydrive/field`.

## NR-Scope's sensing chain on our estimates (2026-11-09)

| Script | What |
|---|---|
| `gs_to_nrsr.py`, `gs_to_nrsr_ep.py` | dumps -> `sensing.record_estimates` for `build/nrscope/src/tests/sensing_offline` (variants raw / proc / csi / trs; env WSH=1 undoes the logged window moves, LOS0 puts each antenna's direct path at range zero) |
| `so_csi.py` | CSI-RS / TRS reference corrections, window-shift (GSW1) and drift helpers |
| `so_light.py` | light network: TRS + every DM-RS except the tablet's, one stream |
| `so_simprec.py`, `so_synth.py` | synthetic estimates on the real occupancy (simulated precoders, timing jumps) |
| `so_score.py` | scores `map2d.csv` against the injected targets (`--inject`) |
| `csi_beamdiv.py` | how much the precoder direction changes per PRG (from the CSI-RS) |

Needs the NR-Scope options `sensing.allow_gaps`, `sensing.align_search_bins`, `sensing.min_snapshots` and the
GoldenSniffer GSW1 record. Per-antenna replays must be aligned (window moves, direct path) before anything
multi-antenna: each replay has its own sync loop.

## Reference, calibration and merging (2026-11-10/11)

All streams: `GCAL=1` calibrates the inter-antenna phase per 4-PRB group from the TRS (one phase per snapshot is not
enough: it swings by turns across the band). Run `sensing_offline` with `sensing.align_search_bins: 0`
(`cfg_*_noalign.yaml`): the converters already undo the window moves and zero the direct path, and the delay search
can lock on a secondary path (every target a constant few bins off).

TRS reference (`so_csi.ref_fix`, `ref_fix_hyb`), per stream: `REF_EP` (`gs_to_nrsr_ep.py`), `TAB_EP` (decoded UE)
and `BLIND_EP` (blind record) in `so_light.py`:

| Mode | Reference numbers | Notes |
|---|---|---|
| `snap` | one per snapshot and 4-PRB group | default; eats the array-mean part of near slow targets (+15 deg read as +22 at 25 m) |
| `epoch` | one per precoder epoch and group | angles right, many detections lost |
| `slot` | one per slot and group | no measurable change vs `snap` |
| `hyb` | `c(t) * D(epoch, group)`: one per snapshot common to all groups x one per epoch and group | with `HYB_MAXEP=0.01` (epochs <= 10 ms): angles right, detection kept |

`BLIND_TRACK=1` splits the blind record (several UEs per symbol) into UE tracks, one epoch key per group, so
`BLIND_EP=hyb` keeps its detections. Proposed: tablet `TAB_EP=hyb BLIND_EP=hyb BLIND_TRACK=1 HYB_MAXEP=0.01`,
no tablet `BLIND_EP=hyb BLIND_TRACK=1 HYB_MAXEP=0.01`.

| Script | What |
|---|---|
| `run_v3.sh`, `run_hyb4.sh`, `run_hyb4_all.sh` | TRS-ref stream: reference variants, delay search on/off, all bursts |
| `run_light4*.sh`, `run_light5_all.sh`, `run_light6_all.sh` | light streams: reference variants, UE tracks, slot pooling |
| `run_minprb.sh`, `run_usemap.sh`, `gen_um7.py` | narrow DM-RS symbols; usage map (env `NR_SENSING_TRACE`) |
| `run_kpi.sh`, `kpi_pfa.py` | no injection: detections per map/cell (Pfa upper bound), clutter suppression |
| `sum_v3.py`, `aoa_v3.py`, `cellpow.py` | main table, 3-target angles, target cell power |
| `trs_consist.py`, `wsh_resid.py` | TRS consistency between occasions; residual slope vs logged window moves |
| `track_real.py` | tracking on real-data maps (stand-in geometry) |
