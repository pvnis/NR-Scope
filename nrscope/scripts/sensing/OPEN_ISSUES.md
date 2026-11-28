# PCI 632 passive sensing — open issues

Cell: Sunrise SA n78 PCI 632, PDSCH DM-RS (type 1, ports 1000+1001, symbols 2+11), PRG = 4 PRBs.
Data: ~/field/20261104_164101_sunrise_632_4ch_bursts_tabletdl_iq (6 bursts x 2.45 s, 4 antennas; DELETED 2026-11-05, results kept in ~/field/ana).
Tablet RNTI 0x7466, second UE 0x7176 (wideband precoder, burst 0 only).

How to use this file: when an item is solved, change its status, write the result and the
script/commit next to it, and add a line to the changelog at the bottom. Never delete items.

Status: OPEN / IN PROGRESS / DONE / DROPPED

---

## Results of 2026-11-05 (recording ~/field/20261105_135044_*, 26 x 2.45 s, + test 20261105_132906_*)
Report: https://claude.ai/artifact/CBCpsHkpSuTW9DcpfeFmdE (local ~/Downloads/pci632_multiue_report_5nov.html).
Raw tables: sens6 ~/field/ana/*.md; scripts ~/field/gs_all.py (all tests, per burst), gs_all_report.py, rec4c.sh (guarded recording).
Note: today antennas 1-2 decode only the tablet's PDCCH; other UEs analysed with CH=0,3 (both days, for comparison).
- 1b: stitching on one-precoder UEs (true jumps 0): bias 1.6 deg / 0.23 dB per edge, ghost unchanged (-32.2 -> -32.3 dB);
  tablet ghost -13.6 -> -30.2 dB (per-snapshot profiles). Stitching is not the limit; old -9..-11 dB used a coherent
  time-averaged metric that includes the channel's own variation.
- 1c: amplitude jumps are real: tablet p90 7.6 dB vs bias p90 0.65 dB; corrected by the complex stitching gain.
- 2c: holes, one-precoder UE: 0/4/8/16/32 PRBs -> 2.3/3.7/5.8/10.8/22 deg. Per-4-PRB UE: any hole fails (hidden jumps).
- 3a: 25 second UEs, 2121 side-by-side DM-RS symbols.
- 3b: per-4-PRB UEs: one number per PRG 49% at 2 ms (ideal 22%); self-stitch + one number 80%. Next: average per epoch.
- 3c: merge floor = how far the UE is from tablet x number (ideal 12-20%), not stitching.
- 3e: edge 47% vs anchor 30% vs ideal 20% (touching, self-stitched); edge = fallback when no anchor.
- NEW gate: keep merge if anchor fit residual < 30%: keeps 52%, error 14% (ideal 12%). No ground truth needed.
- Bridge: today ~10 ms (wideband UEs: 27/41/71% at 2/10/30 ms; tablet own change 4/19/64%); 4 Nov ~30 ms. Lower SNR today
  (two-symbol noise 4-8% vs 0.7-1.5%) and faster channel change.
- 3g: 1 s anomaly not reproduced (100% at 1 s today).
- Coverage: tablet alone 240 PRBs/slot median, 86% slots >=160 PRBs; + all mergeable UEs +1%. Merging adds little band
  while a heavy downloader is present -> 4a deferred.

## Usage map 2026-11-05 (gs_usemap.py, https://claude.ai/artifact/Ma5iGHZvtcguvYiKrCBT8Y)
Per 150 ms window, reference = UE with most PRG-snapshots; same processing as gs_all (bridge 10 ms, gate 30%).
- With the tablet: 85% of on-air DM-RS PRB-slots used; median window band 62% (p10 25%); 8 precoder epochs per window
  (longest 47 ms median) -> Doppler coherence breaks ~every 19 ms (4c still open).
- Without the tablet: 5.5% used; 68% of windows have nothing placed (light network = we are blind today).
- 14% (with tablet) of on-air PCI DM-RS belongs to UEs whose grants we do not place (survey labels them by RNTI, e.g. 0x79b9):
  median 4.9% of the band per window, 42% of windows above 7.1%. -> new item 6e.
- Periodic: TRS (symbols 4+8) in 2 slots every 20 ms (100 slots/s); CSI-RS on symbol 13 in ~283 slots/s (several resources).
- Compared with LaSen (SenSys'26): they use their own monostatic SDR at 5.8 GHz with occupancy masks from real gNBs
  (no precoding, all REs known, no UE placement); range+velocity tracking only, no angle.

## Sparse models / OMP 2026-11-06 (gs_omp.py, gs_omp_t.py; math page https://claude.ai/artifact/Nj3rYQYg1txP8sb958haML)
- NEW: per-subframe phase step: +57..62 deg (std 0.7) at every 1 ms boundary, ~2 deg inside a subframe (~168 Hz =
  0.04 ppm gNB/X410 offset; NR-Scope likely restarts its CFO phase each subframe). Remove y*exp(-j*phi*floor(slot/2))
  before any time processing; fix in NR-Scope (new item 4d).
- Frequency (127 tablet snapshots, known numbers): OMP with one number per piece never beats edge matching
  (no hole 18-29 deg vs 1.2; 8-PRB hole 29-44 vs 11); separated grants fail beyond ~8 PRBs apart either way.
  Even with the true number the sparse delay model fills 4-8 PRB holes (10%), 16 PRBs 30-40%, fails >= 32 PRBs.
  Channel not sparse enough in delay -> 2b, 3f DROPPED (frequency).
- Time (56 controlled split epochs, 286 real epoch changes): sparse Doppler model (4 tones shared by all series, phase scan)
  3.4/5.1/11 deg at 0/5/20 ms gap vs anchor 3.8/19/66; real epoch changes are mostly back-to-back (gap 0.7 ms), anchor and
  Doppler model agree to 10 deg. -> 4c IN PROGRESS: glue epochs with anchor (short gaps) / Doppler model (long gaps).

## 1. One UE across frequency (precoder jump every 4 PRBs)

Works: the jump is one complex number common to all 4 antennas; stitching removes it.
4 antennas: per-edge precision 0.5 deg vs 1.1 deg (1 antenna); far end ~3x more repeatable.
(scripts: gs_edge4.py)

- [OPEN] 1a. Error accumulates along the chain: far end of PRBs 32-191 median 7 deg, p90 36 deg.
- [DONE 2026-11-05] 1b. Better precision did not improve the range profile: PRG ghost (+208 m) raw -4.6 dB,
  1-antenna epoch -11.2 dB, 4-antenna epoch -9.2 dB, 4-antenna snapped to +-90 deg -4.9 dB. Unexplained.
- [DONE 2026-11-05] 1c. Amplitude part of the jump never studied separately (only clipped to 0.1-10).
  Does the precoder change amplitude per PRG? Should we correct it?

## 2. Gaps (holes in frequency inside/between grants)

Works: filling holes from the same PRBs in other slots of the same epoch (gs_ue_sensing_holes.py).

- [OPEN] 2a. Stitching cannot cross a hole in the same slot: phase between the two sides unknown.
- [DROPPED 2026-11-06] 2b. Not tried: model fit (same path delays on both sides, one unknown number per piece).
- [DONE 2026-11-05] 2c. Edge matching across holes (gs_edge_ue.py, controlled test on 0x7176's own grant, 4 antennas,
  fit 2-4 PRBs per side): gap 0 -> 2.5 deg median, 4 PRBs -> 3-4 deg, 8 PRBs -> 4.5 deg (p90 15), 16 PRBs -> 9-10 deg
  (p90 27), 32 PRBs -> 19-23 deg (fails). 1 antenna is ~2x worse. Usable up to ~8 PRBs; not yet plugged into sensing.

## 3. Merging different UEs in frequency (MAIN GOAL)

Works (small sample): second UE = tablet x one complex number per PRG (8% residual, noise 1.7%).
Cross-slot anchoring: number measured from tablet slots on the same PRBs.
Merged band error 11% vs 63% raw join (8 side-by-side slots). Bridge distance test
(gs_merge_window.py, 70 cases): same quality up to ~30 ms (8 deg, 1 dB, ~27% channel error),
breaks at ~100 ms (31 deg, 3.3 dB, 82%). (scripts: gs_merge4.py, gs_merge_window.py)

- [DONE 2026-11-05] 3a. Too few cases: need a recording with a second phone we control downloading next to the tablet.
- [IN PROGRESS] 3b. Only one second UE tested, with a wideband precoder. A second UE with per-PRG precoder
  (one number per 4 PRBs) untested.
- [DONE 2026-11-05] 3c. Error floor ~25% on the second UE's PRBs, mostly from the tablet's own stitching (see 1a/1b).
- [OPEN] 3d. Not tried: chaining through a third UE (A <-> C <-> B).
- [DONE 2026-11-05] 3e. Edge matching between UEs (gs_edge_ue.py, 2026-11-05). 14 real slots with tablet + 0x7176
  side by side (28 DM-RS symbols), 13 touching (gap 0), 1 with a 16-PRB gap.
  Touching: agrees with cross-slot anchoring to 6 deg / 1 dB (median); error on the UE's PRBs vs full-band tablet
  reference: edge 18%, anchor 25%, best possible single number 17% (only 5 slots have a reference).
  Works in 5 slots where anchoring had no tablet snapshot. Fails when one side is thin (UE run 16 PRBs: -28 deg)
  and in 1 unexplained slot (6.381 s: 27-31 deg); the 16-PRB-gap slot fails (200-276%).
  Open: too few real cases (needs the second-phone recording); combine edge + anchor (weighted) when both exist.
- [DROPPED 2026-11-06] 3f. Not tried: model fit (common delays, one number per grant).
- [DONE 2026-11-05] 3g. 1 s row of the bridge test is good again (8.5 deg) although the tablet changed 90%. Unexplained.

## 4. Merging over time (range/Doppler maps)

- [DEFERRED] 4a. UE merge not folded into the 150 ms range/Doppler maps.
- [OPEN] 4b. 2-beam bilinear model (good for time consistency/clutter: static residue -34 vs -19.5 dB,
  bad for range because of per-PRG gauge freedom) not integrated in the pipeline.
- [OPEN, FOUND BY CODE READING 2026-11-14] 4b3. Clutter kernel ignores gaps. nr_ue_map.c clutter_kernel models a
  static path in snapshot m as K_m(u) = sum over l = 0..n-1 of w[l] e^{j2pi(a_m+l)u/N}: Hann over the whole span from
  first to last pilot, every position assumed measured. The snapshot keeps only n_pilots (span incl. gaps) and
  k_first, no mask. Real-data configs run allow_gaps: true + clutter_kernel: true (so/cfg_light.yaml etc.), where gaps
  are zeros in the data but present in the model -> residue alpha * sum_{gaps} w e^{...} after the fit, per snapshot
  different (gap pattern changes) -> modulated in slow time, Doppler stripes / clutter left, scales with the gap
  fraction. Still valid: one alpha per window (needs a common reference = TRS anchoring), Hann, the LS fit itself.
  Fix: keep the pilot mask per snapshot, K_m(u) = sum over measured l of w[l] e^{...} (kernel table no longer dedups
  by (n, a_m): per snapshot, or FFT of w*mask). Test: sim on the real occupancy masks, clutter residue mask-aware vs
  today; then real data cfg_light.

  2026-11-16 FIX IMPLEMENTED (NR-Scope sensing, uncommitted): snapshot keeps its measured mask (bitset), taper
  coherent gain and mask hash; clutter kernel, direct-path matched filter, gain and the TDD detector range window
  (mean of gain x taper x mask over the map) use it; kernels shared by mask hash, 64-kernel cap removed (would skip
  the clutter removal). Build clean, sensing_selftest PASSED (no gaps in it). Real data, run_light4 eval bursts 0-5,
  na (OLD light6all -> NEW): tablet 87/82/86/82 -> 96/89/97/92 % (slow-10/slow-20/fast-10/fast-20), other/map same;
  no tablet 86/81/60/48 -> 95/75/65/45, other/map 6.0/2.1/1.7/0.7 -> 3.4/1.1/1.3/0.4.
  Gain ROC (run_roc.sh, roc_sum.py, Pfa 1e-9/1e-6/1e-3; g0 = no equalisation, g1 = full, g1u = full for clutter then
  undone before Doppler): at matched other/map g0 >= g1u > g1; full equalisation only moves to a stricter point and
  loses weak targets (no tablet -20 dB ~62% vs 70% for g0 at ~2.1 other/map; tablet 89 vs 92 at ~3.3). Total found
  equal at equal false rate. -> proposal: no gain equalisation (g0). Temporary env switches NR_SENSING_GAIN_EXP,
  NR_SENSING_GAIN_UNDO, NR_TDD_PFA in the code (to remove). Open: gap unit test in sensing_selftest, malloc-fail leak
  in remove_path, Python map builders (sim_sense/sim_track/gs_inc, sens2 precode_maps) still gap-unaware, snapshot
  dump lacks the mask.

  2026-11-16 cleanup: experiment env switches removed; gain equalisation is now yaml sensing.gain_equalise (default
  true = as before); new yaml sensing.min_measured (drop a symbol whose MEASURED pilots cover less than this fraction
  of the carrier; default 0 = span test only); malloc-fail leak fixed. Build clean, selftest PASSED.
  Fill sweep (run_fill.sh, roc_sum.py /mnt/mydrive/field/fill, bursts 0-5, Pfa 1e-6), found % (-20 dB) @ other/map:
    no tablet  gain off: m0 74.5 (75.5) @3.09 | m.25 72.3 @3.08 | m.5 66.6 @3.37 | m.75 53.5 @3.49
               gain on : m0 69.4 (59.6) @1.54 | m.25 72.2 (67.3) @2.02 | m.5 69.0 @2.80 | m.75 55.9 @3.34
    tablet     gain off: 94.8 / 94.6 / 93.8 / 93.0 ; gain on: 92.7 / 94.9 / 94.4 / 93.2
  -> dropping symbols does not help: with gain off it only loses targets; with gain on, a threshold just removes
  the sparse snapshots the gain over-amplifies (same effect as gain off). Keep every symbol (min_measured 0).
  Next (LaSen lessons): per-map quality score (masked-atom coherence + occupancy) as output and tracker input;
  mask-aware target sidelobe removal (per-snapshot masked PSF or 2D-OMP) instead of periodogram + CFAR.

- [DONE prototype 2026-11-16] 4b4. Mask-aware 2D-OMP + LaSen ID score (omp_eval3.py, run_omp3.sh, omp_roc.py,
  /mnt/mydrive/field/omp3; bursts 0-5, tablet / no tablet, slow/fast, -10/-20 dB). Same windows (150 ms), same
  preprocessing, same threshold rule for the periodogram (PER, Hann range/time, local maxima) and simultaneous 2D-OMP
  (atoms on measured REs only, off-grid refinement, joint LS, noise fixed from the first map). Self-checks: map =
  explicit masked-atom inner products, residual non-increasing (pass).
  Found on the way: (a) per-snapshot direct-path jumps of 3-5 bins left upstream (b0 window 9: 84/372 snapshots)
  smear the static scene -> per-snapshot direct-path alignment added to the shared preprocessing (likely also hurts
  the C pipeline in na configs: check). (b) OMP noise must be fixed once (a residual median falls as atoms go).
  (c) The injection (sensing_offline --inject and the copy here) adds a copy of the WHOLE channel, so each target
  brings copies of the static scene at its Doppler: OMP resolves them (92% of its "false" with the tablet), PER blurs
  them. Scored without detections at a target Doppler, at equal other/map:
    no tablet  0.5/map PER 33.8% OMP 53.7% (-20 dB 10 -> 23) | 1/map 43.6 -> 64.8 (15 -> 38) | 2/map 62.1 -> 74.4
    tablet     0.5/map 90.0 -> 91.8 | 1/map 93.5 -> 94.7 | 2/map 95.2 -> 96.6
  -> OMP helps where data are sparse (no tablet +12..+21 points), little with dense data. ID score: no predictive
  value (corr with found -0.06..0.16; mu only 0.25-0.30 over windows): failures come from residue/jumps/SNR, not
  from the sampling geometry of a 150 ms window. Not adopted. Next: point-echo injection (direct-path gain x one
  path, no static copies) to score cleanly; port OMP to C if kept.

  2026-11-16 clean check (omp_eval4.py, /mnt/mydrive/field/omp4): targets injected as ONE path each (the direct
  path's own per-snapshot, per-antenna gain and delay moved to the target; no static copies), strict scoring (every
  non-target detection is false). Found at equal other/map 0.5 / 1 / 2 / 3:
    no tablet PER 27.1 / 35.8 / 50.3 / 59.1 %  OMP 49.3 / 56.0 / 64.7 / 69.9 %  (-20 dB at 1/map 8.6 -> 21.3)
    tablet    PER 76.7 / 82.2 / 87.1 / 90.4 %  OMP 81.8 / 87.0 / 91.6 / 93.4 %
  -> gain confirmed: C port as sensing.backend option (offline first), FFT correlation over bins (bins to grow).

  2026-11-16 C PORT (NR-Scope, uncommitted): nr_ue_omp.c/.h, yaml sensing.backend: periodogram|omp, omp_max_paths,
  omp_threshold_db; slow-time blocks carry snap[]/gain[]; OMP fills the TDD detector's target list (markers, AoA,
  trackers unchanged); noise-weighted (data / gain, so gain_equalise does not rule the fit); exact sub-bin refinement
  (iterated parabola on the exact correlation in delay, golden section in Doppler); 50 dB dynamic-range floor;
  exclusion of one cell around taken paths. Test sensing_omp_test (data through the real delay transform): 4 cases
  PASS (full band, gaps, gaps + mean removed, gaps + gains). sensing_offline --inject ...,point=1: one path following
  the direct path (gain + delay from the delay-profile peak; a phase-slope estimate was off by up to 89 m - bug found
  and fixed, the first NR-Scope OMP comparison /mnt/mydrive/field/ompc is INVALID).
  Findings on b0 (single point target, -10 dB, 20/50/100 m): both back ends see the target in 11-15/16 maps, but the
  MIRROR CHECK (mirror_reject: true in our configs, verdict 3) sets most of them aside (found 2-5/16); with it off
  per 10-14/16, omp 12-15/16. The TDD detector also marks some true targets as replicas (verdict 1). -> our earlier C
  evaluations lost targets to mirror_reject too (cause of the strong -v mirrors not understood yet).
  OMP extra detections (4.7-7.4/map at 11.4 dB): 39% at the target range at another speed, 14% at the target speed at
  another range, 20% near the direct path: the target's slow-time / range signature is not one clean atom on real data
  (per-snapshot gain/phase/delay errors) and the direct-path residue. Next: cluster detections like the tracker does
  before counting, ROC at equal false rate with mirror_reject off, then make OMP robust to per-snapshot errors.

- [IN PROGRESS 2026-11-16] 4b5. Precoded simulator (sens2 worktree ~/GIOVANNI_prec, branch precoding-sim,
  tools/precoding_sim/: prec_channel.py PortChannel + Type I codebook, gen_dumps.py, checks; masks from 5 Nov via
  ~/field/masks_export.py -> /mnt/mydrive/field/precsim_masks). Sionna static scene per gNB port once (3GPP 38.901
  element, N1xN2x2 panel, TS 38.214 port order), moving paths per 5 ms trace; real 5 Nov scheduling; per UE/PRG/epoch
  Type I PMI on the UE's own channel (LoS + 4 rays); TRS one element; writes NR-Scope GSD1/GSC1 dumps (real pipeline
  runs unchanged) + truth. New site epfl_bc_wide (OSM radius 300 m, hand edits copied). Geometry: gNB SV roof edge
  (16.3, 105.2, 20.5 m); the old receivers were blocked or blind to the agents; rx B (-20.2, 9.0) has LoS (direct
  path +19.6 dB over the next) and sees 2.1/2.3 agents. Found: rx1 LoS -45 dB (strongest path a reflection +59.6 m)
  -> truth now from the strongest static path (the pipeline zero).
  First result (burst b9 masks, 64T 8x4x2, PRG 4, epochs ~19 ms, rx B, light pipeline MODE=tab, mirror off):
    no precoding: car (-4.7 m/s, 14-40 m) per 14/16 omp 15/16   |   Type I precoding: per 0/16 omp 1/16
  -> precoding destroys detection here; NOT yet trusted: validate the simulator against 5 Nov statistics first
  (anchor residuals tablet 10% / others 61%, epoch durations, TRS stability), score with all target paths (car
  scatterers, wall echoes), drop static targets from truth.

  2026-11-16 GUI: second instance from the worktree, port 8051 (ssh -L 8051:127.0.0.1:8051 sens@sens2), logs in
  ~/GIOVANNI_prec/logs; runs sim_{none,t1}_b9_{per,omp} published (publish_run.py: map_clock utc + time.csv, Sionna's
  own ground_truth*.csv), scenario precoded_bc_wide. Orchestration: scratchpad sim_burst.sh (gen on sens2, so_light +
  sensing_offline on sens8, publish on sens2).
  VALIDATION, same script as on real data (gs_trs.py on the simulated burst, layout rec_<tag>/b0):
                         real 5 Nov (6f)        sim Type I (b9 masks)     sim no precoding
    TRS 1-|rho|^2 20 ms  0.03 (floor .02-.05)   0.006                     0.006
    anchor pass, tablet  10%                    91%                       100%
    anchor pass, others  61%                    32%                       100%
    anchor pass, blind   32%                    23%                       100%
  -> the simulator does NOT match yet: too clean (no CFO / phase noise / timing jumps, SNR), the tablet fits the TRS
  far too well (static scene dominated by the direct path, +19.6 dB over the next path; real multipath richer?),
  other UEs fit too badly (random positions in the sector put rx in their sidelobes). Tune before any result.

  2026-11-16 calibration (trs_profile.py, real_trs_echoes.json, gen_dumps --clutter/--clutter-key, PortChannel
  .add_scatterers: measured echoes as point scatterers on their bistatic ellipse, exact; level set through the TRS
  beam). CORRECTION: averaging the real TRS profile WITHOUT undoing NR-Scope's timing jumps made a fake -1..-5 dB echo
  at 13-19 m; with each symbol aligned on its direct path the real static scene is (b9 ~ b2): 9 m -19.4, 20 m -20.5,
  71 m -20.8, 79 m -22.9, 54 m -28.6, 42 m -30.9 dB, floor -42.3 dB, delay spread 35 ns (Sionna alone: -32..-45 dB,
  7 ns). Calibrated sim (b9 echoes, SNR 17 dB): echoes within ~1 dB, floor -44.3, spread 30 ns, TRS 1-|rho|^2 0.039
  (real 4-symbol 0.025). Precoded statistics (NOT fitted), anchor pass sim / real: tablet 87% / 10%, others 23% / 61%,
  blind 16% / 32% -> still wrong. Tablet-specific? (heavy UE: Type II multi-beam CSI, or rank 2, or its real position);
  others fit too badly (random sector positions / 64T beams too narrow?). Need: tablet position on 5 Nov, gNB array.

  2026-11-16 REAL GEOMETRY (user: 5 Nov gNB on Cubotron ~1 km away, receiver on the BC roof, 32x32 MIMO): Cubotron
  from OSM at local (814, 249) m (851 m from the origin), ground -3.4 m (swisstopo); gNB assumed 20 m over its ground
  (z 16.6), rx BC roof NE corner (15, -20, 19.1), boresight 19 deg; 32T = 8x2x2; users 30-600 m from the gNB;
  tablet next to the receiver (assumed). gen_dumps --absolute --ue-dist; sim_burst.sh GEO_GNB/GEO_RX/GEO_BORE.
  TRS profile sim vs real b9: 8 m -18.6 / 9 m -19.4, 20 m -20.6 / -20.5, 71 -20.8 / -20.8, 79 -22.9 / -22.9 dB,
  spread 34 / 35 ns -> static scene matches. Precoded statistics, same script (gs_trs.py) on real b9 vs sim:
    anchor eps median / pass:  tablet real 0.87 / 3%   sim 0.14 / 77%
                               other  real 0.16 / 74%  sim 0.74 / 15%
                               blind  real 0.64 / 29%  sim 0.74 / 17%
  -> tablet and others are INVERTED. Hypotheses: others: their DMRS reaches our receiver much stronger in reality
  (real sidelobes / scattering near the gNB) - compare DMRS-to-TRS power per kind; tablet: something tablet specific
  (rank 2 / Type II / how its records are despread) - inspect the real tablet DMRS/TRS ratio structure.
  No scenario target is visible in this geometry yet (needed for the experiments, not for the validation).

  2026-11-16 CALIBRATED MODEL (prg_stats.py = power vs TRS, PRG-to-PRG precoder change, same-precoder share per kind):
  real b9 showed reciprocity-like tablet precoding (77 deg / 5% same between adjacent PRGs) and every kind at ~TRS
  power -> model: TRS = fixed beam toward the receiver (SSB-like), users +-5 deg around BC's direction 500-1100 m
  from the gNB, tablet = reciprocity per PRG on its own (receiver's) channel with arbitrary phase per PRG, -8.6 dB
  (rank-2 split / TRS boost), 3 dB per-PRG amplitude; users' rays 100 ns; SNR 19 dB (gen_dumps --trs-beam toward_rx
  --tablet-precoding recip --tablet-power-db -8.6 --tablet-amp-db 3 --ue-ds 100 --ue-az 5 --ue-dist 500,1100).
  Tuned on b9, HELD OUT b2 with the same parameters (scratchpad validate.sh / run_validation.sh):
                      sim b2            real b2           (real b9)
    other  pass/eps   70% / 0.15        72% / 0.13        74% / 0.16
    other  jump/same  17 deg / 29%      17 deg / 30%      14 deg / 42%
    tablet pass/same  0% / 5%           27% / 13%         3% / 5%
    blind  pass/eps   32% / 0.41        21% / 0.65        29% / 0.64
    TRS 1-|rho|^2 20 ms  0.025          0.249 (!)         0.075
  -> codebook users reproduced on held-out data; tablet reproduced within the real burst-to-burst spread (3% b9 vs
  27% b2); blind a bit optimistic; real b2 TRS much less stable (environment/cart?) - not modelled. Accept for the
  experiments with these caveats.

  2026-11-16 EXPERIMENT GEOMETRY (option 3; real 1 km geometry: no ground near BC sees the gNB, even at 40-50 m
  antenna height; nearest spot seeing both 160-225 m away). find_gnb.py: ELG roof edge (198, -14), antenna 17 m above
  ground, 203 m from the area, sees 52% of it; receiver (2, 4, 1.5) clean LoS (next path -13.8 dB at +30 m), array
  facing west (180 deg). Visible area = the east-west corridor between BC/INN and MXH + a NW patch (vis_area.py).
  Scenario precoded_elg_200m (gui/scenarios in the worktree; ~/precsim/scen_elg/{spec,settings,agents}.json,
  scenario.npz): 2 walkways along the corridor (2.5 + 2 /min), lane y=-4 (bikes/cars/vans 1.5 /min), 4 people in the
  NW patch, 300 s: 39 agents, ~11 visible at a time. First burst (b9 masks, t0 200 s, calibrated precoding, users
  50-300 m +-20 deg): runs sim_elg_t1_b9_{per,omp} on the 8051 GUI.

- [DONE 2026-11-16] 4b6. SIMPLE PRECODING LADDER (user: the calibrated geometry was too complicated; scenarios and
  results of 4b5 deleted, simulator code / masks / validation logs kept). Old scenario, old geometry and site, one
  receiver, 1 port, the TABLET's real scheduling of b9 (masks_tablet.py; 60/68 PRGs median), gen_dumps --precoding
  random (per UE/PRG/epoch one random complex gain, same for all paths). L0 full band / L1 real holes / L2 holes +
  random precoding; so_light MODE=tab, sensing_offline mirror off.
  Found and fixed in NR-Scope (uncommitted):
  (1) detector edge: almost every false detection with gaps sat in the last stored bin -> NR_TDD_VERDICT_EDGE (5),
      NR_TDD_EDGE_GUARD_BINS 2;
  (2) static path fit with its delay derivative (remove_path: alpha K + beta dK/du) - did NOT remove the residue;
  (3) ROOT CAUSE of the gap residue: with gaps each snapshot sees the static scene through its own mask, so mean
      removal after the transform cannot reach it and the kernel fits model only the resolved paths -> NEW
      sensing.static_subtract (+ static_weight, default 1/1024): per stream and lattice position a running mean of
      the static scene over the snapshots that measured it, subtracted right after alignment, before the taper and
      the transform (nr_ue_sensing_static_subtract, arrays in the history). Exact for any static scene.
  Results (car / person found, other per map):            periodogram           OMP
    L0 full band                                          14/16 7/9 0.50        13/16 8/9 1.0
    L1 holes (before)                                      3/16 5/9 9.7         10/16 7/9 13.6
    L1 holes, static_subtract                             12/16 7/9 0.38        13/16 8/9 0.75
    L2 holes + random precoding, static_subtract          12/16 7/9 0.38        13/16 8/9 0.31
  -> map residue flat (~1 dB, as L0); holes + per-PRG precoding cost ~1-2/16; OMP's extra detections were the static
  residue, not the method. GUI 8051: sim_tab_L{0,1,2}_*_{per,omp} (before) and *_staticsub (after).
  Next: check static_subtract on REAL data (light eval, injected point targets), slow targets (weight vs speed),
  long runs + trackers; ROC per vs OMP at equal false rate.

  2026-11-16 GUI scoring (user: keep the GUI's own analyze.py, longer runs). Why the short runs scored badly:
  (a) one 2.45 s burst -> no track can confirm; (b) so_light GCAL forces the direct path to 0 deg (needed on the
  X410's unknown cable phases, wrong in simulation): AoA biased -38 deg, positions 11 m off -> GCAL=0 for sims;
  (c) with static_subtract the maps' bin_los was the residual's peak (a target) -> analyze.py excess = (bin - bin_los)
  shifted. FIX in NR-Scope: the static estimate's own delay profile gives the direct path's bin per stream
  (hist->st_bin, refreshed every 64 updates), used as bin_los when static_subtract is on; config with static_subtract
  uses clutter_kernel: false (mean mode; no static paths left to fit). 2-burst check in the GUI: detection 95%, range
  bias -0.1 m (rmse 0.9), angle median 0.6-1.9 deg, position median 0.55-2.0 m, coverage 80%, false 0.31/map.
  Long runs (scratchpad long_level.sh + publish_long.py): 24 consecutive bursts (tablet masks of b0-b5,b9 cycled, t0
  every 2.5 s), GUI runs long_{L0_full,L1_holes,L2_prec}_{per,omp}.

  2026-11-16 LONG RUNS (60 s = 24 bursts, 388 maps, 15 targets >= 10 maps in view; GUI analyze.py scoring; tablet
  occupancy; 1 receiver; impairments OFF; static_subtract on, clutter_kernel off, GCAL=0):
                         detect  range rmse  speed rmse  angle med  pos med  coverage  false/map
    L0 full    per       0.78    1.46 m      0.29        1.3 deg    1.18 m   0.76      1.67
    L0 full    omp       0.90    1.02 m      0.12        1.0 deg    1.07 m   0.82      1.28
    L1 holes   per       0.76    1.51 m      0.29        1.4 deg    1.41 m   0.75      1.05
    L1 holes   omp       0.88    1.07 m      0.13        1.0 deg    1.25 m   0.82      0.87
    L2 prec    per       0.75    1.52 m      0.28        1.3 deg    1.34 m   0.72      1.06
    L2 prec    omp       0.88    1.11 m      0.16        1.3 deg    1.27 m   0.79      0.77
  -> holes + per-PRG precoding cost 2-4 points of coverage; OMP better than the periodogram at every level
  (+10-12 points detection, lower range error, fewer false). Missing: 'precoding ignored' baseline at L2,
  impairments rerun, 10 trackers on these detections, all users' occupancy.

  2026-11-16 night: new switches. so_light PREC_IGNORE=1 (precoding-unaware baseline: one complex number per
  snapshot against the TRS, no per-PRG reference). gen_dumps --impair: the twin's nonideal.py defaults through
  SionnaChannel (array errors/coupling, RCS fluctuation, micro-Doppler, vehicle aspect, wall echoes) + clock per OFDM
  symbol (Wiener phase noise 50 Hz, residual CFO N(0,5 Hz) re-drawn every 20 ms SSB, sample clock 0.1 ppm continuous,
  removed by so_light dedrift). Scripts: scratchpad long_run.sh (one generation, several so_light variants),
  night.sh (8 runs: notab L1/L2 (+ign), tab L2 (+ign), imp tab L0/L1/L2 (+ign), imp notab L1/L2 (+ign); notab with
  MINPRB 0 and 48). Smoke (1 burst, notab + impair + random): OMP det 0.36-0.40, per 0.08, false 5.4/map (omp);
  precoding ignored: nothing detected, angle 41 deg off.

  2026-11-16 TRACKERS on trace-driven detections (long_L1_holes_omp, 1 rx, harness + tune.py --setups 1, 3 min /
  1.4 GB per harness pass; tune.py got a --setups option): GOSPA best (m) / missed / false per map / false tracks /
  cover people-cars: BSense 3.46/0.87/0.31/11 of 39/65-70%; GGIW 3.57; Bauhofer 3.58; Gunes 3.62; GM-PHD 3.65;
  GNN 3.73; OURS 3.89/0.43/0.82/6 of 31/83-83% (7th); JIPDA 4.00; CellSense 4.09; MHT 4.11; BP 4.15.
  [OPEN] ours: fewest missed and false tracks but 0.82 false estimates per map -> check whether its positions are
  > 5 m off (GOSPA c) or coasting tracks; grid edge again (p_d 0.9, r_confirm 0.95): extend the grid.

  WHY ours ranks 7th on GOSPA (long_L1_holes_omp, defaults): missed+false per map ours 0.32+0.99 = 1.31 vs BSense
  1.12+0.12 = 1.24 (about equal); the gap is the localisation term (1.58 vs 0.84 m): ours tracks the hard targets
  BSense drops (per-target coverage ours 60-96%, BSense 0-84%; target 12 81 vs 0%, 16 81 vs 17%). GOSPA rewards
  abstaining on hard targets -> report GOSPA with coverage, missed, false; also GOSPA on targets both cover.

  2026-11-17 night phase 1 done. PRECODING IGNORED vs TRS ANCHORING (tablet, L2 random precoding, same dumps, 60 s):
    OMP anchored det 0.88 cov 0.79 false 0.76 | OMP ignored det 0.61 cov 0.51 false 0.51
    per anchored det 0.75 cov 0.72 false 1.06 | per ignored det 0.26 cov 0.27 false 1.71
  -> anchoring is worth +28 (OMP) / +45 (per) points of coverage with the tablet.

  NO TABLET, L1 (other users + blind pieces, tablet PRGs removed, no impairments, 60 s):
    OMP MINPRB 0 det 0.72 cov 0.63 false 0.50 | MINPRB 48 det 0.73 cov 0.64 false 0.49  (tablet L1 OMP: 0.88/0.82/0.87)
    per MINPRB 0 det 0.50 cov 0.46 false 0.82 | MINPRB 48 det 0.49 cov 0.49 false 1.11
  -> sparse occupancy costs ~18 points of coverage (OMP); dropping symbols < 48 PRBs makes no difference.

  NO TABLET, L2 (random precoding per user/PRG/epoch, no impairments): OMP anchored cov 0.59-0.60 false 2.2-2.3/map
  (det 0.63); per anchored cov 0.37-0.40 false 0.8-1.0; PRECODING IGNORED: OMP cov 0.02 false 10.0, per cov 0.00
  false 10.0, angle 49-56 deg off -> without the tablet, ignoring precoding destroys sensing; anchoring keeps most
  of L1's coverage (0.60 vs 0.63) but OMP false detections rise (2.3 vs 0.5): [OPEN] many users' precoder residue
  read as extra paths by OMP -> check the OMP threshold / epoch length on notab L2.
  IMPAIRMENTS, tablet L0: OMP det 0.91 cov 0.79 false 1.01, angle 2.4 deg (array errors), pos 1.9 m; per cov 0.68
  false 2.38 (wall ghosts?).

  IMPAIRMENTS, tablet: L1 OMP cov 0.78 false 1.66 angle 2.2 | L2 OMP cov 0.77 false 1.59 | L2 ignored OMP cov 0.56
  false 1.13, per cov 0.13 false 2.42 (per anchored 0.69). Ladder with impairments: 0.79 -> 0.78 -> 0.77 (OMP).
  All 8 night runs done 03:49, no failed burst.

  IMPAIRMENTS, NO TABLET L1: OMP cov 0.56-0.60 but FALSE 7.4-7.7/map (0.5 without impairments, 1.66 imp+tablet);
  per cov 0.30-0.34 false 0.9-1.25. ABLATION (4 bursts, OMP false/map, det): all 8.4/0.48; no phase noise 0.33/0.73;
  no residual CFO 8.7/0.48; no wall echoes 8.6/0.52 -> CAUSE = receiver phase noise (50 Hz: ~2.5 rad random walk
  between TRS occasions) not removed on sparse blind snapshots. clutter_los_norm cannot be used (map stage, after
  static_subtract removed the direct path). FIX (NR-Scope, toggle sensing.static_phase_norm, off by default):
  in nr_ue_sensing_static_subtract, the snapshot's common phase against the static estimate (>= 12 shared
  subcarriers) is removed from the whole snapshot before learning/subtracting. Configs c_{per,omp}_stkpn.yaml.
  selftest + omp test PASS. Test running (notab imp / tab imp / notab clean, 4 bursts, stk0 vs stkpn same dumps).

  static_phase_norm test (4 bursts, same dumps): NO EFFECT (notab imp OMP false 8.40 -> 8.29; no harm elsewhere):
  so_light's per-snapshot c(t) already removes the common phase. The false detections are OMP-only (periodogram on the
  same maps 0.65/map): markers at < 20 m excess range spread over ALL speeds (|v| > 10 m/s many) = the direct path
  smeared in Doppler by phase noise; the periodogram's guards/verdicts drop it, OMP's global median noise does not.
  FIX 2 (toggle sensing.omp_local_noise, off by default; config c_omp_stkln.yaml): a cell must also exceed thr x its
  own range bin's Doppler median (first correlation map). Deterministic rerun (stk0 reproduces exactly):
    notab imp  false 8.40 -> 0.48, det 0.48 -> 0.35, cov 0.70 -> 0.27 (2 targets only)
    tab imp    false 0.92 -> 0.22, det 0.88 -> 0.90, cov 0.86 -> 0.86
    notab clean false 0.22 -> 0.08, det 0.65 -> 0.50, cov 0.52 -> 0.52
  [OPEN] trade-off without the tablet: needs ROC at equal false rate on 24-burst runs (stk0 vs stkln), and a
  check of which true targets the local floor loses (near the direct path?). static_phase_norm kept as a toggle.

  2026-11-17 MORNING SUMMARY (60 s runs, 1 rx, GUI analyze.py; no-tablet with MINPRB 48). det / coverage / false per map:
                         OMP               periodogram
    tab   L0            .90/.82/1.28      .78/.76/1.67
    tab   L1            .88/.82/0.87      .76/.75/1.05
    tab   L2            .88/.79/0.76      .75/.72/1.06
    tab   L2 ignored    .61/.51/0.51      .26/.27/1.71
    notab L1            .73/.64/0.49      .49/.49/1.11
    notab L2            .63/.60/2.23      .41/.40/0.97
    notab L2 ignored    .07/.02/10.0      .02/.00/10.0
    imp tab L0          .91/.79/1.01      .79/.68/2.38
    imp tab L1          .85/.78/1.66      .75/.65/1.39
    imp tab L2          .85/.77/1.59      .75/.69/1.38
    imp tab L2 ignored  .64/.56/1.13      .21/.13/2.42
    imp notab L1        .64/.60/7.38      .31/.34/1.25
    imp notab L2        .50/.49/9.35      .22/.15/0.93
    imp notab L2 ign.   .07/.00/10.2      .01/.00/9.66
  TRACKERS (tuned, OMP detections), GOSPA best: tab L1 ours 3.89 (7th, best BSense 3.46); tab L2 ours 3.87 (best
  BSense 3.52); notab L2 ours 4.04 (best BP 3.44); imp tab L2 ours 4.31 (best BSense 3.93); imp notab L2 (9.4 false
  dets/map) ours 3.98 = 2nd (MHT 3.68) with coverage 72/47% (people/cars) while BSense 19/0, GNN-EKF 29/21, Bauhofer
  36/18. Ours: highest or near-highest coverage everywhere; GOSPA favours conservative trackers on clean detections,
  ours is most robust under heavy clutter. Table builder: sens2 tools/precoding_sim/track_table.py.

  2026-11-17 PRECODING DEBUG, NO TABLET (12 bursts, 5 targets, OMP det/cov/false): L1 (no precoding) .61/.54/.28 =
  the CEILING; L2 epochs cut 10 ms .52/.56/1.81, 20 ms .50/.51/2.51, 50 ms .49/.51/3.66, detected only .47/.53/4.48
  (longer epochs mix precoders: worse). Local floor: false -> 0.04-0.06 but det -10..12 points (also at L1).
  -> the limit is the SPARSE OCCUPANCY (69 snapshots/map vs 237 with the tablet), not precoding (~9 points det).
  Levers: longer maps, CSI-RS (283 slots/s; the 8-port conversion worked in sim but failed on real data earlier),
  thresholds tuned by ROC, why cars are covered less (45% vs people 72%).
  NR-Scope: sensing.target_dv (m/s) now a runtime key (was NR_SENSING_TARGET_DV_MS, default 0.3 = 0.146 s maps);
  configs c_{per,omp}_dv15 / _dv10 (0.29 / 0.44 s). long_run.sh: KEEP_EST=1 keeps so_light outputs for re-sensing;
  publish step bug with 5-field variants fixed. Test running: dbg2_notab_L1/L2 x {stk0, dv15, dv10}.

  2026-11-17 REAL PRECODER PERSISTENCE (b9, antenna 0, per-PRG HW/H against the nearest TRS, slot common phase
  removed; similarity |sum v1* v2| / sum|v1||v2|): tablet 20 ms 0.97, 50 ms 0.91, 100 ms 0.78, 200 ms 0.50, 500 ms
  0.23 (random 0.27); blind 0.92 / 0.85 / 0.83 / 0.48 / 0.36 (random 0.40); symbols 2 vs 11 of a slot 1.00.
  -> precoders hold ~50-100 ms, NOT the 19 ms the simulator draws (gen_dumps --epoch-ms 19) nor the 10 ms so_light
  cuts (HYB_MAXEP=0.01): the sim L2 and the processing were both too pessimistic. Test running: dbg3_notab_L2 with
  --epoch-ms 80, HYB_MAXEP 10/50/100 ms/detected only. Also to redo: HYB_MAXEP on the real 5 Nov data.
  so_light: env FIG_DUMP=<npz> saves antenna 0's TRS (H) and DM-RS (HW) as anchoring starts (figure input).
  Figure: ~/Documents/EPFL/pdm/images/make_fig_precoding.py (laptop).

  MAP LENGTH (no tablet, 12 bursts, OMP det/cov/false): L1 0.146 s .61/.54/.28, 0.29 s .74/.73/.73, 0.44 s .80/.76/1.13;
  L2 (sim precoders ~19 ms) .52/.56/1.81, .59/.63/3.44, .54/.59/4.37 -> integration works without precoding; with
  fast-changing precoders longer maps add false detections. Maps per 30 s 189 / 96 / 60 (fewer tracker updates).
  Next: dbg3 (precoders persisting 80 ms) re-sensed at 0.29/0.44 s (scratchpad resense.sh on kept estimates).

  ROOT CAUSE no-tablet L2 = ASSOCIATION of blind pieces to users. Oracle (sim writes owners.npz: true owner x precoder
  epoch per PRG/slot; so_light ORACLE_OWN=1 uses it as reference keys), precoders persisting 80 ms, OMP det/cov/false:
    0.146 s maps: estimated .51/.56/2.49 -> oracle .59/.64/0.28 (L1 ceiling .61/.54/.28)
    0.29 s maps : estimated .65/.61/4.17 -> oracle .72/.72/0.66 (L1 ceiling .74/.73/.73)
  snapshots per reference number: estimated 3.4, oracle 14.1. With correct association precoding costs ~nothing.
  Longer epochs alone (50/100 ms, detected) are worse: they average different users' precoders on a PRG.
  track_keys (so_light) limits: templates expire after 20 ms (dt), >= 4 shared PRGs, + HYB_MAXEP 10 ms cut.
  Now env BT_DT / BT_CMIN / BT_NMIN (defaults 0.02 / 0.9 / 4). long_run.sh KEEP_DUMPS=1 keeps dumps on sens8.
  Test running dbg5 (base, oracle, dt 0.1, + nmin 2, + cmin 0.8; maps 0.146/0.29 s).
  Path to > 80%: association (oracle level) + 0.29-0.44 s maps (L1 at 0.44 s: det .80 cov .76) + tracker.

  2026-11-17 DCI CHECK (5 Nov b9, ch0 golden_sniffer.csv / _grants.csv vs the DM-RS masks): GS found 7675 DCIs from
  119 RNTIs (85% tracked, 15% search); 5897 DL (sizes 54-59); grants with PRBs 2467 (2054 = the tablet). Grants by DL
  size: 54 2054 of 3040 (only the tablet), 57 138/1210, 58 0/714, 56 275/502, 55 9/338, 59 0/93.
  Blind DM-RS PRG-slots not owned by a decoded record 28318: 60% in a slot with a decoded DL DCI without placed PRBs
  (FDRA), 40% with none (DCI missed or k0 > 0). 514 slots have exactly ONE un-placed DL DCI.
  -> association is mostly a DCI problem: (1) FDRA per DCI size AND per UE (others' size-54 never lock: BWP/config
  per UE?); (2) association by elimination with the slot's un-placed DCIs; (3) per-PRG continuity for the rest.

  2026-11-17 FDRA OFFLINE (b9, scratchpad dci_check/fdra_lib.py, fdra_llr.py): parsing validated on the tablet (size
  54: offset 2, type 0, BWP 273 -> 100% of 2054 grants exact). Scored against the slot's DM-RS occupancy (chance-
  corrected LLR): size 56 locks off 2 type 0 (3.3 vs <=1.2); per UE: 0x70bd s56 off2 t0 4.79, 0x789e s57 off3 t0 3.72,
  0x74a2 s56 off2 t0 3.61, 0x7106 s58 off3 t0 2.72 (GS places 0 size-58 grants). Layout: id bit + 1-2 bit BWP
  indicator + 18-bit RBG bitmap (type 0, P 16). Many busy UEs never lock (their DM-RS invisible here: beamed away).
  RULE (type 0, off 2 for 54-56, 3 for 57-59) on all non-tablet DL DCIs: 31% decode fully onto DM-RS (62 UEs) but they
  own only 8% of the blind PRG-slots (2371/28318; k0=1: 4%). -> the '60% next to a DCI' was co-occurrence: the blind
  DM-RS mostly belong to UEs whose DCI is not decoded. FDRA fix = small, exact gain; per-PRG continuity association
  is the main lever for the blind pieces.

  2026-11-17 NEW ASSOCIATION (so_light BLIND_TRACK=2 prg_keys: per-PRG precoder continuity, several templates per
  PRG alive PP_DT 0.1 s, snapshot common factor by consensus, PP_TOL) + ROBUST FIT (so_csi ref_fix_hyb HYB_ROBUST=1:
  members with relative residual > 0.3 down-weighted (0.3/r)^2, refit). ASSOC_EVAL also prints the share of PRG-
  snapshots in MIXED references (all members of a mixed key get a wrong reference): old 1.6-2.1%; PP_TOL 0.10: 7.5-9%
  (detection collapsed to 30%); 0.05: 1.5-3.2% with oracle-level snapshots per reference (13.6/9.8 vs 14.1/10.5).
  12 bursts no tablet L2 (80 ms precoders), OMP det/cov/false at 0.146 s: old .51/.56/2.49; per-PRG 0.03 .49/.48/1.07;
  per-PRG 0.03 + robust .51/.59/0.76; oracle .59/.64/.29. ROC (OMP threshold 11.4/10/9 dB): old det .51/.58/.59 at
  false 2.49/3.11/3.58; per-PRG+robust .51/.56/.59 at 0.76/1.49/2.21 -> +8 points of detection at equal false (~2.3).
  0.29 s maps per-PRG+robust thr 10/9: .61/.64/2.71, .64/.66/3.68. rx2 generated (dbg5r2, same transmissions),
  fused runs fused_notab_L2_{base,base_dv15,orac_dv15,orac_dv10}_omp, trackers 1/2/fused running.

  2026-11-17 FUSED TRACKING, NO TABLET L2 (12 bursts, rx1 + rx2 same transmissions, harness defaults, OURS people /
  cars coverage, rx1 | rx2 | fused): old assoc 0.146 s 82/56 | 100/54 | 94/45; old 0.29 s 79/72 | 98/31 | 96/57;
  oracle 0.29 s 86/79 | 100/67 | 96/67; oracle 0.44 s 88/65 | 96/0 | 96/61.
  -> PEOPLE > 90% fused, fully non-cooperative, already with the old association. CARS 45-67% [OPEN]: fusion can
  LOWER car coverage below rx1 alone (72 -> 57): rx2 noisier (1-2 false/map vs ~0.5) -> check fused gating /
  association of cars. Baselines fused: MHT/BP/GNN similar or lower GOSPA on oracle detections, BSense and Gunes
  collapse (people 20-56%). Small sample (4-5 people, few cars), untuned. To do: per-PRG+robust association in the
  fused runs (rx2 reprocess), tuning, car investigation, longer runs.

  CORRECTION (2026-11-17): the fused 'people 94-96%' above is ONE person (28-68 maps); 'cars' = 2 cars + 3 bikes.
  Per class (OURS, 0.29 s maps, rx1 | rx2 | fused): old assoc person 79/98/96, car 72/31/57, bike 45/35/46; oracle
  person 86/100/96, car 79/67/67, bike 55/60/90. -> no conclusion on 80% from 12 bursts; bikes suggest association
  matters for vehicles and fusion helps with clean detections. Always report target counts per class.
  Running: imp_r1/imp_r2 (24 bursts, impairments, rx1+rx2, base vs per-PRG+robust, stk0/dv15/stkln/dv15ln) -> fused.

  2026-11-18 IMPAIRED FUSED RUN (24 bursts, no tablet, L2, 80 ms precoders, impairments, rx1+rx2; 5 people/5 cars/5
  bikes), OURS coverage people/cars/bikes: rx1 55-68 / 37-53 / 7-31; fused best 59-60 / 62-66 / 27-39. Per-map OMP
  (rx1): 0.146 s det .51-.53 cov .45-.49 false 6.6-9.2; with local floor false 0.3-1.3 but det .28-.48. New
  association: small gains. Fusion helps cars (rx1 ~47 -> fused ~64), not people.
  REAL PHASE NOISE (6 bursts b0,2,4,9,15,21; TRS common-phase increments, constant CFO removed; scratchpad
  phase_noise.py; so_light FIG_DUMP now works with USE_DMRS=0): std 0.046 rad at 0.14 ms, 0.16 at 0.5 ms (slot 5->6,
  likely a resource/beam difference), 0.035 at 20 ms, 0.053 at 40 ms -> ~0.01 Hz equivalent linewidth vs the 50 Hz
  modelled (2.5 rad at 20 ms): the impaired runs were ~70x too pessimistic in phase noise; residual CFO < ~0.3 Hz
  (modelled 5 Hz). Rerun with calibrated impairments (linewidth 0.1 Hz, cfo_res 0.3 Hz, still pessimistic): cal_r1 /
  cal_r2 -> fused_cal_{base,rob}{,_dv15}_omp, running.
  FIELD OF VIEW (fused_notab_L2_base_dv15): vehicles in rx1's +-60 deg 86-89% of the time, rx2's 23-29%, both 12-15%
  -> little overlap for fusion in this scene. FUSED TRACKER BUG [OPEN]: car 12 (100% in rx1 FOV, 0% in rx2) rx1 75%
  -> fused 0%: rx2 frames likely count as misses outside rx2's view (P_D map per receiver not applied?).

  2026-11-18 PRECODER STATISTICS (5 Nov, 26 bursts, 4 antennas, ports 1000/1001; /mnt/mydrive/field/precstats,
  analyse.py; so_light FIG_DUMP + FIG_ALL, gs_all P1001=1). Per-PRG g = r^H y / |r|^2 vs nearest TRS.
  1 adjacent-PRG phase steps (4th-power slope removal): tablet +-90 deg 17%+17%, 0 3%, 180 3% (uniform 5.6% each);
    blind +-90 18%+17%, 0 7%. 2 magnitude change at +-90 steps 0.6 dB vs 1.8 dB at 0, 2.5 dB at 180: quarter turns are
    unit-modulus. 3 port1001/port1000 ratio not constant across PRGs (0.12) and steps by 0 or 180 deg (11%/10% vs 3%):
    rank-2 [v; phi v], [v; -phi v] with phi in {1,j,-1,-j} per PRG (Type-I-like co-phasing; inferred).
  4 time: tablet 1 wideband change in 37604 slot pairs; persistence 10 ms .997, 20 .97, 50 .89, 100 .75, 200 .47;
    temporal SVD flat (26/19/15%): slow continuous drift, not beam switching.
  5 blind pieces next to tablet slots match its precoders (median sim 0.998, 80% >= 0.9): tablet DMRS with a missed
    DCI appear as blind -> association idea: match blind pieces to known UEs' recent precoders.
  6 MODEL FIT (corrects '8-12%', single-antenna): 4-antenna shape residual DMRS vs TRS (< 2 ms) median 0.21-0.27 vs
    noise floors DMRS sym2/sym11 0.02-0.04, TRS sym4/sym8 0.035-0.04 -> one gain per PRG explains ~75-80% of the
    energy; the rest is path-dependent precoding (the 4b2 limitation, now measured).

  2026-11-18 CALIBRATED FUSED RUN (24 bursts, no tablet, L2, impairments calibrated on the real cell, 5/class), OURS
  people/cars/bikes fused: old assoc 82/84/46, new assoc 79/76/59 (GOSPA 4.71). FUSED TRACKER FIX (sens2 sim/track.py,
  EKF_DEFAULTS view_gate, view_gate_m): car 12 drifted 7 -> 24 m in fused mode, dragged by rx2 detections through rx2's
  bistatic model while outside rx2's view. 'hard' (no update from a receiver that cannot see the track): people/cars/
  bikes 76/63/37 (worse: the LOS+FOV grid is too strict); 'soft' (only a detection with a position within 5 m):
  old assoc 91/78/46, miss 0.70, GOSPA 5.26; new assoc 86/77/54. Default now 'soft'.
  18 Nov RECORDING: the first 120 s capture (TX/RX0, gain 38) held no signal (40 dB down, no SSB): the antennas are on
  RX1. Re-recorded on RX1 gain 45 (0 gaps, levels 2600-4760 rms): 20261118_154424 (archived); empty captures deleted
  with the user's OK. rec_cont.sh takes ANT= / GAIN=; rec_try.sh wraps it (rec_direct exit hang after power cycle).
  Night 18-19: 4b2 test (32-port Type I vs scalar, rx1 12 bursts) + scenarios dens3_seed2, dens6_seed1 (rx1+rx2).

- [OPEN, TO TEST] 4b2. Geometric precoder-ratio model (proposed 2026-11-14, replaces nothing yet).
  Problem: today each DM-RS PRG/epoch = ONE complex number x nearest TRS, i.e. the DMRS/TRS beam-gain ratio
  r(theta) = g_D(theta)/g_T(theta) is assumed equal for every path. False when paths leave the gNB at different angles:
  target Delta from the receiver is ~Delta/D off the direct path (D = 200 m planned: 20 m -> 6 deg, 50 m -> 14 deg,
  100 m -> 27 deg; sidelobe ~12-14 deg for 8 columns). Silent failure: the fit follows the direct path, targets and
  clutter keep a precoder-dependent residue (Doppler spread, peak loss, clutter left). Evidence it exists: 6f anchor
  test, tablet DM-RS fits the TRS with one number only 10% (other UEs 61%).
  Model: r(theta) ~ a0 + a1 (theta - theta0) per PRG/epoch, fitted on strong static paths; each range-angle cell
  corrected with r(theta_cell). theta of a cell from bistatic geometry (delay + receiver AoA -> reflection point ->
  departure angle). Unlike 2-beam (4b): no free beams, angles computed, so no per-PRG gauge freedom.
  Assumes: r smooth over the sector (fails near beam nulls), single bounce, far field / azimuth only, TRS beam fixed
  (measured, 6f test 2), precoder constant within PRG+epoch.
  Needs: gNB + receiver positions (few m) and gNB array azimuth; receiver AoA = 4 channels with ONE common timing from
  the IQ (not the per-antenna replays) + calibration; strong static reflectors spread in angle.
  Tests, cheapest first: (1) quantify the 6f tablet/others anchor residual vs geometry (existing data);
  (2) sim twin, 64T array, random UE precoders, targets at Delta/D 5-30 deg, D = 200 m: scalar vs 2-beam vs geometric
  -> target SNR loss, range bias, Doppler spread, clutter residue; is a0+a1 enough?; (3) target absorption by the fit
  with 2-4 PRB PRGs (share a1 across PRGs/epochs if needed); (4) real data: common-timing 4-ch TRS extraction, stable
  AoA of static reflectors, departure-angle spread, then scalar vs geometric on clutter + injected targets.
  Kill criterion: no gain over the scalar in (2) at D = 200 m.
  Related: Nataraja et al. TVT 2025 (only prior precoding algorithm) groups identical-precoder sub-bands; needs
  per-sub-band target detection, not applicable to 2-4 PRB PRGs / ~19 ms epochs; possible simplified baseline.
- [IN PROGRESS] 4c. Precoder epochs last ~19 ms median: joining epochs inside one map is still open.
- [OPEN] 4d. Subframe phase step (~60 deg per ms) in NR-Scope output: correct at the source.

- [DONE 2026-11-06] 4e. Coherent merge vs incoherent sum of range-Doppler maps (gs_inc.py, ana/inc_*.json, 144 windows,
  tablet band cut into pieces with unknown numbers + synthetic per-PRG beam mismatch). Two halves, map similarity to the
  full-band map / moving peaks found / range width: true number 0.98 / 100% / 9.8 m, edge 0.96 / 90% / 9.8 m,
  no correction 0.92 / 80% / 12 m, incoherent 0.89 / 65% / 12 m (40% false peaks), one piece only 0.78 / 30% / 15 m.
  Coherent wins whenever possible (even uncorrected); incoherent only for pieces that cannot be joined.
  Real other UEs (CH=0,3, ana/incr_*.json): only 12 windows usable, inconclusive.

## 5. Validation

- [OPEN] 5a. No ground truth for absolute range: need a known (moving) target, or a check of the
  direct path against the geometry.

- [IN PROGRESS 2026-11-09] 5b. Simulation with known targets on the real 5 Nov occupancy (sim_sense.py, sim/*_m*.json,
  sim_sum.py): 8-port gNB, 4 RX, 20 static scatterers calibrated to the measured stitching error (~2 deg/edge in sim vs
  1.6 deg measured), 2 movers (-25 dB, 15-250 Hz), 150 ms windows, static clutter removed before Doppler; 376 windows.
  Found / false detections per window (precoder mix 0, then 0.2):
    ideal (full band every slot) 100% / 1.3
    without tablet: DM-RS with precoder removed perfectly 99% / 86 (sparse irregular sampling, not precoding, is the
      limit); TRS alone 94% / 21 (Doppler aliases every 50 Hz); TRS + DM-RS via TRS anchor 93-95% / 62-73 (worse);
      TRS + incoherent 58-60% / 14 (targets 23 dB over floor); CSI-RS alone 100% / 27 (aliases)
    with tablet: DM-RS oracle 100% / 4.7; today's pipeline 54% / 34 and 42% / 41; + TRS no change; TRS + incoherent 9-11%.
  Caveats: precoder model (2 beams, mix per PRG) and scene only calibrated on the stitching error; threshold-based
  detection; stitching diagnostic gives nan for mix 0 (bug to check).

## 6. Seeing more UEs (limits everything above)

- [OPEN] 6d. 2026-11-05: antennas 1-2 decode only the tablet's DCIs (other UEs only on 0 and 3). Fix: NR-Scope takes DCIs
  from the best antenna and extracts DM-RS on all four.

- [IN PROGRESS 2026-11-09] 6e. Blind DM-RS snapshots: dump the DM-RS channel of every survey run (PCI-scrambled, no DCI needed) so
  unplaced UEs (14% of on-air DM-RS; the only source in a light network) enter stitching/merging. Needs NR-Scope change.
  Implemented 2026-11-09 (uncommitted, golden_sniffer.cc/.h): golden_sniffer.dump_blind writes GSD1/GSD2 records with
  rnti 0 (strict map of the first DM-RS symbol as mask, symbols 2 + type-A additional positions sharing its PRBs); the
  survey now runs in every slot (also with no converged UE). Check: blind values identical to the grant dump on shared
  PRGs (4088 records, max diff 0). Replay: replay4ch_blind.sh -> b*/replay_blind_ch{0,3}; analysis gs_usemap_v6e.py
  (DIR=replay_blind_ch; blind = one more UE minus placed PRGs, reference only if no placed UE) -> ana/usemap6e.
  Burst 8: on air not placed 139k -> 59k PRB-slots; with tablet the blind PRGs are mostly rejected by the gate
  (anchor 1826, edge 219, rejected 24433, not joinable 7639 PRGs); without tablet DM-RS used 5% -> 21%, windows with
  nothing 13% -> 0%, but median band per window only 0.3% -> 0.7% (no anchor: most blind pieces not joinable).
  All 23 bursts (2026-11-09, ana/usemap6e): with tablet: coherent 86% of DM-RS (blind 1.3%), incoherent only 9.9%,
  still no channel 4.1% (was 14.4%); median band per window 62.8%. Without tablet: coherent 44% (39% blind, mostly
  blind as its own reference, i.e. stitched across UE boundaries inside a slot: optimistic), incoherent only 31%,
  no channel 24%; median coherent band per window 1.6% (was 0%), p90 14%, windows with nothing 4% (was 68%), but
  32 epochs per window, longest 2 ms: blind pieces are short. Usage map page updated (blind + incoherent colours).
  Change needed (2026-11-06): in the survey block of golden_sniffer.cc write a GSD record with rnti 0, mask = strict
  bitmap minus decoded grants, values from gs_dmrs_rb_pairs on symbols 2 and 11 (only DM-RS with N_ID = PCI and n_SCID = 0 is UE-independent;
  UEs with an RRC scramblingID0/1 or n_SCID = 1 stay invisible; cheap extension: also test n_SCID = 1); move the
  survey out of if(!tracked_now.empty()) so it runs with no tracked UE. ~16 MB/s per antenna. Testable on the 5 Nov IQ.
- [IN PROGRESS 2026-11-09] 6f. Periodic signals (TRS every 20 ms, CSI-RS sym 13) as the always-present time anchor between DM-RS bursts.
  2026-11-09 (gs_trs.py, gs_trs_lib.py, run_trs_all.sh; NR-Scope: GSC1 raw dump also for TRS symbols 4/8 under
  dump_blind, uncommitted). TRS: slots 5+6, symbols 4+8, every 20 ms, REs 2/6/10, n_ID = PCI, 96% of PRBs, coherence
  0.98. After removing the delay slope 1-|rho|^2 = 0.03 at 20/40/100 ms (stable reference). NEW: NR-Scope timing
  window jumps by up to ~16 samples between TRS occasions (relative slope, half of the 20 ms pairs) -> must be
  corrected before 150 ms maps and anchors. TRS as anchor (one number + relative slope, gate 0.3, <= 10 ms): other
  placed UEs pass 61% (tablet anchor 48%), blind 32% (tablet anchor 33%), tablet itself 10% (narrow beam != TRS beam).
  Light network per 150 ms window: DM-RS joined via TRS median 0.7% (p90 2.7%), TRS itself 6.6%, both 7.1% (p90
  11%), no empty window. CSI-RS (8 ports, slots 2/10 every 20 ms + 1/11/12): per port 1-|rho|^2 0.15 at 1 slot,
  0.30 at 20 ms (port 0). Tablet DM-RS per PRG in the span of the 8 ports (4 antennas, bursts 2/9/21): residual
  0.06 near, 0.17 with CSI-RS 500 ms away, 0.50 random -> precoder estimation from CSI-RS looks feasible (next).
  Usage map with the TRS (gs_usemap_v6f.py -> ana/usemap6f, page updated with two greens): median band per 150 ms
  window without the tablet: DM-RS alone 1.8% -> with the TRS 8.3% (p10 6.8%, p90 27.6%); with the tablet 63% -> 67%.
  DM-RS joined via the TRS anchor: 4.6% of DM-RS without the tablet, 1.4% with it. No maps built (no ground truth yet).
  CSI-RS span test on every kind (gs_csi.py, ana/csi_others.json, bursts 2/9/21, 4 antennas, values from the blind
  records): residual per PRG with 8 CSI-RS ports x 8 weights: tablet 0.05, other placed UEs 0.05-0.06, unplaced
  (blind) 0.06-0.07, TRS 0.05-0.06; flat over 0-10 ms from the CSI-RS. Stale CSI-RS (500 ms) 0.17-0.25, best single
  port 0.20-0.27, random 8-dim level 0.50. -> every PDSCH DM-RS and the TRS are the same 8-port channel H times
  different weights: a common object to fuse them (not tried yet; needs a simulation with known targets first).

  2026-11-14 TRS test 2, beam stability (trs_stab.py, ana/trs_stab.json, 26 bursts, 4 antennas, one complex number +
  relative delay slope): per antenna 1-|rho|^2 stays at the noise floor (4-symbol lag 0.02-0.05) at 20 ms, 100 ms,
  500 ms, 2 s AND across bursts 15 s-12 min (ch0 0.033, ch3 0.023) -> TRS beam fixed up to one number over minutes;
  w_T static is now measured, not assumed. Joint antennas 0+3 with ONE common number: fine at 4 symbols/1 slot (0.03)
  but 0.37-0.53 at >= 100 ms (p90 0.9 already at 20 ms, bimodal) -> inter-antenna phase jumps between the separate
  per-antenna NR-Scope replays (each run its own phase/timing corrections), not a beam change (the per-antenna shape
  does not move). Spatial signature needs the 4 channels processed with one common timing (open).

- [OPEN] 6a. Very few UEs placed (tablet, 0x7176, a few others occasionally).
- [OPEN] 6b. 58-bit DCIs cannot be locked from air: needs a per-size override in NR-Scope
  (only ~6% of them show DM-RS at our site anyway).
- [OPEN] 6c. FDRA cannot lock in a fully loaded cell.

## 7. Housekeeping

- [DONE 2026-11-09] 7a. Port-1001 dump (NRSCOPE_GS_DUMP_P1001 -> GSD2, gs_fdra.cc/golden_sniffer.cc) and the
  analysis scripts are not committed (never push).
  2026-11-09: committed (c311b0e, d011553 with 6e/6f) as flxinxout and pushed: github pvnis/NR-Scope branch goldenSniffer.
  Analysis scripts copied into the repo (nrscope/scripts/sensing, commit 1100d33, pushed). ~/field keeps working copies:
  copy changes back into the repo before committing.
- [OPEN] 7b. Run scripts/field_prep.sh before any recording or live run (sens6 rebooted 2026-11-05).
- [DONE 2026-11-09] 7c. Disk: 55 GB free after the 5 Nov recording (125 GB). Delete old recordings before the next long one.
  2026-11-05: 4 Nov recording deleted (29 GB) on request -> 73 GB free. 4 Nov results survive only in ~/field/ana (pkl/md).

  2026-11-09: 5 Nov recording (130 GB) moved to /mnt/mydrive/field (1.8 TB SATA SSD), checksums verified, link left
  at the old path; root disk 42 -> 172 GB free. New heavy data goes to /mnt/mydrive/field.
- [DONE 2026-11-09] 5c. Track-level simulation (sim_track.py: continuous street scene, 3 cars per 2.45 s burst, a 150 ms
  map every 100 ms, AoA from 4 antennas, ~/GIOVANNI/sim/track.py run_ekf_tracker, 25 bursts per scenario).
  False tracks per burst / coverage / GOSPA (m), without the tablet:
    ideal 0.7 / 75% / 6.6; DM-RS oracle 99 / 81% / 26.5; TRS 9.2 / 81% / 18.0; CSI-RS 10.6 / 81% / 17.9;
    TRS+DMRS@TRS 170 / 72% / 34.1; TRS+incoherent 9.8 / 73% / 17.0
  with the tablet:
    ideal 0.4 / 73% / 6.4; DM-RS oracle 1.9 / 71% / 8.1; today's pipeline 42 / 42% / 15.2 (+TRS no change);
    TRS 9.4 / 76% / 17.7; CSI-RS 9.4 / 77% / 17.7; TRS+DMRS@TRS 125 / 53% / 27.7; TRS+incoherent 11.6 / 16% / 10.8
  -> merging errors become ghost tracks (one-number TRS anchor worst); with a dense UE, perfect precoder handling
  would be near ideal, today's pipeline is far from it; without a dense UE the DM-RS sampling itself is the limit
  and the TRS/CSI-RS alone are the best source (Doppler aliases give ~10 false tracks per burst).
- [IN PROGRESS 2026-11-09] 5d. NR-Scope's own sensing chain (sensing_offline: alignment, clutter, Doppler, CFAR, AoA) on
  our estimates with injected targets (gs_to_nrsr.py writes record_estimates; so_score.py scores map2d.csv; so/ on sens6).
  NR-Scope changes (local, uncommitted): sensing.allow_gaps (fragmented grants), sensing.align_search_bins (timing jumps
  > 3 bins), alignment ignores gaps, CMake fixes for sensing_offline. Map range is only 50 bins (~120 m).
  Burst 8, targets 40 m +4 m/s and 80 m -3 m/s at -20 dB (found / maps):
    clean synthetic, full band: 32/32, SNR 35 dB.  Simulated precoders on the tablet's real occupancy: oracle 15/16,
    raw 0 maps (alignment fails), our stitching+gluing 15/16, CSI-RS one-number-per-PRG 15/16 (all SNR ~24 dB).
    Real recording: raw 3/14 (range smeared, ~10 ghosts/map), our processing 1/14, CSI-RS 1/8.
  Missing real effect found: the real precoder changes DIRECTION per PRG, not only phase. Per-PRG weights from the
  8-port CSI-RS: share along the snapshot's main beam median 0.38 (10th pct 0.09; random 0.125, pure phase 1); burst's
  main beam holds 29% of the energy (sim 55%). Estimates reliable: symbols 2 vs 11 agree |<w2,w11>|^2 = 0.87 median.
  -> every one-number-per-PRG correction (stitching, CSI scalar, anchors) is wrong on this cell. Next: use the full
  per-PRG weights (convert each PRG to one virtual beam through the 8-port channel), and a sim with this diversity.
  2026-11-09 later: full 8-weight CSI-RS conversion (so_csi.csi_fix8) = oracle in sim (16/16) but fails on real data.
  (a) Sim with the real direction diversity (K=3, MIX=1: share 0.34) still works for every method: diversity is not
  what breaks the real data. (b) NR-Scope timing-window jumps (1..16 samples) added to the sim break EVERY method,
  even the oracle; so_csi.dejump (per-snapshot slope inside PRGs) restores 15-16/16 in sim. On real data still 0-4/14.
  (c) The real CSI-RS is beamformed and changes: slot 2 vs slot 10 correlate 0.15 (0.85 up to an 8x8 port mixing),
  the same slot 20 ms later 0.38, 80 ms 0.11, no cycle; TRS stays 0.97 at 100 ms. So the 8-port H is valid only
  within one CSI-RS occasion; csi_fix8 per epoch leaves 29-50% residual on real data. Precoder itself is stable
  within an epoch (weights agree 0.8 a few slots apart, 0.3 across epochs).
  (d) OPEN: on real data our stitching+gluing puts injected targets at the right speed but 30-45 m too close (raw:
  range smeared around the truth). Next: TRS as ground truth for the stitched DM-RS delay profile, per PRG.
  (e) Not a ramp: after dejump the stitched direct path sits at ~5 m (raw 12 m), stitching gains trend -6 deg/PRG
  (~4 m). TRS check impossible on burst 8 (tablet never scheduled in TRS slots 5/6). Found instead: at each epoch
  change, the per-PRG phase of new vs old snapshot (<= 2 ms apart) after stitching+glue spreads 50 deg on the real
  cell (12 deg inside an epoch); sim with diversity + jumps: 29 / 12 deg and still detects. Per-PRG gluing
  (so_csi.glue_per_prg) makes it worse (sim 12/14 -> 3/8, real 0): errors random-walk per PRG. DROPPED.
  STATE: no method detects injected targets on the real 5 Nov data in NR-Scope's chain; all work in simulation
  (with dejump). Open: what in the real data the sim still lacks beyond the larger epoch-change mismatch; whether
  injection into processed estimates is a fair test; a recording with a real controlled moving target is needed.
  2026-11-09 evening, ROOT CAUSE FOUND: the 4 antennas of the 5 Nov recording come from 4 SEPARATE NR-Scope replays,
  each with its own sync loop and window: their direct paths sat up to ~40 m apart (b9 stitched: -41.5 / -13.4 /
  -6.1 / 0.0 m). NR-Scope applies one alignment per symbol to all chains and averages their maps, so every
  multi-antenna real-data result was broken (sim has aligned antennas: "works in sim, fails on real"). Live
  NR-Scope (one run, 4 chains, one window) does not have this. Also found on the way: my scoring was too strict
  (targets a few m off counted as missed), and range zero was not on the direct path.
  Fixes (gs_to_nrsr.py / so_csi.py): GSW1 record = NR-Scope window shift per slot (golden_sniffer.cc, uncommitted),
  undone per antenna (WSH=1) + straight-line drift; range zero on each antenna's direct path (los_to_zero, average
  profile, LOS0=1). Burst 9, targets at 40/80 m, found (10 m tolerance):
    raw 10/22 (-10 dB), 3/22 (-20), 1/20 (-25); stitching per epoch 3/14, 4/16, 2/16; stitching once per burst 0;
    TRS-referenced (so_csi.ref_fix) 12/12, 6/12, 8/12, all within 2 bins, range error -0.3..-1.6 m.
  Open: TRS variant gives 6 maps per burst vs 11 raw; other bursts need antennas 1-2 replayed with the GSW1/TRS
  dumps; TRS alone through NR-Scope; recheck earlier 4-antenna results for the antenna misalignment.
  Light network (so_light.py: TRS written on the DM-RS pair grid + every DM-RS except the tablet's, each PRG on the
  nearest TRS's reference; window moves undone; range zero per antenna on the TRS direct path). NR-Scope change
  (uncommitted): sensing.min_snapshots (default 32; TRS alone gives ~28 per 133 ms window). Burst 9, found / 30:
    TRS only: slow targets (0.6/-0.8 m/s) 7 (-10 dB) and 6 (-20 dB) with the right speed, but 16-20 more at the
    right range with the speed off by 2/4/6.. m/s (Doppler folding: 20 ms spacing = 50 Hz = 2 m/s); fast (4/-3 m/s) 0.
    TRS + other UEs' DM-RS, no tablet: slow 17/32 and 13/30, fast 6/30 and 2/30.
    TRS + DM-RS + tablet: slow 17/32, 15/30; fast 6/30, 18/30.
  Next: resolve the TRS speed ambiguity with the tracker (multi-hypothesis over 2 m/s copies) and the DM-RS samples.
  ALL BURSTS (run_gsw_all.sh + run_eval_all.sh, so/all, so_sum_all.py; 19 bursts with the tablet, light-network
  converter fails on b6/12/16/17/19/24: no UE decoded on all 4 antennas). Found (right speed, within 10 m):
                            slow -10 dB  slow -20  fast -10  fast -20   other/map
    tablet, raw                  5%        1%       62%       25%      0.4-10.6
    tablet, stitching           29%       14%       40%       27%      0.4-2.2
    tablet, TRS-referenced      72%       38%       75%       71%      0.6-4.5
    no tablet, TRS + DM-RS      63%       44%       22%       15%      0.3-4.7
    no tablet, TRS only         13%       14%        0%        0%      ~5   (+68-76% right range, speed off by 2 m/s steps)
  -> TRS-referenced DM-RS is the method; light network limited by the TRS speed ambiguity (tracker next).
  NIGHT 2026-11-09:
  (1) Range accuracy (all bursts, found targets): TRS-referenced median ~0 m, RMS 2.0-2.3 m, 83-89% within 1 bin
      (2.44 m); stitching RMS ~3.5 m (~60%); raw RMS 4-5 m (30-55%). The "10 m" was only the scoring window.
  (2) Precoder changes are periodic: 294 of 1135 located changes at slot 2 mod 40 (the slot-2 CSI-RS, every 20 ms),
      10x uniform; next peaks slots 28/35/30; most common interval 40 slots. -> gNB updates precoders with the
      CSI-RS/report cycle; ~1/4 of the changes are predictable.
  (3) CSI-RS as extra TRS-beam samples (so_light.py USE_CSI=1): per occasion ONE whole-band 8-port weight vector
      reproduces the nearest TRS with 3.4% residual (167/167 occasions). All 25 bursts, NO tablet, found:
        TRS + CSI-RS + DM-RS: slow 79%/56% (-10/-20 dB), fast 49%/38%, faster (7/-9 m/s) 62%/49%; range median
        -0.5 m, RMS 2.6 m, 72% within 1 bin (before: TRS + DM-RS slow 63/44, fast 22/15). TRS + CSI-RS alone: slow
        26-31% (most at a folded speed), fast ~0 (any 20 ms-periodic pattern keeps copies every 2 m/s).
      CSI-RS-derived samples as DM-RS references too (CSI_REF=1): +1-2 targets without tablet, neutral with: optional.
  (4) Light-stream per-snapshot delay alignment of the TRS to each DM-RS snapshot moved origins (-3.5 m bias): now
      off by default (ALIGN=0), bias -1.4..-2.8 m remains (< 1 bin), unexplained.
  (5) Fair test (targets injected into the RAW data, gs_to_nrsr.py INJ=): b9 fast -20 dB: TRS-referenced 14/16,
      raw 3/22, stitching 2/18. Low-rank gap filling (FILL=1) 9/12, fewer maps: NOT useful (fills the static
      channel; a weak target is not recovered where it was never measured). DROPPED.
  (6) Bug fixed: the window-move undo picked its sign by a near-tie jitter metric and flipped once; sign fixed to -1.
  (7) Range-zero on the earliest peak within 6 dB instead of the strongest: no change (the ~-2 m light-stream bias on
      b9 is not the LoS choice; all-burst median -0.5 m).
  (8) Grouping for the TRS-referenced numbers (gs_to_nrsr.py REF_EP): one number per snapshot gives 11 maps/burst
      vs 6 per detected epoch (same snapshots pushed; NR-Scope completes fewer windows), ~2x found fast targets
      (13/22 vs 6/12 at -20 dB on b9), slow similar; 20 ms cycles cut at slot 2 mod 40 in between. Unexplained why.
  (9) Lighter network (so_light.py KEEP=share of the DM-RS slots, no tablet, b9, -15 dB): slow targets 22/32 (100%),
      21/30 (50%), 15/30 (25%), 14/28 (10%), 12/30 (0%, TRS+CSI only, +17 at a folded speed); fast (4/-3 m/s)
      13/32, 4/30, 0, 0, 1; faster (7/-9) 18/32, 14/30, 8/30, 2/30, 0. Slow targets survive a quiet cell; fast ones
      need >= ~50% of the 5 Nov non-tablet load. 4 m/s folds exactly onto 0 m/s with the 20 ms TRS/CSI pattern
      (copies every 2 m/s) and is removed with the clutter: blind speeds without enough DM-RS.
  (10) Speed ambiguity resolved from the range drift (alias_resolve.py on NR-Scope's markers: chain detections across
      maps by range, fit range vs time, range rate = -2 v, pick the 2 m/s copy closest to it; only when the chain
      residual <= 1.5 m and the disagreement > 1 m/s). b9: TRS + CSI-RS only, slow 12/30 -> 28/30; TRS alone 7 -> 23/30;
      10% DM-RS 14 -> 20/28; full DM-RS slow 22 -> 26/32; fast/faster unchanged (without the gating, faster dropped
      18 -> 11: fast targets move metres per map and mis-chain). Fast targets folded onto 0 m/s are removed by the
      clutter stage and cannot be resolved. Constant-velocity targets only so far.
      ALL 25 BURSTS (so/light maps, NR-Scope injection): TRS + CSI-RS only, slow 26% -> 73% (-10 dB), 31% -> 78% (-20);
      with DM-RS: slow 79 -> 86%, 56 -> 63%, but fast/faster lose 2-8 points. Rule: resolve the copies only when the
      DM-RS are scarce (quiet cell); with enough DM-RS the speeds are already right.
  (11) ROBUSTNESS, all bursts, targets injected into the RAW signals (run_robust_all.sh, sum_robust.py; 3 targets:
      25 m +2.5 m/s -20 dB / 60 m -5 m/s -22 dB / 100 m +8 m/s -25 dB), found:
        tablet raw 13% / 29% / 0%; tablet TRS-referenced 38% / 81% / 14% (69-80% within 1 bin, other 1.1/map);
        tablet + TRS + CSI-RS + DM-RS 20% / 72% / 17%; NO tablet TRS + CSI-RS + DM-RS 10% / 44% / 3% (other 0.2/map).
      The -25 dB target near the 120 m map edge is mostly lost; the 25 m +2.5 m/s one folds to +0.5 m/s with the
      20 ms TRS pattern and sits near the direct path. Harder (and more realistic) than the NR-Scope-injected runs.
  (12) NR-Scope map power averaging (sensing.avg_maps 2): b9 14/16 -> 16/16, but on 6 bursts with the 3-target set
      only +-2-5 points (tablet TRS-ref 35/75/15% -> 37/75/13%; no tablet 1/32/0% -> 1/37/0%). Not worth it.
  (13) Why per-epoch TRS-reference numbers give fewer maps: same snapshots, none empty, but NR-Scope builds no maps in
      the first ~0.6 s and last ~0.5 s of b9 (frames 91-159, 238-287): where epoch detection merges snapshots of
      different precoders, the per-epoch number is wrong for some, NR-Scope's alignment silently drops them
      (corr < 0.5, not 64 in a row) and those windows never complete. 6 bursts, 3-target set: per snapshot 73 maps,
      84 found (21/50/13), 63% within 1 bin, 0.5 other/map vs per epoch 60 maps, 75 found, 57%, 1.0/map.
      -> use REF_EP=snap (one number per snapshot and PRG onto the nearest TRS).
  (14) 2026-11-10: NR-Scope trace (env NR_SENSING_TRACE, off by default; nr_ue_sensing.c / nr_ue_map.c): per symbol
      alignment drop / random drop / pushed, and the snapshots of each map. Usage map rebuilt on it (gen_um7.py,
      run_usemap.sh, data /mnt/mydrive/field/usemap7). Of all candidate PRB-symbols, 64% (tablet) / 63% (no tablet) end
      in a map; random_drop (cfg: one symbol in every other slot) discards 24% / 22%, alignment 5% / 4%, too few pilots
      3% / 4%, no TRS within 10 ms 2% / 3%. Fixed in so_light.py MODE=tab: the tablet's record overwrote the blind
      record of the other UEs in the same symbol (now merged).
      random_drop off (cfg_*_nrd.yaml), 26 bursts, 3-target raw injection: tablet TRS-ref 30/74/17% -> 34/75/18%,
      no tablet 10/44/3% -> 9/48/4%. Small gain: the dropped symbol is mostly the second one of the same slot.
      Receiver: BC building 4th floor, gNB on Cubotron ~1 km away: map range is bistatic EXCESS range.
  (15) 2026-11-10 WINDOW-MOVE UNDO WAS WRONG ABOUT DC (fixed, so_csi.undo_wsh, env WSH_Q0=409.25, 0 = old): the phase
      ramp was taken from the band edge, leaving a constant phase PH*shift*409 per antenna at every move; each antenna
      replay moves at its own times, so the INTER-antenna phase jumped (per antenna alone: coherence 1.000 either way;
      4 antennas jointly across a move: 0.80 -> 1.000). All bursts, 3-target raw injection, random_drop off:
        tablet TRS-ref 34/75/18% -> 79/88/29%, <1 bin 72 -> 79%, other/map 0.6 -> 2.2
        no tablet TRS+CSI+DM-RS 9/48/4% -> 25/42/7%, <1 bin 62 -> 84%, other/map 0.3 -> 0.8
        light stream with the tablet (merge bug fixed): 55/79/26%, <1 bin 91%
      Without injection: 0.0-0.1 detections per map, so the extra "other" are replicas of the injected targets (same
      speed 10+ m off 47%, same range other speed 29%, unrelated 22%), not scene false alarms.
  (16) TRS consistency (trs_consist.py, 6 bursts, DC fix): pairs 20 ms apart, coherence of the per-group phase change
      (1 = only the channel moved): TRS median 1.000, p10 1.000, 0% below 0.9; tablet DM-RS 0.967, p10 0.87, 16% below
      0.9; at 40 ms TRS 1.000 / 0% vs DM-RS 0.910 / 45%. Shape similarity 0.99 for both: the channel is the same, the
      DM-RS precoder changes. -> TRS is a fixed reference (paper validation (ii)).
  (17) Narrow DM-RS symbols (MINPRB in so_light.py): 0 / 24 / 48 PRBs minimum. With the tablet 55/79/26 -> 58/83/30%,
      <1 bin 91 -> 94% at 48; without the tablet no detection change, <1 bin 84 -> 91% at 48. Small; 48 helps range.
  (18) TRACKING on real-data maps (track_real.py: ~/GIOVANNI/sim/track.py EKF; stand-in geometry: rx at origin 15 m up,
      gNB 1 km behind the array, boresight -y; positions not surveyed). First result: coverage 12% (tablet TRS-ref)
      because the AoA was random (injected -25/+15/+30 deg measured -3/-8/-1, IQR +-35 deg).
      Cause: each antenna is a separate NR-Scope replay, so the INTER-ANTENNA phase wanders symbol to symbol (TRS:
      3-11 deg median step, p90 9-37, 30-50 deg per burst). Fix ACAL=1 (so_light.py, gs_to_nrsr_ep.py): per snapshot,
      rotate antennas 1-3 onto antenna 0 on the snapshot's static part (direct path): angles become relative to the
      direct path. AoA then -22.8 / +30.6 deg for -25 / +30 (tight); +15 deg comes out mirrored (-15) about half the
      time (open). Detection unchanged (within 1-4 points). Tracking, 26 bursts, 3 targets:
        light with tablet: coverage 13% -> 41%, false tracks/burst 1.8 -> 0.6, GOSPA 7.3 -> 6.0 m
        no tablet:         coverage  5% -> 16%, false tracks/burst 0.4 -> 0.2
        tablet TRS-ref:    coverage 12% -> 21%, false tracks/burst 2.2 -> 1.2 (fewer maps: 11.7 vs 16.2 per burst)
      Live NR-Scope (one process, 4 chains) should not have the wander, but needs a hardware phase calibration anyway.
  (19) CSI-RS window-move undo in so_light.py had its own copy measured from the band edge: after the DC fix it no
      longer matched the TRS and ~all CSI-RS occasions failed the fit (items 15/17/18 light numbers were without CSI-RS).
      Fixed (residual back to 3.4%). MAIN TABLE v2 (all bursts, NR-Scope --inject, random_drop off, MINPRB 48 for light,
      found %, slow/fast/faster at -10/-20 dB):
        tablet raw          9/3   71/36  46/34
        tablet stitching   67/35  71/68  68/68
        tablet TRS-ref     85/81  85/83  84/85   (within 1 bin 80-91%)
        no tablet TRS+DM   80/70  47/39  63/57
        no tablet +CSI-RS  86/78  61/51  73/66   (within 1 bin 90-96%)
        light with tablet  86/81  85/83  86/84   (within 1 bin 95-99%)
      Tracking (ACAL=1, 3-target raw injection): light with tablet coverage 44%, 0.5 false tracks/burst, GOSPA 5.8 m;
      no tablet 29%, 0.3, 5.6 m. Usage map republished (version 8).
  NEXT (proposed): a real recording with a controlled moving target; maps with overlap for the tracker; longer
  windows / better near-LoS clutter suppression for slow and close targets; PDCCH DM-RS dump as extra time samples;
  understand why per-epoch TRS-reference numbers give fewer NR-Scope maps than per-snapshot ones.
- [OPEN] 7d. sens6 rebooted 2026-11-09 09:35 UTC (no clean shutdown; journal lagging minutes before: memory exhaustion,
  very likely from 8 parallel analysis jobs). Run scripts/field_prep.sh before the next recording. Analysis jobs now
  limited by memory (sim_track.py waits for 3 GB free; at most 3 jobs).

## Suggested order

1. Recording with our own second phone downloading next to the tablet + a known moving target (3a, 3b, 5a).
2. Edge matching between UEs (3e), then fold the merge into the maps (4a).

---

## Changelog

- 2026-11-05: file created from the state after the bridge-distance test.
- 2026-11-05: edge matching between UEs tested (3e, 2c -> IN PROGRESS).
- 2026-11-05: 26-burst recording + gs_all.py on all data: 1b, 1c, 2c, 3a, 3c, 3e, 3g done; 3b in progress; 4a deferred; quality gate added; 6d, 7c added.
- 2026-11-05: 4 Nov recording deleted to free disk (7c).
- 2026-11-05: usage map (gs_usemap.py) + LaSen comparison; 6e, 6f added.
- 2026-11-06: OMP tests (frequency: no gain, 2b/3f dropped; time: epoch gluing works, 4c in progress); subframe phase step found (4d).
- 2026-11-06: incoherent-map test (4e); usage-map before/after view; 6e change specified.
- 2026-11-09: 6e implemented (blind DM-RS dump), burst 8 tested; all bursts replaying.
- 2026-11-09: 6e on all bursts, usage map with blind/incoherent colours.
- 2026-11-09: 6f TRS/CSI-RS first results (TRS anchor, timing-window jumps, CSI-RS port span).
- 2026-11-09: usage map with TRS colours (6f).
- 2026-11-09: CSI-RS span test on all UE kinds (6f); usage map toggle for incoherent-only.
- 2026-11-09: simulation with known targets on the real occupancy (5b).
- 2026-11-09: track-level simulation (5c).
- 2026-11-09: NR-Scope chain with injected targets; CSI-RS test; per-PRG precoder direction diversity (5d).
- 2026-11-09: per-antenna replay misalignment found and fixed; TRS-referenced DM-RS works on real data (5d).
- 2026-11-09: all-burst evaluation in NR-Scope's chain (TRS-referenced best; light network: TRS + DM-RS).

## (20) Mirrored / wrong AoA of injected targets (2026-11-10) -- root causes found, fixes tested on b20
Symptom: 25 m @+15 deg target read at ~-17 deg (tablet TRS-ref 60% of matches, light-tab 33%, no-tablet 15%);
60 m @-25 deg read at +60 deg in tabref. Not a true mirror: the estimate was stuck near a fixed wrong angle.
Two independent causes:
 a) Calibration. los_to_zero (so_csi) applied its per-antenna delay ramp from the band EDGE (same bug as undo_wsh):
    constant per-antenna phase 2 pi b0 409/NF. ACAL (one phase per snapshot over the PRBs of that snapshot) hid it,
    but the inter-antenna phase varies by up to ~2 turns across the 98 MHz band (chain 2 vs 0 on TRS), so the
    calibration depended on which PRBs each DM-RS grant covered -> per-map angle offset (sin offset -0.69 = -44 deg
    for the direct path in ~60% of tabref maps). TRS (always full band) unaffected, hence the light stream was OK.
    Fix: los_to_zero ramp about DC (so_csi.py, backup .bak_los0) + GCAL=1 in gs_to_nrsr_ep.py: inter-antenna phase
    per 4-PRB group from the nearest TRS (replaces ACAL). NR-Scope --inject at -15/0/+15: 100% correct (was 33-42%).
 b) TRS-ref per-group number d (common to all antennas), with REF_EP=snap, absorbs the part of a near slow target
    that is equal on all antennas (array mean): angle pushed outward (+15 -> +23), a target at 0 deg vanishes.
    Fix: REF_EP=epoch (d averaged over the precoder epoch): raw 3-target b20 all three 100% correct
    (+14.4 / -25.9 / +27.9); target 0 detected in 10/16 maps vs 14/16 with snap -> rerun all bursts to see the cost.
Why not seen before: angles were random until ACAL (this morning); main tables inject at 0 deg AFTER processing
(sensing_offline --inject) and score range/speed only. Diagnostics: INJ_AT=late, DUMPA=<pkl> (both converters),
mirror.py, mirror2.py, los_angle.py, xphase_f.py, run_mirror*.sh; outputs /mnt/mydrive/field/mirror.
TODO: port GCAL to so_light.py; rerun all 26 bursts (tables + 3-target angles) with los DC fix + GCAL + REF_EP=epoch.
(20, all bursts, run_v3.sh -> /mnt/mydrive/field/main_v3, aoa_v3; scorers sum_v3.py, aoa_v3.py)
 Main table (found, in brackets within 1 bin): GCAL + DC ramps change nothing for detection (trs, light, lightcsi,
 lighttab within 1-3 points of v2). Epoch-averaged tablet reference costs a lot: trsep slow -20 45% (snap 81%),
 lighttabep 39-74% (snap 81-86%). -> keep REF_EP=snap / TAB_EP=snap for detection.
 Angles, 3 raw-injected targets (det = share of maps detected, ok = within 6 deg), before (ACAL) -> after (GCAL):
   tabref   25m det 70->71 ok 2->25 (med +22) | 60m ok 70->98 | 100m ok 57->95
   lighttab 25m det 56->57 ok 5->34 (med +22) | 60m ok 92->100 | 100m ok 88->100
   notab    25m det 31->31 ok 46->71 (med +19) | 60m ok 98->100 | 100m ok 98->100
   trsep 25m det 47 ok 98 (med +14.2); lighttabep 25m det 22 ok 62.
 Remaining: per-snapshot TRS-ref biases the near slow target outward (+15 -> +22), cause (20b). Epoch fixes the angle
 but loses detections. Option: detect on the snap stream, read the angle from the epoch stream (or TRS/CSI only).

## (21) Near-target angle bias of the per-snapshot reference; epoch reference "losing detections" (2026-11-10 night)
 a) Model of (20b) confirmed: one number per 4-PRB group per snapshot, common to all chains, removes kappa*mean_a(target)
    with kappa = beta^2, beta = sin(G x)/(G sin x), x = pi*scs*range/c, G = 48 subcarriers (kappa 0.95 @25 m, 0.76 @60,
    0.44 @100). On the GCAL-calibrated array the common direction is broadside, so the target steering loses its
    broadside projection: predicted +15 -> +23.6 @25 m (measured +22), 0 deg target -27 dB (vanishes), -25/+30 almost
    unchanged (their array mean is small). proj_pred.py (scratch) reproduces it.
 b) Projected-steering AoA in NR-Scope (nr_ue_aoa.c, opt-in env NR_AOA_PROJ_PRB=4, NR_AOA_PROJ_SCS; backup
    ~/field/nr_ue_aoa.c.bak_proj): Bartlett against a(theta) - kappa(range)*mean(a)*1. All bursts (run_proj.sh ->
    /mnt/mydrive/field/aoa_proj, aoa_proj.py): 25 m ok 25->99% trs (med +14.6), 34->99% lighttab, but over-corrects
    60/100 m (trs 98->77%, 95->66%) and the light 25 m (+9..+11). Real absorption is weaker than the model away from
    25 m and in mixed streams. Kept as an option, NOT a default.
 c) Hybrid reference REF_EP=hyb (so_csi.ref_fix_hyb, backups .bak_hyb): d(t,b) = c(t) * D(epoch,b), alternating LS;
    HYB_MAXEP=<s> also cuts epochs. hybs (+ per-snapshot slope) is broken (moves the direct path), do not use.
 d) KEY: NR-Scope's alignment DELAY SEARCH (sensing.align_search_bins: 8) is what made the epoch/hyb references lose
    detections. It searches around its previous delay and can lock on a secondary path: hyb40 on b20 put every
    target 10 m short (4 bins), 69% found / 0% within 1 bin; with align_search_bins: 0 (phase+amplitude alignment per
    symbol kept, no delay search; so/cfg_sb8_noalign.yaml, cfg_trs_noalign.yaml) 81% / 100%. Our converters already
    undo the window jumps (WSH) and put the direct path at zero, so the search is not needed for these streams.
    NR-Scope's per-symbol phase alignment is itself the c(t) of the hybrid (wideband, does not absorb targets).
    b20+b5 (run_hyb4.sh -> /mnt/mydrive/field/hyb4), main slow-10/slow-20/fast-10/fast-20:
      snap  al 79/81/76/81   na 82/77/79/77     epoch al 24/18/70/42   na 87/71/93/93
      hyb   al 38/33/55/40   na 75/61/70/81     hyb40 na 90/82/85/74
    3-target (det / ok): snap al 25m 81/24, 60m 84/100;  epoch na 25m 60/100, 60m 84/100;  hyb na 68/100, 77/100.
    ALL 26 bursts (hyb4all), main slow-10/slow-20/fast-10/fast-20 found% (within 1 bin%), then 3-target det/ok:
      snap  al  85(87) 81(91) 85(85) 83(86) | 25m 71/25  60m 58/98  100m 24/95   <- current default
      snap  na  84(97) 78(97) 84(98) 83(99) | 25m 50/36  60m 76/99  100m 24/96
      epoch na  81(94) 54(98) 88(98) 83(99) | 25m 58/100 60m 83/97  100m 11/100
      hyb   na  77(98) 51(94) 82(96) 76(98) | 25m 60/98  60m 80/96  100m 10/100
      hyb10 na  82(98) 75(98) 85(98) 81(98) | 25m 60/94  60m 78/96  100m 12/96   <- best compromise
      hyb40 na  83(97) 65(98) 85(96) 80(97) | 25m 66/96  60m 80/93  100m  9/100
    -> delay search off improves range precision for every reference (within 1 bin 85-91% -> 94-99%).
    -> REF_EP=hyb HYB_MAXEP=0.01 + align_search_bins 0: detection -0..-6 points vs default, 25 m angle fixed,
       60 m detected +20 points; 100 m target detection halves (24 -> 12%), open. Epoch alone still weak on slow -20.
    Epoch cuts (b20,b5,b9,b14, hyb5): precoder cycle (HYB_CYCLE=1, slot 2 mod 40) ~ same as HYB_MAXEP=0.01; 5 ms
      brings the 25 m bias back (ok 77-88%, med +18..19). Keep HYB_MAXEP=0.01.
    100 m target with hyb10 na: cell power only ~1 dB below snap al (13.8 vs 14.8 dB over map median): at threshold,
      not a separate mechanism. epoch na makes fewer maps (184 vs 234 over the bursts).
    LIGHT streams, all bursts (light4all; ltab = light+tablet, lcsi = no tablet; ep = TAB_EP=epoch, ep2/lcsiep = also
      BLIND_EP=epoch, new so_light option; backup .bak_blindep): delay search off ~neutral (ltab al 86/81/85/82 -> na
      87/82/86/83); epoch references lose a lot with or without it (ltabep na 71/41/73/62, ltabep2 na 61/25/68/53,
      lcsiep na 72/45/39/25), angles right (25 m ok 95-98%) but detected 11-23%. -> TAB_EP/BLIND_EP=hyb ported to
      so_light (backup .bak_hyb), all bursts queued (run_light5_all.sh -> light5all).

## (22) Sensing-only KPI gaps measured (2026-11-10 night; run_kpi.sh -> /mnt/mydrive/field/kpi, kpi_pfa.py)
 All 26 bursts, nothing injected; Pfa is an UPPER bound (real movers count). Clutter suppression = drop of the
 zero-Doppler clutter-to-noise ratio (noise = median of |f| > f_max/2 cells) between clutter removal off and on.
   stream      maps dets/map Pfa/cell  | CNR at LoS off -> on (suppr.) | all range bins off -> on (suppr.)
   lightcsi     413   0.03   1.6e-06   | 38.0 -> 13.8 dB (24.2 dB)     | 25.7 ->  6.1 dB (19.6 dB)
   lighttab     422   0.07   4.0e-06   | 44.8 -> 18.1 dB (26.7 dB)     | 32.4 -> 11.5 dB (20.9 dB)
   trs (snap)   234   0.15   8.5e-06   | 45.1 -> 20.0 dB (25.2 dB)     | 32.7 -> 13.7 dB (19.0 dB)
   trshyb(na)   231   0.12   7.1e-06   | 43.8 -> 21.2 dB (22.6 dB)     | 31.0 -> 13.6 dB (17.4 dB)
 Residual direct-path clutter stays 14-21 dB above noise (the detector guards those bins).
 Other KPI numbers (from main_v3): CPI ~133 ms -> 0.30 m/s (NR-Scope monostatic convention), 7.5 maps/s, range
 resolution c/B = 3.05 m (bins 2.44 m), range error 1.2-2.2 m rms, speed error 0.019 m/s std (synthetic targets).
    LIGHT + hybrid (light5all, all bursts; so_light TAB_EP/BLIND_EP=hyb HYB_MAXEP=0.01; hyb0 = TAB_ALIGN=0):
      ltab     al (default) 86/81/85/82 | 25m det 57 ok 34 (+21.6) | 60m 77/100 | 100m 39/100
      ltabhyb0 na           87/81/86/81 | 25m det 63 ok 87 (+17.2) | 60m 83/100 | 100m 30/99   <- best overall
      ltabhyb2 na (+blind)  84/78/86/80 | 25m det 52 ok 98 (+16.5) | 60m 72/100 | 100m 32/99
      lcsi     al           86/78/61/50 | 25m det 31 ok 71 (+18.9) | 60m 53/100 | 100m 11/100
      lcsihyb  na           83/74/57/47 | 25m det 25 ok 97 (+16.3) | 60m 47/99  | 100m 11/96
    -> hybrid in the light streams costs no detection with the tablet and fixes most of the 25 m bias; with BLIND_EP
       the bias is gone (98%) for ~5 points of 25 m detection. TAB_ALIGN 0/1 equivalent.
    PROPOSED defaults (not applied): light+tablet TAB_EP=hyb HYB_MAXEP=0.01 (+ BLIND_EP=hyb when angles matter),
    TRS-ref REF_EP=hyb HYB_MAXEP=0.01 with align_search_bins: 0. Open: 100 m (-25 dB) detection 39 -> 30%.

(22b) 2026-11-11 cleanup: removed what did not work. AoA projection (NR_AOA_PROJ_PRB/SCS) taken out of nr_ue_aoa.c
(restored from nr_ue_aoa.c.bak_proj, sensing_offline rebuilt); hybrid slope variant (REF_EP=hybs, slope= arg of
ref_fix_hyb) and HYB_CYCLE (same as REF_EP=both) removed; run_proj.sh, aoa_proj.py, patch_aoaproj.py, patch_hyb.py
deleted. ref_fix_hyb output checked bit-identical to before. Backups: *.bak_clean. Data in /mnt/mydrive/field/aoa_proj kept.

(23) 2026-11-11 merging ideas (run_light6_all.sh -> /mnt/mydrive/field/light6all, all 26 bursts, delay search off "na").
  New in so_light.py: TAB_EP/BLIND_EP=slot (one reference number per slot and PRG, both DM-RS symbols pooled);
  BLIND_TRACK=1 (blind record split into UE tracks: each contiguous PRG run joins the track it overlaps on >= 4 PRGs
  with precoder coherence >= 0.9, epochs() test; per-PRG epoch keys, so_csi._epk). ref_fix_hyb vectorised (global
  ALS, c(t) shared by all epochs of a symbol), bit-identical on old inputs.
  Data (b9): runs that overlap an earlier run within 20 ms match at 94-100%, but only 40% of blind runs have one;
  20% of blind PRG-snapshots are in runs < 4 PRGs. Snapshots per reference number, blind: 2.4 -> 2.1 (10 ms cap).
  Main table (slow -10/-20, fast -10/-20): all variants within 1-3 points of their baselines.
  3-target (det / ok within 6 deg, 25 m | 60 m | 100 m):
    ltab (snap)        59/33 | 82/100 | 40/100      ltabslot     60/34 | 83/100 | 40/99   -> slot pooling: no effect
    ltabhyb0 (blind snap) 63/87 | 83/100 | 30/99    ltabhyb2 (blind hyb) 52/98 | 72/100 | 32/99
    ltabhybtrk         61/98 | 82/100 | 28/98   -> blind hyb WITHOUT the detection loss: best tablet variant
    lcsi (no tablet)   32/74 (+18.9) | 59/100 | 11/100    lcsislot 33/73 | 60/100 | 11/100
    lcsihyb            25/97 | 47/99 | 11/96      lcsitrk 33/98 (+16.8) | 58/100 | 10/100 -> best no-tablet variant
    (lcsitrk main table 86/76/59/48 vs lcsi 85/79/61/49: -3 at slow -20 dB)
  Idea 3 (weight snapshots by reference SNR) not done: NR-Scopes alignment rescales each symbol to its reference, so
  a per-symbol weight from the converter is undone; needs an NR-Scope hook. Slot pooling showing no effect suggests
  reference noise is not what limits, so the expected gain is small.
  Proposed defaults: tablet TAB_EP=hyb BLIND_EP=hyb BLIND_TRACK=1 HYB_MAXEP=0.01; no tablet BLIND_EP=hyb BLIND_TRACK=1
  HYB_MAXEP=0.01; both with align_search_bins 0. Remaining cost: 100 m (-25 dB) detection 40 -> 28% (tablet).

(23b) 2026-11-11 committed and pushed (pvnis/NR-Scope goldenSniffer, 1100d33..5955b28): sensing trace + random_drop
off (6a36a3b); scripts (5955b28). Removed before commit: ACAL (replaced by GCAL), FILL/lowrank_fill (not useful),
one-off scripts (cmp_ref, los_angle, mirror2, ant_phase, xphase_f, run_acal*, run_csifix, run_hyb1/2, run_main_v2,
run_noinj, run_nrd). ~/field copies of so_light/so_csi/gs_to_nrsr*/run_v3 match the repo; backups *.bak_drop.
- 2026-11-11: cleanup of ~/field: removed all *.bak_* backups, patch_*.py, and one-off/superseded scripts (mirror*, los_bad, chain_delay, tap_angle, xslope, so_light_keep/los, ACAL/hyb1-2/csifix/noinj/nrd runs, rec4.sh/rec4b.sh/restart_rec/swap_runner). Everything kept is in the repo (5955b28) except recording/replay helpers (rec4c.sh, run_all*.sh, move_rec.sh, replay*.sh).

## (24) MU-MIMO on PCI 632? (2026-11-12) -- NOT USED in the 5 Nov recording
Scripts: ~/field/mumimo.py (on-air DM-RS of decoded grants + whole-carrier survey), ~/field/mumimo3.py (DCI antenna-ports field).
- DCI 1_1 antenna-ports field (Table 7.3.1.2.2-1) located at bit 47/48/49/50 for sizes 54/55/56/57 (4-bit window matched to on-air ports; validated: value 0 -> port 1001 power 0.02x port 1000; value 9 -> comb-1 coherence 0.72; value 10 -> 0.99).
- 124954 DL DCIs (96 RNTIs, 26 bursts, 4 antennas): only values 0 (1 layer p1000, 8%), 2 (2 layers p1000-1001, 84%), 9 (3 layers, 7%), 10 (4 layers, 2%). All SU-MIMO starting at port 1000; never a port allocation that leaves room for a co-scheduled UE (values 1, 3-8, 11: zero). n_SCID bit = 0 in all but one DCI; on-air n_SCID=1 coherence > 0.5 on 0.5% of grants (floor).
- So "both ports" on 98% of grants = rank-2 SU-MIMO of the same UE, not a co-user on port 1001.
- 1167 slots show two decoded grants overlapping >= 4 PRBs, but 722 pairs both claim ports 1000-1001 (would collide): FDRA misreads / false DCIs, not MU-MIMO.
- Survey (every slot, whole carrier): comb-1 DM-RS excess over floor median +2 PRBs, explained by rank-3/4 grants.
- Caveat: users beamed away from our site are not decoded; a MU partner of a decoded UE would still need a compatible port allocation in that UE's DCI, which never happens. Light cell: MU-MIMO may appear under heavy load (controlled test with several phones, busy hour).
- Consequence for us: the converters read port 1000 only (OCC despread cancels 1001, the UE's 2nd layer); ranks 3-4 put port 1002 on comb 1, not read. Fine for sensing (one layer = one precoder column).

## (25) PRS on PCI 632? (2026-11-12) -- NOT FOUND
- SIB1 (decoded, /mnt/mydrive/field/prs/sib1.json): schedules SIB2/4/5 only, no nonCriticalExtension at all -> no posSIBs / posSI-SchedulingInfo broadcast.
- Sequence search on raw IQ (~/field/prs_corr.py; burst b9, antenna 0, 2.4 s): every DL symbol of every frame vs the PRS Gold sequence for all 4096 ids, combs 2/4/6/12 x every offset, chunk-normalised coherence (1 = none). Frame alignment and positive control from the TRS (same formula for id < 1024).
- Everything above noise (99th pct 1.57) is CSI-RS-family (id < 1024): TRS id 632 (slots 5-6, symbols 4/8, metric 2.2-3.0, ~half the frames), a second tracking set id 1022 (same slots/symbols), a weak one id 783 (slot 7, symbols 4/8, 1.97), CSI-RS ids 633/635 (slot 16, symbol 13). Nothing in the PRS-only id range (>= 1024), nothing in other symbols.
- Caveats: one burst (periods <= 1.28 s covered), PRS sequence assumed indexed from point A; a narrow (24-PRB) comb-12 PRS would be diluted below threshold.
- Side finding: two more tracking-like CSI-RS (ids 1022, 783) besides TRS 632 -> other beams; candidate extra references for sensing.

## (26) Recording on sens8 (2026-11-13) -- continuous 4-ch recording WITHOUT gaps works
- sens8 NVMe (Micron 3500 1 TB): direct writes 6.5 GB/s for the first ~310 GB, then falls to ~0.5 GB/s (SLC cache);
  4 ch x 122.88 Msps sc16 needs 1.97 GB/s -> <= ~150 s per recording (~300 GB), pause between recordings, archive to /mnt/mydrive.
- UHD rx_samples_to_file, 4 ch, 120 s straight to NVMe: lost 0.37% (~300 drop events): it writes through the page cache from
  the receive thread, RAM fills in ~25 s, then the writer blocks and the radio packets are dropped.
- New recorder NR-Scope nrscope/scripts/rec_direct/rec_direct.cpp (uncommitted): RAM ring (--ram-gb), one O_DIRECT writer per
  channel, gaps found from radio timestamps -> gaps.csv + zero fill (timing exact), meta.txt, same iq_ch<c>.sc16 files.
  4 ch, 120 s: 14745600000 samples/ch, 0 overflows, 0 seq errors, 0 gaps, ring max 0.2%.
- Signal level today ~17 counts rms vs ~1000 on 5 Nov (-35 dB): antennas on the sens8 X410 not connected / no coverage. Fix before recording.

## 2026-11-23: old scripts archived
One-off analysis scripts cited above (gs_omp, omp_eval3/4, run_light*, sim_*, trs_*, mumimo, ...) are in ~/field/old_scripts_20261123.tar.gz (tar xzf <archive> <name> to restore one).
