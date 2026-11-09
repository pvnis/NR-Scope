# Blind fallback for tracked UEs after RRC reconfiguration: integration plan

Status: plan for review (2026-11-02). No code written for this plan yet. It builds on
what exists on branch `goldenSniffer` (sens6, commits cf888ab…e38f2e8).

## 0. What exists today (starting point)

- **GoldenSniffer discovery** (`nrscope/src/libs/golden_sniffer/`, C): power-clustered PDCCH
  candidates, PDCCH DM-RS → c_init → N_ID by Gold-sequence inversion, GF(2) solve of
  [payload | RNTI] per DCI size. Runs on its own thread pool off the workers' critical
  path. Validated: GoldenSniffer dataset 100% (conf3, conf17); live sens4 cell, phone,
  3236/3236 DCIs (UL 0_x + DL 1_x), 0 false.
- **Limits of what exists:** hard-decision decode only (needs an error-free codeword);
  N_ID only from the DM-RS; no DCI field interpretation; no PDSCH; no sensing hook; no
  per-RNTI state; runs on every slot regardless of whether a UE is tracked.
- **NR-Scope fixes already in** (e38f2e8): a UE attach no longer collapses the pipeline.

## 1. Current decode path, and where a UE is lost (findings)

```
ue_sync (SSB) → worker → SIB1 (sibs_decoder) → RACH/Msg4 (rach_decoder: TC-RNTI, RRCSetup)
   → DCIDecoder::DCIDecoderandReceptionInit   ← RRC config enters HERE, once per run
   → DecodeandParseDCIfromSlot (per known RNTI): srsran_ue_dl_nr_find_dl_dci_nrscope_dciloop
   → srsran_ra_dl_dci_to_grant_nr → nrscope_sensing_process_grant (PDSCH DM-RS)
```

- RRC config enters only in `DCIDecoderandReceptionInit` (`dci_decoder.cc:63-1092`), from
  the **first** RRCSetup of the run (`task_scheduler.cc:321-324` stores masterCellGroup only
  when `!rach_found`) plus the per-radio yaml (`rrc_recfg_config/*.yaml`: mcs tables, a few
  dci_cfg knobs). CORESET 1, search spaces, DCI sizes (srsRAN derives them from `dci_cfg`),
  PDSCH DM-RS type/positions, TDRA, layers all come from there.
- **Nothing after the RRCSetup is ever parsed** (`WorkState::rrc_recfg` exists but is never
  filled; `rach_decoder.cc:536-574` only notices ciphered SRB PDUs).
- **How a UE is lost:** the reconfiguration changes `dci_cfg` inputs (e.g. DCI 1_1 +5 bits
  on the Sunrise cell; sizes 43/47 here), every CRC fails, `last_seen` stops updating
  (`task_scheduler.cc:389,408`), and the UE is dropped after `rnti_expiry_s` (default 5 s)
  together with its PDCCH DM-RS ID (`:426-441`). Seen live today: the phone's RRCSetup is
  decoded, then the regular decoder finds none of its DCIs and drops the RNTI.

## 2. What the sensing pipeline consumes (finding)

- **PDSCH DM-RS only** for the maps (`nrscope_sensing_process_grant`,
  `nrscope_sensing.h:173`), from **DCI 1_1**, same-slot grants (k0 = 0), **≤ 2 layers**
  (`NR_DMRS_SENSING_MAX_PORTS = 2`; rank 3/4 grants contribute nothing), DM-RS type 1,
  length 1. PUSCH DM-RS (`ul_dmrs`) is optional and feeds a covariance CSV only.
  PDCCH DM-RS is not used by sensing.
- Per UE per slot it needs: full-carrier grid per chain, placement
  `{grid_crb0, bwp_start_crb, reference_crb}`, `srsran_sch_cfg_nr_t` (PRBs, S/L, layers,
  CDM groups without data, n_scid, DM-RS type/length/additional positions/typeA_pos,
  scrambling IDs), the DCI `ports` field, PCI, SFN/slot.
- Today N_ID is in effect **always the PCI and nSCID always 0** (srsRAN `ra_dl_dmrs`
  forces `scrambling_id_present = false`; `n_scid` is never set from the DCI). A network
  that configures `dmrs-ScramblingID0/1` is already sensed with wrong pilots.
- **Conclusion: the heavy case.** We must restore DCI decode **and** interpret it into a
  PDSCH grant (PRBs, symbols, layers, ports), **and** recover the PDSCH DM-RS parameters
  (N_ID, nSCID, positions), not just CORESET localization + PDCCH N_ID.

## 3. srsRAN primitives we build on (already linked)

PDCCH: `srsran_pdcch_nr_encode/decode`, `srsran_dmrs_pdcch_*`, `srsran_dci_nr_size/dl_unpack/ul_unpack`;
FEC: `srsran_polar_code_get`, `srsran_polar_encoder_encode`, `srsran_polar_rm_tx`,
`srsran_polar_interleaver_run_u8`, `srsran_polar_chanalloc_tx`, `srsran_crc_*`;
sequences: `srsran_sequence_apply_bit`, `srsran_sequence_state_*`;
RA: `srsran_ra_dl_dci_to_grant_nr`, `srsran_ra_dl_nr_time/freq`, `srsran_dl_nr_mimo_layers`;
PDSCH DM-RS: `srsran_dmrs_sch_get_symbols_idx`, plus NR-Scope's own `dmrs_pilots.c`
(`nrscope_dmrs_c_init`, `nrscope_dmrs_sequence`) and `nr_ue_dmrs_estimate_symbol`.
The GF(2) generator keeps being built from srsRAN's own encoder chain (bit-exact with what
srsRAN transmits); the only new math is the inverses, the error correction and the
Gold-sequence inversion (srsRAN has no inverse).

## 4. Reference algorithms extracted

### 5GDescrambler (Rust; arXiv 2609.07367 v2, Alg. 1-4; repo kit-ps/5GDescrambler)
- **Unknowns** u = [1²⁴ | a (A bits) | 1 | RNTI(16) | N_ID(16)], U = K+33; G = U×E built from
  CRC + RNTI mask + polar + rate matching + Gold-sequence matrix S (x2 linear in c_init).
  **N_ID is solved algebraically** (no DM-RS needed). Two matrix sets: USS (n_RNTI = RNTI)
  and "shared" (n_RNTI = 0: CSS, and USS without pdcch-DMRS-ScramblingID).
- **Repetition cut:** only the first RML = min(E, N) LLRs are used (E = 864/1728 → 512).
- **Configs:** (K,E) pairs, K = A+24, A = 12..140 (AL1 K ≤ 107): 588 pairs × {USS, shared}
  = 1176 decode matrices in the code (the paper says 1416; not reproducible from the code).
  Per config: T (RML×RML column transform), reduced G, free bits. Built at compile time.
- **Soft error correction (Alg. 2):** take the t least reliable positions (t = 10% of E
  for USS, 5% for CSS, AL1 capped at 97−K), treat them as erasures, and solve ONE GF(2)
  system whose equations are the known bits of u (24 CRC prefix ones, the constant, the
  zero padding rows): (24+1+RML−K−32) equations × (t + free) unknowns. Then syndrome check.
  Errors outside the t positions are not corrected (no retry). Costs 1.5-476 µs/candidate
  (AL1-AL16) vs 0.2-2.3 µs without EC; reaches 99% detection at ≈ 5-6 dB.
- **Active groups (Alg. 3):** per PRB per symbol, ≥ 9-11 of 12 REs above a noise cutoff,
  QPSK phase spread σθ < π/3, amplitude CV < 0.3; phase slope + common phase recovery;
  4 rotations tried in the decoder.
- **Incremental config finding (Alg. 4):** learned "PRB set / symbols → DCI sizes, shared?"
  candidates tried every slot; windows → grouped combinations → bounded C(n,k) search;
  fast / super-fast modes; can stop after N frames. Nothing per RNTI.
- **No PDSCH**, no DCI field parsing beyond the UL/DL bit.

### GoldenSniffer (MATLAB) — PDSCH DM-RS side (not ported yet)
- `PDSCH_DMRS_autodetect`: DM-RS symbols from the even/odd subcarrier energy pattern
  ({2}, {2,11}, {2,7,11}, {2,5,8,11} → additional positions 0..3); N_ID and nSCID from
  `DMRS_process` on the even subcarriers + `c_init = 2^17(14s+l+1)(2N+1) + 2N + nSCID`;
  needs ≥ 16 consecutive DM-RS REs (≥ 3 PRBs), else the encrypted underdetermined solver
  (which we already replaced for PDCCH; same idea extends to PDSCH).
- `PDSCH_symbAlloc_autodetect`, `find_FDA_BWP_hypotheses`: BWP and allocation from energy.
- Single layer only.

## 5. Design

### 5.1 Modules
```
nrscope/src/libs/golden_sniffer/        (C core, existing)
  gs_search.c, gs_gold.c                 discovery (existing)
  gs_gf2.c          → extended: 5GDescrambler decode (N_ID-known and N_ID-unknown systems,
                       USS/shared), repetition cut, soft-LLR error correction
  gs_pdsch_dmrs.c   NEW: PDSCH DM-RS detection, N_ID/nSCID recovery, pilot regeneration
  gs_dci_layout.c   NEW: DCI 1_0/0_0 field parsing; 1_1/0_1 layout inference
nrscope/src/libs/blind_track.cc          NEW: per-RNTI tracking state, stale trigger,
                                          targeted decode in the worker, sensing hook
tools/gs_matgen   NEW: offline generator of the per-config matrices → binary file
```

### 5.2 DCI decode (5GDescrambler on top of what exists)
- Keep G built from srsRAN's chain; add the 5GDescrambler variant with N_ID as 16 more
  unknowns (U = K+33) and the repetition cut (RML).
- **Two operating modes per candidate:**
  1. *N_ID known* (from the PDCCH DM-RS inversion, or cached per RNTI): U = K+17 → 16
     more parity bits → stronger error correction and fewer false decodes. Default.
  2. *N_ID unknown* (DM-RS inversion failed, e.g. low SNR or < 72 subcarriers in one
     symbol): 5GDescrambler's U = K+33 system.
- **Soft error correction:** feed equalized REs as LLRs (scaled by the noise estimate);
  Alg. 2 erasure solve over the t least reliable bits; acceptance by syndrome; confidence
  = 5GDescrambler's weight. Replaces today's hard-only decode; the hard path stays as the
  fast first try (EC only when the hard try fails the syndrome).
- **Precomputed matrices:** generated offline by `gs_matgen` (from srsRAN's encoder, so
  bit-exact) into a file loaded at start (≈ 1176 configs; with RML ≤ 512 the T matrices
  are ≤ 32 KB each, ≈ 25-40 MB total). Runtime decode is one bit-matrix × vector + an
  optional small LES: O(1) per candidate.

### 5.3 Per-RNTI tracking and the fallback trigger (`blind_track.cc`)
- State per tracked RNTI (from RACH, or found blind): `{mode: RRC | BLIND, N_ID_pdcch,
  n_rnti flag, CORESET PRB set + symbols (learned), DCI sizes per direction (learned),
  DCI layout id, PDSCH DM-RS {N_ID, nSCID, type, additional_pos}, confidences, last_seen}`.
- **States per RNTI:** RRC (fast path) → ARMED → STALE (blind search for this RNTI, the
  known RNTI used as 16 extra check bits) → CONVERGED (targeted O(1) decode, regular
  decode for this RNTI skipped, search off) → back to STALE if blind decodes fail. The
  RACH/Msg4 decoder runs on every DL slot in every state (new UEs, re-attach → fresh
  RRCSetup → RRC mode).
- **No silence timeout** (an idle or DRX UE gets no DCIs for seconds without anything
  changing). Triggers:
  1. **RRC reconfiguration seen:** the first ciphered SRB1 PDU after the (integrity-only)
     Security Mode Command, from its clear MAC/RLC/PDCP headers (`rach_decoder.cc:536-574`
     already reports these when the gNB schedules them with DCI 1_0 in the CSS; to be
     confirmed live on the sens4 cell) → ARMED; the next regular-decode miss → STALE.
     Every later SRB1 PDU re-arms.
  2. **Evidence from a sweep:** the blind search sees the RNTI in a slot where its regular
     decode found nothing → STALE (backstop when the SRB PDU was missed).
- **Discovery windows** (whole-slot search, all RNTIs): at launch for `startup_s`
  (UEs already connected: tracked in BLIND mode from the start, no RRCSetup prior) and
  periodically (`scan_period_s` / `scan_window_s`): new RNTIs become blind-tracked UEs;
  tracked ones seen blind but missed by the regular path become STALE. `mode: always`
  keeps today's continuous discovery.
- Config:
  ```yaml
  golden_sniffer:
    startup_s: 10
    scan_period_s: 30
    scan_window_s: 2
    trigger_on_srb: true
    mode: fallback   # | always
  ```
- **Reconfiguration with sync inside the cell** (new C-RNTI inside the ciphered message,
  no PRACH): a sweep finds the new RNTI, but up to `scan_period_s` later, and as a new UE.
  Mitigation: when a blind-tracked UE's decodes stop, start a sweep at once (delay ≈ the
  sweep window), and link a new RNTI found there with the same PDCCH N_ID, CORESET
  location and DCI sizes to the vanished one (reliable with a per-UE
  pdcch-DMRS-ScramblingID, a heuristic with N_ID = PCI).
- **Expiry decoupled:** blind sightings refresh `last_seen`, so neither an idle UE nor a
  stale one is dropped while it is still seen.
- **Targeted decode in the worker** (synchronous, cheap): for each BLIND RNTI, decode only
  its learned candidates with known N_ID + size → O(1) per candidate; then grant → PDSCH
  DM-RS → sensing in the same worker slot (the grid is only valid there). Discovery
  (whole-slot search, config finding) stays async on the searcher pool, opportunistic,
  and can be switched off after convergence.

### 5.4 DCI → PDSCH grant
- **1_0 / 0_0** (fallback, in scope of the papers): fully parsed with srsRAN
  (`srsran_dci_nr_dl_unpack` with the fallback size), BWP from SIB1 or learned.
- **1_1 / 0_1** (what this gNB uses in USS: sizes 43/47): out of the papers' scope.
  Approach: start from the RRCSetup-derived `dci_cfg` (known for UEs tracked from RACH)
  and enumerate the few fields a reconfiguration typically changes (BWP count, TDRA rows,
  antenna-port table, SRS/CSI/ZP-CSI/rate-match bits, DAI, PRB bundling) to get candidate
  layouts with exactly the observed size; **validate each layout on the air**: its FDA/TDRA
  must point at PRBs/symbols that carry energy, and the PDSCH DM-RS recovery (5.5) must
  succeed exactly there. Keep the layout that validates; re-check periodically.
  Blind-found UEs with no RRCSetup at all: layout inference from scratch is a research
  item; v1 handles them for 1_0/0_0 only.

### 5.5 PDSCH DM-RS recovery (GoldenSniffer method, extended)
- From the grant: PRBs, symbols, layers/ports, CDM groups. DM-RS symbols confirmed by the
  even/odd energy pattern.
- **N_ID and nSCID** by Gold-sequence inversion on the DM-RS REs of one CDM group
  (`c_init` with the nSCID term), with the joint multi-symbol solver when a symbol has
  < 16 DM-RS (the phone's ping grants are 2 PRBs = 12 DM-RS per symbol; additional DM-RS
  symbols make the joint solve feasible).
- **Multi-layer:** with FD-OCC, ports 1000/1001 share the REs of CDM group 0 with ±1 on
  alternating DM-RS REs, so the "channel" alternates between H0+H1 and H0−H1: split the
  DM-RS REs into the two k′ sub-sequences, unwrap each separately, same c_init. CDM group
  1 (ports 1002/1003) the same, on the odd subcarriers. This makes N_ID recovery work at
  rank 2-4; the sensing estimator itself is limited to 2 ports today (see §7).
- Output: N_ID, nSCID, positions → an `srsran_sch_cfg_nr_t` with scrambling IDs set, and
  regenerated pilots (`dmrs_pilots.c`, `nr_ue_dmrs_seed`).

### 5.6 Hook into sensing
- Call the existing `nrscope_sensing_process_grant(...)` per chain with the blind grant and
  the recovered DM-RS config (placement from SIB1 as today). Needs one change in the
  sensing side: honour `scrambling_id0/1` + `n_scid` instead of the forced PCI/0.
- Blind-mode estimates are tagged (source = BLIND, confidence) so the maps can weight or
  exclude them.

## 6. Test strategy

- **(a) generator + inverse:** for every config, random u → srsRAN encode → decode
  recovers u (both N_ID modes, USS/shared, with repetition cut); rank/nullity matches
  5GDescrambler's table.
- **EC:** BLER vs SNR on synthetic AWGN slots for hard vs soft decode; target ≈ the
  paper's 99% at 5-6 dB; false-positive rate with blind K.
- **(b) DM-RS → c_init → N_ID:** PDCCH (done) and PDSCH: synthetic grids with srsRAN's
  `srsran_dmrs_sch_put_sf` (types, positions, nSCID, 1-4 layers), recovered vs truth.
- **(c) end-to-end vs ground truth, in this order:**
  1. *Known config:* GoldenSniffer dataset + 5GDescrambler dataset (Zenodo 22659196),
     scored against their gNB logs.
  2. *Our testbed, pre-security:* a UE tracked from RACH with its RRCSetup known: the
     blind path must give the same RNTI, PDCCH N_ID, CORESET location, DCI sizes; grants
     must match the gNB MAC log (`rb=[a..b)`, k1/k2, ri); regenerated PDSCH DM-RS must
     match the legitimate receiver's pilots (per-RE equality, and LS-estimate correlation
     with the regular path on the same grants).
  3. *Post-reconfiguration:* ground truth from the gNB itself: we own the gNB, so its RRC
     messages can be logged in clear (rrc_level debug / RRC pcap) and the MAC log gives
     every grant.
- **(d) sensing hook:** recorded grids (GS_DUMP) replayed through the hook; maps from the
  blind path vs the regular path on the same pre-security segment.
- Robustness: replay with added AWGN (gs_grid_replay + noise) and the AWGN sweep of the
  5GDescrambler evaluation harness.

## 7. Scope violations on our real traffic (flagged)

1. **PDSCH rank 4:** 99.5% of the phone's grants on the sens4 cell are `ri=4`. Both
   papers assume single layer, and the sensing estimator takes ≤ 2 ports. N_ID recovery
   can be extended (5.5); sensing at rank 3-4 is a separate sensing change. For validation
   the gNB can run with `max_rank: 1` (your call: it changes your gNB config).
2. **Non-fallback DCIs (1_1/0_1)** in the USS: outside both papers; needs the layout
   inference of 5.4.
3. **Interleaved CORESET#0, CORESET#0 at point A, CORESET 1 off the 6-RB grid:** handled by
   the existing discovery; 5GDescrambler assumes contiguous non-interleaved (its
   combination search covers interleaving at a cost).
4. **Short PDSCH allocations** (2 PRBs = 12 DM-RS per symbol < 16): need the joint
   multi-symbol solve.
5. **PDSCH DM-RS scrambling IDs:** today's sensing ignores them (always PCI, nSCID 0);
   real networks that set them are sensed with wrong pilots even on the RRC path.
6. **AL1 with large DCIs (K > 76) is ambiguous** (5GDescrambler): decodes there need the
   learned size or are rejected.

## 8. Phases

| Phase | Content | Tests |
|---|---|---|
| P1 | 5GDescrambler decode core in C: N_ID-unknown system, repetition cut, soft EC, offline matrix file | (a), EC BLER |
| P2 | Per-RNTI tracking, stale trigger, targeted decode in the worker, expiry refresh | (c)1-2 for PDCCH |
| P3 | DCI parsing: 1_0/0_0 full; 1_1/0_1 layout inference validated on air | grants vs MAC log |
| P4 | PDSCH DM-RS recovery (N_ID, nSCID, positions, CDM split) + pilot regeneration | (b) PDSCH, pilots vs legit |
| P5 | Sensing hook + scrambling-ID support in sensing | (d), live |

## 9. Decisions (2026-11-02)

1. **Rank: the gNB stays as configured (rank up to 4).** PDSCH N_ID recovery handles
   rank 2-4 (CDM-group split, 5.5). Sensing estimator: extend from 2 to 4 ports (two CDM
   groups × FD-OCC pair, independent) as a sensing-side item in P5.
2. **Ground truth: plaintext RRC logging on the sens4 gNB is approved.** Needs one gNB
   restart (rrc_level debug / RRC pcap) at the start of the post-reconfiguration tests;
   check RU/PTP after it (LiteOn PTP glitch).
3. **Blind search per UE only after it goes stale, triggered by the RRC reconfiguration
   (ciphered SRB1 PDU) or by sweep evidence, never by a silence timeout**; plus discovery
   at launch (`startup_s`) and periodic sweeps (`scan_period_s`/`scan_window_s`) for UEs
   already connected or missed (§5.3). `mode: always` keeps continuous discovery.
4. **DCI interpretation without the RRC config: option b (confirmed)** — partial parse of 1_1/0_1 (frequency and time allocation, the front
   fields, few offset/BWP hypotheses; the RRCSetup layout as prior for stale tracked
   UEs) + layers, ports, N_ID, nSCID measured on the PDSCH DM-RS; hypotheses confirmed
   on air (PDSCH energy at the decoded PRBs, DM-RS recovered there, consistency over
   DCIs). 1_0/0_0 parsed exactly. Alternatives considered: fallback formats only (too
   little traffic after RRC setup), full layout inference (large ambiguous space, mostly
   irrelevant fields), energy only (ambiguous with several UEs per slot).
