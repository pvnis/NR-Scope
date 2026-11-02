/*
 * SPDX-License-Identifier: MIT
 *
 * Self-test of the GoldenSniffer port (nrscope/src/libs/golden_sniffer), built as
 * gs_selftest. No radio, no gNB: srsRAN's own PDCCH encoder and DM-RS mapper put
 * DCIs into a synthetic slot, and the blind search must find them without being
 * told the RNTI, the scrambling ID, the size, the AL or where they are.
 *
 *  1. DM-RS inversion (DMRS_process.m port): N_ID from one symbol of DM-RS.
 *  2. Underdetermined DM-RS solver (replacement of the encrypted .p): N_ID from
 *     2 or 3 short symbols, checked against an exhaustive search.
 *  3. GF(2) DCI decoder: RNTI and payload from noiseless hard bits of
 *     srsran_pdcch_nr_encode, for every AL and both scrambling cases.
 *  4. Blind search of whole slots: non-interleaved UE CORESETs (point A DM-RS
 *     reference, pdcch-DMRS-ScramblingID, n_RNTI = C-RNTI), CORESET#0 (own DM-RS
 *     reference, PCI), interleaved CORESETs; a timing slope, a symbol-to-symbol phase
 *     step and AWGN on every RE.
 *  5. Tracked UE: direct decode of one UE's DCIs from what was learned of it.
 *  6. Soft-information error correction (gs_ec.c): detection vs SNR, hard vs soft;
 *     false alarms on noise and on other UEs' DCIs; time per call.
 *
 * Pass: 1-3 exact; 4 finds >= 95% of the DCIs with no wrong DCI reported, except
 * CORESET#0 with unrestricted sizes (n_RNTI = 0, nested polar codes): <= 1% wrong.
 * 6: see test_soft(). "gs_selftest <seed> <n>" runs test n only.
 */
#include "nrscope/hdr/golden_sniffer/gs_core.h"
#include "nrscope/src/libs/golden_sniffer/gs_ec.h"
#include "srsran/phy/ch_estimation/dmrs_pdcch.h"
#include "srsran/phy/phch/pdcch_nr.h"
#include "srsran/phy/utils/random.h"
#include "srsran/phy/utils/vector.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static srsran_random_t rnd;
static int             failures = 0;

#define CHECK(cond, ...)                                                                                               \
  do {                                                                                                                 \
    if (!(cond)) {                                                                                                     \
      printf("FAIL %s:%d: ", __FILE__, __LINE__);                                                                      \
      printf(__VA_ARGS__);                                                                                             \
      printf("\n");                                                                                                    \
      failures++;                                                                                                      \
    }                                                                                                                  \
  } while (0)

static int rint_(int lo, int hi)
{
  return srsran_random_uniform_int_dist(rnd, lo, hi);
}

static double now_s(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

/* DM-RS of one symbol as the receiver sees it after unwrap: on the axes, one common
 * unknown j^k */
static void dmrs_symbols(uint32_t cinit, uint32_t D, uint32_t n, uint32_t k, cf_t* out)
{
  uint8_t c[2048];
  gs_gold_bits(cinit, D, 2 * n, c);
  cf_t amb = cexpf(-I * (float)M_PI / 4.0f) * cpowf(I, (float)k);
  for (uint32_t m = 0; m < n; m++) {
    cf_t d = ((1.0f - 2.0f * c[2 * m]) + I * (1.0f - 2.0f * c[2 * m + 1])) * (float)M_SQRT1_2;
    out[m] = d * amb;
  }
}

static void test_dmrs_direct(void)
{
  int ok = 0, trials = 500;
  for (int t = 0; t < trials; t++) {
    uint32_t n_id = rint_(0, 65535), slot = rint_(0, 19), l = rint_(0, 2), D = 6 * rint_(0, 300);
    uint32_t n = 18 * (1 << rint_(0, 3)), k = rint_(0, 3);
    cf_t     c[18 * 8];
    dmrs_symbols(gs_pdcch_dmrs_cinit(slot, l, n_id), D, n, k, c);
    uint32_t ci[4], rot[4], got = UINT32_MAX;
    int      nf = gs_dmrs_process(c, n, D, ci, rot);
    for (int i = 0; i < nf; i++) {
      uint32_t id;
      if (gs_pdcch_cinit_check(ci[i], slot, l, &id)) {
        got = id;
        break;
      }
    }
    ok += got == n_id;
  }
  CHECK(ok == trials, "DM-RS direct: %d/%d", ok, trials);
  printf("1. DM-RS direct inversion: %d/%d\n", ok, trials);
}

static void test_dmrs_underdet(void)
{
  /* 3 symbols of 2 RBs (n = 6) give 12 bits each: some N_IDs are truly ambiguous
   * (bits 13..15 of N can vanish from so short a window), and the solver must then
   * refuse rather than guess. So: never a wrong N_ID, and nearly always one. */
  int    ok = 0, wrong = 0, trials = 200;
  double t0 = now_s();
  for (int t = 0; t < trials; t++) {
    uint32_t dur = rint_(2, 3), n = dur == 2 ? 9 : (rint_(0, 1) ? 6 : 12);
    uint32_t n_id = rint_(0, 65535), slot = rint_(0, 19), D = 6 * rint_(0, 300), k = rint_(0, 3);
    cf_t     c[3 * 32];
    for (uint32_t l = 0; l < dur; l++) {
      dmrs_symbols(gs_pdcch_dmrs_cinit(slot, l, n_id), D, n, k, &c[l * n]);
    }
    uint32_t id, rot;
    bool     found = gs_dmrs_solve_underdet(c, n, dur, D, slot, &id, &rot);
    ok += found && id == n_id;
    wrong += found && id != n_id;
  }
  double dt = (now_s() - t0) / trials;
  CHECK(wrong == 0, "underdetermined DM-RS: %d wrong N_IDs", wrong);
  CHECK(ok >= trials * 97 / 100, "underdetermined DM-RS: %d/%d", ok, trials);
  printf("2. Underdetermined DM-RS solver: %d/%d, %d wrong (%.2f ms each)\n", ok, trials, wrong, dt * 1e3);

  /* exhaustive cross-check on a few: every N_ID that explains the bits */
  for (int t = 0; t < 3; t++) {
    uint32_t dur = 2, n = 9, n_id = rint_(0, 65535), slot = rint_(0, 19), D = 6 * rint_(0, 100);
    cf_t     c[2 * 9];
    for (uint32_t l = 0; l < dur; l++) {
      dmrs_symbols(gs_pdcch_dmrs_cinit(slot, l, n_id), D, n, 0, &c[l * n]);
    }
    uint32_t nof = 0;
    for (uint32_t cand = 0; cand < 65536; cand++) {
      for (uint32_t k = 0; k < 4; k++) {
        bool all = true;
        for (uint32_t l = 0; l < dur && all; l++) {
          cf_t r[9];
          dmrs_symbols(gs_pdcch_dmrs_cinit(slot, l, cand), D, n, k, r);
          for (uint32_t m = 0; m < n && all; m++) {
            all = cabsf(r[m] - c[l * n + m]) < 1e-3f;
          }
        }
        nof += all;
      }
    }
    uint32_t id, rot;
    bool     found = gs_dmrs_solve_underdet(c, n, dur, D, slot, &id, &rot);
    CHECK((nof == 1) == (found && id == n_id), "exhaustive: %u solutions, solver %s", nof, found ? "found" : "none");
  }
}

typedef struct {
  srsran_carrier_nr_t carrier;
  srsran_coreset_t    coreset;
  srsran_pdcch_nr_t   pdcch;
  cf_t*               grid;
} tx_t;

static void tx_init(tx_t* tx, uint32_t nof_prb, uint32_t pci)
{
  memset(tx, 0, sizeof(*tx));
  tx->carrier.pci             = pci;
  tx->carrier.nof_prb         = nof_prb;
  tx->carrier.scs             = srsran_subcarrier_spacing_30kHz;
  tx->carrier.max_mimo_layers = 1;
  srsran_pdcch_nr_args_t args = {};
  srsran_pdcch_nr_init_tx(&tx->pdcch, &args);
  tx->grid = srsran_vec_cf_malloc(nof_prb * 12 * 14);
}

static void tx_free(tx_t* tx)
{
  srsran_pdcch_nr_free(&tx->pdcch);
  free(tx->grid);
}

static int tx_put(tx_t* tx, uint32_t slot, uint32_t Ll, uint32_t ncce, srsran_search_space_type_t ss, uint16_t rnti,
                  uint32_t nbits, uint8_t* payload)
{
  srsran_pdcch_nr_set_carrier(&tx->pdcch, &tx->carrier, &tx->coreset);
  srsran_dci_msg_nr_t msg = {};
  msg.ctx.location.L      = Ll;
  msg.ctx.location.ncce   = ncce;
  msg.ctx.ss_type         = ss;
  msg.ctx.rnti            = rnti;
  msg.ctx.coreset_id      = tx->coreset.id;
  msg.nof_bits            = nbits;
  srsran_random_bit_vector(rnd, msg.payload, nbits);
  memcpy(payload, msg.payload, nbits);
  if (srsran_pdcch_nr_encode(&tx->pdcch, &msg, tx->grid) < SRSRAN_SUCCESS) {
    return -1;
  }
  srsran_slot_cfg_t slot_cfg = {.idx = slot};
  return srsran_dmrs_pdcch_put(&tx->carrier, &tx->coreset, &slot_cfg, &msg.ctx.location, tx->grid);
}

/* GF(2) decoder on the noiseless REs of one candidate, read back in mapping order */
static void test_gf2(void)
{
  tx_t tx;
  int  ok = 0, trials = 0;
  tx_init(&tx, 106, 0);
  for (int t = 0; t < 300; t++) {
    uint32_t Ll = rint_(0, 4), nbits = rint_(20, 80), pci = rint_(0, 1007);
    if (nbits + 24 > (108U << Ll)) {
      continue;
    }
    bool ue                       = rint_(0, 1);
    tx.carrier.pci                = pci;
    memset(&tx.coreset, 0, sizeof(tx.coreset));
    tx.coreset.id                 = 1;
    tx.coreset.duration           = 1;
    tx.coreset.mapping_type       = srsran_coreset_mapping_type_non_interleaved;
    tx.coreset.dmrs_scrambling_id_present = ue;
    tx.coreset.dmrs_scrambling_id = rint_(0, 65535);
    for (int r = 0; r < 16; r++) {
      tx.coreset.freq_resources[r] = true; /* 96 RBs = 16 CCEs */
    }
    uint16_t rnti = (uint16_t)rint_(1, 65535);
    uint8_t  payload[128], got[GS_MAX_DCI_BITS];
    srsran_vec_cf_zero(tx.grid, 106 * 12 * 14);
    tx_put(&tx, 0, Ll, 0, ue ? srsran_search_space_type_ue : srsran_search_space_type_common_3, rnti, nbits, payload);
    uint8_t  bits[108 << 4];
    uint32_t E = 0;
    for (uint32_t k = 0; k < (72U << Ll); k++) {
      if (k % 4 == 1) {
        continue;
      }
      bits[E++] = crealf(tx.grid[k]) < 0;
      bits[E++] = cimagf(tx.grid[k]) < 0;
    }
    uint16_t r_got = 0;
    bool     found = gs_dci_try(bits, E, nbits, ue ? tx.coreset.dmrs_scrambling_id : pci, ue, got, &r_got);
    trials++;
    if (found && r_got == rnti && memcmp(got, payload, nbits) == 0) {
      ok++;
    } else {
      printf("   GF(2) miss: AL%u size %u %s found=%d rnti %04x/%04x\n", 1U << Ll, nbits, ue ? "UE" : "CSS", found,
             r_got, rnti);
    }
  }
  tx_free(&tx);
  CHECK(ok == trials, "GF(2) decoder: %d/%d", ok, trials);
  printf("3. GF(2) DCI decoder: %d/%d\n", ok, trials);
}

typedef struct {
  uint16_t rnti;
  uint32_t n_id, L, nbits, ncce;
  uint8_t  payload[128];
  bool     found;
} truth_t;

/* One random slot; returns the number of DCIs placed */
static int make_slot(tx_t* tx, int scenario, uint32_t slot, truth_t* truth, gs_slot_t* s, float snr_db)
{
  const uint32_t nof_prb = tx->carrier.nof_prb, nsc = nof_prb * 12;
  uint32_t       pci     = rint_(0, 1007);
  tx->carrier.pci        = pci;
  memset(&tx->coreset, 0, sizeof(tx->coreset));
  srsran_vec_cf_zero(tx->grid, nsc * 14);

  srsran_search_space_type_t ss;
  uint32_t                   n_freq_res, first_res = 0, dur;
  s->ref_crb0 = 0;
  s->ref_cs0  = GS_NO_REF;
  if (scenario == 0) { /* UE CORESET, non-interleaved, point A reference */
    dur                                    = rint_(1, 2);
    first_res                              = rint_(0, 2);
    n_freq_res                             = 5;
    tx->coreset.id                         = 1;
    tx->coreset.dmrs_scrambling_id_present = true;
    tx->coreset.dmrs_scrambling_id         = rint_(0, 65535);
    ss                                     = srsran_search_space_type_ue;
  } else if (scenario == 1 || scenario == 3) { /* CORESET#0: DM-RS from its first RB, PCI, common */
    dur                   = rint_(1, 2);
    n_freq_res            = 4;
    tx->coreset.id        = 0;
    tx->coreset.offset_rb = rint_(1, 10);
    s->ref_cs0            = (int32_t)tx->coreset.offset_rb * 12;
    ss                    = srsran_search_space_type_common_0;
  } else { /* interleaved UE CORESET, REG bundles of 6, R = 2 */
    dur                                    = 1;
    n_freq_res                             = 8;
    tx->coreset.id                         = 2;
    tx->coreset.mapping_type               = srsran_coreset_mapping_type_interleaved;
    tx->coreset.reg_bundle_size            = srsran_coreset_bundle_size_n6;
    tx->coreset.interleaver_size           = srsran_coreset_bundle_size_n2;
    tx->coreset.shift_index                = 0;
    tx->coreset.dmrs_scrambling_id_present = true;
    tx->coreset.dmrs_scrambling_id         = rint_(0, 65535);
    ss                                     = srsran_search_space_type_ue;
  }
  tx->coreset.duration = dur;
  for (uint32_t r = first_res; r < first_res + n_freq_res; r++) {
    tx->coreset.freq_resources[r] = true;
  }
  const uint32_t n_cce = n_freq_res * 6 * dur / 6;

  /* DCIs on free CCEs, each at its own power */
  bool     busy[64] = {};
  int      nof      = 0;
  uint32_t want     = rint_(1, 3);
  for (uint32_t d = 0; d < want; d++) {
    uint32_t Ll = rint_(0, scenario == 2 ? 2 : 3);
    if (scenario == 2 && Ll == 0) {
      Ll = 1; /* AL1 of a bundle-6 interleaved CORESET is one contiguous REG bundle */
    }
    uint32_t L = 1U << Ll;
    if (L > n_cce) {
      continue;
    }
    uint32_t ncce = L * rint_(0, n_cce / L - 1);
    bool     free_ = true;
    for (uint32_t j = ncce; j < ncce + L; j++) {
      free_ &= !busy[j];
    }
    if (!free_) {
      continue;
    }
    /* a cell uses a handful of sizes; the CORESET#0 scenario draws from three, so the
     * size most often seen settles which of several solving sizes a weak decode is */
    static const uint32_t cs0_sizes[3] = {29, 37, 41};
    uint32_t nbits = scenario == 3 ? 41 : scenario == 1 ? cs0_sizes[rint_(0, 2)] : rint_(20, 70);
    if (nbits + 16 + 8 > 108U * L) {
      continue;
    }
    uint16_t rnti = ss == srsran_search_space_type_common_0 ? (rint_(0, 1) ? 0xffff : rint_(1, 0x4600))
                                                            : (uint16_t)rint_(0x4601, 0xfff0);
    if (tx_put(tx, slot, Ll, ncce, ss, rnti, nbits, truth[nof].payload) < 0) {
      continue;
    }
    for (uint32_t j = ncce; j < ncce + L; j++) {
      busy[j] = true;
    }
    truth[nof].rnti  = rnti;
    truth[nof].n_id  = tx->coreset.dmrs_scrambling_id_present ? tx->coreset.dmrs_scrambling_id : pci;
    truth[nof].L     = L;
    truth[nof].nbits = nbits;
    truth[nof].ncce  = ncce;
    truth[nof].found = false;
    nof++;
  }

  /* PDSCH-like QPSK after the CORESET so the duration detector sees a change */
  for (uint32_t l = dur; l < 3; l++) {
    for (uint32_t k = 0; k < nsc; k++) {
      tx->grid[l * nsc + k] = ((rint_(0, 1) ? 1.0f : -1.0f) + I * (rint_(0, 1) ? 1.0f : -1.0f)) * (float)M_SQRT1_2 *
                              (0.3f + 0.7f * srsran_random_uniform_real_dist(rnd, 0, 1));
    }
  }

  /* channel: timing slope, phase step between symbols, a mild frequency ripple, noise */
  float tau    = srsran_random_uniform_real_dist(rnd, -4.0f, 4.0f) / 1024.0f; /* samples / FFT size */
  float ph0    = srsran_random_uniform_real_dist(rnd, -M_PI, M_PI);
  float dph    = srsran_random_uniform_real_dist(rnd, -0.3f, 0.3f);
  float sigma  = powf(10.0f, -snr_db / 20.0f) * (float)M_SQRT1_2;
  float ripple = srsran_random_uniform_real_dist(rnd, 0.0f, 0.15f);
  for (uint32_t l = 0; l < 3; l++) {
    for (uint32_t k = 0; k < nsc; k++) {
      float mag = 1.0f + ripple * sinf(2.0f * (float)M_PI * k / 300.0f);
      cf_t  h   = mag * cexpf(I * (ph0 + dph * l + 2.0f * (float)M_PI * tau * k));
      cf_t  n   = sigma * (srsran_random_gauss_dist(rnd, 1.0f) + I * srsran_random_gauss_dist(rnd, 1.0f));
      tx->grid[l * nsc + k] = tx->grid[l * nsc + k] * h + n;
    }
  }

  s->grid     = tx->grid;
  s->nof_sc   = nsc;
  s->slot_idx = slot;
  s->pci      = pci;
  return nof;
}

static void test_search(void)
{
  const char* names[4] = {"UE CORESET (point A ref)", "CORESET#0 (own ref, PCI)", "interleaved UE CORESET",
                          "CORESET#0, CSS size known"};
  tx_t        tx;
  tx_init(&tx, 51, 0);
  for (int scenario = 0; scenario < 4; scenario++) {
    gs_args_t args;
    gs_args_default(&args);
    args.max_duration       = 2;
    const uint32_t css41[1] = {41};
    if (scenario == 3) {
      args.css_sizes     = css41;
      args.nof_css_sizes = 1;
    }
    gs_t* gs = gs_create(&args);
    CHECK(gs != NULL, "gs_create");
    if (gs == NULL) {
      return;
    }
    int    placed = 0, found = 0, wrong = 0, slots = 300;
    double t_tot  = 0;
    for (int t = 0; t < slots; t++) {
      truth_t   truth[8];
      gs_slot_t s;
      uint32_t  slot = rint_(0, 19);
      int       nof  = make_slot(&tx, scenario, slot, truth, &s, 25.0f);
      placed += nof;
      gs_dci_t out[GS_MAX_DCI_PER_SLOT];
      double   t0 = now_s();
      int      n  = gs_search_slot(gs, &s, out, GS_MAX_DCI_PER_SLOT);
      t_tot += now_s() - t0;
      for (int i = 0; i < n; i++) {
        bool match = false;
        for (int j = 0; j < nof; j++) {
          if (!truth[j].found && out[i].rnti == truth[j].rnti && out[i].nof_bits == truth[j].nbits &&
              out[i].n_id == truth[j].n_id && out[i].L == truth[j].L &&
              memcmp(out[i].payload, truth[j].payload, truth[j].nbits) == 0) {
            truth[j].found = true;
            match          = true;
            break;
          }
        }
        if (match) {
          found++;
        } else {
          wrong++;
          printf("   wrong DCI: rnti %04x size %u AL%u n_id %u n_rnti %d cs0 %d iv %u:[%d,%d] slot pci %u\n", out[i].rnti,
                 out[i].nof_bits, out[i].L, out[i].n_id, out[i].n_rnti, out[i].coreset0, out[i].cand.nof_iv,
                 out[i].cand.iv[0].first, out[i].cand.iv[0].last, s.pci);
          for (int j = 0; j < nof; j++) {
            printf("      truth: rnti %04x size %u AL%u n_id %u found %d\n", truth[j].rnti, truth[j].nbits, truth[j].L,
                   truth[j].n_id, truth[j].found);
          }
        }
        if (scenario == 1 || scenario == 3) {
          CHECK(out[i].coreset0, "CORESET#0 DCI not flagged as such");
        }
      }
    }
    printf("4.%d %-28s found %d/%d DCIs, %d wrong, %.2f ms/slot\n", scenario, names[scenario], found, placed, wrong,
           1e3 * t_tot / slots);
    CHECK(found >= placed * 95 / 100, "%s: found %d/%d", names[scenario], found, placed);
    /* n_RNTI = 0 decodes with unrestricted sizes keep ~0.4% at a wrong size */
    CHECK(scenario == 1 ? wrong <= placed / 100 : wrong == 0, "%s: %d wrong DCIs", names[scenario], wrong);
    gs_destroy(gs);
  }
  tx_free(&tx);
}


/* 5. Tracked UE: its learned PDCCH (CORESET extent, N_ID, sizes, ALs) decodes its DCIs
 * directly; another UE's DCI in the same CORESET must not be reported. */
static void test_tracked(void)
{
  gs_args_t args;
  gs_args_default(&args);
  gs_t* gs = gs_create(&args);
  tx_t  tx;
  tx_init(&tx, 51, 0);
  const uint16_t rnti = 0x4613, other = 0x4620;
  const uint32_t n_id = 1234, sizes[2] = {43, 47};
  int            placed = 0, found = 0, wrong = 0, others = 0;
  double         t_tot  = 0;
  for (int t = 0; t < 300; t++) {
    uint32_t dur  = 1 + (t & 1), slot = rint_(0, 19), pci = rint_(0, 1007);
    tx.carrier.pci = pci;
    memset(&tx.coreset, 0, sizeof(tx.coreset));
    tx.coreset.id                         = 1;
    tx.coreset.duration                   = dur;
    tx.coreset.dmrs_scrambling_id_present = true;
    tx.coreset.dmrs_scrambling_id         = n_id;
    for (uint32_t r = 1; r < 6; r++) {
      tx.coreset.freq_resources[r] = true; /* RBs 6..35 */
    }
    srsran_vec_cf_zero(tx.grid, 51 * 12 * 14);
    const uint32_t n_cce = 5 * dur;
    uint32_t       Ll    = rint_(0, dur == 1 ? 1 : 2), L = 1U << Ll;
    uint32_t       ncce  = L * rint_(0, n_cce / L - 1);
    uint32_t       nb    = sizes[rint_(0, 1)];
    uint8_t        payload[128], other_payload[128];
    tx_put(&tx, slot, Ll, ncce, srsran_search_space_type_ue, rnti, nb, payload);
    placed++;
    /* another UE on the next free CCE, if any */
    uint32_t oc = ncce + L;
    if (oc < n_cce) {
      tx_put(&tx, slot, 0, oc, srsran_search_space_type_ue, other, 43, other_payload);
    }
    float snr   = getenv("GS_TEST_SNR") ? (float)atof(getenv("GS_TEST_SNR")) : 20.0f;
    float sigma = powf(10.0f, -snr / 20.0f) * (float)M_SQRT1_2;
    float ph = srsran_random_uniform_real_dist(rnd, -M_PI, M_PI), tau = srsran_random_uniform_real_dist(rnd, -3, 3) / 1024.0f;
    for (uint32_t l = 0; l < 3; l++) {
      for (uint32_t k = 0; k < 51 * 12; k++) {
        cf_t h = cexpf(I * (ph + 0.2f * l + 2.0f * (float)M_PI * tau * k));
        cf_t n = sigma * (srsran_random_gauss_dist(rnd, 1.0f) + I * srsran_random_gauss_dist(rnd, 1.0f));
        tx.grid[l * 51 * 12 + k] = tx.grid[l * 51 * 12 + k] * h + n;
      }
    }
    gs_slot_t sl = {.grid = tx.grid, .nof_sc = 51 * 12, .slot_idx = slot, .ref_crb0 = 0, .ref_cs0 = GS_NO_REF, .pci = pci};
    gs_track_t tr = {.rnti = rnti, .n_id = n_id, .n_rnti = true, .coreset0 = false, .sc_lo = 6 * 12,
                     .sc_hi = 36 * 12 - 1, .duration = dur, .al_mask = 0x1f, .sizes = {43, 47}, .nof_sizes = 2};
    gs_dci_t out[8];
    double   t0 = now_s();
    int      n  = gs_decode_tracked(gs, &sl, &tr, out, 8);
    t_tot += now_s() - t0;
    int hit = 0;
    for (int i = 0; i < n; i++) {
      hit |= out[i].rnti == rnti;
    }
    if (!hit && getenv("GS_TEST_DEBUG")) {
      printf("   miss: dur %u L %u ncce %u size %u other %d\n", dur, L, ncce, nb, oc < n_cce);
    }
    for (int i = 0; i < n; i++) {
      if (out[i].rnti != rnti) {
        others++;
      } else if (out[i].nof_bits == nb && out[i].L == L && memcmp(out[i].payload, payload, nb) == 0) {
        found++;
      } else {
        wrong++;
      }
    }
  }
  printf("5. Tracked UE decode at %s dB: %d/%d, %d wrong, %d of another UE (%.3f ms/slot)\n",
         getenv("GS_TEST_SNR") ? getenv("GS_TEST_SNR") : "20", found, placed, wrong, others, 1e3 * t_tot / 300);
  CHECK(found >= placed * 98 / 100, "tracked decode: %d/%d", found, placed);
  CHECK(wrong == 0 && others == 0, "tracked decode: %d wrong, %d of another UE", wrong, others);
  tx_free(&tx);
  gs_destroy(gs);
}

/* 6. Soft-information error correction (gs_ec.c, 5GDescrambler's LGS correction).
 * DCIs from srsran_pdcch_nr_encode, QPSK over AWGN, LLR = 2 sqrt(2) Re/Im(y) / sigma^2
 * (SNR = Es/N0 per RE), decoded with the true size, N_ID and scrambling case, as
 * gs_decode_tracked does:
 *   6.0 the soft path alone (hard decision skipped) decodes noiseless DCIs: its own
 *       tables agree with srsRAN's encoder;
 *   6.1 detection vs SNR, hard decision (gs_dci_try) vs soft (gs_dci_try_soft), per AL;
 *   6.2 false alarms on noise only: per candidate per size (both scrambling cases tried),
 *       at the default acceptance (min_parity 16) and at 8 and 12 against the bound
 *       2^-min_parity per call;
 *   6.3 false alarms on another UE's DCI: every other size, AL, CCE window and
 *       scrambling case over a DCI (and the empty CCEs next to it) at 4 and 12 dB;
 *   6.4 time per soft call (error-correction path, and noise-only inputs).
 * Three soft decoders: "auto" = 5GDescrambler's error budget (10% of E with n_RNTI =
 * RNTI, 5% with n_RNTI = 0), "max" = GS_EC_EFFORT_MAX (t = Q - min_parity with n_RNTI =
 * RNTI), "max+m" = max with metric >= EC_MIN_METRIC; all at min_parity 16 and the
 * default neighbour check.
 * Pass: 6.0 exact; soft never worse than hard and no wrong DCI in 6.1; 6.2 <= 1e-4 per
 * candidate per size and 6.3 <= 1e-4 above the hard decision (see 6.3), for auto and
 * max+m (max sits at its bound, 2 x 2^-16, and is reported only). GS_EC_TRIALS scales
 * 6.2 and 6.3 (default 1e5 noise candidates per AL). */
#define EC_NSIZES 8
static const uint32_t ec_sizes[EC_NSIZES] = {39, 43, 47, 51, 55, 59, 63, 70};

/* hard bits of one DCI at CCE 0, in mapping order */
static uint32_t ec_encode_dci(tx_t* tx, uint32_t Ll, uint32_t nbits, bool ue, uint32_t n_id, uint16_t rnti,
                              uint8_t* payload, uint8_t* bits)
{
  memset(&tx->coreset, 0, sizeof(tx->coreset));
  tx->carrier.pci                        = ue ? 1 : n_id;
  tx->coreset.id                         = 1;
  tx->coreset.duration                   = 1;
  tx->coreset.mapping_type               = srsran_coreset_mapping_type_non_interleaved;
  tx->coreset.dmrs_scrambling_id_present = ue;
  tx->coreset.dmrs_scrambling_id         = ue ? n_id : 0;
  for (int r = 0; r < 16; r++) {
    tx->coreset.freq_resources[r] = true;
  }
  srsran_vec_cf_zero(tx->grid, 106 * 12 * 14);
  tx_put(tx, 0, Ll, 0, ue ? srsran_search_space_type_ue : srsran_search_space_type_common_3, rnti, nbits, payload);
  uint32_t E = 0;
  for (uint32_t k = 0; k < (72U << Ll); k++) {
    if (k % 4 == 1) {
      continue;
    }
    bits[E++] = crealf(tx->grid[k]) < 0;
    bits[E++] = cimagf(tx->grid[k]) < 0;
  }
  return E;
}

/* N(0, 1). Not srsran_random_gauss_dist(): its std::normal_distribution(sigma) takes
 * sigma as the mean, so it returns N(sigma, 1). */
static float ec_gauss(void)
{
  float u1 = srsran_random_uniform_real_dist(rnd, 1e-12f, 1.0f), u2 = srsran_random_uniform_real_dist(rnd, 0.0f, 1.0f);
  return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
}

/* LLRs of QPSK bits through AWGN (bits == NULL: noise only, nothing transmitted) */
static void ec_channel(const uint8_t* bits, uint32_t E, float snr_db, float* llr)
{
  const float n0 = powf(10.0f, -snr_db / 10.0f), sd = sqrtf(n0 / 2.0f);
  for (uint32_t i = 0; i < E; i++) {
    float r = (bits ? (bits[i] ? -1.0f : 1.0f) * (float)M_SQRT1_2 : 0.0f) + sd * ec_gauss();
    llr[i]  = 2.0f * (float)M_SQRT2 * r / n0;
  }
}

static int ec_cmp_double(const void* a, const void* b)
{
  double x = *(const double*)a, y = *(const double*)b;
  return x < y ? -1 : x > y;
}

static void ec_time_report(const char* what, double* t, uint32_t n)
{
  if (n == 0) {
    printf("      %-34s -\n", what);
    return;
  }
  double sum = 0;
  for (uint32_t i = 0; i < n; i++) {
    sum += t[i];
  }
  qsort(t, n, sizeof(double), ec_cmp_double);
  printf("      %-34s mean %6.1f us, median %6.1f us, p99 %6.1f us, max %8.1f us (%u calls)\n", what, 1e6 * sum / n,
         1e6 * t[n / 2], 1e6 * t[(uint32_t)(0.99 * (n - 1))], 1e6 * t[n - 1], n);
}

#define EC_MIN_METRIC 0.10f
static void test_soft(void)
{
  tx_t tx;
  tx_init(&tx, 106, 0);
  gs_ec_args_t def;
  gs_ec_args_default(&def);
  for (uint32_t i = 0; i < EC_NSIZES; i++) {
    gs_ec_prepare(ec_sizes[i]);
  }
  static uint8_t bits[108 << 4], payload[128], got[GS_MAX_DCI_BITS];
  static float   llr[108 << 4];
  uint16_t       r_got;
  /* soft decoders compared: 5GDescrambler's error budget (10% / 5% of E), the largest
   * budget min_parity allows (t = Q - min_parity), and that with a metric threshold */
#define EC_NCFG 3
  const uint32_t    effort[EC_NCFG] = {GS_EC_EFFORT_AUTO, GS_EC_EFFORT_MAX, GS_EC_EFFORT_MAX};
  const char* const ename[EC_NCFG]  = {"auto", "max", "max+m"};
  gs_ec_args_t      cfg[EC_NCFG]    = {def, def, def};
  cfg[0].min_metric = cfg[1].min_metric = -INFINITY;
  cfg[2].min_metric = EC_MIN_METRIC;

  /* 6.0 */
  {
    gs_ec_args_t a = def;
    a.hard_first   = false;
    int ok = 0, n = 0;
    for (int t = 0; t < 200; t++) {
      uint32_t Ll = rint_(0, 4), nb = ec_sizes[rint_(0, EC_NSIZES - 1)];
      bool     ue   = rint_(0, 1);
      uint32_t n_id = ue ? rint_(0, 65535) : rint_(0, 1007);
      uint16_t rnti = (uint16_t)rint_(1, 0xfff0);
      uint32_t E    = ec_encode_dci(&tx, Ll, nb, ue, n_id, rnti, payload, bits);
      ec_channel(bits, E, 40.0f, llr);
      gs_ec_info_t info;
      n++;
      ok += gs_dci_try_soft_ex(&a, llr, E, nb, n_id, ue, GS_EC_EFFORT_AUTO, got, &r_got, &info) && r_got == rnti &&
            memcmp(got, payload, nb) == 0 && !info.hard && info.flips == 0;
    }
    CHECK(ok == n, "soft path alone, noiseless: %d/%d", ok, n);
    printf("6.0 Soft path alone (no hard decision), noiseless: %d/%d\n", ok, n);
  }

  /* 6.1 detection vs SNR, and the time of the correcting calls */
  static const float snr[] = {-2, 0, 2, 4, 6, 8, 10, 12};
  const uint32_t     NSNR  = sizeof(snr) / sizeof(snr[0]);
  static double      t_ec[EC_NCFG][4][8000];
  uint32_t           n_ec[EC_NCFG][4] = {{0}};
  printf("6.1 Detection (%%) vs SNR (Es/N0 per RE), true size/N_ID/scrambling, sizes %u..%u, min_parity 16;\n"
         "    auto = 5GDescrambler budget (10%%/5%% of E), max = t = Q - min_parity (<= %u; CSS: auto),\n"
         "    +m = metric >= %.2f:\n",
         ec_sizes[0], ec_sizes[EC_NSIZES - 1], GS_EC_T_MAX, EC_MIN_METRIC);
  printf("      %-17s", "SNR dB");
  for (uint32_t i = 0; i < NSNR; i++) {
    printf("%6.0f", snr[i]);
  }
  printf("\n");
  int wrong_hard = 0, wrong_soft[EC_NCFG] = {0}, worse = 0;
  for (int ue = 1; ue >= 0; ue--) {
    for (uint32_t Ll = 0; Ll <= 3; Ll++) {
      const int trials = ue ? 400 : 200;
      int       ok_h[8] = {0}, ok_s[EC_NCFG][8] = {{0}};
      for (int t = 0; t < trials; t++) {
        uint32_t nb   = ec_sizes[rint_(0, EC_NSIZES - 1)];
        uint32_t n_id = ue ? rint_(0, 65535) : rint_(0, 1007);
        uint16_t rnti = (uint16_t)rint_(0x4601, 0xfff0);
        uint32_t E    = ec_encode_dci(&tx, Ll, nb, ue, n_id, rnti, payload, bits);
        for (uint32_t si = 0; si < NSNR; si++) {
          ec_channel(bits, E, snr[si], llr);
          uint8_t hb[108 << 4];
          for (uint32_t i = 0; i < E; i++) {
            hb[i] = llr[i] < 0;
          }
          if (gs_dci_try(hb, E, nb, n_id, ue, got, &r_got)) {
            bool good = r_got == rnti && memcmp(got, payload, nb) == 0;
            ok_h[si] += good;
            wrong_hard += !good;
          }
          for (int c = 0; c < EC_NCFG; c++) {
            gs_ec_info_t info;
            double       t0 = now_s();
            bool         s  = gs_dci_try_soft_ex(&cfg[c], llr, E, nb, n_id, ue, effort[c], got, &r_got, &info);
            double       dt = now_s() - t0;
            if (!info.hard && info.t > 0 && n_ec[c][Ll] < 8000) {
              t_ec[c][Ll][n_ec[c][Ll]++] = dt;
            }
            if (s) {
              bool good = r_got == rnti && memcmp(got, payload, nb) == 0;
              ok_s[c][si] += good;
              wrong_soft[c] += !good;
              if (!good) {
                printf("   wrong %s: %s AL%u size %u %.0f dB: t %u rank %u Q %u flips %u eff %.1f metric %.3f rnti %04x/%04x\n",
                       ename[c], ue ? "UE" : "CSS", 1U << Ll, nb, snr[si], info.t, info.rank, info.Q, info.flips,
                       info.eff_parity, info.metric, r_got, rnti);
              }
            }
          }
        }
      }
      printf("      %s AL%u hard  ", ue ? "UE " : "CSS", 1U << Ll);
      for (uint32_t i = 0; i < NSNR; i++) {
        printf("%6.1f", 100.0 * ok_h[i] / trials);
      }
      for (int c = 0; c < EC_NCFG; c++) {
        printf("\n      %s AL%u %-6s", ue ? "UE " : "CSS", 1U << Ll, ename[c]);
        for (uint32_t i = 0; i < NSNR; i++) {
          printf("%6.1f", 100.0 * ok_s[c][i] / trials);
          worse += ok_s[c][i] < ok_h[i];
        }
      }
      printf("\n");
    }
  }
  printf("      wrong DCIs: hard %d, auto %d, max %d, max+m %d\n", wrong_hard, wrong_soft[0], wrong_soft[1],
         wrong_soft[2]);
  CHECK(worse == 0, "soft below hard at %d points", worse);
  CHECK(wrong_soft[0] + wrong_soft[1] + wrong_soft[2] == 0, "soft decoder: %d / %d / %d wrong DCIs", wrong_soft[0],
        wrong_soft[1], wrong_soft[2]);

  /* 6.2 noise only */
  const uint32_t ntr = getenv("GS_EC_TRIALS") ? (uint32_t)atoi(getenv("GS_EC_TRIALS")) : 100000;
  const uint32_t nrl = ntr / 5; /* relaxed acceptances, to check the bound */
  static double  t_noise[EC_NCFG][4][20000], t_hard[4][20000];
  uint32_t       n_noise[EC_NCFG][4] = {{0}}, n_hard[4] = {0};
  double         worst_fa[EC_NCFG]   = {0};
  printf("6.2 Noise only, %u candidates per AL, each tried at every size and both scrambling cases (= 2 calls);\n"
         "    false DCIs per candidate per size at min_parity 16; per call at 12 and 8 with the largest budget:\n",
         ntr);
  for (uint32_t Ll = 0; Ll <= 3; Ll++) {
    const uint32_t E = 108U << Ll;
    uint32_t       fa[EC_NCFG][EC_NSIZES] = {{0}}, fa_h = 0, fa8 = 0, fa12 = 0, calls = 0;
    gs_ec_args_t   a8 = def, a12 = def;
    a8.min_parity     = 8;
    a12.min_parity    = 12;
    a8.min_metric = a12.min_metric = -INFINITY;
    for (uint32_t t = 0; t < ntr; t++) {
      ec_channel(NULL, E, 0.0f, llr);
      uint8_t hb[108 << 4];
      for (uint32_t i = 0; i < E; i++) {
        hb[i] = llr[i] < 0;
      }
      uint32_t n_id = rint_(0, 65535);
      for (uint32_t z = 0; z < EC_NSIZES; z++) {
        for (int h = 0; h < 2; h++) {
          double t0h = now_s();
          fa_h += gs_dci_try(hb, E, ec_sizes[z], n_id, h, got, &r_got);
          if (n_hard[Ll] < 20000) {
            t_hard[Ll][n_hard[Ll]++] = now_s() - t0h;
          }
          for (int c = 0; c < EC_NCFG; c++) {
            double t0 = now_s();
            fa[c][z] += gs_dci_try_soft_ex(&cfg[c], llr, E, ec_sizes[z], n_id, h, effort[c], got, &r_got, NULL);
            double dt = now_s() - t0;
            if (n_noise[c][Ll] < 20000) {
              t_noise[c][Ll][n_noise[c][Ll]++] = dt;
            }
          }
          if (t < nrl) {
            fa8 += gs_dci_try_soft_ex(&a8, llr, E, ec_sizes[z], n_id, h, GS_EC_T_MAX, got, &r_got, NULL);
            fa12 += gs_dci_try_soft_ex(&a12, llr, E, ec_sizes[z], n_id, h, GS_EC_T_MAX, got, &r_got, NULL);
            calls++;
          }
        }
      }
    }
    for (int c = 0; c < EC_NCFG; c++) {
      uint32_t tot = 0, worst = 0;
      printf("      AL%u %-6s per size:", 1U << Ll, ename[c]);
      for (uint32_t z = 0; z < EC_NSIZES; z++) {
        printf(" %u:%u", ec_sizes[z], fa[c][z]);
        tot += fa[c][z];
        worst = fa[c][z] > worst ? fa[c][z] : worst;
      }
      double rate = (double)worst / ntr;
      worst_fa[c] = rate > worst_fa[c] ? rate : worst_fa[c];
      printf("  -> %.1e (worst size %.1e)\n", (double)tot / ((double)ntr * EC_NSIZES), rate);
    }
    printf("      AL%u hard %u; max at min_parity 12: %.2e per call (bound %.2e), 8: %.2e (bound %.2e)\n",
           1U << Ll, fa_h, (double)fa12 / calls, 1.0 / 4096, (double)fa8 / calls, 1.0 / 256);
  }
  CHECK(worst_fa[0] <= 1e-4 && worst_fa[2] <= 1e-4, "noise-only false alarms %.2e / %.2e per candidate per size",
        worst_fa[0], worst_fa[2]);

  /* 6.3 another UE's DCI. The hard decision is the baseline: with n_RNTI = 0 it already
   * solves other sizes in the DCI's own window (nested polar codes, ~2^-8 per size),
   * which gs_search handles with its CSS size list; the soft decoders must not add more
   * than 1e-4 to it (that one scaled by the soft detection gain). */
  const int nc = getenv("GS_EC_TRIALS") ? (int)(ntr / 30) : 3000;
  printf("6.3 Another UE's DCI (CCE 0, empty CCEs after it), %d per AL and scrambling, half at 4 dB, half at 12;\n"
         "    tried at every other size, AL, CCE window and scrambling; false DCIs per candidate per size\n"
         "    (its own window at other sizes | other ALs and windows):\n",
         nc);
  double worst_col = 0;
  for (int ue = 1; ue >= 0; ue--) {
    for (uint32_t Lb = 0; Lb <= 3; Lb++) {
      uint64_t events = 0, events_ss = 0, falses[EC_NCFG + 1] = {0}, falses_ss[EC_NCFG + 1] = {0},
               hits[EC_NCFG + 1] = {0};
      for (int t = 0; t < nc; t++) {
        uint32_t nb   = ec_sizes[rint_(0, EC_NSIZES - 1)];
        uint32_t n_id = ue ? rint_(0, 65535) : rint_(0, 1007);
        uint16_t rnti = (uint16_t)rint_(0x4601, 0xfff0);
        uint32_t Eb   = ec_encode_dci(&tx, Lb, nb, ue, n_id, rnti, payload, bits);
        float    sn   = (t & 1) ? 12.0f : 4.0f;
        ec_channel(bits, Eb, sn, llr);
        ec_channel(NULL, 864 - Eb, sn, llr + Eb); /* empty CCEs: noise at the same level */
        for (uint32_t Ld = 0; Ld <= 3; Ld++) {
          const uint32_t Ld_n = 1U << Ld, Lb_n = 1U << Lb, Ed = 108U << Ld;
          for (uint32_t o = 0; o < (Ld <= Lb ? Lb_n : 1); o += Ld_n) {
            const bool same = Ld == Lb && o == 0;
            uint8_t    hb[108 << 3];
            for (uint32_t i = 0; i < Ed; i++) {
              hb[i] = llr[108 * o + i] < 0;
            }
            for (uint32_t z = 0; z < EC_NSIZES; z++) {
              const uint32_t nz = ec_sizes[z];
              for (int h = 0; h < 2; h++) {
                for (int c = 0; c <= EC_NCFG; c++) {
                  bool ok = c == 0 ? nz + 16 + def.min_parity <= Ed && gs_dci_try(hb, Ed, nz, n_id, h, got, &r_got)
                                   : gs_dci_try_soft_ex(&cfg[c - 1], llr + 108 * o, Ed, nz, n_id, h, effort[c - 1],
                                                        got, &r_got, NULL);
                  if (!ok) {
                    continue;
                  }
                  if (same && nz == nb && h == ue && r_got == rnti && memcmp(got, payload, nb) == 0) {
                    hits[c]++;
                  } else if (same) {
                    falses_ss[c]++;
                  } else {
                    falses[c]++;
                  }
                }
              }
              events_ss += same;
              events += !same;
            }
          }
        }
      }
      printf("      %s AL%u:", ue ? "UE " : "CSS", 1U << Lb);
      for (int c = 0; c <= EC_NCFG; c++) {
        double r1 = (double)falses[c] / events, r2 = (double)falses_ss[c] / events_ss;
        double h1 = (double)falses[0] / events, h2 = (double)falses_ss[0] / events_ss;
        /* n_RNTI = 0, own window: a size-N DCI is a size-N' DCI too (nested codes), so the
         * soft decoder finds those wherever it finds more DCIs: held to the hard rate
         * scaled by the detection gain */
        if (!ue) {
          h2 *= hits[0] > 0 ? (double)hits[c] / (double)hits[0] : 1.0;
        }
        if (c == 1 || c == 3) {
          worst_col = r1 - h1 > worst_col ? r1 - h1 : worst_col;
          worst_col = r2 - h2 > worst_col ? r2 - h2 : worst_col;
        }
        printf("%s %s %lu/%d %.1e|%.1e", c ? "," : "", c ? ename[c - 1] : "hard", (unsigned long)hits[c], nc, r2, r1);
      }
      printf("\n");
    }
  }
  CHECK(worst_col <= 1e-4, "collision false alarms %.2e per candidate per size above the hard decision", worst_col);

  /* 6.4 */
  printf("6.4 Time per call (one size, one scrambling case), hard decision for reference:\n");
  for (uint32_t Ll = 0; Ll <= 3; Ll++) {
    char what[64];
    snprintf(what, sizeof(what), "gs_dci_try AL%u noise only (6.2)", 1U << Ll);
    ec_time_report(what, t_hard[Ll], n_hard[Ll]);
  }
  for (int c = 0; c < EC_NCFG; c++) {
    for (uint32_t Ll = 0; Ll <= 3; Ll++) {
      char what[64];
      snprintf(what, sizeof(what), "%s AL%u correcting (6.1)", ename[c], 1U << Ll);
      ec_time_report(what, t_ec[c][Ll], n_ec[c][Ll]);
      snprintf(what, sizeof(what), "%s AL%u noise only (6.2)", ename[c], 1U << Ll);
      ec_time_report(what, t_noise[c][Ll], n_noise[c][Ll]);
    }
  }
  tx_free(&tx);
}

int main(int argc, char** argv)
{
  uint32_t seed = argc > 1 ? (uint32_t)atoi(argv[1]) : 1234;
  rnd           = srsran_random_init(seed);
  gs_args_t args;
  gs_args_default(&args);
  gs_t* warm = gs_create(&args); /* builds the shared tables once */
  gs_destroy(warm);

  int only = argc > 2 ? atoi(argv[2]) : 0;
  if (only == 0 || only == 1) {
    test_dmrs_direct();
  }
  if (only == 0 || only == 2) {
    test_dmrs_underdet();
  }
  if (only == 0 || only == 3) {
    test_gf2();
  }
  if (only == 0 || only == 4) {
    test_search();
  }
  if (only == 0 || only == 5) {
    test_tracked();
  }
  if (only == 0 || only == 6) {
    test_soft();
  }

  srsran_random_free(rnd);
  printf("%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
