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
 *
 * Pass: 1-3 exact; 4 finds >= 95% of the DCIs with no wrong DCI reported, except
 * CORESET#0 with unrestricted sizes (n_RNTI = 0, nested polar codes): <= 1% wrong.
 */
#include "nrscope/hdr/golden_sniffer/gs_core.h"
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


/* 7. Many tracked UEs at once (gs_decode_tracked_multi): 40 UEs in one CORESET of 102 RBs,
 * 30 of them tracked with a learned extent of one CCE and AL1 only (full_band must find
 * them anywhere, at any AL), the other 10 must come out as untracked sightings; up to 6
 * DCIs per slot, PCI or pdcch-DMRS-ScramblingID scrambling, duration 1-2, 15 dB. Every
 * fourth slot is noise only: nothing may be reported there. Timed with 30 and with 200
 * tracks (the cost must not grow with the number of tracked UEs). */
static void test_tracked_multi(void)
{
  gs_args_t args;
  gs_args_default(&args);
  gs_t*          gs = gs_create(&args);
  tx_t           tx;
  const uint32_t nprb = 106, NU = 40, NT = 30, NX = 200;
  tx_init(&tx, nprb, 0);
  static const uint32_t szs[4] = {43, 47, 51, 55};
  uint16_t              rnti[40];
  uint32_t              dl[40], ul[40];
  for (uint32_t u = 0; u < NU; u++) {
    rnti[u] = (uint16_t)(0x4601 + 7 * u);
    dl[u]   = szs[(u / 2) % 2 ? 1 : 3];
    ul[u]   = szs[(u / 2) % 2 ? 0 : 2];
  }
  int    placed_t = 0, found_t = 0, placed_u = 0, found_u = 0, wrong = 0, false_noise = 0, nslots = 0;
  double t30 = 0, t200 = 0;
  for (int t = 0; t < 400; t++) {
    const bool     scr = t & 1;
    const uint32_t dur = 1 + ((t >> 1) & 1), slot = rint_(0, 19), pci = rint_(0, 1007);
    const uint32_t n_id = scr ? (uint32_t)rint_(0, 65535) : pci;
    tx.carrier.pci      = pci;
    memset(&tx.coreset, 0, sizeof(tx.coreset));
    tx.coreset.id                         = 1;
    tx.coreset.duration                   = dur;
    tx.coreset.dmrs_scrambling_id_present = scr;
    tx.coreset.dmrs_scrambling_id         = n_id;
    for (uint32_t r = 0; r < 17; r++) {
      tx.coreset.freq_resources[r] = true; /* RBs 0..101 */
    }
    srsran_vec_cf_zero(tx.grid, nprb * 12 * 14);
    const uint32_t n_cce = 17 * dur;
    bool           busy[64] = {false};
    struct {
      uint32_t u, L, ncce, nb;
      uint8_t  pay[128];
      bool     hit;
    } ex[6];
    uint32_t   nex   = 0;
    const bool noise = (t % 4) == 3;
    bool       used_u[40] = {false};
    for (uint32_t k = 0, want = noise ? 0 : (uint32_t)rint_(1, 6); k < want; k++) {
      uint32_t u = (uint32_t)rint_(0, NU - 1);
      if (used_u[u]) {
        continue;
      }
      uint32_t Ll = (uint32_t)rint_(0, dur == 1 ? 3 : 4), L = 1U << Ll;
      if (L > n_cce) {
        continue;
      }
      uint32_t ncce = L * (uint32_t)rint_(0, n_cce / L - 1);
      bool     free_ = true;
      for (uint32_t c = ncce; c < ncce + L; c++) {
        free_ = free_ && !busy[c];
      }
      if (!free_) {
        continue;
      }
      uint32_t nb = rint_(0, 1) ? dl[u] : ul[u];
      if (tx_put(&tx, slot, Ll, ncce, srsran_search_space_type_ue, rnti[u], nb, ex[nex].pay) < 0) {
        continue;
      }
      for (uint32_t c = ncce; c < ncce + L; c++) {
        busy[c] = true;
      }
      used_u[u]    = true;
      ex[nex].u    = u;
      ex[nex].L    = L;
      ex[nex].ncce = ncce;
      ex[nex].nb   = nb;
      ex[nex].hit  = false;
      nex++;
      if (u < NT) {
        placed_t++;
      } else {
        placed_u++;
      }
    }
    float sigma = powf(10.0f, -15.0f / 20.0f) * (float)M_SQRT1_2;
    float ph = srsran_random_uniform_real_dist(rnd, -M_PI, M_PI), tau = srsran_random_uniform_real_dist(rnd, -3, 3) / 1024.0f;
    for (uint32_t l = 0; l < 3; l++) {
      for (uint32_t k = 0; k < nprb * 12; k++) {
        cf_t h = cexpf(I * (ph + 0.2f * l + 2.0f * (float)M_PI * tau * k));
        cf_t n = sigma * (srsran_random_gauss_dist(rnd, 1.0f) + I * srsran_random_gauss_dist(rnd, 1.0f));
        tx.grid[l * nprb * 12 + k] = tx.grid[l * nprb * 12 + k] * h + n;
      }
    }
    gs_slot_t  sl = {.grid = tx.grid, .nof_sc = nprb * 12, .slot_idx = slot, .ref_crb0 = 0, .ref_cs0 = GS_NO_REF, .pci = pci};
    gs_track_t tr[200];
    for (uint32_t i = 0; i < NX; i++) {
      const uint32_t u = i < NT ? i : 0;
      tr[i]            = (gs_track_t){.rnti = i < NT ? rnti[u] : (uint16_t)(0x9000 + i), .n_id = n_id, .n_rnti = scr,
                                      .coreset0 = false, .sc_lo = 0, .sc_hi = (int32_t)(72 / dur) - 1, .duration = dur,
                                      .al_mask = 1, .sizes = {dl[u], ul[u]}, .nof_sizes = 2};
    }
    gs_dci_t out[16];
    uint8_t  trk[16];
    double   t0 = now_s();
    int      n  = gs_decode_tracked_multi(gs, &sl, tr, NT, true, true, out, trk, 16);
    t30 += now_s() - t0;
    gs_dci_t out2[16];
    uint8_t  trk2[16];
    t0 = now_s();
    gs_decode_tracked_multi(gs, &sl, tr, NX, true, true, out2, trk2, 16);
    t200 += now_s() - t0;
    nslots++;
    for (int i = 0; i < n; i++) {
      if (noise) {
        false_noise++;
        continue;
      }
      bool ok = false;
      for (uint32_t k = 0; k < nex && !ok; k++) {
        const uint32_t u = ex[k].u;
        if (ex[k].hit || out[i].rnti != rnti[u] || out[i].nof_bits != ex[k].nb || out[i].L != ex[k].L ||
            out[i].cand.iv[0].first != (int32_t)(ex[k].ncce * 72 / dur) ||
            memcmp(out[i].payload, ex[k].pay, ex[k].nb) != 0 || trk[i] != (u < NT)) {
          continue;
        }
        ex[k].hit = ok = true;
        if (u < NT) {
          found_t++;
        } else {
          found_u++;
        }
      }
      wrong += !ok;
    }
  }
  printf("7. Tracked multi-UE: tracked %d/%d, untracked sightings %d/%d, %d wrong, %d on noise-only slots; "
         "%.3f ms/slot with 30 tracks, %.3f with 200\n",
         found_t, placed_t, found_u, placed_u, wrong, false_noise, 1e3 * t30 / nslots, 1e3 * t200 / nslots);
  CHECK(found_t >= placed_t * 98 / 100, "tracked multi: %d/%d", found_t, placed_t);
  CHECK(found_u >= placed_u * 95 / 100, "tracked multi, untracked sightings: %d/%d", found_u, placed_u);
  CHECK(wrong == 0 && false_noise == 0, "tracked multi: %d wrong, %d on noise", wrong, false_noise);
  tx_free(&tx);
  gs_destroy(gs);
}

int main(int argc, char** argv)
{
  uint32_t seed = argc > 1 ? (uint32_t)atoi(argv[1]) : 1234;
  rnd           = srsran_random_init(seed);
  gs_args_t args;
  gs_args_default(&args);
  gs_t* warm = gs_create(&args); /* builds the shared tables once */
  gs_destroy(warm);

  test_dmrs_direct();
  test_dmrs_underdet();
  test_gf2();
  test_search();
  test_tracked();
  test_tracked_multi();

  srsran_random_free(rnd);
  printf("%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
