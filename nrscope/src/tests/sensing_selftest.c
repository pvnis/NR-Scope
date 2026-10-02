/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * \brief Ground-truth test of the sensing pipeline, with no radio and no gNB.
 *
 * The pipeline came from OpenAirInterface5G, where it could be validated against
 * rfsimulator with a target at a known range and Doppler. Nothing equivalent
 * exists here, and this is the replacement. It is also the narrower instrument:
 * rather than simulating a PHY and inferring the channel from it, it writes the
 * channel down directly, so the expected answer is exact rather than
 * approximate.
 *
 * A channel is synthesised from a list of paths, each with a delay expressed in
 * range bins, a speed and a complex gain:
 *
 *     H(k, t) = sum_p  g_p exp(-j 2 pi k df tau_p) exp(+j 2 pi f_D,p t)
 *
 * with tau_p chosen so the path lands on an exact bin of the delay axis, and
 * f_D,p = 2 v_p / lambda so it lands on a known speed. That grid is fed through
 * the same entry the receive path uses, nr_ue_sensing_slot_profile(), which
 * exercises lattice extraction, the Hann taper, the IDFT and the history ring;
 * then nr_ue_sensing_range_doppler() builds the map. The test asserts the peak
 * is in the cell the paths were placed in.
 *
 * What this is for: catching the day a change to the transform, the windowing,
 * the ring or the clutter filter moves a target off its true cell. It says
 * nothing about whether the DM-RS were extracted from a real grid correctly,
 * which is the job of the stage above it.
 */

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_dmrs_despread.h"
#include "nrscope/hdr/sensing/nr_ue_map.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"
#include "nrscope/hdr/sensing/nr_ue_sensing_align.h"
#include "nrscope/hdr/sensing/nrscope_sensing.h"
#include "nrscope/hdr/sensing/sensing_defs.h"
#include "srsran/phy/ch_estimation/dmrs_sch.h"
#include "srsran/phy/common/sequence.h"
#include "srsran/phy/phch/phch_cfg_nr.h"

/* Carrier the X410 configuration targets: the 100 MHz Sunrise cell, 273 RBs at
30 kHz. The numbers matter because the delay axis scale and the Doppler-to-speed
conversion both come from them. */
#define TEST_N_RB 273
#define TEST_N_SC (TEST_N_RB * 12)
#define TEST_SCS_HZ 30000
#define TEST_OFDM_SYMBOL_SIZE 4096
#define TEST_CARRIER_HZ 3711360000ULL
#define TEST_SYMBOLS_PER_SLOT 14

/// Sample clock implied by the numerology: 4096 * 30 kHz = 122.88 Msps
#define TEST_FS_HZ ((double)TEST_OFDM_SYMBOL_SIZE * TEST_SCS_HZ)
/// 0.5 ms at 30 kHz
#define TEST_SAMPLES_PER_SLOT ((uint64_t)(TEST_FS_HZ / 2000.0))

/* DM-RS type 1 at rank 1: every second subcarrier, comb offset 0. This is the
densest routinely available pilot and the case the pipeline is tuned for. */
#define TEST_K_STEP 2
#define TEST_K_OFFSET 0
#define TEST_IDFT_SIZE NR_SENSING_IDFT_SIZE(TEST_K_STEP)

/// Three DM-RS symbols per slot, at positions a type A mapping typically gives
static const int test_dmrs_symbols[] = {2, 7, 11};
#define TEST_N_DMRS_SYM ((int)(sizeof(test_dmrs_symbols) / sizeof(test_dmrs_symbols[0])))

#define C_LIGHT 299792458.0

typedef struct {
  /// delay as a position on the range axis, in bins of m_per_bin
  double range_bin;
  /// radial speed in m/s; the map converts Doppler back to this
  double speed_ms;
  double amplitude;
} test_path_t;

static int failures = 0;

static void check(bool ok, const char* what, double got, double want, double tol)
{
  printf("  %-38s got %9.3f  want %9.3f  +-%.3f   %s\n", what, got, want, tol, ok ? "ok" : "FAIL");
  if (!ok) {
    failures++;
  }
}

/* Metres of path length per delay bin. Pinned by idft_size * k_step * scs, which
is why it does not move when the grant does. */
static double test_m_per_bin(void)
{
  return C_LIGHT / ((double)TEST_IDFT_SIZE * TEST_K_STEP * TEST_SCS_HZ);
}

static double test_speed_to_doppler_hz(double speed_ms)
{
  const double lambda = C_LIGHT / (double)TEST_CARRIER_HZ;
  return 2.0 * speed_ms / lambda;
}

/* Write one OFDM symbol of the synthetic channel onto the pilot lattice.

   Only the lattice positions are filled and marked valid, exactly as the DM-RS
   estimator would leave them: the data REs between pilots are never measured. */
static void test_fill_symbol(cf_t* H_row, bool* valid_row, const test_path_t* paths, int n_paths, double t_seconds)
{
  const double m_per_bin = test_m_per_bin();

  for (int k = TEST_K_OFFSET; k < TEST_N_SC; k += TEST_K_STEP) {
    cf_t acc = 0.0f;
    for (int p = 0; p < n_paths; p++) {
      const double tau = paths[p].range_bin * m_per_bin / C_LIGHT;
      /* Delay is a linear phase across subcarriers; Doppler is a phase that
      advances with absolute time. Both are written in continuous form so the
      path sits where it was asked to sit rather than on a rounded grid. */
      const double ph_delay = -2.0 * M_PI * (double)k * TEST_SCS_HZ * tau;
      const double ph_dopp  = 2.0 * M_PI * test_speed_to_doppler_hz(paths[p].speed_ms) * t_seconds;
      const double ph       = ph_delay + ph_dopp;
      acc += (cf_t)(paths[p].amplitude * cos(ph)) + I * (cf_t)(paths[p].amplitude * sin(ph));
    }
    H_row[k]    = acc;
    valid_row[k] = true;
  }
}

/* Run one scene end to end and report where the strongest moving cell landed.

   n_slots sets the observation span, which is what fixes the Doppler
   resolution; the pipeline needs the span, not the sample count. */
static bool test_run_scene(const char* name,
                           const test_path_t* paths,
                           int n_paths,
                           int n_slots,
                           nr_sensing_clutter_t clutter,
                           double expect_range_m,
                           double expect_speed_ms)
{
  printf("\n%s\n", name);

  nr_sensing_history_t hist;
  if (!nr_ue_sensing_history_init(&hist, NR_SENSING_HISTORY_DEPTH, TEST_SCS_HZ, TEST_OFDM_SYMBOL_SIZE, TEST_CARRIER_HZ)) {
    printf("  history init failed\n");
    failures++;
    return false;
  }

  /* One slot of channel at a time, so the grids stay on the stack and the test
  mirrors how the receive path actually delivers symbols. */
  static cf_t H_est[TEST_SYMBOLS_PER_SLOT][TEST_N_SC];
  static bool H_valid[TEST_SYMBOLS_PER_SLOT][TEST_N_SC];
  uint64_t    t_sample[TEST_SYMBOLS_PER_SLOT];

  const nr_sensing_stream_t stream = {.ports = 1, .layer = 0, .k_step = TEST_K_STEP, .k_offset = TEST_K_OFFSET};

  int pushed = 0;
  for (int slot = 0; slot < n_slots; slot++) {
    memset(H_est, 0, sizeof(H_est));
    memset(H_valid, 0, sizeof(H_valid));

    for (int i = 0; i < TEST_N_DMRS_SYM; i++) {
      const int m = test_dmrs_symbols[i];
      /* Absolute instant of this symbol on the sample clock. The cyclic prefix
      is ignored: it shifts every symbol by the same small amount, which the
      Doppler transform cannot see. */
      t_sample[m] = (uint64_t)slot * TEST_SAMPLES_PER_SLOT + (uint64_t)m * TEST_OFDM_SYMBOL_SIZE;
      test_fill_symbol(H_est[m], H_valid[m], paths, n_paths, (double)t_sample[m] / TEST_FS_HZ);
    }

    pilot_lattice_t lat        = {0};
    uint16_t        used_mask  = 0;
    const int       idft_size  = nr_ue_sensing_slot_profile(TEST_SYMBOLS_PER_SLOT,
                                                            TEST_N_SC,
                                                            H_est,
                                                            H_valid,
                                                            t_sample,
                                                            &hist,
                                                            stream,
                                                            (uint64_t)slot,
                                                            false, // no random drop: the test wants every sample
                                                            0,
                                                            &lat,
                                                            &used_mask);
    if (idft_size <= 0) {
      printf("  slot %d produced no profile\n", slot);
      failures++;
      nr_ue_sensing_history_free(&hist);
      return false;
    }
    if (slot == 0) {
      check(idft_size == TEST_IDFT_SIZE, "idft size", idft_size, TEST_IDFT_SIZE, 0);
      check(lat.k_step == TEST_K_STEP, "lattice k_step", lat.k_step, TEST_K_STEP, 0);
    }
    pushed += TEST_N_DMRS_SYM;
  }

  nr_sensing_map_t* map = calloc_or_fail(1, sizeof(*map));
  const int         n   = nr_ue_sensing_range_doppler(&hist,
                                                      &stream,
                                                      pushed < NR_SENSING_HISTORY_DEPTH ? pushed : NR_SENSING_HISTORY_DEPTH,
                                                      20.0, // max speed of the grid, m/s; see the span criteria below
                                                      clutter,
                                                      NR_CLUTTER_MAX_PATHS,
                                                      0,
                                                      map,
                                                      NULL);
  if (n <= 0) {
    printf("  range-doppler produced no map (%d snapshots)\n", n);
    failures++;
    free(map);
    nr_ue_sensing_history_free(&hist);
    return false;
  }

  /* A shorter window sizes its own Doppler grid differently, as a second Rx chain one
  snapshot short of the first does, and its map can then not be summed with this one
  cell by cell. Built on this map's grid (n_freq_fixed) it has to land on it exactly. */
  {
    nr_sensing_map_t* other = calloc_or_fail(1, sizeof(*other));
    const int         n_short = n - n / 20;
    const int n_own = nr_ue_sensing_range_doppler(&hist, &stream, n_short, 20.0, clutter, NR_CLUTTER_MAX_PATHS, 0, other, NULL);
    const int own_freq = other->n_freq;
    const int n_fix =
        nr_ue_sensing_range_doppler(&hist, &stream, n_short, 20.0, clutter, NR_CLUTTER_MAX_PATHS, map->n_freq, other, NULL);
    check(n_own > 0 && own_freq != map->n_freq, "shorter window sizes its own grid differently", own_freq,
          map->n_freq, 0);
    check(n_fix > 0 && other->n_freq == map->n_freq && other->f_max_hz == map->f_max_hz,
          "shorter window on the first map's grid", other->n_freq, map->n_freq, 0);
    free(other);
  }

  /* Strongest cell away from zero Doppler. The direct path is static and, with
  clutter removal off, dominates the map; the target is what moves. */
  int   best_b = -1, best_f = -1;
  float best_p = -1.0f;
  for (int b = 0; b < map->n_bins; b++) {
    for (int f = 0; f < map->n_freq; f++) {
      const double speed = ((double)f - map->n_freq / 2.0) * (2.0 * map->f_max_hz / map->n_freq)
                           * (C_LIGHT / (double)TEST_CARRIER_HZ) / 2.0;
      if (fabs(speed) < 0.5) {
        continue; // skip the static ridge
      }
      const float p = map->power[(size_t)b * map->n_freq + f];
      if (p > best_p) {
        best_p = p;
        best_b = b;
        best_f = f;
      }
    }
  }

  const double got_range = best_b * map->m_per_bin;
  const double got_speed = ((double)best_f - map->n_freq / 2.0) * (2.0 * map->f_max_hz / map->n_freq)
                           * (C_LIGHT / (double)TEST_CARRIER_HZ) / 2.0;

  printf("  snapshots %d, span %.1f ms, %.3f m/bin, f_max %.1f Hz, %d x %d cells\n",
         map->n_snapshots,
         map->t_span_s * 1e3,
         map->m_per_bin,
         map->f_max_hz,
         map->n_bins,
         map->n_freq);

  /* Tolerances come from resolution, not from grid spacing. The Doppler grid
  deliberately holds two points per resolution cell, so asserting to within a
  grid cell would demand better than the observation span can deliver and fail
  on a correct answer. Speed resolution is (1/t_span) * lambda/2; range
  resolution is c/B, which at a full-band grant is a little coarser than
  m_per_bin. Allowing one resolution cell on each axis is the honest limit of
  where a peak can be said to be. */
  const double lambda        = C_LIGHT / (double)TEST_CARRIER_HZ;
  const double speed_res_ms  = (1.0 / map->t_span_s) * lambda / 2.0;
  const double range_res_m   = C_LIGHT / ((double)TEST_N_SC * TEST_SCS_HZ);
  const double range_tol     = fmax(range_res_m, map->m_per_bin);

  check(fabs(got_range - expect_range_m) <= range_tol, "target range (m)", got_range, expect_range_m, range_tol);
  check(fabs(got_speed - expect_speed_ms) <= speed_res_ms, "target speed (m/s)", got_speed, expect_speed_ms,
        speed_res_ms);

  free(map);
  nr_ue_sensing_history_free(&hist);
  return true;
}


/* Round trip through the DM-RS stage: place real pilots on a grid behind a known
   channel, estimate, and check the channel comes back.
   
   This is the one stage whose failure mode is silent. The pilot sequence has no
   random access: it must be walked forward and stepped over PRBs the grant does
   not use, and a skip that is wrong by a single PRB yields pilots that are
   perfectly valid values in the wrong places. Every estimate then becomes noise
   and nothing reports an error.

   The transmitter here therefore walks the sequence a different way from the
   receiver. It generates the pilots of the whole carrier in one call and indexes
   them absolutely, pilot (p, n', k') at index p*n_pilot_rb + 2n' + k', which is
   the property that makes a gNB's grid independent of who is scheduled. The
   receiver streams and skips. If the two agree on a gapped allocation, the skip
   arithmetic is right; if they do not, the error is total rather than subtle. */
static void test_dmrs_roundtrip(void)
{
  printf("\nDM-RS round trip: gapped allocation, type 1, rank 1\n");

  /* A small carrier keeps the test quick; the sequence walk does not care how
  wide it is, only that it is stepped correctly. */
  const uint32_t nof_prb  = 51;
  const int      n_sc     = (int)nof_prb * SRSRAN_NRE;
  const uint32_t slot_idx = 7;
  const uint32_t symbol   = 2;

  srsran_carrier_nr_t carrier = {};
  carrier.pci                 = 632; // the Sunrise cell; N_ID falls back to this
  carrier.nof_prb             = nof_prb;

  srsran_dmrs_sch_cfg_t dmrs_cfg = {};
  dmrs_cfg.type                  = srsran_dmrs_sch_type_1;
  // scrambling_id0/1 left absent on purpose: that is the case a sniffer sees

  srsran_sch_grant_nr_t grant = {};
  grant.n_scid                = false;
  grant.beta_dmrs             = 0.0f; // unset reads as 1

  /* Two disjoint runs with a gap between them and a gap before the first. This
  is what exercises the advance; a contiguous allocation from PRB 0 would pass
  even with the skip logic removed entirely. */
  const uint32_t runs[][2] = {{4, 12}, {20, 31}};
  uint32_t       n_alloc   = 0;
  for (uint32_t r = 0; r < sizeof(runs) / sizeof(runs[0]); r++) {
    for (uint32_t p = runs[r][0]; p < runs[r][1]; p++) {
      grant.prb_idx[p] = true;
      n_alloc++;
    }
  }

  nr_dmrs_layout_t lay;
  if (!nr_ue_dmrs_layout(0, srsran_dmrs_sch_type_1, &lay)) {
    printf("  layout rejected\n");
    failures++;
    return;
  }
  check(lay.k_step == 2, "layout k_step (rank 1, type 1)", lay.k_step, 2, 0);
  check(!lay.despread, "layout despread flag", lay.despread, 0, 0);

  const uint32_t cinit = nr_ue_dmrs_seed(&carrier, &dmrs_cfg, &grant, slot_idx, symbol);

  /* The channel the test hides behind the pilots: two delay paths, so H varies
  across subcarriers and a constant would not pass. */
  static cf_t H_true[51 * SRSRAN_NRE];
  for (int k = 0; k < n_sc; k++) {
    const double ph1 = -2.0 * M_PI * (double)k * 0.013;
    const double ph2 = -2.0 * M_PI * (double)k * 0.041;
    H_true[k]        = (cf_t)(cos(ph1) + 0.4 * cos(ph2)) + I * (cf_t)(sin(ph1) + 0.4 * sin(ph2));
  }

  /* Transmitter: the whole carrier's pilots in one call, placed by absolute
  index at the REs 38.211 gives, and only on the allocated PRBs. */
  static cf_t tx_pilots[51 * 6];
  srsran_sequence_state_t tx_state = {};
  srsran_sequence_state_init(&tx_state, cinit);
  srsran_sequence_state_gen_f(&tx_state, M_SQRT1_2, (float*)tx_pilots, nof_prb * 6 * 2);

  static cf_t rxF[51 * SRSRAN_NRE];
  memset(rxF, 0, sizeof(rxF));
  int n_placed = 0;
  for (uint32_t p = 0; p < nof_prb; p++) {
    if (!grant.prb_idx[p]) {
      continue;
    }
    for (int np = 0; np < 3; np++) {
      for (int kp = 0; kp < 2; kp++) {
        const int k = (int)p * SRSRAN_NRE + 4 * np + 2 * kp; // delta = 0 for port 0
        const int i = (int)p * 6 + 2 * np + kp;
        rxF[k]      = H_true[k] * tx_pilots[i];
        n_placed++;
      }
    }
  }

  /* Receiver: the estimator under test, on the shared generator (dmrs_pilots.h).
  The transmitter above uses srsRAN's sequence module instead, so this also
  checks the two generators agree. */
  static cf_t H_est[51 * SRSRAN_NRE];
  static bool H_valid[51 * SRSRAN_NRE];
  memset(H_est, 0, sizeof(H_est));
  memset(H_valid, 0, sizeof(H_valid));

  const nr_dmrs_placement_t place = {.grid_crb0 = 0, .bwp_start_crb = 0, .reference_crb = dmrs_cfg.reference_point_k_rb};
  const int n_written = nr_ue_dmrs_estimate_symbol(&lay, 0, cinit, &grant, &place, rxF, n_sc, H_est, H_valid);

  check(n_written == n_placed, "REs estimated vs placed", n_written, n_placed, 0);

  /* Worst error over every pilot the estimator claims to have written. A
  mis-stepped sequence does not degrade this gracefully: the pilots are unit
  modulus, so a wrong one leaves an error of order |H| rather than a small bias. */
  double worst = 0.0;
  int    n_cmp = 0;
  for (int k = 0; k < n_sc; k++) {
    if (!H_valid[k]) {
      continue;
    }
    const double e = cabs(H_est[k] - H_true[k]);
    if (e > worst) {
      worst = e;
    }
    n_cmp++;
  }
  check(n_cmp == n_placed, "REs marked valid", n_cmp, n_placed, 0);
  check(worst < 1e-4, "worst |H_est - H_true|", worst, 0.0, 1e-4);

  /* Nothing outside the allocation may be touched: the lattice stage reads
  valid_row to decide what was measured, and a stray true there would feed an
  unmeasured RE into the transform as if it were data. */
  int n_outside = 0;
  for (uint32_t p = 0; p < nof_prb; p++) {
    if (grant.prb_idx[p]) {
      continue;
    }
    for (int k = (int)p * SRSRAN_NRE; k < (int)(p + 1) * SRSRAN_NRE; k++) {
      if (H_valid[k]) {
        n_outside++;
      }
    }
  }
  check(n_outside == 0, "REs written outside the grant", n_outside, 0, 0);
}

/* Two-layer despreading against srsRAN's own DM-RS mapper.
 *
 * srsran_dmrs_sch_put_sf() writes the port-1000 pilots of a 273-PRB carrier; port
 * 1001 is the same pilots with the frequency cover w_f = [+1 -1] on alternate
 * pilots, as in 38.211 table 7.4.1.1.2-1. Each layer goes through its own
 * channel, a slow delay ramp with its own gain, and the grid holds their sum,
 * which is what a two-layer PDSCH puts on each DM-RS RE. Despreading layer j must
 * return layer j's channel with the other one cancelled.
 *
 * The receiver takes its pilots from the shared generator (dmrs_pilots.h) and
 * looks them up by absolute position, the transmitter here is srsRAN's walk over
 * the allocation, so a gapped grant checks the two agree across the gaps. The
 * shifted case hands the despreader a grid and a BWP starting shift CRBs above
 * point A, with the grant renumbered from the BWP, and must find the same pilots:
 * that checks the placement arithmetic the live grid will depend on. */
static void test_despread_two_layers(const char* name, const uint32_t runs[][2], int nof_runs, uint32_t shift)
{
  printf("\nDM-RS despread, 2 layers on srsRAN's mapper: %s\n", name);
  const uint32_t nof_prb = 273, n_sc = nof_prb * SRSRAN_NRE, slot_idx = 11;

  srsran_carrier_nr_t carrier = SRSRAN_DEFAULT_CARRIER_NR;
  carrier.pci                 = 1;
  carrier.nof_prb             = nof_prb;
  carrier.scs                 = srsran_subcarrier_spacing_30kHz;

  srsran_sch_cfg_nr_t cfg                 = {};
  cfg.dmrs.type                           = srsran_dmrs_sch_type_1;
  cfg.dmrs.typeA_pos                      = srsran_dmrs_sch_typeA_pos_2;
  cfg.dmrs.additional_pos                 = srsran_dmrs_sch_add_pos_2;
  cfg.dmrs.length                         = srsran_dmrs_sch_len_1;
  srsran_sch_grant_nr_t* g                = &cfg.grant;
  g->mapping                              = srsran_sch_mapping_type_A;
  g->S                                    = 1;
  g->L                                    = 13;
  g->nof_layers                           = 1; // srsRAN writes port 1000; port 1001 is built below
  g->nof_dmrs_cdm_groups_without_data     = 1;
  g->beta_dmrs                            = 1.0f;
  uint32_t n_alloc                        = 0;
  for (int r = 0; r < nof_runs; r++) {
    for (uint32_t p = runs[r][0]; p < runs[r][1]; p++) {
      g->prb_idx[p] = true;
      n_alloc++;
    }
  }
  g->nof_prb = n_alloc;

  static cf_t tx[273 * SRSRAN_NRE * SRSRAN_NSYMB_PER_SLOT_NR];
  memset(tx, 0, sizeof(tx));
  srsran_dmrs_sch_t dmrs = {};
  srsran_slot_cfg_t slot = {.idx = slot_idx};
  if (srsran_dmrs_sch_init(&dmrs, false) < SRSRAN_SUCCESS || srsran_dmrs_sch_set_carrier(&dmrs, &carrier) < SRSRAN_SUCCESS ||
      srsran_dmrs_sch_put_sf(&dmrs, &slot, &cfg, g, tx) < SRSRAN_SUCCESS) {
    printf("  srsRAN could not place the DM-RS\n");
    failures++;
    return;
  }
  srsran_dmrs_sch_free(&dmrs);
  uint32_t  symbols[SRSRAN_DMRS_SCH_MAX_SYMBOLS];
  const int nof_symbols = srsran_dmrs_sch_get_symbols_idx(&cfg.dmrs, g, symbols);

  // Two layers, each with its own gain and delay ramp, summed on every DM-RS RE
  const cf_t   a[2]   = {1.0f + 0.3f * I, -0.5f + 0.7f * I};
  const double tau[2] = {0.0015, 0.0025}; // cycles per subcarrier: under a degree across a pilot pair
  static cf_t  rx[273 * SRSRAN_NRE * SRSRAN_NSYMB_PER_SLOT_NR];
  memset(rx, 0, sizeof(rx));
  for (int s = 0; s < nof_symbols; s++) {
    for (uint32_t k = 0; k < n_sc; k++) {
      const cf_t x = tx[symbols[s] * n_sc + k];
      if (x == 0) {
        continue;
      }
      const float w1 = ((k / 2) % 2 == 1) ? -1.0f : 1.0f; // port 1001 cover, k' = (k/2) % 2 on comb 0
      cf_t        y  = 0;
      for (int j = 0; j < 2; j++) {
        const cf_t h = a[j] * cexpf(-I * (float)(2.0 * M_PI * tau[j] * k));
        y += h * x * (j == 1 ? w1 : 1.0f);
      }
      rx[symbols[s] * n_sc + k] = y;
    }
  }

  nr_dmrs_layout_t lay;
  if (!nr_ue_dmrs_layout(0x3, srsran_dmrs_sch_type_1, &lay)) { // ports 1000 and 1001
    printf("  layout rejected\n");
    failures++;
    return;
  }
  check(lay.despread && lay.k_step == 4, "layout: despread, k_step 4", lay.k_step, 4, 0);

  // The receiver's view: grid and BWP start `shift` CRBs above point A
  srsran_sch_grant_nr_t rx_grant = *g;
  memset(rx_grant.prb_idx, 0, sizeof(rx_grant.prb_idx));
  for (uint32_t p = shift; p < nof_prb; p++) {
    rx_grant.prb_idx[p - shift] = g->prb_idx[p];
  }
  const nr_dmrs_placement_t place = {.grid_crb0 = shift, .bwp_start_crb = shift, .reference_crb = 0};
  const int                 n_grid = (int)((nof_prb - shift) * SRSRAN_NRE);

  static cf_t H[273 * SRSRAN_NRE];
  static bool V[273 * SRSRAN_NRE];
  double      worst[2] = {0, 0};
  int         written[2] = {0, 0};
  for (int s = 0; s < nof_symbols; s++) {
    const uint32_t cinit = nr_ue_dmrs_seed(&carrier, &cfg.dmrs, g, slot_idx, symbols[s]);
    for (int j = 0; j < 2; j++) {
      memset(V, 0, sizeof(V));
      const int n = nr_ue_dmrs_estimate_symbol(&lay, j, cinit, &rx_grant, &place,
                                               rx + symbols[s] * n_sc + shift * SRSRAN_NRE, n_grid, H, V);
      written[j] += n;
      for (int kg = 0; kg < n_grid; kg++) {
        if (!V[kg]) {
          continue;
        }
        // Despread output sits at the k' = 0 RE of its pair; compare with the pair's mean channel
        const uint32_t k   = (uint32_t)kg + shift * SRSRAN_NRE;
        const cf_t     h0  = a[j] * cexpf(-I * (float)(2.0 * M_PI * tau[j] * k));
        const cf_t     h1  = a[j] * cexpf(-I * (float)(2.0 * M_PI * tau[j] * (k + 2)));
        const double   err = cabsf(H[kg] - 0.5f * (h0 + h1)) / cabsf(a[j]);
        worst[j]           = fmax(worst[j], err);
      }
    }
  }
  const int want = (int)n_alloc * 3 * nof_symbols; // one estimate per pilot pair, 3 pairs per PRB
  check(written[0] == want, "layer 0: estimates written", written[0], want, 0);
  check(written[1] == want, "layer 1: estimates written", written[1], want, 0);
  // A wrong pilot leaves an error of order 1; the cross-layer leak of a slowly varying channel is ~1e-2
  check(worst[0] < 0.02, "layer 0: worst |H - h0| / |a0|", worst[0], 0, 0.02);
  check(worst[1] < 0.02, "layer 1: worst |H - h1| / |a1|", worst[1], 0, 0.02);
}

/* The glue end to end: slots of a TDD pattern (3 full downlink slots, a special
 * slot with a 5-symbol PDSCH, an uplink slot) go through
 * nrscope_sensing_process_grant() exactly as the DCI decoder hands them over,
 * with grids built by srsRAN's DM-RS mapper on two layers. The scene is a static
 * direct path and a weaker moving target at known delays. The per-slot stage
 * aligns every snapshot on the static scene (nr_ue_sensing_align), so the direct
 * path must come out at 0 Hz and the target at its own Doppler, read from its
 * phase across the snapshots' sample times. Both must peak at their delay bins.
 * The receive window also moves one sample later every 30 slots, as the timing
 * tracking does, with the ramp centred on DC as on air: the glue must undo the
 * moves it is told about, leaving the direct path on its bin in every snapshot.
 * Every DM-RS symbol must reach the history (the count is exact), which guards
 * the second-symbol rejection seen on the first port. That covers placement,
 * lattice, symbol timing, window shifts, alignment and the dump format together. */
static void test_glue(void)
{
  printf("\nSensing glue: TDD slots through nrscope_sensing_process_grant, dump read back\n");
  const uint32_t nof_prb = 273, n_sc = nof_prb * SRSRAN_NRE;
  const double   f_d_hz  = 120.0; // Doppler of the target
  const int      d_bins  = 23;    // delay of the direct path, in bins of 1 / (4096 * 30 kHz)
  const int      t_bins  = 41;    // delay of the target
  const float    t_amp   = 0.3f;  // target amplitude against the direct path
  const char*    dir     = "sensing_selftest_out";

  /* the Benetel config's options. symbols is only the cost cap, left at the history
    depth so the window alone decides what a map holds: at 0.30 m/s that is 0.145 s,
    which at 4400 reference symbols/s is about 640 of them. Short windows do not work
    here: with 7 D, 1 S, 2 U the TDD replicas are 200 Hz apart, and a window of a few
    TDD periods is too short for the detector to tell the target from them. */
  nrscope_sensing_default_args(&nrscope_sensing_args);
  nrscope_sensing_args.enable         = true;
  nrscope_sensing_args.symbols        = NR_SENSING_HISTORY_DEPTH;
  nrscope_sensing_args.max_speed_ms   = 10.0; // 0.30 m/s needs 0.145 s, feasible up to ~10 m/s
  nrscope_sensing_args.clutter_kernel = true;
  nrscope_sensing_args.tdd_detect     = true;
  nrscope_sensing_args.antenna_avg    = true; // two chains, one map per layer averaged over them
  snprintf(nrscope_sensing_args.dump, sizeof(nrscope_sensing_args.dump), "%s/map2d.csv", dir);
  nr_sensing_tdd_period_slots = 10; // 7 D, 1 S, 2 U: the Benetel cell, as SIB1 sets it
  char rm[128];
  snprintf(rm, sizeof(rm), "rm -rf %s", dir);
  if (system(rm) != 0) {
    printf("  could not clear %s\n", dir);
  }
  nrscope_sensing_t*         s  = nrscope_sensing_get(3450000000ULL, 122.88e6, 4096, 30000, 2);
  nrscope_sensing_scratch_t* sc = nrscope_sensing_scratch_alloc(n_sc);
  if (s == NULL || sc == NULL) {
    printf("  could not create the context\n");
    failures++;
    return;
  }

  const int warnings_before = nr_sensing_log_warnings();
  srsran_carrier_nr_t carrier = SRSRAN_DEFAULT_CARRIER_NR;
  carrier.pci = 1, carrier.nof_prb = nof_prb, carrier.scs = srsran_subcarrier_spacing_30kHz;
  srsran_dmrs_sch_t dmrs = {};
  srsran_dmrs_sch_init(&dmrs, false);
  srsran_dmrs_sch_set_carrier(&dmrs, &carrier);

  static cf_t tx[273 * SRSRAN_NRE * SRSRAN_NSYMB_PER_SLOT_NR], rx[273 * SRSRAN_NRE * SRSRAN_NSYMB_PER_SLOT_NR];
  /* Chain 1 sees the same scene through a fixed 50 degree offset, as a cable or LO
    path would give it: same delays and Dopplers, so its map matches chain 0's. */
  const cf_t        cable   = cexpf(I * (float)(50.0 * M_PI / 180.0));
  int               pushed1 = 0;
  /* Chain 1 of a slot is handed over LAG downlink slots after its chain 0, as a dozen
    parallel workers deliver them: when a slot's last chain triggers a map, chain 0's
    ring already holds several newer slots than chain 1's. One slot late is not
    enough to show it, as both windows then often hold the same count. */
  enum { LAG = 7 };
  static cf_t rx1q[LAG][273 * SRSRAN_NRE * SRSRAN_NSYMB_PER_SLOT_NR];
  struct {
    srsran_sch_cfg_nr_t cfg;
    uint32_t            sfn, slot_idx;
    int64_t             shift;
  } pend[LAG];
  int pend_n = 0, pend_head = 0;
  const nr_dmrs_placement_t place = {0, 0, 0};
  const cf_t a[2] = {1.0f, 0.6f - 0.4f * I};
  int pushed = 0;
  for (int sl = 0; sl < 1800; sl++) {
    const int phase = sl % 10; // D D D D D D D S U U
    if (phase >= 8) {
      continue;
    }
    srsran_sch_cfg_nr_t cfg = {};
    cfg.dmrs.type = srsran_dmrs_sch_type_1, cfg.dmrs.typeA_pos = srsran_dmrs_sch_typeA_pos_2;
    cfg.dmrs.additional_pos = srsran_dmrs_sch_add_pos_2, cfg.dmrs.length = srsran_dmrs_sch_len_1;
    cfg.grant.mapping = srsran_sch_mapping_type_A, cfg.grant.S = 1, cfg.grant.L = (phase == 7) ? 5 : 13;
    cfg.grant.nof_layers = 2, cfg.grant.nof_dmrs_cdm_groups_without_data = 1, cfg.grant.beta_dmrs = 1.0f;
    cfg.grant.nof_prb = nof_prb;
    for (uint32_t p = 0; p < nof_prb; p++) {
      cfg.grant.prb_idx[p] = true;
    }
    const uint32_t sfn = 100 + sl / 20, slot_idx = sl % 20;
    srsran_slot_cfg_t slot = {.idx = slot_idx};
    memset(tx, 0, sizeof(tx));
    srsran_sch_cfg_nr_t tx_cfg = cfg;
    tx_cfg.grant.nof_layers    = 1; // srsRAN writes port 1000; 1001 is its covered copy
    srsran_dmrs_sch_put_sf(&dmrs, &slot, &tx_cfg, &tx_cfg.grant, tx);
    memset(rx, 0, sizeof(rx));
    const int64_t shift = sl / 30; // window moved one sample later every 30 slots
    for (uint32_t l = 0; l < SRSRAN_NSYMB_PER_SLOT_NR; l++) {
      // same symbol clock as the glue: long CP on symbol 0, normal on the rest
      const uint32_t t_in = (l == 0) ? 352 : 352 + 4096 + (l - 1) * 4384 + 288;
      const double   t    = ((double)(100 * 20 + sl) * 61440 + t_in) / 122.88e6;
      const cf_t     dop  = t_amp * cexpf(I * (float)(2.0 * M_PI * f_d_hz * t));
      for (uint32_t k = 0; k < n_sc; k++) {
        const cf_t x = tx[l * n_sc + k];
        if (x == 0) {
          continue;
        }
        // paths appear `shift` samples earlier, as a ramp about DC (the carrier centre)
        const cf_t  win = cexpf(I * (float)(2.0 * M_PI * shift * ((double)k - n_sc / 2.0) / 4096.0));
        const cf_t  h   = win * (cexpf(-I * (float)(2.0 * M_PI * d_bins * k / 4096.0)) +
                                 dop * cexpf(-I * (float)(2.0 * M_PI * t_bins * k / 4096.0)));
        const float w1 = ((k / 2) % 2 == 1) ? -1.0f : 1.0f;
        rx[l * n_sc + k] = h * x * (a[0] + a[1] * w1);
      }
    }
    const int n = nrscope_sensing_process_grant(s, sc, 0, rx, n_sc, &place, &cfg, 2, 1, sfn, slot_idx, shift);
    pushed += n > 0 ? n : 0;
    if (pend_n == LAG) { // the oldest pending slot's chain 1
      const int n1 = nrscope_sensing_process_grant(s, sc, 1, rx1q[pend_head], n_sc, &place, &pend[pend_head].cfg, 2, 1,
                                                   pend[pend_head].sfn, pend[pend_head].slot_idx, pend[pend_head].shift);
      pushed1 += n1 > 0 ? n1 : 0;
      pend_head = (pend_head + 1) % LAG;
      pend_n--;
    }
    const int q = (pend_head + pend_n) % LAG;
    for (uint32_t i = 0; i < n_sc * SRSRAN_NSYMB_PER_SLOT_NR; i++) {
      rx1q[q][i] = cable * rx[i];
    }
    pend[q].cfg = cfg, pend[q].sfn = sfn, pend[q].slot_idx = slot_idx, pend[q].shift = shift;
    pend_n++;
  }
  for (; pend_n > 0; pend_n--, pend_head = (pend_head + 1) % LAG) {
    const int n1 = nrscope_sensing_process_grant(s, sc, 1, rx1q[pend_head], n_sc, &place, &pend[pend_head].cfg, 2, 1,
                                                 pend[pend_head].sfn, pend[pend_head].slot_idx, pend[pend_head].shift);
    pushed1 += n1 > 0 ? n1 : 0;
  }
  srsran_dmrs_sch_free(&dmrs);
  /* 1800 slots, 0.9 s: 1260 full downlink (3 DM-RS symbols), 180 special (1), 360
    uplink; 2 layers. 22 symbols per 5 ms TDD period, i.e. 4400 per second per layer. */
  check(pushed == 7920, "snapshots pushed, every DM-RS symbol", pushed, 7920, 0);
  check(pushed1 == 7920, "snapshots pushed on chain 1", pushed1, 7920, 0);

  nrscope_sensing_wait_maps(s);

  /* OAI's one-shot snapshot dump, written with the 10th map (NR_SENSING_SNAP_DUMP_MAP):
    one layer's snapshots, direct path and target on their own */
  char path[512];
  snprintf(path, sizeof(path), "%s/map2d.csv.snap.csv", dir);
  FILE* f = fopen(path, "r");
  if (f == NULL) {
    printf("  no snapshot dump written\n");
    failures++;
    return;
  }
  char  line[65536];
  int    n_rows = 0, peak_ok = 0;
  double t_prev = 0, ph_prev[2] = {0, 0}, turn[2] = {0, 0}, dt = 0;
  fgets(line, sizeof(line), f); // header
  while (fgets(line, sizeof(line), f)) {
    double v[11];
    char*  p = line;
    for (int i = 0; i < 11; i++) {
      v[i] = strtod(p, &p);
      p++;
    }
    const int n_bins = (int)v[7];
    double    re[NR_SENSING_MAP_MAX_BINS_RANGE], im[NR_SENSING_MAP_MAX_BINS_RANGE];
    int       best_b = -1;
    double    best   = -1;
    for (int b = 0; b < n_bins; b++) {
      re[b] = strtod(p, &p), p++;
      im[b] = strtod(p, &p), p++;
      if (re[b] * re[b] + im[b] * im[b] > best) {
        best = re[b] * re[b] + im[b] * im[b], best_b = b;
      }
    }
    peak_ok += (best_b == d_bins);
    // phase advance of the direct path and of the target between consecutive snapshots
    const int bins[2] = {d_bins, t_bins};
    for (int j = 0; j < 2; j++) {
      const double ph = atan2(im[bins[j]], re[bins[j]]);
      if (n_rows > 0) {
        double d = ph - ph_prev[j];
        while (d > M_PI) d -= 2 * M_PI;
        while (d < -M_PI) d += 2 * M_PI;
        turn[j] += d;
      }
      ph_prev[j] = ph;
    }
    if (n_rows > 0) {
      dt += v[2] - t_prev;
    }
    t_prev = v[2];
    n_rows++;
  }
  fclose(f);
  const double f_los = dt > 0 ? turn[0] / dt / (2 * M_PI) : 0, f_tgt = dt > 0 ? turn[1] / dt / (2 * M_PI) : 0;
  check(fabs(n_rows - 638) <= 8, "snapshots in the dump", n_rows, 638, 8);
  check(peak_ok == n_rows, "snapshots peaking at the direct path", peak_ok, n_rows, 0);
  check(fabs(f_los) < 1.0, "direct path Doppler (Hz)", f_los, 0, 1.0);
  check(fabs(f_tgt - f_d_hz) < 2.0, "target Doppler (Hz)", f_tgt, f_d_hz, 2.0);

  /* The maps: with the static scene removed, the strongest cell of each is the target,
    at its delay and its Doppler */
  snprintf(path, sizeof(path), "%s/map2d.csv", dir);
  f = fopen(path, "r");
  int n_maps = 0, maps_ok = 0, snap_min = 1 << 30, maps_avg = 0, maps_aoa = 0;
  static char mline[4 << 20];
  if (f != NULL) {
    fgets(mline, sizeof(mline), f); // header
    while (fgets(mline, sizeof(mline), f)) {
      double v[21];
      char*  p = mline;
      for (int i = 0; i < 21; i++) {
        v[i] = strtod(p, &p);
        p++;
      }
      const int    n_bins = (int)v[12], n_freq = (int)v[13];
      snap_min            = (int)v[7] < snap_min ? (int)v[7] : snap_min;
      maps_avg += ((int)v[2] == -1); // aarx -1: averaged over the chains
      maps_aoa += ((int)v[15] > 0);  // n_aoa: the AoA ran on equal windows of both chains
      const double f_max  = v[10];
      int          best_i = 0;
      double       best   = -1;
      for (int i = 0; i < n_bins * n_freq; i++) {
        const double pw = strtod(p, &p);
        p++;
        if (pw > best) {
          best = pw, best_i = i;
        }
      }
      const int    b    = best_i / n_freq;
      const double f_hz = -f_max + (best_i % n_freq) * 2.0 * f_max / n_freq;
      const double cell = 2.0 * f_max / n_freq;
      if (n_maps == 0) {
        printf("  map 0: peak at bin %d, %.1f Hz (cell %.1f Hz), %.1f above the floor\n", b, f_hz, cell, best);
      }
      maps_ok += (abs(b - t_bins) <= 1 && fabs(f_hz - f_d_hz) <= 2 * cell);
      n_maps++;
    }
    fclose(f);
  }
  /* 0.9 s of slots over a 0.145 s window: 6 maps per layer, each one window long,
    back to back. The trigger is the window, so neither the cadence nor the span
    depends on how many symbols the traffic offered. */
  check(n_maps == 12, "maps built (6 per layer)", n_maps, 12, 0);
  // one map per layer per window, each over both chains, not one per chain
  check(maps_avg == n_maps, "maps averaged over the 2 chains", maps_avg, n_maps, 0);
  check(maps_aoa == n_maps, "maps with angles of arrival", maps_aoa, n_maps, 0);
  check(maps_ok == n_maps, "maps peaking at the target", maps_ok, n_maps, 0);
  /* The scene is exactly what the pipeline models, so nothing should warn. A map
    averaged from fewer chains than captured, for one, shows up only as a warning:
    two identical chains give the same map whether one or both went into it. */
  check(nr_sensing_log_warnings() == warnings_before, "warnings logged", nr_sensing_log_warnings() - warnings_before,
        0, 0);
  /* One history per layer, so each map holds everything its own layer collected in
    the window (about 640), not half of a ring shared with the other layer. */
  check(fabs(snap_min - 638) <= 8, "snapshots per map, fewest", snap_min, 638, 8);
  nrscope_sensing_scratch_free(sc);
}

/* Two receive chains of the same symbols get the same alignment, in the order the
 * sniffer feeds them: all of chain 0's DM-RS symbols of a slot, then chain 1's. The
 * correction is common to the chains (one LO, one sample clock); giving chain 1 its
 * own would remove the phase between the chains, which is what the AoA reads. Chain
 * 1 carries a fixed 50 degree offset, as a cable would, so a correction computed on
 * its own symbols would differ from chain 0's by exactly that. OAI remembered only the
 * last symbol, which here would leave two of every three to be recomputed. */
static void test_align_chains(void)
{
  printf("\nAlignment across two chains, chain 0's symbols first\n");
  enum { NP = 819, KS = 4 }; // rank 2 lattice over the full carrier
  const pilot_lattice_t     lat   = {.n = NP, .n_lattice = NP, .k_first = 0, .k_step = KS};
  const nr_sensing_stream_t st    = {.ports = 0x8, .layer = 0}; // a key no other test uses
  const cf_t                cable = cexpf(I * (float)(50.0 * M_PI / 180.0));
  static cf_t               p[NP];
  int                       same = 0, total = 0;
  for (int slot = 0; slot < 8; slot++) {
    nr_sensing_align_t r0[3];
    uint64_t           t[3];
    for (int c = 0; c < 2; c++) {
      for (int m = 0; m < 3; m++) {
        t[m]              = 1000000 + (uint64_t)slot * 61440 + (uint64_t)m * 4 * 4384;
        const float drift = 0.2f * (float)(slot * 3 + m); // common phase drift, as the LO gives
        for (int k = 0; k < NP; k++) {
          p[k] = (c ? cable : 1.0f) * cexpf(I * (drift - (float)(2.0 * M_PI * 23 * k / (4096 / KS))));
        }
        nr_sensing_align_t r;
        nr_ue_sensing_align_symbol(t[m], st, &lat, p, &r);
        if (c == 0) {
          r0[m] = r;
        } else {
          same += (r.phase_rad == r0[m].phase_rad && r.delay_bins == r0[m].delay_bins);
          total++;
        }
      }
    }
  }
  check(same == total, "chain 1 symbols given chain 0's correction", same, total, 0);
}

/* The delay transform must be safe to call from several workers at once, as the
 * sniffer's decoders do. 64 fixed inputs are transformed once on this thread as
 * the reference; 8 threads then transform them 400 times each, interleaved, and
 * every output must match its reference exactly. Before the fix, srsRAN's
 * srsran_dft_run_c() shared the plan's buffers between callers and outputs were
 * mixes of several inputs. */
#define IDFT_T_SIZE 1024
#define IDFT_T_INPUTS 64
static cf_t idft_in[IDFT_T_INPUTS][IDFT_T_SIZE] __attribute__((aligned(32)));
static cf_t idft_ref[IDFT_T_INPUTS][IDFT_T_SIZE] __attribute__((aligned(32)));

static void* idft_worker(void* arg)
{
  const long id  = (long)arg;
  long       bad = 0;
  static __thread cf_t out[IDFT_T_SIZE] __attribute__((aligned(32)));
  for (int it = 0; it < 400; it++) {
    const int i = (int)((id * 7 + it) % IDFT_T_INPUTS);
    nr_ue_sensing_idft(IDFT_T_SIZE, idft_in[i], out);
    if (memcmp(out, idft_ref[i], sizeof(out)) != 0) {
      bad++;
    }
  }
  return (void*)bad;
}

static void test_idft_threads(void)
{
  printf("\nDelay transform from 8 threads at once\n");
  srand(3);
  for (int i = 0; i < IDFT_T_INPUTS; i++) {
    for (int k = 0; k < IDFT_T_SIZE; k++) {
      idft_in[i][k] = (float)rand() / RAND_MAX - 0.5f + I * ((float)rand() / RAND_MAX - 0.5f);
    }
    nr_ue_sensing_idft(IDFT_T_SIZE, idft_in[i], idft_ref[i]);
  }
  pthread_t th[8];
  for (long t = 0; t < 8; t++) {
    pthread_create(&th[t], NULL, idft_worker, (void*)t);
  }
  long bad = 0;
  for (int t = 0; t < 8; t++) {
    void* r;
    pthread_join(th[t], &r);
    bad += (long)r;
  }
  check(bad == 0, "outputs differing from the reference", bad, 0, 0);
}

int main(void)
{
  printf("sensing self-test: %d RB, %d kHz SCS, %.2f MHz carrier, %.3f m/bin\n",
         TEST_N_RB,
         TEST_SCS_HZ / 1000,
         TEST_CARRIER_HZ / 1e6,
         test_m_per_bin());

  const double m_per_bin = test_m_per_bin();

  /* A direct path at zero delay plus one mover. 40 bins is about 98 m of path
  length, well clear of the direct path's mainlobe, and 8 m/s is a walking-pace
  vehicle: fast enough to leave the static ridge, slow enough to need a real
  observation span. */
  const test_path_t scene[] = {
      {.range_bin = 0.0, .speed_ms = 0.0, .amplitude = 1.0},
      {.range_bin = 40.0, .speed_ms = 8.0, .amplitude = 0.05},
  };

  /* 400 slots is 200 ms, giving a Doppler resolution of 5 Hz, which at this
  carrier is about 0.2 m/s. The history holds 512 snapshots and three arrive per
  slot, so the ring wraps; the map takes the most recent window, which is what
  the live path does too. */
  test_run_scene("scene 1: direct path + mover at 40 bins, 8 m/s, clutter removal by mean",
                 scene,
                 2,
                 400,
                 NR_CLUTTER_MEAN,
                 40.0 * m_per_bin,
                 8.0);

  /* The same scene receding, to confirm the Doppler axis is not folded or
  mirrored -- a sign error here would place every target on the wrong side and
  would be invisible in a test that only looked at magnitude. */
  const test_path_t scene_neg[] = {
      {.range_bin = 0.0, .speed_ms = 0.0, .amplitude = 1.0},
      {.range_bin = 25.0, .speed_ms = -12.0, .amplitude = 0.05},
  };
  test_run_scene("scene 2: mover at 25 bins receding at 12 m/s (sign of the Doppler axis)",
                 scene_neg,
                 2,
                 400,
                 NR_CLUTTER_MEAN,
                 25.0 * m_per_bin,
                 -12.0);

  test_dmrs_roundtrip();

  {
    const uint32_t full[][2]   = {{0, 273}};
    const uint32_t gapped[][2] = {{3, 40}, {60, 120}, {200, 273}};
    const uint32_t upper[][2]  = {{33, 273}}; // the allocation above the SSB, seen live at 600 Mb/s
    test_despread_two_layers("full carrier", full, 1, 0);
    test_despread_two_layers("gapped allocation", gapped, 3, 0);
    test_despread_two_layers("PRBs 33-272, grid and BWP from CRB 5", upper, 1, 5);
  }

  test_glue();
  test_align_chains();
  test_idft_threads();

  nr_ue_sensing_idft_free();

  printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
