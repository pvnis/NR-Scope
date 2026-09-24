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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_dmrs_despread.h"
#include "nrscope/hdr/sensing/nr_ue_map.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"
#include "nrscope/hdr/sensing/sensing_defs.h"
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
                                                      map,
                                                      NULL);
  if (n <= 0) {
    printf("  range-doppler produced no map (%d snapshots)\n", n);
    failures++;
    free(map);
    nr_ue_sensing_history_free(&hist);
    return false;
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

  // Receiver: the streaming walk under test
  static cf_t H_est[51 * SRSRAN_NRE];
  static bool H_valid[51 * SRSRAN_NRE];
  memset(H_est, 0, sizeof(H_est));
  memset(H_valid, 0, sizeof(H_valid));

  const int n_written =
      nr_ue_dmrs_estimate_symbol(&lay, 0, cinit, &grant, nof_prb, dmrs_cfg.reference_point_k_rb, rxF, n_sc, H_est,
                                 H_valid);

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

  nr_ue_sensing_idft_free();

  printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
