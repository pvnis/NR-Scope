/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
#include "nrscope/hdr/sensing/nrscope_sensing.h"

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "nrscope/hdr/sensing/nr_ue_aoa.h"
#include "nrscope/hdr/sensing/nr_ue_localize.h"
#include "nrscope/hdr/sensing/nr_ue_map.h"
#include "nrscope/hdr/sensing/nr_ue_music.h"
#include "nrscope/hdr/sensing/nr_ue_sensing.h"
#include "nrscope/hdr/sensing/nr_ue_sensing_align.h"
#include "nrscope/hdr/sensing/nr_ue_tdd_detect.h"
#include "srsran/phy/ch_estimation/dmrs_sch.h"
#include "srsran/phy/utils/vector.h"

#define NSYMB SRSRAN_NSYMB_PER_SLOT_NR
#define SLOTS_PER_FRAME_30KHZ 20
#define SLOTS_PER_HYPERFRAME (1024 * SLOTS_PER_FRAME_30KHZ)
#define MAP_THREADS 2
#define MAP_QUEUE 16

nrscope_sensing_args_t nrscope_sensing_args = {0};

void nrscope_sensing_default_args(nrscope_sensing_args_t* args)
{
  // OAI's defaults (nr-uesoftmodem.h), except symbols, which is 0 there as it doubles as the switch
  memset(args, 0, sizeof(*args));
  args->enable          = false;
  args->symbols         = NR_SENSING_HISTORY_DEPTH;
  args->max_speed_ms    = 30.0;
  args->clutter_removal = true;
}

struct nrscope_sensing_scratch_s {
  uint32_t n_sc;
  cf_t*    H;     // [NSYMB][n_sc]
  bool*    valid; // [NSYMB][n_sc]
};

nrscope_sensing_scratch_t* nrscope_sensing_scratch_alloc(uint32_t n_sc_grid)
{
  nrscope_sensing_scratch_t* sc = calloc(1, sizeof(*sc));
  if (sc == NULL) {
    return NULL;
  }
  sc->n_sc  = n_sc_grid;
  sc->H     = srsran_vec_cf_malloc(NSYMB * n_sc_grid);
  sc->valid = calloc((size_t)NSYMB * n_sc_grid, sizeof(bool));
  if (sc->H == NULL || sc->valid == NULL) {
    nrscope_sensing_scratch_free(sc);
    return NULL;
  }
  return sc;
}

void nrscope_sensing_scratch_free(nrscope_sensing_scratch_t* sc)
{
  if (sc != NULL) {
    free(sc->H);
    free(sc->valid);
    free(sc);
  }
}

struct nrscope_sensing_s {
  nrscope_sensing_args_t args;
  /* One history per stream, i.e. per (DM-RS ports, layer) on chain 0, each
    NR_SENSING_HISTORY_DEPTH deep: at rank 2 each layer holds that many snapshots of
    its own. The comb is fixed by the ports, so (ports, layer) is the
    stream */
  struct {
    bool                 used;
    uint16_t             ports;
    uint8_t              layer;
    nr_sensing_history_t h;
    /* Where this stream's last map ended, on the OFDM sample clock, and whether a
      window has been opened at all. The map trigger is this time plus one window,
      never a snapshot count: the resolution a map reaches is lambda / (2 t_span),
      so the window is the requirement and the count is whatever the scheduler
      granted inside it. */
    bool                 armed;
    uint64_t             t_map_end;
  } hist[NR_AOA_MAX_ANT][NR_SENSING_MAX_STREAMS];
  pthread_mutex_t        hist_lock; // creation only; each history has its own lock
  /* Receive chains feeding the rings, 1 .. NR_AOA_MAX_ANT. The chain indexes the
    table above and is not part of the stream key: the map task gathers one
    unchanged key from every chain's ring, which is what lets the per-chain maps be
    averaged cell by cell and the AoA compare the same beam across the array. */
  uint32_t               nof_antennas;
  uint64_t               carrier_hz;
  uint32_t               scs_hz;
  uint32_t               samples_per_slot;
  uint32_t               ofdm_size;

  /* Absolute slot count across the SFN wrap. Slots arrive from several workers,
    slightly out of order, so each one is placed in the hyperframe that puts it
    nearest the latest slot seen, which is unambiguous for any disorder under
    10.24 s. */
  pthread_mutex_t lock;
  int64_t         last_slot_abs;

  /* Maps are built off the workers, by MAP_THREADS threads fed from a queue: OAI
    pushes the same task to its thread pool. A full queue drops the new map. */
  pthread_t       map_thr[MAP_THREADS];
  pthread_mutex_t q_lock;
  pthread_cond_t  q_cond;
  void*           q[MAP_QUEUE]; // nr_sensing_map_task_t*
  int             q_head, q_len;
  int             q_busy; // tasks being built
  pthread_cond_t  q_idle;
};

static void* map_thread(void* arg);

/* True when this stream's window is complete, i.e. its newest snapshot is a full
   window past where the last map ended, and claims the window so that only one worker
   builds that map. The window is nr_ue_sensing_window_s(), the same call the transform
   gathers with, so the cadence and the window cannot drift apart.

   Workers decode slots in parallel and so arrive slightly out of order; a snapshot
   older than the current window's start is simply early for the next one. */
static bool map_due(nrscope_sensing_t* s, const nr_sensing_history_t* hist, int aarx, uint16_t ports, int layer,
                    uint64_t t_now)
{
  // Time window of the map t_span required from the different parameters
  // Initially values 0,0 for the last two arguments because we did not compute already
  // the grant size and the comb
  const double window_s = nr_ue_sensing_window_s(hist, s->args.max_speed_ms, 0, 0);

  const double fs       = (double)s->ofdm_size * s->scs_hz;
  if (!(window_s > 0.0) || !(fs > 0.0)) {
    return false;
  }

  const int64_t window_samples = (int64_t)(window_s * fs);

  bool due = false;
  pthread_mutex_lock(&s->hist_lock);
  for (int i = 0; i < NR_SENSING_MAX_STREAMS; i++) {
    if (!s->hist[aarx][i].used || s->hist[aarx][i].ports != ports || s->hist[aarx][i].layer != layer) {
      continue;
    }
    if (!s->hist[aarx][i].armed) {
      // the first snapshot of this stream opens the first window
      s->hist[aarx][i].armed     = true;
      s->hist[aarx][i].t_map_end = t_now;
    } else if ((int64_t)(t_now - s->hist[aarx][i].t_map_end) >= window_samples) {
      s->hist[aarx][i].t_map_end = t_now;
      due                        = true;
    }
    break;
  }
  pthread_mutex_unlock(&s->hist_lock);
  return due;
}

/* The history of one stream, created on first use. NULL when every slot is taken by
   other streams, which needs more distinct DM-RS port sets than the scheduler uses. */
static nr_sensing_history_t* hist_of(nrscope_sensing_t* s, int aarx, uint16_t ports, int layer)
{
  if (aarx < 0 || aarx >= (int)s->nof_antennas) {
    return NULL;
  }
  nr_sensing_history_t* h = NULL;
  pthread_mutex_lock(&s->hist_lock);
  int free_i = -1;
  for (int i = 0; i < NR_SENSING_MAX_STREAMS && h == NULL; i++) {
    if (!s->hist[aarx][i].used) {
      if (free_i < 0)
        free_i = i;
    } else if (s->hist[aarx][i].ports == ports && s->hist[aarx][i].layer == layer) {
      h = &s->hist[aarx][i].h;
    }
  }
  if (h == NULL && free_i >= 0 &&
      nr_ue_sensing_history_init(&s->hist[aarx][free_i].h, NR_SENSING_HISTORY_DEPTH, (int)s->scs_hz,
                                 (int)s->ofdm_size, s->carrier_hz)) {
    s->hist[aarx][free_i].used  = true;
    s->hist[aarx][free_i].ports = ports;
    s->hist[aarx][free_i].layer = (uint8_t)layer;
    h                           = &s->hist[aarx][free_i].h;
  }
  pthread_mutex_unlock(&s->hist_lock);
  return h;
}

static nrscope_sensing_t* g_ctx      = NULL;
static pthread_mutex_t    g_ctx_lock = PTHREAD_MUTEX_INITIALIZER;

nrscope_sensing_t* nrscope_sensing_get(uint64_t carrier_hz, double srate_hz, uint32_t ofdm_size, uint32_t scs_hz,
                                       uint32_t nof_antennas)
{
  if (!nrscope_sensing_args.enable) {
    return NULL;
  }
  pthread_mutex_lock(&g_ctx_lock);
  if (g_ctx == NULL) {
    nrscope_sensing_t* s = calloc(1, sizeof(*s));
    if (s != NULL && nrscope_sensing_args.symbols > NR_SENSING_HISTORY_DEPTH) {
      // OAI asserts the same in nr-uesoftmodem.c: a map cannot hold more than a history
      LOG_E(NR_PHY, "sensing-symbols %d exceeds the history depth %d\n", nrscope_sensing_args.symbols,
            NR_SENSING_HISTORY_DEPTH);
      free(s);
      s = NULL;
    }
    if (s != NULL) {
      s->args             = nrscope_sensing_args;
      s->nof_antennas     = nof_antennas < 1 ? 1 : (nof_antennas > NR_AOA_MAX_ANT ? NR_AOA_MAX_ANT : nof_antennas);
      s->carrier_hz       = carrier_hz;
      s->scs_hz           = scs_hz;
      pthread_mutex_init(&s->hist_lock, NULL);
      s->samples_per_slot = (uint32_t)(srate_hz / 1000.0 / (scs_hz / 15000));
      s->ofdm_size        = ofdm_size;
      s->last_slot_abs    = -1;
      pthread_mutex_init(&s->lock, NULL);
      if (s->args.dump[0] != 0) {
        // the dump's directory, one level, as OAI's /tmp needed none
        char dir[sizeof(s->args.dump)];
        snprintf(dir, sizeof(dir), "%s", s->args.dump);
        char* slash = strrchr(dir, '/');
        if (slash != NULL && slash != dir) {
          *slash = 0;
          mkdir(dir, 0755);
        }
      }
      pthread_mutex_init(&s->q_lock, NULL);
      pthread_cond_init(&s->q_cond, NULL);
      pthread_cond_init(&s->q_idle, NULL);
      for (int i = 0; i < MAP_THREADS; i++) {
        pthread_create(&s->map_thr[i], NULL, map_thread, s);
        pthread_setname_np(s->map_thr[i], "sensing_map");
        pthread_detach(s->map_thr[i]);
      }
      g_ctx = s;
    } else {
      free(s);
      LOG_E(NR_PHY, "sensing: could not create the context; sensing stays off\n");
    }
  }
  pthread_mutex_unlock(&g_ctx_lock);
  return g_ctx;
}

uint16_t nrscope_dmrs_ports_type1_len1(uint32_t value, uint32_t* cdm_groups)
{
  // TS 38.212 table 7.3.1.2.2-1: {CDM groups without data, DM-RS ports as bitmap}
  static const struct {
    uint8_t  cdm;
    uint16_t ports;
  } t[16] = {{1, 0x1}, {1, 0x2}, {1, 0x3}, {2, 0x1}, {2, 0x2}, {2, 0x4}, {2, 0x8}, {2, 0x3},
             {2, 0xc}, {2, 0x7}, {2, 0xf}, {2, 0x5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
  if (value >= 16 || t[value].ports == 0) {
    return 0;
  }
  if (cdm_groups != NULL) {
    *cdm_groups = t[value].cdm;
  }
  return t[value].ports;
}

/* Absolute slot index, unwrapped across the 10.24 s SFN period. */
static uint64_t slot_abs_of(nrscope_sensing_t* s, uint32_t sfn, uint32_t slot_idx)
{
  const int64_t in_period = (int64_t)sfn * SLOTS_PER_FRAME_30KHZ + slot_idx;
  pthread_mutex_lock(&s->lock);
  int64_t abs_slot = in_period;
  if (s->last_slot_abs >= 0) {
    // floor division: the nearest hyperframe, also for a slot from just before a wrap
    const int64_t num = s->last_slot_abs - in_period + SLOTS_PER_HYPERFRAME / 2;
    const int64_t hf  = (num >= 0) ? num / SLOTS_PER_HYPERFRAME : -((-num + SLOTS_PER_HYPERFRAME - 1) / SLOTS_PER_HYPERFRAME);
    abs_slot         = in_period + hf * SLOTS_PER_HYPERFRAME;
  }
  if (abs_slot > s->last_slot_abs) {
    s->last_slot_abs = abs_slot;
  }
  pthread_mutex_unlock(&s->lock);
  return (uint64_t)abs_slot;
}

/* Sample time of the start of OFDM symbol l's FFT window, on the capture clock.
  At 30 kHz the first symbol of each slot has the long cyclic prefix (a slot is
  half a subframe); the others the normal one. Only differences matter to the
  Doppler transform, so the origin is the capture's first slot. */
static uint64_t symbol_time(uint64_t slot_abs, uint32_t l, uint32_t samples_per_slot, uint32_t N)
{
  const uint32_t cp      = N * 144 / 2048;                // normal CP
  const uint32_t cp0     = cp + (samples_per_slot - 14 * (N + cp)); // the long one takes the remainder
  const uint32_t t_in    = (l == 0) ? cp0 : cp0 + N + (l - 1) * (N + cp) + cp;
  return slot_abs * samples_per_slot + t_in;
}

typedef struct {
  const nrscope_sensing_args_t *args;
  /// Rx antenna the map belongs to, -1 when it is the average of all of them
  int aarx;
  int frame;
  int slot;
  /// the measurement stream every history in snap[] is to be gathered for
  nr_sensing_stream_t stream;
  /// histories in snap[], 1 unless the antennas are averaged together
  int n_ant;
  /// layers of the allocation summed into this map
  int n_layers;
  nr_sensing_map_t map;
  /* The history keeps being written by the DL actors while this runs, so the task
  owns a copy taken at trigger time rather than a pointer into the ring. */
  nr_sensing_history_t snap[];
} nr_sensing_map_task_t;

/* Copy the samples of one slow-time set into another whose buffers are already
allocated, the buffers themselves staying where they are. */
static void nr_ue_sensing_slowtime_copy(nr_sensing_slowtime_t *dst, const nr_sensing_slowtime_t *src)
{
  const size_t n_snap = src->n_snap > 0 ? (size_t)src->n_snap : 0;
  const size_t n_val = n_snap * (src->n_bins > 0 ? (size_t)src->n_bins : 0);
  memcpy(dst->t_s, src->t_s, n_snap * sizeof(*dst->t_s));
  memcpy(dst->h_re, src->h_re, n_val * sizeof(*dst->h_re));
  memcpy(dst->h_im, src->h_im, n_val * sizeof(*dst->h_im));
  double *t_s = dst->t_s;
  float *h_re = dst->h_re;
  float *h_im = dst->h_im;
  *dst = *src;
  dst->t_s = t_s;
  dst->h_re = h_re;
  dst->h_im = h_im;
}

/* The map of one task: either antenna 0 layer 0, averaged with the other antennas and layers
when the task holds several. slow[a] receives the layer 0 slow-time samples of antenna
a for a < n_slow (slow may be NULL when n_slow is 0). max_paths is passed on to
nr_ue_sensing_range_doppler(). n_combined, when not NULL, receives how many maps went
into the average. Returns the snapshots of the first map, 0 if too few.

obs, when not NULL, receives the slow-time samples of every (antenna a, layer j) at
obs[a * n_layers + j], for the MUSIC map; an entry whose map could not be built is left
with n_snap 0. It is filled from the same calls that build the DFT map, so the two maps
see exactly the same conditioned samples.

Split out so --sensing-clutter-compare can build the same map twice from the same
snapshots, once per max_paths. */
static int nr_ue_sensing_task_map(const nr_sensing_map_task_t *t,
                                  nr_sensing_clutter_t clutter,
                                  int max_paths,
                                  nr_sensing_map_t *map,
                                  nr_sensing_slowtime_t *slow,
                                  int n_slow,
                                  int *n_combined_out,
                                  nr_sensing_slowtime_t *obs)
{
  if (obs != NULL)
    for (int k = 0; k < t->n_ant * t->n_layers; k++)
      obs[k].n_snap = 0;

  int n = nr_ue_sensing_range_doppler(&t->snap[0],
                                      &t->stream,
                                      t->args->symbols,
                                      t->args->max_speed_ms,
                                      clutter,
                                      max_paths,
                                      map,
                                      obs != NULL ? &obs[0] : (n_slow > 0 ? &slow[0] : NULL));
  if (obs != NULL && n > 0 && n_slow > 0)
    nr_ue_sensing_slowtime_copy(&slow[0], &obs[0]);
  if (obs != NULL && n <= 0)
    obs[0].n_snap = 0;

  int n_combined = n > 0 ? 1 : 0;
  if (n > 0 && (t->n_ant > 1 || t->n_layers > 1)) {
    nr_sensing_map_t *other = malloc_or_fail(sizeof(*other));

    // We compute the average map for all layers for this particular run
    for (int a = 0; a < t->n_ant; a++) {
      for (int j = 0; j < t->n_layers; j++) {
        if (a == 0 && j == 0) continue; // done for antenna 0 by default above

        /* The layers this map is built from: all of them for a layer-averaged map,
        else the map's own. OAI set j here in both cases, so without layer
        averaging a layer 1 map took layer 0 from every chain but the first: an
        average of two different beams, and an AoA comparing H w_1 on rx0 with
        H w_0 on the others. */
        nr_sensing_stream_t st = t->stream;
        st.layer = t->n_layers > 1 ? (uint8_t)j : t->stream.layer;

        /* slow[] feeds the AoA, which builds a spatial covariance across the array
        and so assumes one transmit beam per entry. Only layer 0 fills it. */
        nr_sensing_slowtime_t *sl = (j == 0 && a < n_slow) ? &slow[a] : NULL;
        nr_sensing_slowtime_t *ob = obs != NULL ? &obs[a * t->n_layers + j] : NULL;

        /* Here we handle the fact that we do not necessarily use all layers. In fact,
        the scheduler can maybe not schedule layer 1 in the window. */
        if (nr_ue_sensing_range_doppler(&t->snap[a], &st, t->args->symbols, t->args->max_speed_ms,
                                        clutter, max_paths, other, ob != NULL ? ob : sl) <= 0) {
          if (ob != NULL)
            ob->n_snap = 0;
          continue;
        }
        if (ob != NULL && sl != NULL)
          nr_ue_sensing_slowtime_copy(sl, ob);
        /* Each chain sizes its own grid from the window it gathered, so a chain whose
        window came out a snapshot shorter can land on a different n_freq (379 against
        383). Summed cell by cell, its rows would then be read with the wrong stride,
        and past its end when it is the smaller one: the last range bin filled with
        whatever the buffer held. Such a map is left out of the average. */
        if (other->n_bins != map->n_bins || other->n_freq != map->n_freq || other->f_max_hz != map->f_max_hz) {
          LOG_W(NR_PHY,
                "sensing: rx%d layer %d map is %d x %d against %d x %d on the first, left out of the average\n",
                a,
                j,
                other->n_bins,
                other->n_freq,
                map->n_bins,
                map->n_freq);
          continue;
        }
        for (int i = 0; i < map->n_bins * map->n_freq; i++)
          map->power[i] += other->power[i];
        n_combined++;
      }
    }
    free(other);
    /* Mean rather than sum, so the dB scale of a map does not depend on how many
    antennas fed it and captures stay comparable. */
    if (n_combined > 1) {
      for (int i = 0; i < map->n_bins * map->n_freq; i++)
        map->power[i] /= n_combined;
    }

    /* The map no longer describes one beam, so it must not keep claiming to.
    range_doppler() set stream from the first snapshot it gathered, i.e. layer 0, the
    same way aarx would read as antenna 0 without the -1 that marks an average. */
    if (t->n_layers > 1)
      map->stream.layer = NR_SENSING_LAYER_AVG;
  }
  if (n_combined_out != NULL)
    *n_combined_out = n_combined;
  return n;
}

static void nr_ue_sensing_map_task(void *arg)
{
  nr_sensing_map_task_t *t = (nr_sensing_map_task_t *)arg;

  /* Runtime parameter: clutter removal for static objects. --sensing-clutter-removal
  turns it on and --sensing-clutter-kernel picks the mode, so the kernel flag alone
  does nothing and removal stays a single switch. */
  const nr_sensing_clutter_t clutter = !t->args->clutter_removal ? NR_CLUTTER_NONE
                                       : t->args->clutter_kernel ? NR_CLUTTER_KERNEL
                                                                       : NR_CLUTTER_MEAN;

  /* Runtime parameter --sensing-tdd-detect: run the TDD detector (the TDD paper).
  It finds the targets and their true speed, even beyond the map's Doppler window, by
  telling the real peaks from the TDD replicas. It works on rx0 only.
  When on, its targets are the ones the AoA measures. When off, the AoA finds its own
  peaks on the averaged map, and speeds are only right inside the map's window. */
  const bool detect = t->args->tdd_detect;

  /* Runtime parameter --sensing-antenna-avg: average the maps of all Rx antennas and
  estimate the angle of arrival. The AoA needs at least 2 antennas. */
  const bool aoa = t->args->antenna_avg && t->n_ant > 1;

  /* Slow-time samples: rx0 alone feeds the detector, every antenna feeds the AoA. */
  int n_slow = aoa ? t->n_ant : (detect ? 1 : 0);

  // Prepare the buffers for storing the slow-time axis responses
  nr_sensing_slowtime_t slow[NR_AOA_MAX_ANT] = {0};
  for (int a = 0; a < n_slow; a++) {
    slow[a].t_s = malloc_or_fail((size_t)t->snap[a].depth * sizeof(*slow[a].t_s));
    slow[a].h_re = malloc_or_fail((size_t)t->snap[a].depth * NR_SENSING_MAP_MAX_BINS_RANGE * sizeof(*slow[a].h_re));
    slow[a].h_im = malloc_or_fail((size_t)t->snap[a].depth * NR_SENSING_MAP_MAX_BINS_RANGE * sizeof(*slow[a].h_im));
  }

  // Computes the sensing map for the first antenna (see t->snap[0]) independently of the 
  // runtime parameters.

  /* One-shot dump of the raw slow-time samples this map was built from, for checking
  the delay-response model offline. Counted per map rather than per slot so it fires
  once a run, on a window that is already full. The path is derived from the map dump
  so it needs no option of its own. 
  
  This is used for model checking, using the model_checker.py script. Will be removed afterwards
  */
  static int map_count = 0;
  if (__sync_add_and_fetch(&map_count, 1) == NR_SENSING_SNAP_DUMP_MAP && t->args->dump[0] != 0) {
    char snap_path[512];
    snprintf(snap_path, sizeof(snap_path), "%s.snap.csv", t->args->dump);
    nr_ue_sensing_dump_snapshots(snap_path, &t->snap[0], &t->stream, t->args->symbols,
                                 t->args->max_speed_ms);
  }

  /*
  Per antenna and per layer, we compute the 2d map. We average the maps power, but only if
  nr_ue_sensing_range_doppler() returns a positive number. If it returns 0, it was because
  it had not enough snapshots to compute the map. We use only layer 0 for the AoA because mixing
  multiple layers for AoA is wrong because of precoding. Mixing layer 0 at rx0 with layer 1 at rx1
  would give you H·w₀ and H·w₁ — the receive steering a_R(θ_p) is identical, but each 
  carries its own a_T^H(φ_p)·w_j factor, so the phase difference between entries is 
  steering plus an arbitrary precoder rotation
  */

  /* Runtime parameter --sensing-music: a second map of the same window, with Doppler MUSIC
  in place of the slow-time DFT, dumped to <sensing-dump>.music.csv. It needs the
  conditioned samples of every antenna and layer, which are its snapshots. */
  const bool music = t->args->music && t->args->dump[0] != 0;
  int n_obs = music ? t->n_ant * t->n_layers : 0;
  if (n_obs > NR_MUSIC_MAX_OBS) {
    LOG_W(NR_PHY, "sensing music: %d antenna x layer observations, only %d used\n", n_obs, NR_MUSIC_MAX_OBS);
    n_obs = NR_MUSIC_MAX_OBS;
  }
  // task_map indexes obs[a * n_layers + j] over every pair, so the buffers cover them all
  const int n_obs_alloc = music ? t->n_ant * t->n_layers : 0;
  nr_sensing_slowtime_t *obs = music ? calloc_or_fail((size_t)n_obs_alloc, sizeof(*obs)) : NULL;
  for (int k = 0; k < n_obs_alloc; k++) {
    const int a = k / t->n_layers;
    obs[k].t_s = malloc_or_fail((size_t)t->snap[a].depth * sizeof(*obs[k].t_s));
    obs[k].h_re = malloc_or_fail((size_t)t->snap[a].depth * NR_SENSING_MAP_MAX_BINS_RANGE * sizeof(*obs[k].h_re));
    obs[k].h_im = malloc_or_fail((size_t)t->snap[a].depth * NR_SENSING_MAP_MAX_BINS_RANGE * sizeof(*obs[k].h_im));
  }

  int n_combined = 0;
  int n = nr_ue_sensing_task_map(t, clutter, NR_CLUTTER_MAX_PATHS, &t->map, slow, n_slow, &n_combined, obs);

  nr_sensing_marker_t markers[NR_TDD_MAX_TARGETS + NR_TDD_MAX_REJECTED];
  int n_markers = 0;

  // Targets of the TDD detector, handed to the AoA below when the detector is on
  nr_tdd_target_t tgt[NR_TDD_MAX_TARGETS];
  int n_tgt = 0;

  // -----------------------------------------------------------------------------

  /* TDD detector, on the same Doppler grid as the map. Same grid on purpose: a verdict
  then lands on a map cell without any conversion, and there is only one grid in the
  pipeline to reason about. The map's axis already carries one TDD period of margin
  past the requested max speed, which is what the detector needs to look at where the
  replicas of a candidate would be.
  
  The reason why they are replicas is the following: Because of TDD pattern, RS arrive
  only at time t = n * T_TDD + delta_r, where delta indexes where the symbols are in the slot.
  In our setup, supposing all symbols in the ring buffer that are used have the same n_m and c_m,
  then its h_m[b] = \sum_l beta_l(b) e^{j 2 pi f_l (n * T_TDD + delta_r)}, where t_m = n * T_TDD + delta_r
  Now by taking the DFT we get (see math) 2 sums where one of them creates replicas every
  f = k / T_TDD, which is at 5ms of period is every 200 Hz.
  */
  if (detect && n > 0 && slow[0].n_snap >= 8) {
    const double slot_dur_s = 1e-3 / ((double)t->snap[0].scs_hz / 15000.0);
    const double t_tdd_s = NR_SENSING_TDD_PERIOD_SLOTS * slot_dur_s;
    const double t_span = slow[0].t_s[slow[0].n_snap - 1] - slow[0].t_s[0];

    float *hann = malloc_or_fail((size_t)slow[0].win_len * sizeof(*hann));
    if (!nr_ue_sensing_hann_coeffs(slow[0].win_len, hann))
      for (int i = 0; i < slow[0].win_len; i++)
        hann[i] = 1.0f; // degenerate lattice: the transform was left rectangular too

    const double lambda_m = C_M_PER_S / t->map.carrier_hz;
    const double f_max_det = t->map.f_max_hz;
    const int n_freq_det = t->map.n_freq;

    /* The subspace the slow-trend filter removed from these samples, so the detector
    models a target the way the data holds it. Same basis, from the same function, as
    the filter itself; laid out [n_q][n_snap] for the detector. */
    double q_basis[NR_CLUTTER_SLOW_TREND_MAX + 1][NR_SENSING_HISTORY_DEPTH];
    const int trend_n_q = nr_ue_sensing_slow_basis(slow[0].t_s, slow[0].n_snap, slow[0].trend_degree, q_basis);
    double *trend_q = trend_n_q > 0 ? malloc_or_fail((size_t)trend_n_q * slow[0].n_snap * sizeof(*trend_q)) : NULL;
    for (int k = 0; k < trend_n_q; k++)
      memcpy(&trend_q[(size_t)k * slow[0].n_snap], q_basis[k], (size_t)slow[0].n_snap * sizeof(*trend_q));

    const nr_tdd_obs_t obs = {
        .n_snap = slow[0].n_snap,
        .t_s = slow[0].t_s,
        .h_re = slow[0].h_re,
        .h_im = slow[0].h_im,
        .n_bins = slow[0].n_bins,
        .idft_size = slow[0].idft_size,
        .win_len = slow[0].win_len,
        .win_start = slow[0].win_start,
        .win = hann,
        .m_per_bin = slow[0].m_per_bin,
        .n_freq = n_freq_det,
        .f_max_hz = f_max_det,
        .lambda_m = lambda_m,
        .t_tdd_s = t_tdd_s,
        .trend_q = trend_q,
        .trend_n_q = trend_n_q,
    };

    nr_tdd_cfg_t cfg;
    nr_tdd_cfg_default(&cfg);
    nr_tdd_target_t rej[NR_TDD_MAX_REJECTED];
    int n_rej = 0;
    n_tgt = nr_tdd_detect(&obs, &cfg, tgt, NR_TDD_MAX_TARGETS, rej, &n_rej);

    /* Set aside the targets whose mirror at -f is about as strong, see
    NR_TDD_VERDICT_UNCERTAIN. Read on the averaged map, the same grid the detector ran
    on. Each side is the strongest cell within one bin and one Doppler cell, so an
    off-grid peak is not compared against the skirt of its own mirror. A target within a
    cell of zero Doppler is its own mirror and is kept. Set aside rather than dropped:
    they travel with the map as markers, but the AoA and the localisation never see them. */
    {
      const int nf = t->map.n_freq;
      const int nb = t->map.n_bins;
      const double df = nf > 1 ? 2.0 * t->map.f_max_hz / (nf - 1) : 0.0;
      int n_kept = 0;
      for (int i = 0; i < n_tgt; i++) {
        const int b0 = (int)lround(tgt[i].bin);
        const int f0 = df > 0.0 ? (int)lround((tgt[i].f_hz + t->map.f_max_hz) / df) : -1;
        const int fm = nf - 1 - f0;
        bool uncertain = false;
        if (b0 >= 0 && b0 < nb && f0 >= 0 && f0 < nf && abs(f0 - fm) > 2) {
          double p_own = 0.0, p_mir = 0.0;
          for (int b = b0 - 1; b <= b0 + 1; b++) {
            if (b < 0 || b >= nb)
              continue;
            for (int d = -1; d <= 1; d++) {
              if (f0 + d >= 0 && f0 + d < nf && t->map.power[b * nf + f0 + d] > p_own)
                p_own = t->map.power[b * nf + f0 + d];
              if (fm + d >= 0 && fm + d < nf && t->map.power[b * nf + fm + d] > p_mir)
                p_mir = t->map.power[b * nf + fm + d];
            }
          }
          uncertain = p_mir > 0.0 && p_own < p_mir * pow(10.0, NR_TDD_MIRROR_DB / 10.0);
        }
        if (uncertain) {
          if (n_rej < NR_TDD_MAX_REJECTED) { // its marker only; set aside either way
            rej[n_rej] = tgt[i];
            rej[n_rej++].verdict = NR_TDD_VERDICT_UNCERTAIN;
          }
        } else {
          tgt[n_kept++] = tgt[i];
        }
      }
      n_tgt = n_kept;
    }

    // LOG_I(NR_PHY,
    //       "SENSING DETECT %d.%d rx%d ports 0x%x layer %d: %d target(s) from %d snapshots, "
    //       "grid +-%.0f Hz x %d, replicas every %.0f Hz (%.1f m/s)\n",
    //       t->frame,
    //       t->slot,
    //       t->aarx,
    //       t->stream.ports,
    //       t->stream.layer,
    //       n_tgt,
    //       slow[0].n_snap,
    //       f_max_det,
    //       n_freq_det,
    //       1.0 / t_tdd_s,
    //       (1.0 / t_tdd_s) * obs.lambda_m / 2.0);
    // for (int i = 0; i < n_tgt; i++)
    //   LOG_I(NR_PHY,
    //         "SENSING TARGET %d.%d rx%d: %.1f m (bin %.2f), %+.2f m/s (%+.1f Hz), %.1f dB, iter %d\n",
    //         t->frame,
    //         t->slot,
    //         t->aarx,
    //         tgt[i].range_m,
    //         tgt[i].bin,
    //         tgt[i].speed_ms,
    //         tgt[i].f_hz,
    //         tgt[i].snr_dB,
    //         tgt[i].iteration);

    /* Both verdicts travel with the map. A plot showing only the survivors looks the
    same whether the detector ran or not; the rejections are what make it visible. */
    for (int i = 0; i < n_tgt && n_markers < (int)(sizeof(markers) / sizeof(markers[0])); i++)
      markers[n_markers++] = (nr_sensing_marker_t){.range_m = (float)tgt[i].range_m,
                                                   .speed_ms = (float)tgt[i].speed_ms,
                                                   .snr_dB = (float)tgt[i].snr_dB,
                                                   .verdict = tgt[i].verdict};
    for (int i = 0; i < n_rej && i < NR_TDD_MAX_REJECTED && n_markers < (int)(sizeof(markers) / sizeof(markers[0]));
         i++)
      markers[n_markers++] = (nr_sensing_marker_t){.range_m = (float)rej[i].range_m,
                                                   .speed_ms = (float)rej[i].speed_ms,
                                                   .snr_dB = (float)rej[i].snr_dB,
                                                   .verdict = rej[i].verdict};
    free(hann);
    free(trend_q);
  }

  // ----------------------------------------------------------------------------------

  /* Angle of arrival at the strongest cells of the map.

  Last in the pipeline on purpose: by this point the delay and Doppler transforms have
  integrated each cell non-coherently and the clutter mean has taken the direct path out
  of it, so the four numbers the array contributes are at the best SNR they ever
  reach and are not dominated by the LoS steering 

  nr_ue_aoa_process() also reports the speed of each target. With the TDD detector on,
  that is the true speed. Without it, it is only right inside the map's Doppler window.
  */

  nr_sensing_aoa_t aoa_out[NR_SENSING_AOA_MAX];
  int n_aoa = 0;
  if (aoa && n > 0) {
    /* Where the AoA takes its targets from: the TDD detector's targets when it is on.
    When it is off we pass NULL, and the AoA searches the averaged map itself with CFAR. */
    nr_aoa_cell_t cells[NR_TDD_MAX_TARGETS];
    for (int i = 0; i < n_tgt; i++)
      cells[i] = (nr_aoa_cell_t){.bin = tgt[i].bin, .f_hz = tgt[i].f_hz, .snr_dB = tgt[i].snr_dB};

    n_aoa = nr_ue_aoa_process(&t->map, slow, n_slow, detect ? cells : NULL, n_tgt, aoa_out);
    // for (int i = 0; i < n_aoa; i++)
    //   LOG_I(NR_PHY,
    //         "SENSING AOA %d.%d ports 0x%x layer %d: %.1f m, %+.2f m/s, %+.1f deg, SNR %.1f dB, quality %.1f dB\n",
    //         t->frame,
    //         t->slot,
    //         t->stream.ports,
    //         t->stream.layer,
    //         aoa_out[i].range_m,
    //         aoa_out[i].speed_ms,
    //         aoa_out[i].angle_deg,
    //         aoa_out[i].power_dB,
    //         aoa_out[i].quality_dB);
  }

  /* Bistatic localisation of the cells that got an angle. Guarded rather than
  returned from: the map dump below and the frees at the end of this task run
  whether or not any cell cleared the AoA threshold. */
  if (n_aoa > 0) {
    /* TODO placeholder geometry. Until these come from the UE configuration the
    positions below are made up and so is every coordinate derived from them. */
    const float gnb[2] = {0.0f, 0.0f};
    const float ue[2] = {0.0f, 0.0f};

    // This is the angle of the ULA with respect to the reference frame of the map
    const float boresight_deg = 0.0f;

    nr_sensing_target_t targets[NR_SENSING_AOA_MAX];
    const int n_tgt_pos = nr_ue_target_position(aoa_out, n_aoa, gnb, ue, boresight_deg, targets);
    // for (int i = 0; i < n_tgt_pos; i++)
    //   LOG_I(NR_PHY,
    //         "SENSING POS %d.%d ports 0x%x layer %d: (%.1f, %.1f) m, %.1f m from the UE, "
    //         "from %.1f m %+.1f deg\n",
    //         t->frame,
    //         t->slot,
    //         t->stream.ports,
    //         t->stream.layer,
    //         targets[i].pos_x,
    //         targets[i].pos_y,
    //         targets[i].rho,
    //         aoa_out[targets[i].aoa_index].range_m,
    //         aoa_out[targets[i].aoa_index].angle_deg);
  }

  /* Dump the all data for plotting everything */

  if (n > 0) {
    char rx_label[8];
    if (t->aarx < 0)
      snprintf(rx_label, sizeof(rx_label), "avg");
    else
      snprintf(rx_label, sizeof(rx_label), "%d", t->aarx);
    nr_ue_sensing_dump_map(t->args->dump[0] ? t->args->dump : NULL, t->frame, t->slot, t->aarx, &t->map, markers, n_markers, aoa_out, n_aoa);

    /* The MUSIC map of the same window, on the same grid, written to its own file one line
    per map in the same order as the DFT dump, so line i of both files is the same window.
    The markers of the TDD detector are copied along: they were found on the samples, not
    on the DFT map, and show where the real targets are on either picture. Always written,
    also when MUSIC could not run, so the line numbers of the two files never drift apart:
    the header of the DFT map is kept and power[] set to the MUSIC floor. */
    if (music) {
      nr_sensing_map_t *map_music = malloc_or_fail(sizeof(*map_music));
      memcpy(map_music, &t->map, sizeof(*map_music));
      nr_music_stats_t ms = {0};
      if (nr_ue_music_doppler(obs, n_obs, map_music, &ms) > 0) {
        LOG_I(NR_PHY,
              "SENSING MUSIC %d.%d ports 0x%x layer %d: %d obs x %d cols per bin, noise ref %.3g, "
              "%d of %d bins with signal, order up to %d\n",
              t->frame,
              t->slot,
              t->stream.ports,
              t->map.stream.layer,
              ms.n_obs,
              ms.n_cols,
              ms.noise_ref,
              ms.bins_with_signal,
              map_music->n_bins,
              ms.order_max);
      } else {
        LOG_W(NR_PHY, "sensing music: no usable observation for map %d.%d, floor written\n", t->frame, t->slot);
        for (int i = 0; i < map_music->n_bins * map_music->n_freq; i++)
          map_music->power[i] = 1.0f;
      }
      char music_path[512];
      snprintf(music_path, sizeof(music_path), "%s.music.csv", t->args->dump);
      nr_ue_sensing_dump_map(music_path, t->frame, t->slot, t->aarx, map_music, markers, n_markers, NULL, 0);
      free(map_music);
    }

    /* Runtime parameter --sensing-clutter-compare: the same map with the kernel fit
    limited to the direct path, built from the same snapshots, to <sensing-dump>.L1.csv.
    Line i of both files is the same window, so --index i compares them. No detector
    or AoA runs on it. */
    if (t->args->clutter_compare && clutter == NR_CLUTTER_KERNEL && t->args->dump[0] != 0) {
      nr_sensing_map_t *map_l1 = malloc_or_fail(sizeof(*map_l1));
      if (nr_ue_sensing_task_map(t, clutter, 1, map_l1, NULL, 0, NULL, NULL) > 0) {
        char l1_path[512];
        snprintf(l1_path, sizeof(l1_path), "%s.L1.csv", t->args->dump);
        nr_ue_sensing_dump_map(l1_path, t->frame, t->slot, t->aarx, map_l1, NULL, 0, NULL, 0);
      }
      free(map_l1);
    }

    /* Slots the window spans. Not the slots that carried DM-RS: PDSCH is only scheduled
    in some of them, so this is always the larger of the two. A slot is 1 ms / 2^mu, and
    2^mu is scs / 15 kHz, so 0.5 ms at 30 kHz. */
    const double slot_dur_s = 1e-3 / ((double)t->snap[0].scs_hz / 15000.0);
    const int n_slots_elapsed = (int)(t->map.t_span_s / slot_dur_s + 0.5);

    LOG_I(NR_PHY,
          "SENSING MAP %d.%d rx%s ports 0x%x layer %d k_step %d: %d snapshots over %.1f ms (# slots %d), "
          "%d bins x %d doppler, %.2f m/bin, +-%.0f Hz, resolution %.2f m/s\n",
          t->frame,
          t->slot,
          rx_label,
          t->map.stream.ports,
          t->map.stream.layer,
          t->map.stream.k_step,
          t->map.n_snapshots,
          t->map.t_span_s * 1e3,
          n_slots_elapsed,
          t->map.n_bins,
          t->map.n_freq,
          t->map.m_per_bin,
          t->map.f_max_hz,
          /* Doppler resolution is 1/t_span Hz; as a speed that is lambda/2 per Hz,
          so it depends on the carrier the map was taken at. */
          (C_M_PER_S / t->map.carrier_hz) / 2.0 / t->map.t_span_s);
    /* One map per (chain, layer) goes into the average. OAI compared the count with
      the chains alone, so with the layers averaged as well every map warned
      "4 of 2". */
    if ((t->aarx < 0 || t->n_layers > 1) && n_combined != t->n_ant * t->n_layers)
      LOG_W(NR_PHY,
            "sensing: averaged map built from %d of %d maps (%d Rx antennas x %d layers)\n",
            n_combined,
            t->n_ant * t->n_layers,
            t->n_ant,
            t->n_layers);
  }
  for (int a = 0; a < n_slow; a++) {
    free(slow[a].t_s);
    free(slow[a].h_re);
    free(slow[a].h_im);
  }
  for (int k = 0; k < n_obs_alloc; k++) {
    free(obs[k].t_s);
    free(obs[k].h_re);
    free(obs[k].h_im);
  }
  free(obs);

  for (int a = 0; a < t->n_ant; a++)
    free(t->snap[a].ring);
  free(t);
}


static void* map_thread(void* arg)
{
  nrscope_sensing_t* s = (nrscope_sensing_t*)arg;
  for (;;) {
    pthread_mutex_lock(&s->q_lock);
    while (s->q_len == 0) {
      pthread_cond_wait(&s->q_cond, &s->q_lock);
    }
    void* t   = s->q[s->q_head];
    s->q_head = (s->q_head + 1) % MAP_QUEUE;
    s->q_len--;
    s->q_busy++;
    pthread_mutex_unlock(&s->q_lock);
    nr_ue_sensing_map_task(t);
    pthread_mutex_lock(&s->q_lock);
    if (--s->q_busy == 0 && s->q_len == 0) {
      pthread_cond_broadcast(&s->q_idle);
    }
    pthread_mutex_unlock(&s->q_lock);
  }
  return NULL;
}

static void map_push(nrscope_sensing_t* s, nr_sensing_map_task_t* t)
{
  pthread_mutex_lock(&s->q_lock);
  if (s->q_len == MAP_QUEUE) {
    pthread_mutex_unlock(&s->q_lock);
    LOG_W(NR_PHY, "sensing: map threads behind, map %d.%d dropped\n", t->frame, t->slot);
    for (int a = 0; a < t->n_ant; a++)
      free(t->snap[a].ring);
    free(t);
    return;
  }
  s->q[(s->q_head + s->q_len) % MAP_QUEUE] = t;
  s->q_len++;
  pthread_cond_signal(&s->q_cond);
  pthread_mutex_unlock(&s->q_lock);
}

void nrscope_sensing_wait_maps(nrscope_sensing_t* s)
{
  if (s == NULL)
    return;
  pthread_mutex_lock(&s->q_lock);
  while (s->q_len > 0 || s->q_busy > 0) {
    pthread_cond_wait(&s->q_idle, &s->q_lock);
  }
  pthread_mutex_unlock(&s->q_lock);
}

/* One history holding the snapshots of several, for OAI's layer-averaged map task,
   which reads every layer from t->snap[0] and tells them apart by stream. Each part
   is laid out oldest to newest, one after the other, and the ring is exactly full
   with its head at 0: a walk back from the head then visits each part newest first,
   which is all nr_ue_sensing_gather() needs, as it filters by stream and then sorts
   by time. Frees the parts' rings. */
static void hist_merge(nr_sensing_history_t* parts, int n_parts, nr_sensing_history_t* out)
{
  int total = 0;
  for (int j = 0; j < n_parts; j++)
    total += parts[j].count;
  *out      = parts[0];
  out->ring = malloc_or_fail((size_t)(total > 0 ? total : 1) * sizeof(*out->ring));
  int k     = 0;
  for (int j = 0; j < n_parts; j++) {
    const nr_sensing_history_t* p     = &parts[j];
    const int                   first = (p->head - p->count + p->depth) % p->depth;
    for (int i = 0; i < p->count; i++)
      out->ring[k++] = p->ring[(first + i) % p->depth];
    free(p->ring);
    pthread_mutex_destroy(&parts[j].lock);
  }
  out->depth = total > 0 ? total : 1;
  out->count = total;
  out->head  = 0;
  pthread_mutex_init(&out->lock, NULL);
}

/* What one chain's rings hold of a stream, copied into one history for the map
   task: the stream's own ring, or, when the map averages layers, every layer's merged
   (hist_merge). Taken unconditionally, without a minimum count: the caller has already
   claimed the window, and a window holding too few snapshots to condition is dropped
   by nr_ue_sensing_range_doppler() on its own (NR_SENSING_MIN_SNAPSHOTS), where the
   count it gathered is the one that matters. */
static bool take_chain(nrscope_sensing_t* s, int aarx, uint16_t ports, const nr_sensing_stream_t* st, int n_layers,
                       bool avg_lay, nr_sensing_history_t* out)
{
  if (!avg_lay) {
    nr_sensing_history_t* h = hist_of(s, aarx, ports, st->layer);
    return h != NULL && nr_ue_sensing_history_take(h, st, 1, 0, out);
  }
  nr_sensing_history_t parts[2];
  if (n_layers > 2) {
    return false;
  }
  for (int j = 0; j < n_layers; j++) {
    nr_sensing_stream_t   sj = *st;
    nr_sensing_history_t* hj = hist_of(s, aarx, ports, j);
    sj.layer                 = (uint8_t)j;
    if (hj == NULL || !nr_ue_sensing_history_take(hj, &sj, 1, 0, &parts[j])) {
      for (int i = 0; i < j; i++) { // the ones already taken
        free(parts[i].ring);
        pthread_mutex_destroy(&parts[i].lock);
      }
      return false;
    }
  }
  hist_merge(parts, n_layers, out);
  return true;
}

static int cmp_u64(const void* a, const void* b)
{
  const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
  return (x > y) - (x < y);
}

/* The keys of one chain's copy: time and layer of every snapshot of this measurement
   up to the trigger, sorted. Time and layer identify a snapshot on every chain alike
   (the chains are sample-aligned, so a symbol has the same t_sample on all of them).
   Returns the count; *keys is the caller's to free. */
static int chain_keys(const nr_sensing_history_t* h, const nr_sensing_stream_t* st, uint64_t t_end, uint64_t** keys)
{
  *keys = malloc_or_fail((size_t)(h->count > 0 ? h->count : 1) * sizeof(**keys));
  int n = 0;
  for (int k = 1; k <= h->count; k++) {
    const nr_sensing_snapshot_t* e = &h->ring[(h->head - k + h->depth) % h->depth];
    if (e->stream.ports == st->ports && e->stream.k_step == st->k_step && e->stream.k_offset == st->k_offset &&
        e->t_sample <= t_end) {
      (*keys)[n++] = (e->t_sample << 8) | e->stream.layer;
    }
  }
  qsort(*keys, (size_t)n, sizeof(**keys), cmp_u64);
  return n;
}

/* Keep, in every chain's copy, only the snapshots all the chains hold, up to the
   slot that triggered the map. Workers decode slots in parallel, so when the last
   chain of slot N triggers, chain 0's ring may already hold slot N+1 from another
   worker, or still miss slot N-1. The transform measures each window back from that
   copy's newest snapshot, so the chains would cover different instants: the AoA then
   refuses them ("rx1 gathered 519 snapshots against 539 on rx0"), and an average
   mixes two windows. The others are marked with ports 0, which no measurement has,
   so the gather passes over them. */
static void chains_common(nr_sensing_history_t* h, int n_ant, const nr_sensing_stream_t* st, uint64_t t_end)
{
  uint64_t* keys[NR_AOA_MAX_ANT] = {NULL};
  int       n_keys[NR_AOA_MAX_ANT] = {0};
  for (int a = 0; a < n_ant; a++) {
    n_keys[a] = chain_keys(&h[a], st, t_end, &keys[a]);
  }
  // keys present on chain 0 and on every other chain
  uint64_t* common   = malloc_or_fail((size_t)(n_keys[0] > 0 ? n_keys[0] : 1) * sizeof(*common));
  int       n_common = 0;
  for (int i = 0; i < n_keys[0]; i++) {
    bool everywhere = true;
    for (int a = 1; a < n_ant && everywhere; a++) {
      everywhere = bsearch(&keys[0][i], keys[a], (size_t)n_keys[a], sizeof(uint64_t), cmp_u64) != NULL;
    }
    if (everywhere) {
      common[n_common++] = keys[0][i];
    }
  }
  for (int a = 0; a < n_ant; a++) {
    for (int k = 0; k < h[a].count; k++) {
      nr_sensing_snapshot_t* e = &h[a].ring[k];
      if (e->stream.ports != st->ports || e->stream.k_step != st->k_step || e->stream.k_offset != st->k_offset) {
        continue;
      }
      const uint64_t key = (e->t_sample << 8) | e->stream.layer;
      if (e->t_sample > t_end || bsearch(&key, common, (size_t)n_common, sizeof(uint64_t), cmp_u64) == NULL) {
        e->stream.ports = 0;
      }
    }
    free(keys[a]);
  }
  free(common);
}

/* Queue a map for this chain and layer when its window is complete. OAI triggered on
   a snapshot count (args.symbols) in pdsch_processing(), which let the traffic decide
   the window and so the resolution; here args.symbols is only the most snapshots one
   map may use.

   Every chain and layer is fed the same reference symbols, so a map that averages
   them waits for the last one: only then does every ring hold this slot. That is
   OAI's my_turn. It relies on a slot's chains being processed in chain order, which
   the caller guarantees (nrscope_sensing_process_grant). */
static void maybe_map(nrscope_sensing_t* s, int aarx, uint16_t ports, int layer, int n_layers,
                      const pilot_lattice_t* lat, uint32_t sfn, uint32_t slot_idx, uint64_t t_now)
{
  const nr_sensing_stream_t st      = nr_ue_sensing_stream(ports, layer, lat);
  const int                 n_rx    = (int)s->nof_antennas;
  const bool                avg_ant = s->args.antenna_avg && n_rx > 1;
  const bool                avg_lay = s->args.layer_avg && n_layers > 1;
  const bool                my_turn = (!avg_ant || aarx == n_rx - 1) && (!avg_lay || layer == n_layers - 1);
  if (!my_turn) {
    return;
  }
  nr_sensing_history_t* hist = hist_of(s, aarx, ports, layer);
  if (hist == NULL || !map_due(s, hist, aarx, ports, layer, t_now)) {
    return;
  }

  const int              n_ant = avg_ant ? n_rx : 1;
  nr_sensing_map_task_t* t     = malloc_or_fail(sizeof(*t) + n_ant * sizeof(t->snap[0]));
  t->args     = &s->args;
  t->aarx     = avg_ant ? -1 : aarx;
  t->frame    = (int)sfn;
  t->slot     = (int)slot_idx;
  t->stream   = st;
  t->n_ant    = n_ant;
  t->n_layers = avg_lay ? n_layers : 1;
  if (avg_lay) {
    t->stream.layer = 0; // the seed call and the snapshot dump use layer 0
  }
  for (int a = 0; a < n_ant; a++) {
    if (!take_chain(s, avg_ant ? a : aarx, ports, &st, n_layers, avg_lay, &t->snap[a])) {
      for (int i = 0; i < a; i++) {
        free(t->snap[i].ring);
        pthread_mutex_destroy(&t->snap[i].lock);
      }
      free(t);
      return;
    }
  }
  if (n_ant > 1) {
    chains_common(t->snap, n_ant, &st, t_now);
  }
  map_push(s, t);
}

int nrscope_sensing_process_grant(nrscope_sensing_t*          s,
                                  nrscope_sensing_scratch_t*  sc,
                                  uint32_t                    aarx,
                                  const cf_t*                 grid,
                                  uint32_t                    n_sc_grid,
                                  const nr_dmrs_placement_t*  place,
                                  const srsran_sch_cfg_nr_t*  cfg,
                                  int                         dci_ports,
                                  uint32_t                    pci,
                                  uint32_t                    sfn,
                                  uint32_t                    slot_idx,
                                  int64_t                     window_shift)
{
  if (s == NULL || sc == NULL || grid == NULL || cfg == NULL || place == NULL || n_sc_grid != sc->n_sc) {
    return s == NULL ? 0 : -1;
  }

  // DM-RS ports from the DCI; only configuration type 1, single symbol, is modelled
  uint16_t ports = 0x1; // DCI 1_0: port 1000
  if (dci_ports >= 0) {
    uint32_t cdm = 0;
    ports        = (cfg->dmrs.type == srsran_dmrs_sch_type_1 && cfg->dmrs.length == srsran_dmrs_sch_len_1)
                       ? nrscope_dmrs_ports_type1_len1((uint32_t)dci_ports, &cdm)
                       : 0;
    if (ports == 0 || cdm != cfg->grant.nof_dmrs_cdm_groups_without_data ||
        (uint32_t)__builtin_popcount(ports) != cfg->grant.nof_layers) {
      return 0;
    }
  }
  nr_dmrs_layout_t lay;
  if (!nr_ue_dmrs_layout(ports, cfg->dmrs.type, &lay)) {
    return 0;
  }

  uint32_t  dmrs_symbols[SRSRAN_DMRS_SCH_MAX_SYMBOLS];
  const int n_dmrs = srsran_dmrs_sch_get_symbols_idx(&cfg->dmrs, &cfg->grant, dmrs_symbols);
  if (n_dmrs <= 0) {
    return 0;
  }
  uint16_t dmrs_mask = 0;
  for (int i = 0; i < n_dmrs; i++) {
    dmrs_mask |= (uint16_t)(1u << dmrs_symbols[i]);
  }

  const uint64_t      slot_abs = slot_abs_of(s, sfn, slot_idx);
  uint64_t            t_sample[NSYMB] = {0};
  srsran_carrier_nr_t carrier         = {.pci = pci};
  for (int i = 0; i < n_dmrs; i++) {
    t_sample[dmrs_symbols[i]] = symbol_time(slot_abs, dmrs_symbols[i], s->samples_per_slot, s->ofdm_size);
  }
  // this slot's newest reference symbol, which is where the map trigger reads the time
  const uint64_t t_newest = t_sample[dmrs_symbols[n_dmrs - 1]];

  int             n_pushed = 0;
  pilot_lattice_t lats[2]  = {{0}};
  for (int layer = 0; layer < lay.n_ports && layer < 2; layer++) {
    memset(sc->valid, 0, (size_t)NSYMB * n_sc_grid * sizeof(bool));
    for (int i = 0; i < n_dmrs; i++) {
      const uint32_t l     = dmrs_symbols[i];
      const uint32_t cinit = nr_ue_dmrs_seed(&carrier, &cfg->dmrs, &cfg->grant, slot_idx, l);
      nr_ue_dmrs_estimate_symbol(&lay, layer, cinit, &cfg->grant, place, &grid[(size_t)l * n_sc_grid], (int)n_sc_grid,
                                 &sc->H[(size_t)l * n_sc_grid], &sc->valid[(size_t)l * n_sc_grid]);
      /* Undo the window's moves. The window is window_shift samples later than at
        start, so every path sits window_shift bins earlier in the profile (one
        sample is one bin: the delay axis spans the 4096-point symbol). A delay d
        is exp(-j 2 pi k d / N) across subcarriers k, so multiplying by
        exp(-j 2 pi k shift / N) puts the paths back. The ramp's origin only adds
        a phase common to the symbol, which the alignment removes anyway. */
      if (s->args.compensate_window_shifts && window_shift != 0) {
        const double step = -2.0 * M_PI * (double)(window_shift % (int64_t)s->ofdm_size) / s->ofdm_size;
        cf_t*        H    = &sc->H[(size_t)l * n_sc_grid];
        const bool*  V    = &sc->valid[(size_t)l * n_sc_grid];
        for (uint32_t k = 0; k < n_sc_grid; k++) {
          if (V[k]) {
            H[k] *= cexpf(I * (float)fmod(step * k, 2.0 * M_PI));
          }
        }
      }
    }
    /* OAI's per-slot stage, unchanged: the DM-RS symbols of this layer become
      delay responses pushed into the history (grants under half the carrier are
      dropped there). The other symbols are masked out. */
    const cf_t(*H_grid)[n_sc_grid] = (const cf_t(*)[n_sc_grid])sc->H;
    const bool(*V_grid)[n_sc_grid] = (const bool(*)[n_sc_grid])sc->valid;
    uint16_t used = 0;
    nr_sensing_history_t* hist = hist_of(s, (int)aarx, ports, layer);
    const int idft = nr_ue_sensing_slot_profile(NSYMB, (int)n_sc_grid, H_grid, V_grid, t_sample, hist,
                                                (nr_sensing_stream_t){.ports = ports, .layer = (uint8_t)layer},
                                                slot_abs, s->args.random_drop, (uint16_t)~dmrs_mask, &lats[layer], &used);
    if (idft > 0) {
      n_pushed += __builtin_popcount(used);
      maybe_map(s, (int)aarx, ports, layer, lay.n_ports, &lats[layer], sfn, slot_idx, t_newest);
    }
  }

  return n_pushed;
}
