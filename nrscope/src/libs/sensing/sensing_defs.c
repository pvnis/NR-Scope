/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "nrscope/hdr/sensing/sensing_defs.h"

/* Cached IDFT plans, one per distinct transform size.
 *
 * Four is more than the pipeline can ask for: sizes come from
 * NR_SENSING_IDFT_SIZE(k_step) and only a couple of combs exist. The table is
 * linear because it is searched a handful of times per slot against at most a
 * few entries, which is cheaper than anything with structure. */
#define SENSING_MAX_IDFT_PLANS 4

/* Smallest transform accepted. Below this a "lattice" is not a measurement, and
NR_SENSING_DELAY_SPAN/k_step cannot reach here for any comb NR defines. */
#define SENSING_MIN_IDFT 64

typedef struct {
  int               size;
  srsran_dft_plan_t plan;
  bool              valid;
} sensing_idft_plan_t;

static sensing_idft_plan_t plans[SENSING_MAX_IDFT_PLANS];
static int                 nof_plans   = 0;
static pthread_mutex_t     plans_mutex = PTHREAD_MUTEX_INITIALIZER;

/* True when a size is worth planning. Powers of two only: FFTW will happily
build a plan for any size, but a delay transform sitting on the receive path
must not inherit the runtime of an awkward factorisation. The delay axis is
pinned by idft_size * k_step, so a rejected size means that comb is unusable
rather than merely slow, which is the honest outcome. */
static bool sensing_idft_size_ok(int size)
{
  return size >= SENSING_MIN_IDFT && (size & (size - 1)) == 0;
}

/// Find an existing plan, or build one. Caller holds plans_mutex.
static srsran_dft_plan_t* sensing_idft_get_plan(int size)
{
  for (int i = 0; i < nof_plans; i++) {
    if (plans[i].valid && plans[i].size == size) {
      return &plans[i].plan;
    }
  }

  if (nof_plans >= SENSING_MAX_IDFT_PLANS) {
    LOG_E(NR_PHY, "sensing: no room for an IDFT plan of size %d, %d already cached\n", size, nof_plans);
    return NULL;
  }

  sensing_idft_plan_t* slot = &plans[nof_plans];
  if (srsran_dft_plan_c(&slot->plan, size, SRSRAN_DFT_BACKWARD) != SRSRAN_SUCCESS) {
    LOG_E(NR_PHY, "sensing: failed to create an IDFT plan of size %d\n", size);
    return NULL;
  }

  /* Unnormalised, matching what the pipeline was written against: the snapshot
  normalisation in nr_ue_sensing_range_doppler() divides by the pilot count
  itself, and a 1/N here would be applied twice. */
  srsran_dft_plan_set_norm(&slot->plan, false);

  slot->size  = size;
  slot->valid = true;
  nof_plans++;

  LOG_D(NR_PHY, "sensing: created IDFT plan of size %d (%d cached)\n", size, nof_plans);
  return &slot->plan;
}

int nr_ue_sensing_idft(int size, const cf_t* in, cf_t* out)
{
  if (in == NULL || out == NULL) {
    return -1;
  }
  if (!sensing_idft_size_ok(size)) {
    LOG_E(NR_PHY, "sensing: IDFT size %d is not a power of two, refusing to plan it\n", size);
    return -1;
  }

  pthread_mutex_lock(&plans_mutex);
  srsran_dft_plan_t* plan = sensing_idft_get_plan(size);
  pthread_mutex_unlock(&plans_mutex);

  if (plan == NULL) {
    return -1;
  }

  /* The transform runs on the caller's buffers, through FFTW's new-array execute
  (srsran_dft_run_c_zerocopy), which is thread safe for a shared plan.

  srsran_dft_run_c() is not: it copies the input into the plan's own buffer,
  transforms there and reads the plan's output buffer, so two workers
  transforming at once overwrite each other mid-transform. That was the case
  here, with the sniffer's workers estimating in parallel, and every snapshot
  came out as a mix of several: a flat floor about 11 dB under the direct path
  and peaks scattered over the first 20 bins.

  The zero-copy call skips srsran_dft_run_c's pre/post copies and its
  normalisation. Both are plain copies for this plan (backward, no mirror, no DC
  offset) and it is unnormalised, so the result is identical. FFTW requires the
  new arrays to share the plan's SIMD alignment; a caller whose buffers do not
  falls back to the plan's buffers under the lock. */
  if ((((uintptr_t)in | (uintptr_t)out) & 31) == 0) {
    srsran_dft_run_c_zerocopy(plan, in, out);
  } else {
    pthread_mutex_lock(&plans_mutex);
    srsran_dft_run_c(plan, in, out);
    pthread_mutex_unlock(&plans_mutex);
  }
  return 0;
}

void nr_ue_sensing_idft_free(void)
{
  pthread_mutex_lock(&plans_mutex);
  for (int i = 0; i < nof_plans; i++) {
    if (plans[i].valid) {
      srsran_dft_plan_free(&plans[i].plan);
      plans[i].valid = false;
    }
  }
  nof_plans = 0;
  pthread_mutex_unlock(&plans_mutex);
}

static int n_warnings = 0;

int nr_sensing_log_warnings(void)
{
  return __atomic_load_n(&n_warnings, __ATOMIC_RELAXED);
}

void nr_sensing_log(const char* colour, const char* comp, const char* fmt, ...)
{
  if (strcmp(colour, "\033[93m") == 0) { // LOG_W
    __atomic_add_fetch(&n_warnings, 1, __ATOMIC_RELAXED);
  }
  /* Formatted whole, then written in one call, so lines from workers logging at
  the same time do not interleave. */
  char    msg[1024];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(msg, sizeof(msg), fmt, ap);
  va_end(ap);
  char head[16];
  snprintf(head, sizeof(head), "[%s]", comp);
  printf("%s%-8s %s%s", colour, head, msg, colour[0] ? "\033[0m" : "");
}
