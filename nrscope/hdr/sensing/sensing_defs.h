/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * \brief Primitives the sensing pipeline needs from its host.
 *
 * The pipeline was written inside OpenAirInterface5G and moved here whole. This
 * header is what replaces the handful of OAI facilities it used: a complex
 * sample type, an IDFT, logging, allocation and an assert. It is not a
 * portability layer -- nothing here maps back to OAI, and the OAI copy is a
 * frozen reference -- it only spells out, in one place, which srsRAN facility
 * stands in for each of them.
 *
 * The one substantive change made during the move is the sample type. OAI works
 * in Q15 fixed point (c16_t, a pair of int16_t); srsRAN works in complex float
 * (cf_t) from the radio down. Carrying the fixed point representation across
 * would have meant choosing a scale at the estimator boundary, where the least
 * significant bits are exactly what survives clutter removal, so the pipeline
 * was converted to cf_t instead. Arithmetic that needed helper functions and
 * explicit shifts in Q15 is now native.
 */

#ifndef NRSCOPE_SENSING_DEFS_H
#define NRSCOPE_SENSING_DEFS_H

#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "srsran/config.h"
#include "srsran/phy/dft/dft.h"
#include "srsran/phy/utils/debug.h"
#include "srsran/phy/utils/vector.h"

#ifdef __cplusplus
extern "C" {
#endif

/* srsRAN offers DEBUG, INFO and ERROR but no warning level. The pipeline uses
warnings for configurations it does not model and skips rather than estimates
wrongly, which is neither an error nor mere information: the reader needs to see
it, but it does not mean the run is broken. */
#ifndef WARNING
#define WARNING(_fmt, ...) fprintf(stderr, "[WARNING]: " _fmt "\n", ##__VA_ARGS__)
#endif

/* Allocation that cannot fail silently. The pipeline allocates its history
rings and map buffers once per stream and has no path to recover from running
out of memory mid-slot, so it stops instead of returning a partial map. */
#define malloc_or_fail(SIZE)                                                                                           \
  ({                                                                                                                   \
    void* ptr_ = malloc(SIZE);                                                                                         \
    if (ptr_ == NULL) {                                                                                                \
      ERROR("sensing: out of memory allocating %zu bytes", (size_t)(SIZE));                                            \
      abort();                                                                                                         \
    }                                                                                                                  \
    ptr_;                                                                                                              \
  })

#define calloc_or_fail(N, SIZE)                                                                                        \
  ({                                                                                                                   \
    void* ptr_ = calloc((N), (SIZE));                                                                                  \
    if (ptr_ == NULL) {                                                                                                \
      ERROR("sensing: out of memory allocating %zu x %zu bytes", (size_t)(N), (size_t)(SIZE));                         \
      abort();                                                                                                         \
    }                                                                                                                  \
    ptr_;                                                                                                              \
  })

/* srsran_assert.h is C++ only, so the pipeline's assertions get a C equivalent
here. These guard invariants the caller cannot recover from, such as a lattice
that does not fit its transform. */
#define AssertFatal(COND, ...)                                                                                         \
  do {                                                                                                                 \
    if (!(COND)) {                                                                                                     \
      ERROR(__VA_ARGS__);                                                                                              \
      abort();                                                                                                         \
    }                                                                                                                  \
  } while (0)

/* Inverse DFT over a pilot lattice.
 *
 * OAI looked a transform up by size from a static table; srsRAN needs a plan
 * object per size, created once and reused. Plans are therefore built on first
 * use and kept for the life of the process, under a lock, because the map tasks
 * that call this run on the worker pool.
 *
 * The sizes actually requested are NR_SENSING_IDFT_SIZE(k_step) = 4096/k_step:
 * 2048 for DM-RS type 1 at rank 1, and 1024 for the rank 2 despread. Type 2 at
 * rank 2 would ask for 682, whose factorisation (2 x 11 x 31) makes it far
 * slower than its size suggests, so it is rejected rather than silently
 * accepted on a real-time path.
 *
 * Returns 0 on success and -1 if no plan could be made for that size, in which
 * case out[] is left untouched.
 */
int nr_ue_sensing_idft(int size, const cf_t* in, cf_t* out);

/// Release every cached IDFT plan. For shutdown and for leak checking.
void nr_ue_sensing_idft_free(void);

#ifdef __cplusplus
}
#endif

#endif // NRSCOPE_SENSING_DEFS_H
