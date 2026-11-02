/*
 * SPDX-License-Identifier: MIT
 *
 * GoldenSniffer port: internals shared by the gs_*.c files.
 */
#ifndef NRSCOPE_GS_INTERNAL_H
#define NRSCOPE_GS_INTERNAL_H

#include "nrscope/hdr/golden_sniffer/gs_core.h"
#include <complex.h>

#define GS_NSC_REG 12
#define GS_NSC_CCE 72 /* 6 REGs x 12 subcarriers */
#define GS_MAX_RUNS 256
#define GS_MAX_CAND 1024

void gs_gold_init(void);

/* GF(2) tables (gs_gf2.c), built once for all threads */
int  gs_gf2_init(void);
bool gs_gf2_size_valid(uint32_t nof_bits, uint32_t L_log2);
void gs_gf2_prepare(uint32_t nof_bits); /* builds the tables of one size for all ALs */

typedef struct {
  int32_t first, last;
  float   pow;
} gs_run_t;

struct gs_s {
  gs_args_t args;
  uint32_t  sizes[GS_MAX_DCI_BITS + 1];
  uint32_t  nof_sizes;
  uint32_t  size_count[GS_MAX_DCI_BITS + 1]; /* DCI_count of the MATLAB, statistics only */
  bool      css_ok[GS_MAX_DCI_BITS + 1];     /* sizes accepted with n_RNTI = 0 */

  /* interleaved CORESETs: distance between the two halves of a candidate, learned
   * from decoded DCIs, to cut candidates out of runs that merged with neighbours */
  int32_t  il_spacing[4];
  uint32_t il_hits[4];
  uint32_t nof_il;

  float     amp[GS_MAX_SC];
  cf_t      h_amb[GS_MAX_DURATION][GS_MAX_SC]; /* channel up to a j^k per run */
  gs_run_t  runs[GS_MAX_RUNS];
  uint32_t  nof_runs;
  gs_cand_t cand[GS_MAX_CAND];
  uint32_t  nof_cand;

  /* scratch */
  cf_t    tmp[GS_MAX_SC];
  cf_t    cdm[GS_MAX_DURATION][GS_MAX_SC / 4];                 /* DM-RS REs / ambiguous channel */
  cf_t    yeq_iv[GS_MAX_IV][GS_MAX_DURATION * GS_MAX_SC * 3 / 4]; /* equalized data REs per interval */
  uint8_t bits[108 << GS_MAX_L_LOG2];
};

#endif /* NRSCOPE_GS_INTERNAL_H */
