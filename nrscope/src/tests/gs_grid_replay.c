/*
 * SPDX-License-Identifier: MIT
 *
 * Replays grids dumped with GS_DUMP=<file> (golden_sniffer.cc) through the GoldenSniffer
 * search, for debugging on real captures. Run with GS_DEBUG=1 to see every stage.
 * usage: gs_grid_replay <dump> <pci> <ref_crb0> <ref_cs0> [sfn.slot ...]
 * Environment: GS_CSS=<size,size> (CSS size list), GS_CS0_SC=<CORESET#0 width in subcarriers>,
 * GS_MIN_PARITY, GS_SIZE_MAX.
 */
#include "nrscope/hdr/golden_sniffer/gs_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv)
{
  if (argc < 5) {
    fprintf(stderr, "usage: %s <dump> <pci> <ref_crb0> <ref_cs0> [sfn.slot ...]\n", argv[0]);
    return 1;
  }
  FILE* f = fopen(argv[1], "rb");
  if (f == NULL) {
    perror(argv[1]);
    return 1;
  }
  gs_args_t args;
  gs_args_default(&args);
  if (getenv("GS_MIN_PARITY")) {
    args.min_parity = (uint32_t)atoi(getenv("GS_MIN_PARITY"));
  }
  if (getenv("GS_SIZE_MAX")) {
    args.size_max = (uint32_t)atoi(getenv("GS_SIZE_MAX"));
  }
  uint32_t css[8], ncss = 0;
  if (getenv("GS_CSS")) {
    char* save = NULL;
    char  buf[128];
    snprintf(buf, sizeof(buf), "%s", getenv("GS_CSS"));
    for (char* t = strtok_r(buf, ",", &save); t && ncss < 8; t = strtok_r(NULL, ",", &save)) {
      css[ncss++] = (uint32_t)atoi(t);
    }
    args.css_sizes     = css;
    args.nof_css_sizes = ncss;
  }
  gs_t*     gs = gs_create(&args);
  gs_slot_t s  = {};
  s.pci        = (uint32_t)atoi(argv[2]);
  s.ref_crb0   = atoi(argv[3]);
  s.ref_cs0    = atoi(argv[4]);
  s.cs0_nof_sc = getenv("GS_CS0_SC") ? (uint32_t)atoi(getenv("GS_CS0_SC")) : 0;
  uint32_t hdr[4];
  cf_t*    grid = NULL;
  size_t   cap  = 0;
  long     nslots = 0, ndci = 0;
  while (fread(hdr, sizeof(hdr), 1, f) == 1) {
    size_t n = (size_t)hdr[2] * hdr[3];
    if (n > cap) {
      grid = realloc(grid, n * sizeof(cf_t));
      cap  = n;
    }
    if (fread(grid, sizeof(cf_t), n, f) != n) {
      break;
    }
    if (argc > 5) {
      char want[32];
      snprintf(want, sizeof(want), "%u.%u", hdr[0], hdr[1]);
      int hit = 0;
      for (int i = 5; i < argc; i++) {
        hit |= strcmp(argv[i], want) == 0;
      }
      if (!hit) {
        continue;
      }
      printf("== %s\n", want);
    }
    s.grid     = grid;
    s.nof_sc   = hdr[2];
    s.slot_idx = hdr[1];
    gs_dci_t out[GS_MAX_DCI_PER_SLOT];
    int      k = gs_search_slot(gs, &s, out, GS_MAX_DCI_PER_SLOT);
    nslots++;
    for (int i = 0; i < k; i++) {
      ndci++;
      printf("%u.%u rnti 0x%04x n_id %u n_rnti %d L%u size %u first_sc %d evm %.2f\n", hdr[0], hdr[1], out[i].rnti,
             out[i].n_id, out[i].n_rnti, out[i].L, out[i].nof_bits, out[i].cand.iv[0].first, out[i].evm);
    }
  }
  printf("%ld slots, %ld DCIs\n", nslots, ndci);
  gs_destroy(gs);
  return 0;
}
