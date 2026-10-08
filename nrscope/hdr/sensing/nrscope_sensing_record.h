/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */
/*
 * \brief File format of sensing.record_estimates: the DM-RS channel estimates of every
 *        grant, as nrscope_sensing_push_estimates() receives them.
 *
 * The cut is between the receiver and the sensing pipeline. Everything before it
 * (capture, sniffing, DM-RS estimation and despreading) is frozen in the file with
 * whatever the air, the gNB and the radio did to it; everything after it (random
 * drop, alignment, history, maps, detector, AoA) is replayed by sensing_offline
 * through the same code the live receiver runs, so it can be changed and rerun on
 * identical data.
 *
 * Layout, native byte order (both ends are this machine):
 *   nrsr_file_hdr_t, once
 *   per (grant, chain, layer), in the order the live pipeline received them:
 *     nrsr_rec_hdr_t
 *     per DM-RS symbol of dmrs_mask, lowest symbol first:
 *       uint32_t n_valid
 *       uint8_t  valid bitmap, (n_sc_grid + 7) / 8 bytes, bit k of byte k/8 = subcarrier k
 *       cf_t     the n_valid estimates, in subcarrier order
 *
 * File order is arrival order: a slot's chains come in chain order, chain 0 first,
 * which the alignment and the map trigger rely on, exactly as live.
 */
#ifndef NRSCOPE_SENSING_RECORD_H
#define NRSCOPE_SENSING_RECORD_H

#include <stdint.h>

#define NRSR_MAGIC "NRSEST1"
#define NRSR_VERSION 1
#define NRSR_REC_MAGIC 0x4E525352u // "RSRN" read as bytes

typedef struct {
  char     magic[8]; // NRSR_MAGIC, zero padded
  uint32_t version;
  uint32_t nof_antennas;
  uint64_t carrier_hz;
  double   srate_hz;
  uint32_t ofdm_size;
  uint32_t scs_hz;
  uint32_t n_sc_grid;
  uint32_t reserved;
} nrsr_file_hdr_t;

typedef struct {
  uint32_t magic; // NRSR_REC_MAGIC
  uint16_t aarx;
  uint16_t ports; // DM-RS port bitmap, bit i = port 1000 + i
  uint8_t  layer;
  uint8_t  n_layers;
  /* NR_SENSING_TDD_PERIOD_SLOTS when this grant was received: set from SIB1 at run
  time, and read by the map stage (comb removal, TDD detector). */
  uint8_t  tdd_period_slots;
  uint8_t  reserved;
  uint16_t dmrs_mask; // bit l = OFDM symbol l of the slot carries DM-RS
  uint16_t reserved2;
  uint32_t sfn;
  uint32_t slot_idx;
  uint64_t slot_abs; // absolute slot index, unwrapped across the SFN period
} nrsr_rec_hdr_t;

#endif
