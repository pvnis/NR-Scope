#ifndef NRSCOPE_UL_DMRS_H
#define NRSCOPE_UL_DMRS_H

/* Uplink DM-RS of the cell's UEs (sensing.ul_dmrs): which way each RNTI's phone is.

  In an uplink slot the gNB is silent, so what the receiver gets there is the UEs' own
  transmissions, with no direct path from the gNB to mask them. A DCI 0_1 says which RNTI
  sends a PUSCH, in which slot (n + K2) and on which PRBs; its DM-RS is regenerated and
  correlated there (nr_ue_dmrs_estimate_symbol, the same estimate as the PDSCH's: with
  transform precoding off the PUSCH DM-RS is the same Gold sequence), on every chain, and the
  chains' 4 x 4 spatial covariance over the grant's REs is written out, one row per grant
  (<dump>_ul.csv, beside the maps). The angle is computed from it offline.

  ul_dmrs_add() is called by the DCI decoders when they find an uplink grant;
  ul_dmrs_process() by a worker on an uplink-only slot. Several workers run at once: the
  grants wait in a shared store until their slot comes. */

#include "nrscope/hdr/dci_decoder.h"

/* A grant found in slot (sfn, slot_idx), for the PUSCH k = cfg.grant.k slots later.
  dci_ports: the DCI's antenna ports field. */
void ul_dmrs_add(uint32_t sfn, uint32_t slot_idx, uint32_t slots_per_frame, uint16_t rnti,
                 const srsran_sch_cfg_nr_t& cfg, uint32_t dci_ports);

/* The grants due in this uplink slot, estimated on d's grid (its input is this slot's samples). */
void ul_dmrs_process(DCIDecoder& d, uint32_t sfn, srsran_slot_cfg_t* slot, uint32_t slots_per_frame);

/* Opens the CSV next to the sensing dump; prints the counts at exit. */
void ul_dmrs_open(const char* map_dump);

#endif // NRSCOPE_UL_DMRS_H
