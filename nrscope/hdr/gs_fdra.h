/*
 * GoldenSniffer: the frequency-domain resource assignment (FDRA) of a DCI 1_1 whose RRC
 * configuration is not known (ciphered reconfiguration), found on air.
 *
 * The FDRA field follows the format identifier, the carrier indicator (0 or 3 bits) and
 * the BWP indicator (0 to 2 bits), TS 38.212 7.3.1.2.2. Its width and meaning depend on
 * the BWP (size, start), the allocation type (0: RBG bitmap, 1: RIV, or a dynamic switch)
 * and the RBG configuration, all from RRC. Every combination that fits is a hypothesis;
 * each DL DCI of a UE is decoded under each of them, and the PRBs a hypothesis predicts
 * are checked against the received PDSCH DM-RS of that slot (TS 38.211 7.4.1.1, N_ID =
 * PCI, n_SCID 0, CDM group 0, the first DM-RS symbol from the MIB): a correct hypothesis
 * puts the allocation where DM-RS is, every time; a wrong one lands on empty PRBs some of
 * the time. The hypothesis with the best mean is kept once it stands out.
 *
 * Limits: in a fully loaded cell every PRB carries some UE's DM-RS and the hypotheses
 * cannot be told apart this way (they then never lock); PDSCH in a later slot (k0 > 0)
 * is not where it is checked; DM-RS scrambled with a scramblingID0 is not seen.
 */
#ifndef NRSCOPE_GS_FDRA_H
#define NRSCOPE_GS_FDRA_H

#include "srsran/phy/common/phy_common_nr.h"
#include <cstdint>

#define GS_FDRA_MAX_HYP 64

struct GsFdraHyp {
  uint8_t  off;       // first bit of the field (after identifier, carrier and BWP indicators)
  uint8_t  type;      // 0: type 0 (RBG bitmap), 1: type 1 (RIV), 2: dynamic switch (MSB)
  uint8_t  rbg_cfg;   // 1 or 2 (type 0 and dynamic)
  uint8_t  width;     // bits of the field
  uint16_t nbwp;      // BWP size, PRBs
  uint16_t bwp_start; // BWP start, CRB
};

/* The hypotheses of a cell: the carrier and the initial BWP (SIB1) as BWP candidates */
uint32_t gs_fdra_hypotheses(uint32_t carrier_start_crb, uint32_t carrier_nprb, uint32_t init_bwp_start_crb,
                            uint32_t init_bwp_nprb, GsFdraHyp* out, uint32_t max_out);

/* The PRBs (BWP-relative, prb[0..nbwp)) a DCI's payload assigns under one hypothesis;
   false if the field does not fit the payload or decodes to nothing valid */
bool gs_fdra_parse(const GsFdraHyp& h, const uint8_t* payload, uint32_t nof_bits, bool* prb, uint32_t* nof_prb);

/* Per-CRB DM-RS coherence of one PDSCH symbol of the grid (CRB first_crb at grid
   subcarrier 0): s[c] (complex) and d[c] (power) of the pilot pairs inside CRB
   first_crb + c, for type 1 CDM group 0 with N_ID n_id. coherence of a set of CRBs =
   |sum s| / sum d; a CRB with |s|/d above ~0.5 carries that DM-RS. */
void gs_fdra_rb_coherence(const cf_t* sym, uint32_t first_crb, uint32_t nof_prb, uint32_t slot_idx, uint32_t l,
                          uint32_t n_id, cf_t* s, float* d);

/* Score of one hypothesis for one DCI: the fraction of its PRBs that carry DM-RS (-1 if
   the hypothesis does not decode) */
float gs_fdra_score(const GsFdraHyp& h, const uint8_t* payload, uint32_t nof_bits, uint32_t first_crb,
                    uint32_t nof_grid_prb, const cf_t* s, const float* d, float* coherence);

#endif
