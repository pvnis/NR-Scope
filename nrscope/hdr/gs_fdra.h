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
 * Each hypothesis is tried with the PDSCH DM-RS scrambled by the PCI and by the UE's own
 * PDCCH DM-RS N_ID (a cell that gives every UE a pdcch-DMRS-ScramblingID often gives it
 * the same scramblingID0): then only the UE's own PRBs carry that DM-RS, loaded cell or not.
 *
 * Limits: in a fully loaded cell with PCI-scrambled DM-RS every PRB carries some UE's
 * DM-RS and the hypotheses cannot be told apart this way (they then never lock); PDSCH in
 * a later slot (k0 > 0) is not where it is checked; any other scramblingID0 is not seen.
 */
#ifndef NRSCOPE_GS_FDRA_H
#define NRSCOPE_GS_FDRA_H

#include "srsran/phy/common/phy_common_nr.h"
#include <cstdint>

#define GS_FDRA_MAX_HYP 128

struct GsFdraHyp {
  uint8_t  off;       // first bit of the field (after identifier, carrier and BWP indicators)
  uint8_t  type;      // 0: type 0 (RBG bitmap), 1: type 1 (RIV), 2: dynamic switch (MSB)
  uint8_t  rbg_cfg;   // 1 or 2 (type 0 and dynamic)
  uint8_t  width;     // bits of the field
  uint16_t nbwp;      // BWP size, PRBs
  uint16_t bwp_start; // BWP start, CRB
  uint8_t  nid_ue;    // PDSCH DM-RS N_ID: 0 the PCI, 1 the UE's PDCCH DM-RS N_ID (a scramblingID0)
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

/* The DM-RS of one grant as it is on air (type 1, single symbol), for a UE whose
   antenna-ports field cannot be read (its RRC is ciphered). On symbol l, over its PRBs
   (prb[c], c from first_crb on the grid):
   occ   : |z0 - z1|^2 / |z0 + z1|^2 over the CDM pairs (k' = 0, 1) of comb 0, after the
           PCI sequence: ~0 for port 1000 alone, >> 1 for port 1001 alone, in between for
           both (rank 2 on one CDM group), the case a port-1000 estimate gets wrong;
   e1    : power on comb 1 against comb 0: ~0 when comb 1 is left empty (CDM groups
           without data 2: DM-RS boosted 3 dB), ~1 when it carries PDSCH;
   coh1  : DM-RS coherence of comb 1 (ports 1002/1003: rank 3-4, CDM groups 2);
   coh0  : DM-RS coherence of comb 0, the same as gs_fdra_score's. */
struct GsDmrsShape {
  float occ, e1, coh1, coh0;
};
void gs_dmrs_shape(const cf_t* sym, uint32_t first_crb, const bool* prb, uint32_t nof_grid_prb, uint32_t slot_idx,
                   uint32_t l, uint32_t n_id, GsDmrsShape* out);

/* Per CRB of one symbol (CRB first_crb + c at grid subcarrier 12 c), DM-RS type 1 with
   N_ID n_id: coherence of comb 0 and comb 1 (|sum z_j z*_j+2| / power, ~1 where that DM-RS
   is), the CDM-pair ratio of comb 0 (gs_dmrs_shape's occ) and the power of each comb. */
void gs_dmrs_rb_map(const cf_t* sym, uint32_t first_crb, uint32_t nof_prb, uint32_t slot_idx, uint32_t l, uint32_t n_id,
                    float* coh0, float* coh1, float* occ, float* p0, float* p1);

/* Port 1000's channel per CDM pair of comb 0 (type 1: pairs at k = 4n, 4n + 2 of each CRB,
   (z0 + z1) / 2 after the sequence), 3 per CRB into p[3 c + n]; and the type 2 DM-RS
   coherence per CRB (CDM group 0 at k = 0, 1, 6, 7) into coh2, for telling the types apart. */
void gs_dmrs_rb_pairs(const cf_t* sym, uint32_t first_crb, uint32_t nof_prb, uint32_t slot_idx, uint32_t l,
                      uint32_t n_id, cf_t* p, float* coh2);

/* What the shape says: the ports seen on comb 0 (bit i = port 1000 + i; 0x1 when only
   1000 is, 0x2 only 1001, 0x3 both), CDM groups without data 2, and the antenna-ports
   value (TS 38.212 table 7.3.1.2.2-1) the grant is despread with: always ports 1000+1001
   on CDM group 0, so that layer 0 is port 1000's channel whatever the rank and every
   grant of a UE lands in the same stream. */
void gs_dmrs_decide(const GsDmrsShape& s, uint16_t* ports, bool* cdm2, int* dci_ports);

#endif
