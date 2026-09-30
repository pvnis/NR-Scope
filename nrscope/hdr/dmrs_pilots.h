#ifndef NRSCOPE_DMRS_PILOTS_H
#define NRSCOPE_DMRS_PILOTS_H

/* PDSCH DM-RS pilots, generated from TS 38.211 directly: the one source of
 * pilot values in NR-Scope.
 *
 * Both the on-air DM-RS check (dmrs_check.h) and the sensing estimator
 * (sensing/nr_ue_dmrs_despread.h) take their pilots from here. The check measures
 * every decoded grant's pilots against the received grid on air, so a pilot
 * wrong in value or position shows up there; sharing the generator means what
 * sensing divides by is exactly what that check has passed.
 *
 * Written independently of srsRAN's DM-RS code (dmrs_sch.c), so a
 * misunderstanding the two shared would not hide itself; dmrs_check_test
 * cross-checks them on srsRAN's own mapper.
 *
 * The sequence is generated for the whole carrier from its reference point and
 * indexed by absolute position, never walked. A generator stepped PRB by PRB
 * over the gaps of an allocation fails silently when a skip is off by one: it
 * returns valid pilots in the wrong places, every estimate becomes noise, and
 * nothing reports an error.
 */

#include <stdint.h>

#include "srsran/config.h"
#include "srsran/phy/phch/phch_cfg_nr.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Largest pilot index needed: 275 PRBs of type 1 carry 6 pilots each. */
#define NRSCOPE_DMRS_MAX_PILOTS (275 * 6)

/* N_ID^{n_SCID} for a PDSCH, TS 38.211 7.4.1.1.1: scramblingID0/1 when
 * DMRS-DownlinkConfig sets one for this n_SCID, the physical cell identity
 * otherwise. A passive listener normally only has the PCI, as that IE travels
 * over encrypted RRC; the fallback is what deployments mostly use. */
uint32_t nrscope_dmrs_n_id(uint32_t pci, const srsran_dmrs_sch_cfg_t* dmrs_cfg, bool n_scid);

/* c_init of the DM-RS sequence of one OFDM symbol, TS 38.211 7.4.1.1.1.
 * slot_idx: slot within the frame; symbol_idx: OFDM symbol within the slot.
 * The lambda-bar term is omitted: it is zero for CDM groups 0 and 1 of type 1
 * and for groups 0 and 1 of type 2, and only group 2 of type 2 (ports
 * 1004/1005/1010/1011) would need it; the sensing path refuses those. */
uint32_t nrscope_dmrs_c_init(uint32_t slot_idx, uint32_t symbol_idx, uint32_t n_id, bool n_scid);

/* Pilots r(0) .. r(n-1) of one symbol, TS 38.211 7.4.1.1.1:
 *   r(m) = amplitude * ((1 - 2 c(2m)) + j (1 - 2 c(2m+1)))
 * with c the Gold sequence of 5.2.1 for c_init. Pass amplitude M_SQRT1_2 for the
 * transmitted pilots (unit modulus), or M_SQRT1_2 / beta_dmrs to correlate a
 * received RE into H rather than beta_dmrs * H.
 * n is at most NRSCOPE_DMRS_MAX_PILOTS. */
void nrscope_dmrs_sequence(uint32_t c_init, uint32_t n, float amplitude, cf_t* r);

/* Where pilot index m of one symbol sits, TS 38.211 7.4.1.1.2, for configuration
 * type 1 or 2 and a CDM group's comb shift delta:
 *   type 1: k = 4n + 2k' + delta,  m = 2n + k'
 *   type 2: k = 6n +  k' + delta,  m = 2n + k'
 * k counted from subcarrier 0 of the sequence's reference CRB (point A for a
 * C-RNTI PDSCH). Returns k for m, and nrscope_dmrs_index_of_k() the inverse (or
 * -1 when k carries no pilot of that group). */
uint32_t nrscope_dmrs_k_of_index(uint32_t m, srsran_dmrs_sch_type_t type, uint32_t delta);
int      nrscope_dmrs_index_of_k(uint32_t k, srsran_dmrs_sch_type_t type, uint32_t delta);

#ifdef __cplusplus
}
#endif

#endif
