#ifndef NRSCOPE_DMRS_CHECK_H
#define NRSCOPE_DMRS_CHECK_H

#include <stdint.h>

#include "srsran/srsran.h"

/* Checks that the PDSCH DM-RS of a decoded grant is where, and what, NR-Scope
 * thinks it is, by regenerating the pilots and testing them against the
 * received resource grid.
 *
 * Everything sensing does starts from dividing the received DM-RS REs by the
 * regenerated pilots. If the symbols, subcarriers, sequence or scrambling ID
 * are wrong, that division yields noise and every range-Doppler map built on it
 * is garbage, with no error anywhere. So this is measured on every grant rather
 * than assumed.
 *
 * The pilots are generated here from TS 38.211 section 7.4.1.1 directly (the
 * Gold sequence of 5.2.1, c_init of 7.4.1.1.1, the mapping of 7.4.1.1.2),
 * independently of srsRAN's own DM-RS code, so a shared misunderstanding would
 * not hide itself. Only the configuration (symbols, scrambling ID, CDM groups)
 * comes from what NR-Scope decoded.
 *
 * The metric is pilot coherence. For each DM-RS symbol, z_m = y_m * conj(r_m)
 * over the grant's pilots in CDM group 0, where y is the received RE and r the
 * regenerated pilot. With the right pilots z_m is the channel, which varies
 * slowly across subcarriers, so neighbouring values agree:
 *
 *   coherence = |sum z_m conj(z_{m+2})| / sum |z_m|^2   (pairs inside the grant)
 *
 * Lag 2 rather than 1: ports 1000 and 1001 share CDM group 0 and differ by the
 * frequency OCC [+1 -1] on alternate pilots, so with two layers z alternates
 * between h0+h1 and h0-h1; every second pilot has the same OCC weight. A
 * constant timing offset only rotates the sum, it does not shrink it.
 *
 * Right pilots give a coherence close to 1 at good SNR. Wrong pilots give the
 * magnitude of a random sum, about 1/sqrt(pairs). Two wrong hypotheses are
 * measured the same way, as controls that must come out low:
 *   - the same sequence on the PDSCH symbols that carry data, not DM-RS;
 *   - the claimed symbols with scrambling ID N_ID + 1.
 */
struct DmrsCheckResult {
  bool     valid = false;       // false when the grant could not be checked
  uint32_t nof_pilots = 0;      // pilots used across the DM-RS symbols
  float    coherence = 0;       // claimed configuration
  float    coherence_data_symbol = 0;  // control: the data symbols of the same PDSCH
  float    coherence_wrong_nid   = 0;  // control: scrambling ID + 1
  float    pilot_epre_db = 0;   // mean pilot RE power, dB (grid units)
  /* Pilot SNR from the coherence, assuming the channel is flat over 4
    subcarriers; delay spread lowers it, so read it as a lower bound. */
  float    snr_db = 0;
};

/* grid: one slot of the full carrier, nof_carrier_prb * 12 subcarriers per
 * symbol, symbol-major, starting at CRB crb_offset (offsetToCarrier). The grant's
 * PRB indices are BWP-relative; bwp_start_crb places them on the carrier.
 * dmrs_symbols / nof_dmrs_symbols and n_id are what NR-Scope would use to
 * regenerate the pilots. */
DmrsCheckResult dmrs_check_pdsch(const cf_t*                  grid,
                                 uint32_t                     nof_carrier_prb,
                                 uint32_t                     crb_offset,
                                 uint32_t                     bwp_start_crb,
                                 const srsran_sch_grant_nr_t& grant,
                                 const uint32_t*              dmrs_symbols,
                                 uint32_t                     nof_dmrs_symbols,
                                 uint32_t                     n_id,
                                 uint32_t                     slot_idx_in_frame);

#endif
