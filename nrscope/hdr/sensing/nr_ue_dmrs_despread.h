/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

/*
 * \brief PDSCH DM-RS least-squares estimation for the sensing path, including the
 *        frequency-domain CDM despreading that rank 2 needs.
 *
 * This is the only stage that knows the 38.211 RE mapping. Everything
 * downstream works on a (symbol, subcarrier) grid and a uniform pilot lattice,
 * and never asks which port produced a sample -- which is what let the rest of
 * the pipeline move here from OpenAirInterface5G untouched while this file was
 * rewritten against srsRAN.
 *
 * Two things differ from the OAI original beyond the sample type:
 *
 * The pilot sequence. OAI handed back a random-access array of the whole
 * symbol's DM-RS. srsRAN's generator is a stream: it must be walked forward and
 * advanced explicitly over PRBs the grant does not use. Getting that bookkeeping
 * wrong by a single PRB yields pilots that are valid values in the wrong places,
 * so every estimate becomes noise and nothing reports an error.
 *
 * The grid indexing. OAI's rxdataF is indexed by FFT bin, so every access needed
 * (first_carrier_offset + k) % ofdm_symbol_size to undo the DC wrap. srsRAN's
 * demodulator already returns a grid indexed by subcarrier, so that term is gone
 * rather than merely hidden.
 *
 * Why despreading is needed at all. Above rank 1 the DM-RS ports of the
 * scheduled layers share their resource elements: they are told apart by a
 * frequency shift (the CDM group) and, inside a group, by a +-1 cover w_f across
 * the pilot pair k' = 0, 1. When both ports land in one group, every pilot RE
 * carries h0 +- h1 and a per-RE correlation returns neither channel. On the
 * k_step 2 lattice that superposition reads as h0[i] + (-1)^i h1[i], whose
 * transform is h0's delay profile plus h1's shifted by half the delay axis:
 * layer 1 appears as a mirrored ghost 2.5 km up the range axis rather than as
 * noise.
 *
 * Undoing the cover consumes both REs of a pair and returns one estimate per
 * layer, so each layer ends up on a lattice of step 4 (type 1) or 6 (type 2)
 * instead of 2 or a non-uniform comb. That halves the unambiguous delay window
 * and leaves the delay resolution, which the occupied bandwidth sets, untouched.
 *
 * What the caller gets per layer is still H_eff = H w_j, one column of the
 * precoded channel. The layers are different beams and must never be summed
 * coherently.
 */

#ifndef NRSCOPE_SENSING_DMRS_DESPREAD_H
#define NRSCOPE_SENSING_DMRS_DESPREAD_H

#include <stdbool.h>
#include <stdint.h>

#include "nrscope/hdr/sensing/sensing_defs.h"
#include "srsran/phy/phch/phch_cfg_nr.h"

#ifdef __cplusplus
extern "C" {
#endif

/* DM-RS ports of one allocation the sensing path handles. Rank 2 is the ceiling
because despreading a group of two is the only cover geometry implemented here;
a CDM group never holds more than two single-symbol ports anyway, so going higher
means handling several groups at once, not a wider despread. */
#define NR_DMRS_SENSING_MAX_PORTS 2

/* How the DM-RS ports of one PDSCH allocation sit on the resource grid, reduced
to what the sensing estimator needs. Built once per allocation by
nr_ue_dmrs_layout() and read per symbol. */
typedef struct {
  /// DM-RS ports in the allocation, i.e. the rank
  int n_ports;
  /// port numbers p, where the antenna port is 1000 + p, in layer order
  int port[NR_DMRS_SENSING_MAX_PORTS];
  /// comb shift of each port; equal shifts mean one CDM group
  int delta[NR_DMRS_SENSING_MAX_PORTS];
  /* Second element of the layer's frequency cover, w_f(k' = 1); w_f(0) is +1 for
  every port. Only read when despread is set, because srsRAN's generator emits
  the bare sequence and never folds a cover into it. */
  int wf1[NR_DMRS_SENSING_MAX_PORTS];
  /// srsran_dmrs_sch_type_1 or srsran_dmrs_sch_type_2
  srsran_dmrs_sch_type_t config_type;
  /// pilots one CDM group carries per PRB: 6 for type 1, 4 for type 2
  int n_pilot_rb;
  /// true when the ports share a CDM group and the cover has to be undone
  bool despread;
  /// subcarrier spacing of the lattice each layer lands on, for the caller to log
  int k_step;
} nr_dmrs_layout_t;

/* Describe the DM-RS port structure of one PDSCH allocation.

   dmrs_ports  : the port bitmap the DCI carries, bit i meaning antenna port
                 1000 + i. 0 is the DCI 1_0 case and is read as port 1000 alone.
   config_type : srsran_dmrs_sch_type_1 or srsran_dmrs_sch_type_2

   Returns false, having logged why, for an allocation the sensing path does not
   model: more than NR_DMRS_SENSING_MAX_PORTS ports; two ports in one CDM group
   separated by the time cover rather than the frequency one, i.e. a
   double-symbol DM-RS, which is not despread in time here; and a single type 2
   port, whose four REs per PRB sit at k = 0,1,6,7 and are not a uniform
   lattice. */
bool nr_ue_dmrs_layout(uint16_t dmrs_ports, srsran_dmrs_sch_type_t config_type, nr_dmrs_layout_t* out);

/* Scrambling seed for one DM-RS symbol, per 38.211 7.4.1.1.1.

   Reimplemented rather than called because srsRAN keeps its own copy static
   inside dmrs_sch.c. It reads scrambling_id0/1 when the RRC configured them and
   falls back to the physical cell identity otherwise, which is the case that
   matters for a passive listener: DMRS-DownlinkConfig travels over encrypted
   RRC, so a sniffer normally only has the PCI, and the default is what most
   deployments leave in place.

   slot_idx is the slot within the frame; symbol_idx the OFDM symbol within the
   slot. */
uint32_t nr_ue_dmrs_seed(const srsran_carrier_nr_t*   carrier,
                         const srsran_dmrs_sch_cfg_t* dmrs_cfg,
                         const srsran_sch_grant_nr_t* grant,
                         uint32_t                     slot_idx,
                         uint32_t                     symbol_idx);

/* Fill one OFDM symbol's row of the sensing channel grid with the least-squares
   estimate of one layer.

   The estimate is the raw correlation of the received RE with the regenerated
   conjugate pilot: no interpolation, no smoothing and no delay
   pre-compensation. What the demodulation path produces is tuned for equalising
   data and biases any delay estimate built on it.

   Where despreading applies, the two REs of a pair produce one value, written at
   the k' = 0 RE of that pair. Labelling it there rather than at the pair
   midpoint keeps the placement the same for both configuration types, and costs
   only a phase exp(-j 2 pi df tau) that is common to every snapshot and so
   cancels in the Doppler transform. The delay bin of a resolved path is
   unaffected.

   lay        : from nr_ue_dmrs_layout()
   layer      : which layer to estimate, 0 .. lay->n_ports - 1
   cinit      : from nr_ue_dmrs_seed(), for this symbol
   grant      : read for prb_idx[] and beta_dmrs
   nof_prb    : PRBs of the carrier grid, i.e. how far prb_idx[] is walked
   reference_point_k_rb : DM-RS reference point offset, from srsran_dmrs_sch_cfg_t
   rxF        : the received grid row for this symbol, indexed by subcarrier
   n_sc_grid  : subcarriers of the carrier grid, 12 * nof_prb
   H_row      : receives the estimates, indexed by grid subcarrier
   valid_row  : marks the REs written; may be NULL to leave it untouched, which
                the caller wants when several antennas fill the same mask

   Returns the number of REs written, or a negative value if the layer or the
   layout is not usable. */
int nr_ue_dmrs_estimate_symbol(const nr_dmrs_layout_t*      lay,
                               int                          layer,
                               uint32_t                     cinit,
                               const srsran_sch_grant_nr_t* grant,
                               uint32_t                     nof_prb,
                               uint32_t                     reference_point_k_rb,
                               const cf_t*                  rxF,
                               int                          n_sc_grid,
                               cf_t*                        H_row,
                               bool*                        valid_row);

#ifdef __cplusplus
}
#endif

#endif // NRSCOPE_SENSING_DMRS_DESPREAD_H
