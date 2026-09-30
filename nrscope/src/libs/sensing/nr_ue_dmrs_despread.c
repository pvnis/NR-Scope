/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#include <math.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_dmrs_despread.h"
#include "nrscope/hdr/dmrs_pilots.h"
#include "srsran/phy/utils/vector.h"

/* Every rejection below describes a configuration the estimator does not model,
so it repeats for as long as the scheduler keeps sending it, i.e. once a slot.
What the reader needs is the reason, once. */
#define WARN_ONCE(...)           \
  do {                           \
    static bool warned_ = false; \
    if (!warned_) {              \
      warned_ = true;            \
      WARNING(__VA_ARGS__);      \
    }                            \
  } while (0)

/* Pilots one CDM group carries per PRB. Type 1 is a comb of 2, so 6 of the 12
REs; type 2 is two adjacent pairs, so 4. */
#define DMRS_PILOTS_PER_PRB(type) ((type) == srsran_dmrs_sch_type_1 ? 6 : 4)

/* Comb shift of a DM-RS port, 38.211 tables 7.4.1.1.2-1 and -2.
 *
 * Type 1 puts ports 0,1 on comb 0 and 2,3 on comb 1, repeating for the
 * time-covered ports 4..7, which is (p >> 1) & 1. Type 2 has three CDM groups at
 * shifts 0, 2 and 4, so 2 * ((p >> 1) % 3). Written as formulas rather than
 * tables because that is all the tables are. */
static int dmrs_delta(int p, srsran_dmrs_sch_type_t type)
{
  if (type == srsran_dmrs_sch_type_1) {
    return (p >> 1) & 1;
  }
  return 2 * ((p >> 1) % 3);
}

/* Second element of the frequency cover, w_f(k'=1). w_f(0) is +1 for every port
in both tables; the odd port of each CDM pair carries -1. */
static int dmrs_wf1(int p)
{
  return (p & 1) ? -1 : 1;
}

bool nr_ue_dmrs_layout(uint16_t dmrs_ports, srsran_dmrs_sch_type_t config_type, nr_dmrs_layout_t* out)
{
  memset(out, 0, sizeof(*out));

  out->config_type = config_type;
  out->n_pilot_rb  = DMRS_PILOTS_PER_PRB(config_type);

  /* DCI 1_0 carries no antenna-port field and the MAC leaves the bitmap at 0,
  which is read as port 1000 alone. */
  if (dmrs_ports == 0) {
    out->n_ports  = 1;
    out->port[0]  = 0;
    out->delta[0] = 0;
    out->wf1[0]   = 1;
  } else {
    for (int i = 0; i < 12; i++) {
      if (!((dmrs_ports >> i) & 0x01)) {
        continue;
      }
      if (out->n_ports == NR_DMRS_SENSING_MAX_PORTS) {
        WARN_ONCE("sensing: DM-RS port bitmap 0x%x holds more than %d ports, not modelled",
                  dmrs_ports,
                  NR_DMRS_SENSING_MAX_PORTS);
        return false;
      }
      out->port[out->n_ports]  = i;
      out->delta[out->n_ports] = dmrs_delta(i, config_type);
      out->wf1[out->n_ports]   = dmrs_wf1(i);
      out->n_ports++;
    }
  }

  // Both ports on the same CDM group, e.g. rank 2, type 1, ports 1000 and 1001
  if (out->n_ports == 2 && out->delta[0] == out->delta[1]) {
    /* Same comb, so the two ports are superimposed on every pilot RE. Inside a
    CDM group the pair the frequency cover separates is (p, p+1) with p even; a
    pair sharing a comb without that relation is separated by the time cover
    instead, i.e. a double-symbol DM-RS. The sequence is regenerated per symbol
    at l' = 0, so the time cover is never undone and such an allocation cannot be
    estimated here. */
    if ((out->port[0] & 1) != 0 || out->port[1] != out->port[0] + 1) {
      WARN_ONCE("sensing: DM-RS ports %d and %d share comb %d but are not a w_f pair, "
                "double-symbol DM-RS is not handled",
                out->port[0],
                out->port[1],
                out->delta[0]);
      return false;
    }
    out->despread = true;
  }

  /* Subcarriers between consecutive estimates of one layer. Despreading emits
  one value per pilot pair, so the step is the pair pitch: 4 for type 1, 6 for
  type 2. Otherwise every pilot RE of the group is its own estimate, which is
  uniform at step 2 for type 1 and not uniform at all for type 2, whose REs sit
  at 0,1,6,7 of the PRB. */
  if (out->despread) {
    out->k_step = (config_type == srsran_dmrs_sch_type_1) ? 4 : 6;
  } else if (config_type == srsran_dmrs_sch_type_1) {
    out->k_step = 2;
  } else {
    WARN_ONCE("sensing: a single type 2 DM-RS port gives REs at 0,1,6,7 of the PRB, "
              "which is not a uniform lattice; allocation skipped");
    return false;
  }

  return true;
}

uint32_t nr_ue_dmrs_seed(const srsran_carrier_nr_t*   carrier,
                         const srsran_dmrs_sch_cfg_t* dmrs_cfg,
                         const srsran_sch_grant_nr_t* grant,
                         uint32_t                     slot_idx,
                         uint32_t                     symbol_idx)
{
  /* The shared generator's N_ID and c_init (dmrs_pilots.h), the same ones the
  on-air DM-RS check uses: scramblingID0/1 when configured, the PCI otherwise. */
  const uint32_t n_id = nrscope_dmrs_n_id(carrier->pci, dmrs_cfg, grant->n_scid);
  return nrscope_dmrs_c_init(slot_idx, symbol_idx, n_id, grant->n_scid);
}

int nr_ue_dmrs_estimate_symbol(const nr_dmrs_layout_t*      lay,
                               int                          layer,
                               uint32_t                     cinit,
                               const srsran_sch_grant_nr_t* grant,
                               const nr_dmrs_placement_t*   place,
                               const cf_t*                  rxF,
                               int                          n_sc_grid,
                               cf_t*                        H_row,
                               bool*                        valid_row)
{
  if (lay == NULL || grant == NULL || place == NULL || rxF == NULL || H_row == NULL) {
    return -1;
  }
  if (layer < 0 || layer >= lay->n_ports) {
    return -1;
  }

  /* Amplitude the pilots are regenerated at. The transmitter sends them at
  M_SQRT1_2 * beta_dmrs, so correlating against M_SQRT1_2 / beta_dmrs returns H
  itself rather than H scaled by the DM-RS power offset. Matching
  srsran_dmrs_sch_estimate(), beta_dmrs of zero means unset and reads as one. */
  float amplitude = M_SQRT1_2;
  if (isnormal(grant->beta_dmrs)) {
    amplitude /= grant->beta_dmrs;
  }

  /* The allocation's last CRB sets how much of the sequence is needed. The
  sequence starts at the reference CRB (point A for a C-RNTI PDSCH) and is
  indexed by absolute position below, never walked: every allocated pilot is
  looked up at its own index, so a gap in the allocation cannot shift anything. */
  const uint32_t grid_prb = (uint32_t)n_sc_grid / SRSRAN_NRE;
  int            crb_last = -1;
  for (uint32_t p = 0; p < SRSRAN_MAX_PRB_NR; p++) {
    if (grant->prb_idx[p]) {
      crb_last = (int)(place->bwp_start_crb + p);
    }
  }
  if (crb_last < 0 || (uint32_t)crb_last < place->reference_crb) {
    return 0;
  }
  const uint32_t n_seq = ((uint32_t)crb_last + 1 - place->reference_crb) * (uint32_t)lay->n_pilot_rb;
  if (n_seq > NRSCOPE_DMRS_MAX_PILOTS) {
    ERROR("sensing: DM-RS sequence of %u pilots exceeds the generator's %d", n_seq, NRSCOPE_DMRS_MAX_PILOTS);
    return -1;
  }
  cf_t pilots[NRSCOPE_DMRS_MAX_PILOTS];
  nrscope_dmrs_sequence(cinit, n_seq, amplitude, pilots);

  const int delta   = lay->delta[layer];
  const int n_pairs = lay->n_pilot_rb / 2; // pairs (k' = 0, 1) of one CDM group per PRB
  int       n_written = 0;

  for (uint32_t p = 0; p < SRSRAN_MAX_PRB_NR; p++) {
    if (!grant->prb_idx[p]) {
      continue;
    }
    const uint32_t crb = place->bwp_start_crb + p;
    if (crb < place->reference_crb || crb < place->grid_crb0 || crb - place->grid_crb0 >= grid_prb) {
      continue; // before the sequence reference or outside the grid
    }
    const uint32_t k_seq  = (crb - place->reference_crb) * SRSRAN_NRE; // PRB start, sequence frame
    const int      k_grid = (int)(crb - place->grid_crb0) * SRSRAN_NRE; // PRB start, grid frame
    for (int np = 0; np < n_pairs; np++) {
      /* The two REs of pair np, k = 4n' + 2k' + delta (type 1) or 6n' + k' + delta (type 2). */
      const int kin0 = (lay->config_type == srsran_dmrs_sch_type_1 ? 4 * np : 6 * np) + delta;
      const int kin1 = kin0 + (lay->config_type == srsran_dmrs_sch_type_1 ? 2 : 1);
      const int m0   = nrscope_dmrs_index_of_k(k_seq + (uint32_t)kin0, lay->config_type, (uint32_t)delta);
      const int m1   = nrscope_dmrs_index_of_k(k_seq + (uint32_t)kin1, lay->config_type, (uint32_t)delta);
      const int k0   = k_grid + kin0;
      const int k1   = k_grid + kin1;
      if (m0 < 0 || m1 < 0 || k1 >= n_sc_grid) {
        continue;
      }

      /* Least squares: the received RE against the conjugate of the known
      pilot. The pilot has unit modulus once the amplitude above is folded in,
      so this is a division without the divide. */
      const cf_t z0 = rxF[k0] * conjf(pilots[m0]);
      const cf_t z1 = rxF[k1] * conjf(pilots[m1]);
      if (!lay->despread) {
        /* One port on this comb, so each RE is an estimate on its own and both
        REs of the pair are kept. */
        H_row[k0] = z0;
        H_row[k1] = z1;
        if (valid_row) {
          valid_row[k0] = true;
          valid_row[k1] = true;
        }
        n_written += 2;
        continue;
      }

      /* Two ports on this comb. With the bare sequence, z0 = h0 + h1 and
      z1 = h0 - h1, so the matched filter of the cover is (z0 + w_f(1) z1) / 2.
      The 2x2 mixing it inverts is exactly unitary, which is why the residual
      taper and the cross-layer leak it leaves behind cancel identically once the
      layers are combined as a power sum, and only bite when a layer is used on
      its own. */
      H_row[k0] = (z0 + (float)lay->wf1[layer] * z1) * 0.5f;
      if (valid_row) {
        valid_row[k0] = true;
      }
      n_written++;
    }
  }
  return n_written;
}
