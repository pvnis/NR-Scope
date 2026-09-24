/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#include <math.h>
#include <string.h>

#include "nrscope/hdr/sensing/nr_ue_dmrs_despread.h"
#include "srsran/phy/common/sequence.h"
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
  /* N_ID defaults to the physical cell identity. scrambling_id0/1 override it
  when DMRS-DownlinkConfig set them, which a passive listener usually cannot
  read: that IE travels over encrypted RRC. The default is therefore the case
  that matters here, and it is also what most deployments leave in place. */
  uint32_t n_id   = carrier->pci;
  uint32_t n_scid = (grant->n_scid) ? 1 : 0;
  if (!grant->n_scid && dmrs_cfg->scrambling_id0_present) {
    n_id = dmrs_cfg->scrambling_id0;
  } else if (grant->n_scid && dmrs_cfg->scrambling_id1_present) {
    n_id = dmrs_cfg->scrambling_id1;
  }

  /* 38.211 7.4.1.1.1. The floor(lambda/2) term of the specification is omitted,
  matching srsran_dmrs_sch_seed(): it is zero for CDM groups 0 and 1, and only
  group 2 of type 2 -- ports 4,5,10,11, which need a rank this path already
  refuses -- would see it. */
  return SRSRAN_SEQUENCE_MOD(
      (((uint64_t)(SRSRAN_NSYMB_PER_SLOT_NR * slot_idx + symbol_idx + 1UL) * (2UL * n_id + 1UL)) << 17UL)
      + (2UL * n_id + n_scid));
}

/* One contiguous run of allocated PRBs, estimated in one pass.
 *
 * The pilots of the run are generated in a single call, because the generator is
 * a stream and the caller has already advanced it over whatever was skipped.
 * Within a PRB the sequence order is (n', k') with k' fastest, which is the
 * order srsran_dmrs_get_pilots_type{1,2}() walk, so the two orders stay in step
 * without either having to know about the other. */
static int estimate_run(const nr_dmrs_layout_t*  lay,
                        int                      layer,
                        srsran_sequence_state_t* state,
                        uint32_t                 prb_start,
                        uint32_t                 prb_count,
                        float                    amplitude,
                        const cf_t*              rxF,
                        int                      n_sc_grid,
                        cf_t*                    H_row,
                        bool*                    valid_row)
{
  const int delta      = lay->delta[layer];
  const int n_pilot_rb = lay->n_pilot_rb;
  /* Pairs of one CDM group per PRB: the range of n' in k = 4n' + 2k' + delta
  (type 1) or 6n' + k' + delta (type 2). */
  const int n_pairs  = n_pilot_rb / 2;
  const int n_pilots = (int)prb_count * n_pilot_rb;

  /* Worst case is the whole carrier in one run: 275 PRB x 6 pilots. Stack rather
  than a scratch buffer on the layout, because several antennas estimate the same
  symbol concurrently. */
  cf_t pilots[275 * 6];
  if (n_pilots > (int)(sizeof(pilots) / sizeof(pilots[0]))) {
    ERROR("sensing: %d pilots in one run exceeds the scratch buffer", n_pilots);
    return -1;
  }

  /* gen_f writes one float per bit and a QPSK pilot is two bits, so the length
  is twice the pilot count. The amplitude is applied here, which is what makes
  the correlation below an estimate of H rather than of H scaled by beta. */
  srsran_sequence_state_gen_f(state, amplitude, (float*)pilots, (uint32_t)n_pilots * 2);

  int n_written = 0;

  for (uint32_t prb = prb_start; prb < prb_start + prb_count; prb++) {
    const int rb_sc   = (int)prb * SRSRAN_NRE;
    const int seq_prb = (int)(prb - prb_start) * n_pilot_rb;

    for (int np = 0; np < n_pairs; np++) {
      /* RE of each element of the pair, and where its pilot sits in the run. */
      const int k0 = rb_sc + (lay->config_type == srsran_dmrs_sch_type_1 ? 4 * np + 0 : 6 * np + 0) + delta;
      const int k1 = rb_sc + (lay->config_type == srsran_dmrs_sch_type_1 ? 4 * np + 2 : 6 * np + 1) + delta;
      const int i0 = seq_prb + 2 * np;
      const int i1 = i0 + 1;

      if (k0 >= n_sc_grid || k1 >= n_sc_grid) {
        continue;
      }

      /* Least squares: the received RE against the conjugate of the known
      pilot. The pilot has unit modulus once the amplitude above is folded in,
      so this is a division without the divide. */
      const cf_t z0 = rxF[k0] * conjf(pilots[i0]);

      if (!lay->despread) {
        /* One port on this comb, so each RE is an estimate on its own and both
        REs of the pair are kept. */
        const cf_t z1 = rxF[k1] * conjf(pilots[i1]);
        H_row[k0]     = z0;
        H_row[k1]     = z1;
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
      const cf_t z1 = rxF[k1] * conjf(pilots[i1]);
      H_row[k0]     = (z0 + (float)lay->wf1[layer] * z1) * 0.5f;
      if (valid_row) {
        valid_row[k0] = true;
      }
      n_written++;
    }
  }

  return n_written;
}

int nr_ue_dmrs_estimate_symbol(const nr_dmrs_layout_t*      lay,
                               int                          layer,
                               uint32_t                     cinit,
                               const srsran_sch_grant_nr_t* grant,
                               uint32_t                     nof_prb,
                               uint32_t                     reference_point_k_rb,
                               const cf_t*                  rxF,
                               int                          n_sc_grid,
                               cf_t*                        H_row,
                               bool*                        valid_row)
{
  if (lay == NULL || grant == NULL || rxF == NULL || H_row == NULL) {
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

  srsran_sequence_state_t state = {};
  srsran_sequence_state_init(&state, cinit);

  /* Walk the carrier, estimating each contiguous run of allocated PRBs and
  advancing the generator across the gaps. The generator has no random access, so
  a PRB that is skipped must still be stepped over: its pilots exist in the
  sequence whether or not this grant uses them.

  reference_point_k_rb shifts where the sequence is taken to start, which a PDSCH
  carrying SIB1 needs; it is subtracted from the first skip and then spent. */
  uint32_t prb_count = 0;
  uint32_t prb_start = 0;
  uint32_t prb_skip  = 0;
  int      n_written = 0;

  for (uint32_t prb_idx = 0; prb_idx < nof_prb; prb_idx++) {
    if (grant->prb_idx[prb_idx]) {
      if (prb_count == 0) {
        prb_start = prb_idx;

        const uint32_t skip = (prb_skip > reference_point_k_rb) ? prb_skip - reference_point_k_rb : 0;
        srsran_sequence_state_advance(&state, skip * (uint32_t)lay->n_pilot_rb * 2);
        prb_skip = 0;
      }
      prb_count++;
      continue;
    }

    prb_skip++;
    if (prb_count == 0) {
      continue;
    }

    const int n = estimate_run(lay, layer, &state, prb_start, prb_count, amplitude, rxF, n_sc_grid, H_row, valid_row);
    if (n < 0) {
      return -1;
    }
    n_written += n;
    prb_count = 0;
  }

  if (prb_count > 0) {
    const int n = estimate_run(lay, layer, &state, prb_start, prb_count, amplitude, rxF, n_sc_grid, H_row, valid_row);
    if (n < 0) {
      return -1;
    }
    n_written += n;
  }

  return n_written;
}
