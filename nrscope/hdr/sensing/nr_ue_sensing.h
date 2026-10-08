/*
 * SPDX-License-Identifier: LicenseRef-CSSL-1.0
 */

#ifndef __NR_UE_SENSING__H__
#define __NR_UE_SENSING__H__

#include <stdbool.h>
#include <stdint.h>
#include "nrscope/hdr/sensing/sensing_defs.h"

/// Max number of CDM groups in csi_mapping_parms_t (see get_csi_mapping_parms())
#define NR_SENSING_MAX_CSI_CDM 16

/// Max delay taps reported per profile by nr_ue_sensing_find_peaks()
#define NR_SENSING_MAX_PEAKS 8

/* The delay axis is pinned by the product idft_size * k_step, not by the IDFT size
alone.

Bin b is a delay of b / (idft_size * k_step * delta_f), so holding that product
constant makes one bin mean the same delay for every reference signal, whatever its
pilot spacing. DMRS type 1 (k_step 2) then uses a 2048 point IDFT and TRS row 1
(k_step 4) a 1024 point one, and the two profiles land on the same axis with no
resampling. Deriving the size from the pilot count instead would move the axis every
time the scheduler changed the allocation, and the Doppler transform averages bins
across slots.

4096 is the smallest power of two above the largest NR carrier, 3276 subcarriers.
That is what makes idft_size >= n_lattice hold for every comb without a per-signal
check: n_lattice is at most n_sc / k_step, and SPAN / k_step >= n_sc / k_step. */
#define NR_SENSING_DELAY_SPAN 4096

/// IDFT size for a lattice of the given pilot spacing
/// This means that for the bin spacing in meters the formula
/// idft_size * k_step * delta_f is always 4096 * 30 khz, so we keep 2.44m/bin.
#define NR_SENSING_IDFT_SIZE(k_step) (NR_SENSING_DELAY_SPAN / (k_step))

/* Bound on NR_SENSING_IDFT_SIZE(), so profiles can be plain arrays rather than VLAs.
The densest comb NR defines for a single port is DMRS type 1 at k_step 2 (type 2 and
the CSI-RS rows are all sparser), so no lattice can ask for more than SPAN / 2 bins
and every smaller k_step is excluded by construction. */
#define NR_SENSING_MAX_IDFT (NR_SENSING_DELAY_SPAN / 2)

/* A symbol is only used when its grant covers at least this fraction of the
carrier, expressed as a numerator over NR_SENSING_MIN_BAND_DEN. */
#define NR_SENSING_MIN_BAND_NUM 1
#define NR_SENSING_MIN_BAND_DEN 2

/* Random drop: one slot in NR_SENSING_DROP_ONE_IN, drawn independently per slot, loses
one of its reference symbols, itself drawn at random among the slot's DM-RS symbols.
The other slots keep all of theirs.

PDSCH DM-RS sits at the same symbols of every downlink slot (3 per slot here), so the
slow-time sampling repeats a fixed 3-symbol pattern every slot, on top of the TDD
period. Dropping the same symbol everywhere would only move that pattern; drawing
both the slot and the symbol breaks it. At 2 about half the slots lose one symbol,
which keeps 5/6 of the snapshots for 3 DM-RS symbols per slot. Larger keeps more. */
#define NR_SENSING_DROP_ONE_IN 2

/* Slots in one TDD period.

10 slots is 5 ms at 30 kHz subcarrier spacing. Match the gNB configuration.

NR-Scope: a runtime value rather than OAI's constant 10, set from the cell's SIB1
(tdd-UL-DL-ConfigurationCommon) once it is decoded, so the TDD detector uses
the pattern of whatever cell is sniffed (10 slots on the Benetel
cell, 7 D, 1 S, 2 U). 10 until then. */
extern int nr_sensing_tdd_period_slots;
#define NR_SENSING_TDD_PERIOD_SLOTS nr_sensing_tdd_period_slots

typedef struct nr_sensing_history_s nr_sensing_history_t;

/* Describes one CSI-RS resource whose LS channel estimates were exported by
nr_ue_csi_rs_procedures(). The RE mapping is computed inside that function from the
FAPI PDU and is not otherwise visible to
the caller, so it is reported here together with the port stride of the exported
buffer. */
typedef struct {
  /// false when no LS estimate was produced (ZP CSI-RS, or RSRP-only measurement)
  bool valid;
  /// 0: TRS, 1: NZP CSI-RS, 2: ZP CSI-RS
  int csi_type;
  /// number of CSI-RS ports; this is the port stride of the exported LS buffer
  int ports;
  /// number of valid entries in koverline[] / loverline[]
  int size;
  /// max k' : extra subcarriers per CDM group
  int kprime;
  /// max l' : extra OFDM symbols per CDM group
  int lprime;
  /// subcarrier offset within the RB, per CDM group
  int koverline[NR_SENSING_MAX_CSI_CDM];
  /// OFDM symbol index within the slot, per CDM group
  int loverline[NR_SENSING_MAX_CSI_CDM];
  /// first PRB of the resource, on the carrier grid
  int start_rb;
  /// number of PRBs of the resource
  int nr_of_rbs;
  /// 0: dot5 (even RB), 1: dot5 (odd RB), 2: one, 3: three
  int freq_density;
} csi_est_meta_t;

/* The RS's REs of one OFDM symbol sit on a uniform subcarrier lattice: step 2 for
PDSCH DMRS type 1, step 4 for TRS row 1. The lattice is a property of the reference
signal, not of the grant: it is every subcarrier of the carrier congruent to
k_first mod k_step, and a smaller grant only leaves more of its positions
unmeasured. n_lattice counts the positions, n counts the ones actually measured. */
typedef struct {
  /// number of pilots collected
  int n;
  /// number of lattice positions spanned by the carrier grid; n <= n_lattice
  int n_lattice;
  /// grid subcarrier of the first pilot (allocation start plus the DMRS comb offset)
  int k_first;
  /// spacing between consecutive pilots
  int k_step;
} pilot_lattice_t;

/* Receive chains the sensing pipeline handles, which every per-chain array here is
sized by. Four, because that is what the X410 exposes and what the sniffer captures
(NRSCOPE_MAX_RX_ANTENNAS in nrscope_def.h, clamped against this one in
nrscope_sensing_get()). The sensing side is the binding limit of the two: a chain the
capture delivers but this does not cover is simply never gathered.

Defined here, in the lowest header of the sensing library, because both the map
(nr_ue_map.h, the spatial null's weights) and the spatial stage (nr_ue_aoa.h, where it
is NR_AOA_MAX_ANT) size arrays by it, and the map header is included by the other. */
#define NR_SENSING_MAX_RX 4

/* Distinct measurement streams one history is expected to hold at a time. Four is
the realistic set today, rank 1 DM-RS plus the two layers of rank 2 and a spare;
eight leaves room for a rank change to be seen before the older stream ages out.
This is the number of distinct (ports, layer, k_step, k_offset) counters
*/
#define NR_SENSING_MAX_STREAMS 8

/* Which measurement process a snapshot belongs to.

Snapshots may only be transformed together when they are samples of one channel
observed the same way, and the pilot spacing alone does not establish that. Across a
rank change the DM-RS lattice may not move at all, yet layer 0 goes from H w^(1) to
H w_0^(2): a step in the per-path gain, which is indistinguishable from a Doppler
event and comes out as the true line convolved with the spectrum of whatever the
scheduler was doing. The same applies to a port change at fixed rank, and to any
other reference signal that happens to land on the same comb, since an unprecoded
one measures H rather than H w_j.

Keying on the port bitmap and the layer catches all of those: the key changes exactly
when the beam does, which is the property the slow-time record needs. k_step and the
comb offset are carried too, so a stream is self-describing and a lattice that moves
without the ports moving still separates. */
typedef struct {
  /// DM-RS port bitmap of the allocation, as fapi_nr_dl_config_dlsch_pdu_rel15_t carries it
  uint16_t ports;
  /// which of the allocation's layers, 0 for a single-port allocation
  uint8_t layer;
  uint8_t k_step;
  /// comb offset, k_first % k_step
  uint8_t k_offset;
} nr_sensing_stream_t;

/// stream.layer of a map summed over every layer, mirroring aarx = -1 for antennas
#define NR_SENSING_LAYER_AVG 0xff

/// Two streams are the same measurement process when every field agrees
static inline bool nr_ue_sensing_stream_eq(const nr_sensing_stream_t *a, const nr_sensing_stream_t *b)
{
  return a->ports == b->ports && a->layer == b->layer && a->k_step == b->k_step && a->k_offset == b->k_offset;
}

/* Complete a stream key from the lattice its snapshots landed on. The caller knows
which allocation and layer it is estimating; the comb comes from the data. */
static inline nr_sensing_stream_t nr_ue_sensing_stream(uint16_t ports, int layer, const pilot_lattice_t *lat)
{
  return (nr_sensing_stream_t){.ports = ports,
                               .layer = (uint8_t)layer,
                               .k_step = (uint8_t)lat->k_step,
                               .k_offset = (uint8_t)(lat->k_step > 0 ? lat->k_first % lat->k_step : 0)};
}

/* Place the measured REs of one OFDM symbol of H_est on their lattice.

   out[] is indexed by lattice position, (k - k_first % k_step) / k_step, and
   zeroed where nothing was measured. Packing the pilots from index 0 instead would
   leave the response carrying an exp(j*2*pi*k_first*tau) factor, which changes
   whenever the scheduler moves the grant and which the Doppler transform then reads
   as motion.

   Zeroing the unscheduled positions is not the aliasing hazard that zeroing the
   data REs of the slot grid would be: those sit between lattice positions and would
   halve the unambiguous delay by breaking the k_step periodicity, whereas these are
   lattice positions themselves, so they only window the profile.

   n_sc       : number of subcarriers of the carrier grid (12 * N_RB_DL)
   H_row      : H_est[aarx][symbol], indexed by grid subcarrier
   valid_row  : H_valid[symbol], marks the REs that carry a measurement
   max_out    : capacity of out[], must be at least n_sc / k_step
   Returns false if fewer than 2 pilots are present, if they are not on a
   uniform lattice, or if out[] is too small. */

   /* claude: cplusplus needed because of the parameters imscope when we record samples */
   // #ifndef __cplusplus
bool nr_ue_sensing_extract_lattice(int n_sc,
                                   const cf_t H_row[n_sc],
                                   const bool valid_row[n_sc],
                                   int max_out,
                                   cf_t out[max_out],
                                   pilot_lattice_t *lattice);
//#endif

/* Taper the measured span of a lattice vector with a Hann window, in place.

   Without this the IDFT sees a rectangular window: the pilots are 1 over the grant
   and 0 outside it, and multiplying by that rectangle in frequency convolves the
   delay profile with its transform, the Dirichlet kernel. A single path then does
   not give a single peak but a mainlobe plus sidelobes only 13 dB down, decaying as
   1/tau. With the direct path 40 dB above a target, its sidelobes sit 27 dB above
   that target everywhere on the axis, and no amount of averaging removes them:
   the leakage is deterministic, not noise.

   The sidelobes come from the discontinuity at the edges of the rectangle. A window
   going smoothly to zero has none, so its transform decays far faster. Hann puts the
   first sidelobe at -31 dB and the rolloff at -18 dB/octave instead of -6.

   The price is 1.5x the mainlobe width, so the delay resolution reported from
   lattice->n is optimistic by that factor, and 1.76 dB of SNR. The coherent gain is
   0.5, i.e. one bit of the fixed point range, which matters because the response is
   stored as cf_t. That 0.5 is independent of how many pilots were measured, so it
   composes with the pilot count normalisation done per snapshot in
   nr_ue_sensing_range_doppler(): both leave a resolved path at A * const.

   Only the measured span is tapered. The lattice positions left at zero by
   nr_ue_sensing_extract_lattice() are not part of the record and stay zero.

   lattice : as returned by nr_ue_sensing_extract_lattice(), read for k_first,
             k_step and n
   out     : the lattice vector to taper, at least n_lattice entries */
/* Hann coefficients for a record of L points, w[i] = 0.5*(1 - cos(2*pi*i/(L-1))).

   Split out of nr_ue_sensing_apply_hann() because a detector modelling the range
   response needs the same taper the transform was given, and applying it in place
   loses it. Keeping one generator means the model cannot drift from the data.

   w must hold L floats. L < 8 leaves w untouched and returns false: at L = 2 both
   endpoints are zero and the record would disappear. */
bool nr_ue_sensing_hann_coeffs(int L, float *w);

void nr_ue_sensing_apply_hann(const pilot_lattice_t *lattice, cf_t *out);

/* Complex delay response of one symbol: the IDFT over the fast (frequency) axis of
   the lattice-indexed vector from nr_ue_sensing_extract_lattice(), zero-padded from
   n_lattice to idft_size, which only interpolates the delay axis.
   The result is complex on purpose. Its magnitude is the range profile, while its
   phase across symbols is the Doppler, so squaring here would lose half the point. */

   //#ifndef __cplusplus
void nr_ue_sensing_delay_response(int n_lattice, const cf_t H_lattice[n_lattice], int idft_size, cf_t time_out[idft_size]);
//#endif

/* Absolute time of one OFDM symbol on the sample clock, matching the offset
   nr_slot_fep() reads from. Needed because cyclic prefix lengths differ: the
   first symbol of every half subframe is longer. */
uint64_t nr_ue_sensing_symbol_time(int hfn,
                                   int frame,
                                   int slot,
                                   int symbol,
                                   int symbols_per_slot,
                                   int numerology,
                                   uint32_t samples_per_frame,
                                   uint32_t slot_timestamp,
                                   int ofdm_symbol_size,
                                   int nb_prefix_samples,
                                   int nb_prefix_samples0);

/* Non-coherently average the delay profiles of one slot, for one Rx antenna.
   Valid for a static scene only.
   Symbols are grouped by pilot lattice, because DMRS and TRS for example
   produce delay axes with different scales that cannot be summed. One call handles
   the lattice of the first unskipped symbol and averages every symbol matching it.
   Two symbols share a lattice when k_step and the comb offset k_first % k_step
   match; the grant sizes need not, since the lattice indexing makes profiles of
   different allocations share one delay axis.

   Symbols whose grant covers less than NR_SENSING_MIN_BAND_NUM/DEN of the carrier
   are dropped, never averaged and never pushed to hist. Their delay resolution,
   c/(n * k_step * scs), is more than twice the full-band one and their window
   sidelobes are correspondingly worse, so averaging one in degrades the whole
   profile.

   The complex response of each symbol is computed once here and used twice: pushed
   to hist as a slow-time sample for the Doppler axis, and squared into this slot's
   range profile.

   H_est_ant  : H_est[aarx], indexed [symbol][subcarrier]
   t_sample   : absolute time of each symbol, from nr_ue_sensing_symbol_time()
   hist       : history to accumulate into, NULL to only build the range profile
   stream     : which measurement process the grid holds, so the snapshots can be
                told apart later. Only ports and layer are read; k_step and the comb
                offset are filled in from the lattice this pass finds.
   slot_abs   : absolute slot index, read by random_drop and by the grant emulation
                (NR_SENSING_EMULATE_GRANTS in nr_ue_sensing.c, simulation only)
   random_drop : break up the periodicity of the slow-time sampling: in one slot in
                NR_SENSING_DROP_ONE_IN, drawn at random, one DM-RS symbol drawn at
                random is left out of hist. The range profile still averages every
                symbol.
   skip_mask  : bit m set = ignore symbol m on this pass
   lattice    : receives the lattice the group shares, with n set to the smallest
                pilot count averaged, so the reported resolution is the coarsest one
   used_mask  : receives the bitmap of symbols averaged
   Returns NR_SENSING_IDFT_SIZE(k_step) for the lattice handled, or 0 when no
   unskipped symbol carries usable estimates. */
//#ifndef __cplusplus
int nr_ue_sensing_slot_profile(int n_sym,
                               int n_sc,
                               const cf_t H_est_ant[n_sym][n_sc],
                               const bool H_valid[n_sym][n_sc],
                               const uint64_t t_sample[n_sym],
                               nr_sensing_history_t *hist,
                               nr_sensing_stream_t stream,
                               uint64_t slot_abs,
                               bool random_drop,
                               uint16_t skip_mask,
                               pilot_lattice_t *lattice,
                               uint16_t *used_mask);
//#endif
#endif
