#ifndef NRSCOPE_RUN_RECORDER_H
#define NRSCOPE_RUN_RECORDER_H

#include <stdint.h>
#include <string>

#include "srsran/phy/phch/dci_nr.h"
#include "srsran/phy/phch/phch_cfg_nr.h"

/* Per-run CSV recordings for offline statistics, switched on by
  log_config.recording_mode in the yaml. Each run writes, under the project
  root rather than the build tree:

    msg4/msg4_<date>_<time>_pci<N>.csv  one row per RRCSetup decoded
    DCIs/dci_<date>_<time>_pci<N>.csv   one row per DCI found for a known UE

  <date>_<time> is when the run started, so the two files of a run pair up. A
  file is only created with its first row: a run that never sees an RRCSetup,
  or never a UE's DCI, leaves no file of that kind.

  While recording, the terminal keeps only the "Found DCI" and "rrc_setup,
  hooray" lines; the grant dumps and the JSON go to the files instead. */
namespace RunRecorder {

void init(bool enable);
bool enabled();

/* Flush and close the files. Rows are also flushed periodically, since runs
  usually end by being killed. */
void close();

void record_rrc_setup(uint32_t           pci,
                      uint32_t           sfn,
                      uint32_t           slot_idx,
                      uint16_t           tc_rnti,
                      uint16_t           c_rnti,
                      uint32_t           rrc_transaction_id,
                      const char*        dci_str,
                      uint32_t           rrc_offset,
                      const uint8_t*     msg4_bytes,
                      uint32_t           nof_bytes,
                      const std::string& master_cell_group_json);

/* sch_cfg is the PDSCH/PUSCH configuration derived from the DCI, or NULL when
  the decoder did not derive one (downlink formats other than 1_1). dci_str is
  the line printed as "Found DCI", kept whole for fields without a column. */
void record_dci(uint32_t                   pci,
                uint32_t                   sfn,
                uint32_t                   slot_idx,
                bool                       downlink,
                const srsran_dci_ctx_t&    ctx,
                uint32_t                   freq_alloc,
                uint32_t                   time_alloc,
                uint32_t                   mcs,
                uint32_t                   ndi,
                uint32_t                   rv,
                uint32_t                   harq_id,
                uint32_t                   tpc,
                uint32_t                   ports,
                uint32_t                   dmrs_id,
                uint32_t                   srs_request,
                const srsran_sch_cfg_nr_t* sch_cfg,
                const char*                dci_str);

} // namespace RunRecorder

#endif
